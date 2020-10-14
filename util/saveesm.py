
import subprocess
import os
import time
import numpy as np
import xarray

from geometric_features import GeometricFeatures

from scipy.spatial import distance
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from mpas_tools.mesh.conversion import convert, mask, cull
from mpas_tools.io import write_netcdf
from mpas_tools.mesh.creation.jigsaw_to_netcdf import jigsaw_to_netcdf
from mpas_tools.ocean.inject_bathymetry import inject_bathymetry
from mpas_tools.ocean.inject_preserve_floodplain import \
    inject_preserve_floodplain
from mpas_tools.ocean.coastline_alteration import \
    widen_transect_edge_masks, add_critical_land_blockages, \
    add_land_locked_cells_to_mask

HERE = os.path.abspath(os.path.dirname(__file__))



#!! these should move to mpas_tools.ocean.coastline_alteration

def subtract_critical_passages(dsMask, dsPassages):
    '''
    Parameters
    ----------
    dsMask : `xarray.Dataset`
        The mask to which critical passages should be added
    dsPassages : `xarray.Dataset`
        The transect masks defining critical passages to be kept open for
        ocean flow

    Returns
    -------
    dsMask : `xarray.Dataset`
        The mask with critical passages included
    '''

    dsMask = dsMask.copy()

    nTransects = dsPassages.sizes["nTransects"]
    for transectIndex in range(nTransects):
        index = dsPassages.transectCellMasks[:, transectIndex] > 0
        dsMask.regionCellMasks[index, 0] = 0

    return dsMask

def mask_reachable_ocean(dsMesh, dsMask, fcSeed):
    '''
    Return a new land mask that ensures all ocean cells are "reachable" from
    (at least) one of the ocean points. Isolated patches of ocean cells will 
    be added to the land mask.

    Parameters
    ----------
    dsMesh : `xarray.Dataset`
        The unculled base mesh object to which the land mask is applied.

    dsMask : `xarray.Dataset`
        The current land mask.

    fcSeed : `geometric_features.FeatureCollection`
        A set of "seed" points associated with the ocean domain. Used to
        determine which cells in the ocean mesh are "unreachable" wrt. the 
        land mask via a flood-fill. 

    Returns
    -------
    dsMask : `xarray.Dataset`
        The updated land mask, with the set of unreachable ocean cells added.
    '''

    dsMask = dsMask.copy()

    cellMasks = np.array(dsMask.regionCellMasks[:, 0])

    nCells = dsMesh.sizes["nCells"]

    # a sparse matrix representation of masked cell-to-cell connectivity
    iNext = +0
    iPtrs = np.zeros(nCells + 1, dtype=int)
    xData = np.zeros(
        nCells * dsMesh.sizes["maxEdges"], dtype=int)
    iCols = np.zeros(
        nCells * dsMesh.sizes["maxEdges"], dtype=int)
    
    countOnCell = np.array(dsMesh.nEdgesOnCell)
    cellsOnCell = np.array(dsMesh.cellsOnCell)

    # add graph edges between adjacent cells as long as they share common 
    # mask entries
    for iCell in range(nCells):
        for iEdge in range(countOnCell[iCell]):
            iPair = cellsOnCell[iCell, iEdge] - 1
            if (cellMasks[iCell] == cellMasks[iPair]):
                xData[iNext] = +1
                iCols[iNext] = iPair
                iNext = iNext + 1

        iPtrs[iCell + 1] = iNext

    xData.resize((iPtrs[nCells]))
    iCols.resize((iPtrs[nCells]))

    spMesh = csr_matrix((xData, iCols, iPtrs), dtype=int)

    # find "connected" parts of masked mesh, as the connected components of 
    # the masked cell-to-cell matrix
    nLabel, labels = connected_components(
        csgraph=spMesh, directed=False, return_labels=True)

    # set seedMasks = +1 for the cells closest to any "seed" points defined 
    # for the ocean
    seedMasks = np.zeros(nCells, dtype=int)  
    cellCoord = \
        np.vstack((dsMesh.lonCell, dsMesh.latCell)).T

    for feature in fcSeed.features:
        point = feature["geometry"]["coordinates"]
        
        point[0] = point[0] * np.pi / +180.
        point[1] = point[1] * np.pi / +180.

        index = distance.cdist([point], cellCoord).argmin()

        if (cellMasks[index] == 0):
            seedMasks[index] = +1

    # reset the land mask: mark all cells in any component that contains an
    # ocean point as 0
    dsMask.regionCellMasks[:, 0] = +1

    for iLabel in range(nLabel):
        labelMasks = seedMasks[labels == iLabel]
        if (np.any(labelMasks > +0)):
            dsMask.regionCellMasks[labels == iLabel, 0] = 0
    
    return dsMask




def saveesm(path, geom, mesh,
            preserve_floodplain=False,
            floodplain_elevation=20.0,
            do_inject_elevation=False,
            elevation_filename="",
            with_cavities=False,
            lat_threshold=43.00,
            with_critical_passages=True):
    """
    SAVEESM: export a jigsaw mesh obj. to MPAS-style output.

    1. Writes "mesh_triangles.nc" and "base_mesh.nc" files.
    2. (Optionally) injects elevation + floddplain data.
    3. (Optionally) export compatible FVCOM and ATS outputs.
    4. Calls MPAS-Tools + Geometric-Data to cull mesh into 
       ocean/land partitions.
    5. Writes "culled_mesh.nc" (ocean) and "invert_mesh.nc"
       (land) MPAS-spec. output files.

    """
    # Authors: Darren Engwirda
    
    ttic = time.time()

    print("")
    print("Running MPAS mesh-tools...")

    # adapted from BUILD_MESH.py

    if (geom.mshID.lower() == "ellipsoid-mesh"):
        print("Forming mesh_triangles.nc")
        jigsaw_to_netcdf(
            on_sphere=True,
            sphere_radius=np.mean(geom.radii),
            msh_filename=os.path.join(
                path, "tmp", "mesh.msh"),
            output_name=os.path.join(
                path, "tmp", "mesh_triangles.nc"))

    if (geom.mshID.lower() == "euclidean-mesh"):
        print("Forming mesh_triangles.nc")
        jigsaw_to_netcdf(
            on_sphere=False,
            msh_filename=os.path.join(
                path, "tmp", "mesh.msh"),
            output_name=os.path.join(
                path, "tmp", "mesh_triangles.nc"))

    print("Forming base_mesh.nc")
    write_netcdf(
        convert(xarray.open_dataset(
            os.path.join(
                path, "tmp", "mesh_triangles.nc"))),
        fileName=os.path.join(
            path, "out", "base_mesh.nc"))

    if do_inject_elevation:
        print("Injecting cell elevations")
        inject_bathymetry(
            mesh_file="base_mesh.nc", 
            bathymetry_file=elevation_filename)

    if preserve_floodplain:
        print("Injecting floodplain flag")
        inject_preserve_floodplain(
            mesh_file="base_mesh.nc",
            floodplain_elevation=floodplain_elevation)

    args = ["paraview_vtk_field_extractor.py",
            "--ignore_time",
            "-l",
            "-d", "maxEdges=0",
            "-v", "allOnCells",
            "-f", os.path.join(
                path, "out", "base_mesh.nc"),
            "-o", os.path.join(
                path, "out", "base_mesh_vtk")]
    print("")
    print("running:", " ".join(args))
    subprocess.check_call(args, env=os.environ.copy())

    # adapted from CULL_MESH.py

    # required for compatibility with MPAS
    netcdfFormat = "NETCDF3_64BIT"

    gf = GeometricFeatures(
        cacheLocation='{}'.format(os.path.join(
            HERE, "..", "data", "geometric_data")))

    # start with the land coverage from Natural Earth
    fcLandCoverage = gf.read(
        componentName="natural_earth", objectType="region",
        featureNames=["Land Coverage"])

    # remove the region south of 60S so we can replace 
    # it based on ice-sheet topography
    fcSouthMask = gf.read(
        componentName="ocean", objectType="region",
        featureNames=["Global Ocean 90S to 60S"])

    fcLandCoverage = \
        fcLandCoverage.difference(fcSouthMask)

    # add land coverage from either the full ice sheet 
    # or just the grounded part
    if with_cavities:
        fcAntarcticLand = gf.read(
            componentName="bedmap2", objectType="region",
            featureNames=["AntarcticGroundedIceCoverage"])
    else:
        fcAntarcticLand = gf.read(
            componentName="bedmap2", objectType="region",
            featureNames=["AntarcticIceCoverage"])

    fcLandCoverage.merge(fcAntarcticLand)

    # save the feature collection to a geojson file
    fcLandCoverage.to_geojson(
        os.path.join(
            path, "tmp", "land_coverage.geojson"))

    # Create the land mask based on the land coverage, 
    # i.e. coastline data.
    dsBaseMesh = xarray.open_dataset(
        os.path.join(path, "out", "base_mesh.nc"))
    dsLandMask = mask(dsBaseMesh, fcMask=fcLandCoverage)

    dsLandMask = add_land_locked_cells_to_mask(
        dsLandMask, dsBaseMesh, 
        latitude_threshold=lat_threshold, nSweeps=20)

    if with_critical_passages:
        # merge transects for critical passages into 
        # critical_passages.geojson
        fcCritPassages = gf.read(
            componentName="ocean", objectType="transect",
            tags=["Critical_Passage"])

        # create masks from the transects
        dsCritPassMask = \
            mask(dsBaseMesh, fcMask=fcCritPassages)

        # Alter critical passages to be at least two 
        # cells wide, to avoid sea ice blockage.
        dsCritPassMask = widen_transect_edge_masks(
            dsCritPassMask, dsBaseMesh, 
            latitude_threshold=lat_threshold)

        dsLandMask = subtract_critical_passages(
            dsLandMask, dsCritPassMask)

        # merge transects for critical land blockages 
        # into critical_land_blockages.geojson
        fcCritBlockages = gf.read(
            componentName="ocean", objectType="transect",
            tags=["Critical_Land_Blockage"])

        # create masks from the transects for critical 
        # land blockages
        dsCritBlockMask = \
            mask(dsBaseMesh, fcMask=fcCritBlockages)

        dsLandMask = add_critical_land_blockages(
            dsLandMask, dsCritBlockMask)

    # create seed points for a flood fill of the ocean
    # use all points in the ocean directory, on the 
    # assumption that they are, in fact *in* the ocean
    fcSeed = gf.read(
        componentName="ocean",
        objectType="point", tags=["seed_point"])

    # update the land mask to ensure all ocean cells really 
    # are "reachable" from the rest of the global ocean 
    dsLandMask = mask_reachable_ocean(
        dsMesh=dsBaseMesh, 
        dsMask=dsLandMask, fcSeed=fcSeed)

    # cull the (ocean) mesh based on the land mask, and a 
    # cull the (land) mesh using the inverse mask 

    if preserve_floodplain:
    # with "preserve_floodplains", the (ocean) mesh will 
    # contain overlap with the (land) mesh, otherwise the 
    # two are "perfectly" disjoint
        dsCulledMesh = cull(
            dsBaseMesh, dsMask=dsLandMask, 
            dsPreserve=dsBaseMesh,
            graphInfoFileName=os.path.join(
                path, "out", "culled_graph.info"))

        dsInvertMesh = cull(
            dsBaseMesh, dsInverse=dsLandMask,
            graphInfoFileName=os.path.join(
                path, "out", "invert_graph.info"))

    else:
        dsCulledMesh = cull(
            dsBaseMesh, dsMask=dsLandMask, 
            graphInfoFileName=os.path.join(
                path, "out", "culled_graph.info"))

        dsInvertMesh = cull(
            dsBaseMesh, dsInverse=dsLandMask,
            graphInfoFileName=os.path.join(
                path, "out", "invert_graph.info"))

    write_netcdf(
        dsCulledMesh, os.path.join(
            path, "out", "culled_mesh.nc"), netcdfFormat)

    write_netcdf(
        dsInvertMesh, os.path.join(
            path, "out", "invert_mesh.nc"), netcdfFormat)

    args = ["paraview_vtk_field_extractor.py",
            "--ignore_time",
            "-d", "maxEdges=",
            "-v", "allOnCells",
            "-f", os.path.join(
                path, "out", "culled_mesh.nc"),
            "-o", os.path.join(
                path, "out", "culled_mesh_vtk")]
    print("")
    print("running", " ".join(args))
    subprocess.check_call(args, env=os.environ.copy())

    args = ["paraview_vtk_field_extractor.py",
            "--ignore_time",
            "-d", "maxEdges=",
            "-v", "allOnCells",
            "-f", os.path.join(
                path, "out", "invert_mesh.nc"),
            "-o", os.path.join(
                path, "out", "invert_mesh_vtk")]
    print("running", " ".join(args))
    subprocess.check_call(args, env=os.environ.copy())

    ttoc = time.time()

    print("CPUSEC =", (ttoc - ttic))

    return
