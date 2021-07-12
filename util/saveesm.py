
from distutils.util import strtobool

import subprocess
import os
import time
import numpy as np
import xarray
import argparse

from geometric_features import GeometricFeatures

from util.mpasmsh import jigsaw_to_netcdf, \
    subtract_critical_passages, mask_reachable_ocean

from mpas_tools.mesh.conversion import convert, mask, cull
from mpas_tools.io import write_netcdf
from mpas_tools.ocean.inject_preserve_floodplain import \
    inject_preserve_floodplain
from mpas_tools.ocean.coastline_alteration import \
    widen_transect_edge_masks, add_critical_land_blockages, \
    add_land_locked_cells_to_mask

HERE = os.path.abspath(os.path.dirname(__file__))


def saveesm(msh_file, out_path="",
            on_a_sphere=True,
            sphere_radius=1.0,
            preserve_floodplain=False,
            floodplain_elevation=20.0,
            inject_elevation=True,
            with_cavities=False,
            lat_threshold=43.00,
            with_critical_passages=True):
    """
    SAVEESM: write a jigsaw mesh to "unified" MPAS-type output.

    1. Writes "mesh_triangles.nc" and "base_mesh.nc" files.
    2. (Optionally) injects elevation + floodplain data.
    3. Calls MPAS-Tools + Geometric-Data to cull mesh into 
       'unified' ocean/land partitions.
    4. Writes "ocn_cull_mesh.nc" (ocean) and "lnd_cull_mesh.nc"
       (land) MPAS-spec. output + graph files.

    """
    # Authors: Darren Engwirda
    
    ttic = time.time()

    print("")
    print("Running MPAS mesh-tools...")

    raise Exception()

    if (on_a_sphere):
        print("Forming mesh_triangles.nc")
        jigsaw_to_netcdf(
            msh_file=msh_file,
            on_sphere=True,
            sphere_radius=sphere_radius,
            output_name=os.path.join(
                out_path, "mesh_triangles.nc"))
    else:
        raise Exception("Planar meshes not supported.")

    print("Forming base_mesh.nc")
    write_netcdf(
        convert(xarray.open_dataset(
            os.path.join(
                out_path, "mesh_triangles.nc"))),
        fileName=os.path.join(out_path, "base_mesh.nc"))

    """
    if inject_elevation:
        print("Injecting cell elevations")
        inject_elevation_data(
            mesh_file=os.path.join(
                out_path, "base_mesh.nc"), elevation_file)
    """

    if preserve_floodplain:
        print("Injecting floodplain flag")
        inject_preserve_floodplain(
            mesh_file=os.path.join(
                out_path, "base_mesh.nc"),
            floodplain_elevation=floodplain_elevation)

    args = ["paraview_vtk_field_extractor.py",
            "--ignore_time",
            "-l",
            "-d", "maxEdges=0",
            "-v", "allOnCells",
            "-f", os.path.join(out_path, "base_mesh.nc"),
            "-o", os.path.join(out_path, "base_mesh_vtk")]
    print("")
    print("running:", " ".join(args))
    subprocess.check_call(args, env=os.environ.copy())

    # adapted from CULL_MESH.py

    # required for compatibility with MPAS
    netcdfFormat = "NETCDF3_64BIT"

    gf = GeometricFeatures(
        cacheLocation="{}".format(os.path.join(
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
        os.path.join(out_path, "land_coverage.geojson"))

    # Create the land mask based on the land coverage, 
    # i.e. coastline data.
    dsBaseMesh = xarray.open_dataset(
        os.path.join(out_path, "base_mesh.nc"))
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
        ocnGraphName = \
            os.path.join(out_path, "ocn_graph.info")
        lndGraphName = \
            os.path.join(out_path, "lnd_graph.info")

        dsOcnMesh = cull(
            dsBaseMesh, dsMask=dsLandMask, 
            dsPreserve=dsBaseMesh,
            graphInfoFileName=ocnGraphName)

        dsLndMesh = cull(
            dsBaseMesh, dsInverse=dsLandMask,            
            graphInfoFileName=lndGraphName)

    else:
        ocnGraphName = \
            os.path.join(out_path, "ocn_graph.info")
        lndGraphName = \
            os.path.join(out_path, "lnd_graph.info")

        dsOcnMesh = cull(
            dsBaseMesh, dsMask=dsLandMask,
            graphInfoFileName=ocnGraphName)

        dsLndMesh = cull(
            dsBaseMesh, dsInverse=dsLandMask,
            graphInfoFileName=lndGraphName)

    write_netcdf(
        dsOcnMesh, os.path.join(
            out_path, "ocn_cull_mesh.nc"), netcdfFormat)

    write_netcdf(
        dsLndMesh, os.path.join(
            out_path, "lnd_cull_mesh.nc"), netcdfFormat)

    args = ["paraview_vtk_field_extractor.py",
            "--ignore_time",
            "-d", "maxEdges=",
            "-v", "allOnCells",
            "-f", os.path.join(out_path, "ocn_cull_mesh.nc"),
            "-o", os.path.join(out_path, "ocn_cull_mesh_vtk")]
    print("")
    print("running", " ".join(args))
    subprocess.check_call(args, env=os.environ.copy())

    args = ["paraview_vtk_field_extractor.py",
            "--ignore_time",
            "-d", "maxEdges=",
            "-v", "allOnCells",
            "-f", os.path.join(out_path, "lnd_cull_mesh.nc"),
            "-o", os.path.join(out_path, "lnd_cull_mesh_vtk")]
    print("running", " ".join(args))
    subprocess.check_call(args, env=os.environ.copy())

    ttoc = time.time()

    print("CPUSEC =", (ttoc - ttic))

    return


if (__name__ == "__main__"):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "--msh-file", dest="msh_file", type=str,
        required=True, help="Path+name to jigsaw mesh file.")

    parser.add_argument(
        "--out-path", dest="out_path", type=str,
        default="",
        required=False, help="Output path of the MPAS mesh.")

    parser.add_argument(
        "--on-a-sphere", dest="on_a_sphere", 
        type=lambda x: bool(strtobool(x)), default=True,
        required=False, help="True if mesh on a sphere.")

    parser.add_argument(
        "--sphere-radius", dest="sphere_radius", type=float,
        default=1.0,
        required=False, help="Radius (m) of spherical mesh.")

    parser.add_argument(
        "--inject-elevation", dest="inject_elevation", 
        type=lambda x: bool(strtobool(x)), default=True,
        required=False, help="True to build elevation data.")

    parser.add_argument(
        "--preserve-floodplain", dest="preserve_floodplain", 
        type=lambda x: bool(strtobool(x)), default=False,
        required=False, help="True to inc. floodplain mask.")

    parser.add_argument(
        "--floodplain-elevation", dest="floodplain_elevation", 
        type=float, default=20.0,
        required=False, help="Elevation of floodplain mask.")

    parser.add_argument(
        "--with-cavities", dest="with_cavities", 
        type=lambda x: bool(strtobool(x)), default=False,
        required=False, help="True to build ice-shelf mask.")

    args = parser.parse_args()

    saveesm(msh_file=args.msh_file, out_path=args.out_path,
            on_a_sphere=args.on_a_sphere,
            sphere_radius=args.sphere_radius,
            preserve_floodplain=args.preserve_floodplain,
            floodplain_elevation=args.floodplain_elevation,
            inject_elevation=args.inject_elevation,
            with_cavities=args.with_cavities)
