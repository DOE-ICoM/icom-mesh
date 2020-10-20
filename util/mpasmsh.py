
import numpy as np
from jigsawpy import jigsaw_msh_t
from netCDF4 import Dataset as NetCDFFile
from mpas_tools.mesh.creation.util import circumcenter


def inject_edge_tags(mesh):
    """
    Injects "edge" ID-tags onto nodes, flagging MPAS cells associated
    with (Delaunay) edge constraints. 

    Parameters
    ----------
    mesh : jigsaw_msh_t
        A JIGSAW mesh object, modified in-place.    

    """
    # Authors: Darren Engwirda

    if (mesh.edge2 is not None):
        for epos in range(mesh.edge2.size):

            inod = mesh.edge2["index"][epos, 0]
            jnod = mesh.edge2["index"][epos, 1]
            itag = mesh.edge2["IDtag"][epos]

            mesh.point["IDtag"][inod] = max(
                itag, mesh.point["IDtag"][inod])
            
            mesh.point["IDtag"][jnod] = max(
                itag, mesh.point["IDtag"][jnod])

    return


def jigsaw_mesh_to_netcdf(mesh, output_name, on_sphere, sphere_radius=None):
    """
    Converts mesh data defined in JIGSAW format to MPAS NetCDF format

    Parameters
    ----------
    mesh : jigsaw_msh_t
        A JIGSAW mesh object
    output_name: str
        The name of the output file
    on_sphere : bool
        Whether the mesh is spherical or planar
    sphere_radius : float, optional
        The radius of the sphere in meters.  If ``on_sphere=True`` this argument
        is required, otherwise it is ignored.
    """
    # Authors: Phillip J. Wolfram, Matthew Hoffman, Xylar Asay-Davis
    #          and Darren Engwirda

    grid = NetCDFFile(output_name, 'w', format='NETCDF3_CLASSIC')

    # Get dimensions
    # Get nCells
    nCells = mesh.point.shape[0]

    # Get vertexDegree and nVertices
    vertexDegree = 3  # always triangles
    nVertices = mesh.tria3.shape[0]

    if vertexDegree != 3:
        ValueError('This script can only compute vertices with triangular '
                   'dual meshes currently.')

    grid.createDimension('nCells', nCells)
    grid.createDimension('nVertices', nVertices)
    grid.createDimension('vertexDegree', vertexDegree)

    # Create cell variables and sphere_radius
    if mesh.vert3.size > 0:
        xCell_full = mesh.vert3['coord'][:, 0]
        yCell_full = mesh.vert3['coord'][:, 1]
        zCell_full = mesh.vert3['coord'][:, 2]
    else:
        xCell_full = mesh.vert2['coord'][:, 0]
        yCell_full = mesh.vert2['coord'][:, 1]
        zCell_full = np.zeros(nCells, dtype=float)

    for cells in [xCell_full, yCell_full, zCell_full]:
        assert cells.shape[0] == nCells, 'Number of anticipated nodes is ' \
                                         'not correct!'
    if on_sphere:
        grid.on_a_sphere = "YES"
        grid.sphere_radius = sphere_radius
        # convert from km to meters
        xCell_full *= 1e3
        yCell_full *= 1e3
        zCell_full *= 1e3
    else:
        grid.on_a_sphere = "NO"
        grid.sphere_radius = 0.0

    # Create cellsOnVertex
    cellsOnVertex_full = mesh.tria3['index'] + 1
    assert cellsOnVertex_full.shape == (nVertices, vertexDegree), \
        'cellsOnVertex_full is not the right shape!'

    # Create vertex variables
    xVertex_full = np.zeros((nVertices,))
    yVertex_full = np.zeros((nVertices,))
    zVertex_full = np.zeros((nVertices,))

    for iVertex in range(nVertices):
        cell1 = cellsOnVertex_full[iVertex, 0]
        cell2 = cellsOnVertex_full[iVertex, 1]
        cell3 = cellsOnVertex_full[iVertex, 2]

        x1 = xCell_full[cell1 - 1]
        y1 = yCell_full[cell1 - 1]
        z1 = zCell_full[cell1 - 1]
        x2 = xCell_full[cell2 - 1]
        y2 = yCell_full[cell2 - 1]
        z2 = zCell_full[cell2 - 1]
        x3 = xCell_full[cell3 - 1]
        y3 = yCell_full[cell3 - 1]
        z3 = zCell_full[cell3 - 1]

        pv = circumcenter(on_sphere, x1, y1, z1, x2, y2, z2, x3, y3, z3)
        xVertex_full[iVertex] = pv.x
        yVertex_full[iVertex] = pv.y
        zVertex_full[iVertex] = pv.z

    var = grid.createVariable('xCell', 'f8', ('nCells',))
    var[:] = xCell_full
    var = grid.createVariable('yCell', 'f8', ('nCells',))
    var[:] = yCell_full
    var = grid.createVariable('zCell', 'f8', ('nCells',))
    var[:] = zCell_full
    var = grid.createVariable('featureTagCell', 'i4', ('nCells',))
    var[:] = mesh.point['IDtag']
    var = grid.createVariable('xVertex', 'f8', ('nVertices',))
    var[:] = xVertex_full
    var = grid.createVariable('yVertex', 'f8', ('nVertices',))
    var[:] = yVertex_full
    var = grid.createVariable('zVertex', 'f8', ('nVertices',))
    var[:] = zVertex_full
    var = grid.createVariable('featureTagVertex', 'i4', ('nVertices',))
    var[:] = mesh.tria3['IDtag']
    var = grid.createVariable(
        'cellsOnVertex', 'i4', ('nVertices', 'vertexDegree',))
    var[:] = cellsOnVertex_full

    grid.close()
