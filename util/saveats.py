
import subprocess
import os
import time
import xarray

from util.mpasmsh import jigsaw_mesh_to_netcdf, inject_edge_tags

from mpas_tools.mesh.conversion import convert
from mpas_tools.io import write_netcdf

from jigsawpy import savevtk


def saveats(path, mprj):

    ttic = time.time()

    print("")
    print("Running MPAS mesh-tools...")

    inject_edge_tags(mprj)

    # adapted from BUILD_MESH.py

    if (mprj.mshID.lower() == "euclidean-mesh"):
        print("Forming mesh_triangles_ats.nc")
        jigsaw_mesh_to_netcdf(
            mesh=mprj,
            on_sphere=False,
            output_name=os.path.join(
                path, "tmp", "mesh_triangles_ats.nc"))

    print("Forming base_mesh_ats.nc")
    write_netcdf(
        convert(xarray.open_dataset(
            os.path.join(
                path, "tmp", "mesh_triangles_ats.nc"))),
        fileName=os.path.join(
            path, "out", "base_mesh_ats.nc"))

    args = ["paraview_vtk_field_extractor.py",
            "--ignore_time",
            "-d", "maxEdges=0",
            "-v", "allOnCells",
            "-f", os.path.join(
                path, "out", "base_mesh_ats.nc"),
            "-o", os.path.join(
                path, "out", "base_mesh_ats_vtk")]
    print("")
    print("running:", " ".join(args))
    subprocess.check_call(args, env=os.environ.copy())

    # export jigsaw msh_t direct to *.vtk as well

    savevtk(os.path.join(
        path, "out", "base_mesh_ats.vtk"), mprj)

    return
