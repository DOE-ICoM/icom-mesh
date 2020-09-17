
import os
import math
import time
from importlib import import_module

import numpy as np
import argparse
import jigsawpy


def meshify(mesh_path="mesh/vanilla_100", 
            write_esm=True, write_fvc=True, write_ats=True,
            make_geom=True, make_spac=True, make_init=True, 
            make_opts=True):
    """
    MESHIFY: main call to the ICoM mesh-gen. infrastructure.

    1. Call user-defined MESH-PATH/COMPOSE.py to assemble
       mesh geometry, spacing pattern, initial conditions
       and mesh-gen. parameters.
    2. Call JIGSAW to build the triangulation.
    3. Call MPAS meshtools to make the MPAS mesh data files.
    4. (Optionally) export compatible FVCOM and ATS outputs.

    Authors: Darren Engwirda

    """

#---------- call JIGSAW to build the initial triangular mesh

    class obj(object): pass
    
    mesh_args = obj()
    mesh_args.sphere_radius = +6.371E+003

    make_bool = obj()
    make_bool.geom = make_geom
    make_bool.init = make_init
    make_bool.spac = make_spac
    make_bool.opts = make_opts
    
    runjgsw(mesh_path, mesh_args, make_bool)

#---------- run MPAS meshtools to build MPAS data-structures

   #jigsaw_to_MPAS.jigsaw_to_netcdf(
   #    msh_filename=os.path.join(mesh_path, "tmp", "mesh.msh"),
   #    output_name="mesh_triangles.nc", on_sphere=True)


    return


def runjgsw(mesh_path, mesh_args, make_bool):
    """
    RUNJGSW: main call to JIGSAW to build the triangulation.

    MESH-PATH should point to a user-defined mesh directory,
    containing the COMPOSE.py template. 

    Firstly, MESH-PATH/COMPOSE.py is called to build 
    user-defined geometry, initial conditions, mesh spacing
    and mesh configuration information.      

    The boolean flags MAKE-BOOL control whether mesh
    information is built from scratch, or if an exitsing an
    existing file is to be used. For example, setting 
    MAKE-BOOL.SPAC = FALSE relies on an existing spacing
    pattern be available in MESH-PATH/tmp/.

    This information is written to MESH-PATH/tmp/ to be 
    accessed by subsequent calls to JIGSAW.

    Finally, JIGSAW is run to build the triangular mesh, 
    calling either the multi-level (TETRIS) or single-level 
    (JIGSAW) algorithms. 

    Authors: Darren Engwirda

    """
    
    mesh = jigsawpy.jigsaw_msh_t()

#------------------------------------ setup via user COMPOSE
    
    base = mesh_path.replace(os.path.sep, ".")
    
    if (make_bool.geom): 
        geom = getattr(import_module(
            base + ".compose"), "setgeom")(mesh_args)

    if (make_bool.spac): 
        spac = getattr(import_module(
            base + ".compose"), "setspac")(mesh_args)

    if (make_bool.init): 
        init = getattr(import_module(
            base + ".compose"), "setinit")(mesh_args)

    if (make_bool.opts): 
        opts = getattr(import_module(
            base + ".compose"), "setopts")(mesh_args)

#------------------------------------ setup files for JIGSAW

    opts.geom_file = os.path.join(
        mesh_path, "tmp", "geom.msh")

    opts.jcfg_file = os.path.join(
        mesh_path, "tmp", "opts.jig")

    opts.init_file = os.path.join(
        mesh_path, "tmp", "init.msh")

    opts.hfun_file = os.path.join(
        mesh_path, "tmp", "spac.msh")

    opts.mesh_file = os.path.join(
        mesh_path, "tmp", "mesh.msh")

    opts._vtk_file = os.path.join(      # vis. in paraview
        mesh_path, "out", "mesh.vtk")

    jigsawpy.savemsh(opts.geom_file, geom)
    jigsawpy.savemsh(opts.hfun_file, spac)
    jigsawpy.savemsh(opts.init_file, init)

#------------------------------------ make mesh using JIGSAW

    if (not hasattr(opts, "bisection")):

        opts.bisection = +0

    if (opts.bisection < +0):           # bisect heuristic

        rbar = np.mean(geom.radii)
        hbar = np.mean(spac.value)

        nlev = round(math.log2(
            rbar / math.sin(.4 * math.pi) / hbar)
        )

        nlev = nlev - 1

        ttic = time.time()

        jigsawpy.cmd.tetris(opts, nlev - 0, mesh)

        ttoc = time.time()

        print("CPUSEC =", (ttoc - ttic))
        print("BISECT =", +nlev)

    elif (opts.bisection > +0):         # bisect specified

        nlev = opts.bisection

        ttic = time.time()

        jigsawpy.cmd.tetris(opts, nlev - 0, mesh)

        ttoc = time.time()

        print("CPUSEC =", (ttoc - ttic))
        print("BISECT =", +nlev)

    else:                               # do non-recursive

        ttic = time.time()

        jigsawpy.cmd.jigsaw(opts, mesh)

        ttoc = time.time()

        print("CPUSEC =", (ttoc - ttic))
        
#------------------------------------ check mesh for quality

    cost = jigsawpy.triscr2(            # quality metrics!
        mesh.point["coord"],
        mesh.tria3["index"])

    print("TRISCR =", np.min(cost), np.mean(cost))

    cost = jigsawpy.pwrscr2(
        mesh.point["coord"],
        mesh.power,
        mesh.tria3["index"])

    print("PWRSCR =", np.min(cost), np.mean(cost))

    tbad = jigsawpy.centre2(
        mesh.point["coord"],
        mesh.power,
        mesh.tria3["index"])

    print("OBTUSE =",
          +np.count_nonzero(np.logical_not(tbad)))

    ndeg = jigsawpy.trideg2(
        mesh.point["coord"],
        mesh.tria3["index"])

    print("TOPOL. =",
          +np.count_nonzero(ndeg==+6) / ndeg.size)

#------------------------------------ save mesh for Paraview

    jigsawpy.savevtk(opts._vtk_file, mesh)

    return


if (__name__ == "__main__"): 
    parser = argparse.ArgumentParser(
        description=__doc__, 
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "--mesh_path", dest="mesh_path", type=str,
        required=False, help="Path to user-def. mesh-dir.",
        default=os.path.join("mesh", "vanilla_100"))

    parser.add_argument(
        "--write_esm", dest="write_esm", type=bool,
        required=False, help="TRUE to write ESM data file",
        default=True)

    parser.add_argument(
        "--write_fvc", dest="write_fvc", type=bool,
        required=False, help="TRUE to write FVC data file",
        default=True)

    parser.add_argument(
        "--write_ats", dest="write_ats", type=bool,
        required=False, help="TRUE to write ATS data file",
        default=True)

    parser.add_argument(
        "--make_geom", dest="make_geom", type=bool,
        required=False, help="TRUE to re-build GEOM. data",
        default=True)

    parser.add_argument(
        "--make_spac", dest="make_spac", type=bool,
        required=False, help="TRUE to re-build SPAC. data",
        default=True)

    parser.add_argument(
        "--make_init", dest="make_init", type=bool,
        required=False, help="TRUE to re-build INIT. data",
        default=True)

    parser.add_argument(
        "--make_opts", dest="make_opts", type=bool,
        required=False, help="TRUE to re-build OPTS. data",
        default=True)

    args = parser.parse_args()

    meshify(mesh_path=args.mesh_path,
            write_esm=args.write_esm,
            write_fvc=args.write_fvc,
            write_ats=args.write_ats,
            make_geom=args.make_geom,
            make_spac=args.make_spac,
            make_init=args.make_init,
            make_opts=args.make_opts)
