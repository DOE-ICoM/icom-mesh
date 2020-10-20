
import os
import math
import time
import copy
from importlib import import_module
import numpy as np

import argparse
import jigsawpy

from util.inpoly2 import inpoly2
from util.savefvc import savefvc
from util.saveats import saveats
from util.saveesm import saveesm


def meshify(mesh_path="mesh/vanilla_100",
            idtag_esm=[-1], idtag_fvc=[-1], idtag_ats=[-1],
            projector=[0.0, 0.0],
            make_geom=True, make_spac=True, make_init=True,
            make_opts=True):
    """
    MESHIFY: main call to the ICoM mesh-gen. infrastructure.

    1. Call user-defined MESH-PATH/COMPOSE.py to assemble
       mesh geometry, spacing pattern, initial conditions
       and mesh-gen. parameters.
    2. Call JIGSAW to build the triangulation.
    3. (Optionally) export compatible FVCOM and ATS outputs.
    4. Call MPAS meshtools to make the MPAS mesh data files.

    """
    # Authors: Darren Engwirda

    class obj(object): pass                 # dummy object

    make_bool = obj()
    make_bool.geom = make_geom
    make_bool.init = make_init
    make_bool.spac = make_spac
    make_bool.opts = make_opts

#---------- call JIGSAW to build the initial triangular mesh

    geom, gprj, mesh, mprj = \
        runjgsw(mesh_path, make_bool, projector)

#---------- call utilities to write output for ATS and FVCOM

    if (isinstance(idtag_ats, list) and len(idtag_ats) > 0
            and idtag_ats[0] >= +0):

#-------------------------------------- write output for ATS
        mout = zipmesh(mprj, idtag_ats)

        saveats(mesh_path, mout)

    if (isinstance(idtag_fvc, list) and len(idtag_fvc) > 0
            and idtag_fvc[0] >= +0):

#-------------------------------------- write output for FVC
        mout = zipmesh(mprj, idtag_fvc)

    #   savefvc(mesh_path, mout)

#---------- run MPAS meshtools to build MPAS data-structures

    if (isinstance(idtag_esm, list) and len(idtag_esm) > 0
            and idtag_esm[0] >= +0):

#-------------------------------------- write output for ESM
        saveesm(mesh_path, geom, mesh)

    return


def runjgsw(mesh_path, make_bool, projector):
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
    (JIGSAW) algorithms. Cells are assigned ID-tags via
    the polygon/regions defined in GEOM.BOUNDS.

    Returns full-dimensional and 2d-projected msh_t objects.

    """
    # Authors: Darren Engwirda

    mesh = jigsawpy.jigsaw_msh_t()
    mprj = jigsawpy.jigsaw_msh_t()
    gprj = jigsawpy.jigsaw_msh_t()

#------------------------------------ setup via user COMPOSE

    base = \
        mesh_path.replace(os.path.sep, ".")

    if (make_bool.geom):
        geom = getattr(import_module(
            base + ".compose"), "setgeom")()

    if (make_bool.spac):
        spac = getattr(import_module(
            base + ".compose"), "setspac")()

    if (make_bool.init):
        init = getattr(import_module(
            base + ".compose"), "setinit")()

    if (make_bool.opts):
        opts = getattr(import_module(
            base + ".compose"), "setopts")()

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

    opts.hfun_tags = "precision = 9"    # less float prec.

    jigsawpy.savemsh(opts.geom_file, geom,
                     opts.geom_tags)

    jigsawpy.savemsh(opts.hfun_file, spac,
                     opts.hfun_tags)

    jigsawpy.savemsh(opts.init_file, init,
                     opts.init_tags)

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

#------------------------------------ form local projections

    gprj = copy.deepcopy(geom)          # local 2d objects
    mprj = copy.deepcopy(mesh)

    if (mesh.vert3.size > +0):
        project(geom, mesh, gprj, mprj, projector)

#------------------------------------ assign IDtag's to cell

    if (geom.bound is not None and
            geom.bound.size > +0):      # tags per polygon

        imin = np.amin(geom.bound["IDtag"])
        imax = np.amax(geom.bound["IDtag"])

        for itag in range(
                imin + 0, imax + 1):

            tagcell(geom, mesh, gprj, mprj, itag)

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

    jigsawpy.savevtk(os.path.join(
        mesh_path, "out", "geom.vtk"), geom)
    jigsawpy.savevtk(os.path.join(
        mesh_path, "out", "spac.vtk"), spac)
    jigsawpy.savevtk(os.path.join(
        mesh_path, "out", "init.vtk"), init)
    jigsawpy.savevtk(os.path.join(
        mesh_path, "out", "mesh.vtk"), mesh)

    jigsawpy.savevtk(os.path.join(
        mesh_path, "out", "geom_prj.vtk"), gprj)
    jigsawpy.savevtk(os.path.join(
        mesh_path, "out", "mesh_prj.vtk"), mprj)

    return geom, gprj, mesh, mprj


def project(geom, mesh, gprj, mprj, pmid):
    """
    PROJECT: projection of GEOM/MESH objects to a 2d. plane.

    Modifies GPRJ, MPRJ objects "inplace".

    """
    # Authors: Darren Engwirda

    proj = jigsawpy.jigsaw_prj_t()

    mprj.point = np.full(
        mesh.point.size, +0, dtype=mprj.VERT2_t)

    mprj.point["coord"] = \
        jigsawpy.R3toS2(
            geom.radii, mesh.point["coord"])

    proj.prjID = "stereographic"
    proj.radii = np.mean(geom.radii)
    proj.xbase = pmid[0] * np.pi / +180.
    proj.ybase = pmid[1] * np.pi / +180.

    jigsawpy.project(gprj, proj, "fwd")
    jigsawpy.project(mprj, proj, "fwd")

    return


def tagcell(geom, mesh, gprj, mprj, itag):
    """
    TAGCELL: assign ID-tags to mesh cells based-on polygons
    defined in GEOM.BOUND.

    Modifies MESH, MPRJ objects "inplace".

    """
    # Authors: Darren Engwirda

    this = geom.bound["IDtag"] == itag
    cell = geom.bound["index"][this]

    loop = geom.edge2["index"][cell, :]

#------------------------------------ compute cell "centres"
    ball = jigsawpy.tribal2(
        mprj.point["coord"], mesh.tria3["index"])

#------------------------------------ calc. in-polygon tests
    tmsk, _ = inpoly2(
        ball[:, :-1], gprj.point["coord"], loop)

    vmsk, _ = inpoly2(mprj.point["coord"],
                      gprj.point["coord"], loop)

    vbnd = np.full(
        mesh.point.size, False, dtype=np.bool)

    vbnd[mesh.edge2["index"].flatten()] = True

    vmsk[vbnd] = False

    tbnd = np.logical_and.reduce((      # all vert. on bnd
        vbnd[mesh.tria3["index"][:, 0]],
        vbnd[mesh.tria3["index"][:, 1]],
        vbnd[mesh.tria3["index"][:, 2]]))

    tmsk[np.logical_not(tbnd)] = False

    mask = np.logical_or.reduce((       # vert. or cell in
        vmsk[mesh.tria3["index"][:, 0]],
        vmsk[mesh.tria3["index"][:, 1]],
        vmsk[mesh.tria3["index"][:, 2]],
        tmsk))

#------------------------------------ assign IDtags to cells
    mesh.tria3["IDtag"][mask == 1] = itag
    mprj.tria3["IDtag"][mask == 1] = itag

    return


def zipmesh(mesh, tags):
    """
    ZIPMESH: "zip" a mesh down to just the verts/edges/cells
    associated with the list of ID's in TAGS.

    Returns a new "zipped" msh_t object.

    """
    # Authors: Darren Engwirda

    mout = jigsawpy.jigsaw_msh_t()

    keep_cell = np.full(
        mesh.tria3.size, False, dtype=np.bool)
    keep_edge = np.full(
        mesh.edge2.size, False, dtype=np.bool)
    keep_vert = np.full(
        mesh.point.size, False, dtype=np.bool)

#------------------------------------ "flag" tagged entities
    for itag in tags:

        keep_cell[mesh.tria3["IDtag"] == itag] = True

    keep_vert[mesh.tria3[
        "index"][keep_cell].flatten()] = True

    keep_edge = np.logical_and.reduce((
        keep_vert[mesh.edge2["index"][:, 0]],
        keep_vert[mesh.edge2["index"][:, 1]])
    )

    mout.point = mesh.point[keep_vert]
    mout.edge2 = mesh.edge2[keep_edge]
    mout.tria3 = mesh.tria3[keep_cell]

#------------------------------------ update vertex indexing
    redo = \
        np.zeros(mesh.point.size, dtype=np.int32)

    redo[keep_vert] = \
        np.arange(0, np.count_nonzero(keep_vert))

    mout.edge2["index"] = redo[mout.edge2["index"]]
    mout.tria3["index"] = redo[mout.tria3["index"]]

    return mout


if (__name__ == "__main__"):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "--mesh-path", dest="mesh_path", type=str,
        required=False, help="Path to user-def. mesh-dir.",
        default=os.path.join("mesh", "vanilla_100"))

    parser.add_argument(
        "--IDtag-ESM", dest="idtag_esm",
        type=lambda s: [int(it) for it in s.split(",")],
        required=False, help="Mesh ID tags for ESM output",
        default="-1")

    parser.add_argument(
        "--IDtag-FVC", dest="idtag_fvc",
        type=lambda s: [int(it) for it in s.split(",")],
        required=False, help="Mesh ID tags for FVC output",
        default="-1")

    parser.add_argument(
        "--IDtag-ATS", dest="idtag_ats",
        type=lambda s: [int(it) for it in s.split(",")],
        required=False, help="Mesh ID tags for ATS output",
        default="-1")

    parser.add_argument(
        "--projector", dest="projector",
        type=lambda s: [float(it) for it in s.split(",")],
        required=False, help="[lon, lat] for 2-projection",
        default="-75.00, +39.00")

    parser.add_argument(
        "--make-geom", dest="make_geom", type=bool,
        required=False, help="TRUE to re-build GEOM. data",
        default=True)

    parser.add_argument(
        "--make-spac", dest="make_spac", type=bool,
        required=False, help="TRUE to re-build SPAC. data",
        default=True)

    parser.add_argument(
        "--make-init", dest="make_init", type=bool,
        required=False, help="TRUE to re-build INIT. data",
        default=True)

    parser.add_argument(
        "--make-opts", dest="make_opts", type=bool,
        required=False, help="TRUE to re-build OPTS. data",
        default=True)

    args = parser.parse_args()

    meshify(mesh_path=args.mesh_path,
            idtag_esm=args.idtag_esm,
            idtag_fvc=args.idtag_fvc,
            idtag_ats=args.idtag_ats,
            projector=args.projector,
            make_geom=args.make_geom,
            make_spac=args.make_spac,
            make_init=args.make_init,
            make_opts=args.make_opts)
