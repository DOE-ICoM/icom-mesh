
import os
import time
import copy
import numpy as np
import netCDF4 as nc
import jigsawpy
import mpas_tools.mesh.creation.mesh_definition_tools as mdt

from util.inpoly2 import inpoly2
from util.loadshp import loadshp
from util.utility import addpoly, addline, innerto


"""
DELAWARE-LOCAL-2D: boundary-aligned mesh of Delaware-centric 
watersheds, in a local 2d stereographic plane:

"""

HERE = os.path.abspath(os.path.dirname(__file__))
TDIR = os.path.join(HERE, "tmp")
ODIR = os.path.join(HERE, "out")

GEOM = [jigsawpy.jigsaw_msh_t()]
SPAC = [jigsawpy.jigsaw_msh_t()]

FULL_SPHERE_RADIUS = +6.371E+003

PROJ_CENTRE = [-75.2316, 39.1269]


def filt_narivs(feat):

    return feat["properties"]["UP_CELLS"] > 12500


def setgeom():

    geom = jigsawpy.jigsaw_msh_t()
    poly = jigsawpy.jigsaw_msh_t()

#------------------------------------ define JIGSAW geometry

    geom.mshID = "euclidean-mesh"
    
#------------------------------------ set watershed boundary

    print("BUILDING MESH GEOM.")

    filename = os.path.join(
        "data",
        "NHD_H_0204_HU4_Shape", "Shape", "WBDHU4.shp")

    loadshp(filename, poly)

    poly.point["coord"] *= np.pi / +180.

    addpoly(geom, poly, +1)

    filename = os.path.join(
        "data",
        "NHD_H_0205_HU4_Shape", "Shape", "WBDHU4.shp")

    loadshp(filename, poly)

    poly.point["coord"] *= np.pi / +180.

    addpoly(geom, poly, +2)
    
#------------------------------------ approx. stream network
    filename = os.path.join(
        "data", "namerica_rivers", "narivs.shp")

    loadshp(filename, poly, filt_narivs)

    poly.point["coord"] *= np.pi / +180.

    itag = innerto(poly.vert2["coord"], geom)

    keep = np.logical_and.reduce((
        itag[poly.edge2["index"][:, 0]] > +0,
        itag[poly.edge2["index"][:, 1]] > +0
    ))
    poly.edge2 = poly.edge2[keep]

    addline(geom, poly, +0)

#------------------------------------ proj. to local 2-plane
    proj = jigsawpy.jigsaw_prj_t()
    proj.prjID = "stereographic"
    proj.radii = FULL_SPHERE_RADIUS
    proj.xbase = PROJ_CENTRE[0] * np.pi / +180.
    proj.ybase = PROJ_CENTRE[1] * np.pi / +180.

    jigsawpy.project(geom, proj, "fwd")

    GEOM[0] = geom                      # save a "pointer"

    return geom


def setinit():

    init = jigsawpy.jigsaw_msh_t()

#------------------------------------ set initial conditions

    # just an empty msh_t object!

    return init


def setopts():

    opts = jigsawpy.jigsaw_jig_t()

#------------------------------------ define JIGSAW useropts

    opts.hfun_scal = "absolute"
    opts.hfun_hmax = float("inf")       # null spacing lim
    opts.hfun_hmin = float(+0.00)

    opts.mesh_dims = +2                 # 2-dim. simplexes

    opts.bisection = +0                 # just single-lev.

    opts.mesh_rad2 = +1.20
    opts.mesh_eps1 = +0.67

#   opts.optm_kern = "cvt+dqdx"

    opts.optm_qlim = +9.5E-01           # tighter opt. tol
    opts.optm_iter = +64
    opts.optm_qtol = +1.0E-05

    return opts


def setspac():

    spac_wbd = 5.                       # watershed km
    spac_pfz = 2.                       # PFZ km
    elev_pfz = 25.                      # PFZ elev. thresh
    dhdx_lim = 0.0625                   # |dH/dx| thresh

    bbox = [-80., 35., -70., 45.]

    spac = jigsawpy.jigsaw_msh_t()

    opts = jigsawpy.jigsaw_jig_t()
    poly = jigsawpy.jigsaw_msh_t()

    geom = GEOM[0]

    opts.jcfg_file = os.path.join(TDIR, "opts.jig")
    opts.hfun_file = os.path.join(TDIR, "spac.msh")

#------------------------------------ define spacing pattern

    print("BUILDING MESH SPAC.")

    print("Loading elevation assets...")
    
    data = nc.Dataset(os.path.join(
        "data",
        "etopo_gebco", "etopo_gebco_tiled.nc"), "r")

    zlev = np.array(data.variables["z"])

    spac.mshID = "ellipsoid-grid"       # use elev. grid
    spac.radii = np.full(
        +3, FULL_SPHERE_RADIUS, dtype=spac.REALS_t)

    spac.xgrid = np.array(
        data.variables["x"][:], dtype=spac.REALS_t)

    spac.ygrid = np.array(
        data.variables["y"][:], dtype=spac.REALS_t)

    xmsk = np.logical_and.reduce((
        spac.xgrid >= bbox[0] - 1.,
        spac.xgrid <= bbox[2] + 1.
    ))

    ymsk = np.logical_and.reduce((
        spac.ygrid >= bbox[1] - 1.,
        spac.ygrid <= bbox[3] + 1.
    ))

    spac.xgrid = spac.xgrid[xmsk]
    spac.ygrid = spac.ygrid[ymsk]
    zlev = zlev[:, xmsk]
    zlev = zlev[ymsk, :]
    
    spac.xgrid *= np.pi / +180.
    spac.ygrid *= np.pi / +180.

    xgrd, ygrd = np.meshgrid(spac.xgrid, spac.ygrid)

    spac.value = +10. * np.ones(
        (spac.ygrid.size, spac.xgrid.size))

    grid = np.concatenate((             # to [x, y] list
        xgrd.reshape((xgrd.size, +1)),
        ygrd.reshape((ygrd.size, +1))), axis=+1)

#------------------------------------ push watershed(s) = 5.

    print("Compute watersheds h(x)...")

    shed = np.full(
        (grid.shape[0]), False, dtype=bool)

    filename = os.path.join(
        "data",
        "NHD_H_0204_HU4_Shape", "Shape", "WBDHU4.shp")

    loadshp(filename, poly)

    poly.point["coord"] *= np.pi / +180.

    mask, _ = inpoly2(
        grid, poly.point["coord"], poly.edge2["index"])

    shed[mask] = True

    filename = os.path.join(
        "data",
        "NHD_H_0205_HU4_Shape", "Shape", "WBDHU4.shp")

    loadshp(filename, poly)

    poly.point["coord"] *= np.pi / +180.

    mask, _ = inpoly2(
        grid, poly.point["coord"], poly.edge2["index"])

    shed[mask] = True

    shed = np.reshape(shed, spac.value.shape)

    spac.value[shed] = \
        np.minimum(spac_wbd, spac.value[shed])

#------------------------------------ partially flooded zone

    mask = np.logical_and(shed, zlev < elev_pfz)

    spac.value[mask] = \
        np.minimum(spac_pfz, spac.value[mask])

#------------------------------------ push |DH/DX| threshold

    print("Impose |DH/DX| threshold...")

    spac.slope = np.full(
        spac.value.shape, dhdx_lim, dtype=spac.REALS_t)

    jigsawpy.savemsh(opts.hfun_file, spac,
                     "precision = 9")

    jigsawpy.cmd.marche(opts, spac)

    spac.slope = \
        np.empty((+0), dtype=spac.REALS_t)

#------------------------------------ proj. to local 2-plane
    proj = jigsawpy.jigsaw_prj_t()
    proj.prjID = "stereographic"
    proj.radii = FULL_SPHERE_RADIUS
    proj.xbase = PROJ_CENTRE[0] * np.pi / +180.
    proj.ybase = PROJ_CENTRE[1] * np.pi / +180.

    jigsawpy.project(spac, proj, "fwd")

    SPAC[0] = spac                      # save a "pointer"

    return spac
