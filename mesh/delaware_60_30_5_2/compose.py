
import os
import time
import numpy as np
from scipy.ndimage import median_filter
import netCDF4 as nc
import jigsawpy
import mpas_tools.mesh.creation.mesh_definition_tools as mdt

from util.loadshp import loadshp
from util.loadgeo import loadgeo
from util.inpoly2 import inpoly2
from util.spacing import sphdist, blender


"""
DELAWARE-60-30-5-2: a variable-res. Delaware config., inc.:
    - ECC-60-to-30 (global ocean)
    - 45km (global land)
    - 30km (North Atlantic ocean)
    - 5.km (Delaware, etc watersheds)
    - 2.km (Delaware coastline + PFZ)

Authors: Darren Engwirda

"""

HERE = os.path.abspath(os.path.dirname(__file__))
TDIR = os.path.join(HERE, "tmp")
ODIR = os.path.join(HERE, "out")

FULL_SPHERE_RADIUS = +6.371E+003


def setgeom():

    geom = jigsawpy.jigsaw_msh_t()

#------------------------------------ define JIGSAW geometry

    print("BUILDING MESH GEOM.")

    geom.mshID = "ellipsoid-mesh"
    geom.radii = np.full(
        +3, FULL_SPHERE_RADIUS, dtype=geom.REALS_t)

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

    opts.bisection = -1                 # call heutristic!

    opts.optm_kern = "cvt+dqdx"

    opts.optm_qlim = +9.5E-01           # tighter opt. tol
    opts.optm_iter = +256
    opts.optm_qtol = +1.0E-05

    return opts


def setspac():

    spac_ocn = 30.                      # regional ocn. km
    spac_1_m = 2.                       # sqrt(H) ocn. km
    spac_lnd = 45.                      # global land km
    spac_wbd = 5.                       # watershed km
    spac_pfz = 2.                       # PFZ km
    elev_pfz = 25.                      # PFZ elev. thresh
    dhdx_lim = 0.0625                   # |dH/dx| thresh

    fade_pos = [-75.2316, 39.1269]
    fade_len = 700.
    fade_gap = 350.

    spac = jigsawpy.jigsaw_msh_t()

    opts = jigsawpy.jigsaw_jig_t()
    poly = jigsawpy.jigsaw_msh_t()

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

    spac.xgrid *= np.pi / +180.
    spac.ygrid *= np.pi / +180.

    xgrd, ygrd = np.meshgrid(spac.xgrid, spac.ygrid)

    grid = np.concatenate((             # to [x, y] list
        xgrd.reshape((xgrd.size, +1)),
        ygrd.reshape((ygrd.size, +1))), axis=+1)

#------------------------------------ global ocn ec-60-to-30

    print("Compute global ocean h(x)...")

    vals = \
        mdt.EC_CellWidthVsLat(spac.ygrid * 180. / np.pi)

    vals = np.reshape(vals, (spac.ygrid.size, 1))

    spac.value = np.array(np.tile(
        vals, (1, spac.xgrid.size)), dtype=spac.REALS_t)

#------------------------------------ region ocn "eddy" halo

    filename = os.path.join(
        "data",
        "na_ocean_halo", "na_ocean_halo.geojson")

    loadgeo(filename, poly)

    poly.point["coord"] *= np.pi / +180.

    mask, _ = inpoly2(
        grid, poly.point["coord"], poly.edge2["index"])

    mask = np.reshape(mask, spac.value.shape)

    mask = np.logical_and(mask, zlev <= +0.0)

    spac.value[mask] = \
        np.minimum(spac_ocn, spac.value[mask])

#------------------------------------ coastal ocn heuristics

    mask = zlev <= +0.0

    hval = np.sqrt(np.maximum(-zlev, +0.0))
    hval = np.maximum(
        spac_1_m, hval / np.sqrt(+1.0) / spac_1_m)

    dist = sphdist(
        FULL_SPHERE_RADIUS, grid[:, 0], grid[:, 1],
        fade_pos[0] * np.pi / 180.,
        fade_pos[1] * np.pi / 180.)
    dist = np.reshape(dist, spac.value.shape)

    hval = blender(
        hval, spac.value, dist, fade_len, fade_gap)

    spac.value[mask] = \
        np.minimum(hval[mask], spac.value[mask])

#------------------------------------ global lnd const. = 45

    print("Compute global land h(x)...")

    halo = +9                           # filter islands

    zmed = median_filter(zlev, size=halo, mode="wrap")

    spac.value[zmed >= 0.0] = spac_lnd

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

    filename = os.path.join(
        "data",
        "NHD_H_0206_HU4_Shape", "Shape", "WBDHU4.shp")

    loadshp(filename, poly)

    poly.point["coord"] *= np.pi / +180.

    mask, _ = inpoly2(
        grid, poly.point["coord"], poly.edge2["index"])

    shed[mask] = True

    filename = os.path.join(
        "data",
        "NHD_H_0207_HU4_Shape", "Shape", "WBDHU4.shp")

    loadshp(filename, poly)

    poly.point["coord"] *= np.pi / +180.

    mask, _ = inpoly2(
        grid, poly.point["coord"], poly.edge2["index"])

    shed[mask] = True

    filename = os.path.join(
        "data",
        "NHD_H_0208_HU4_Shape", "Shape", "WBDHU4.shp")

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

    return spac
