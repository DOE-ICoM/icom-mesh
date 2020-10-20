
import os
import time
import copy
import numpy as np
import jigsawpy
from mpas_tools.cime.constants import constants
import mpas_tools.mesh.creation.mesh_definition_tools as mdt

from util.loadshp import loadshp
from util.spacing import zipnear
from util.utility import addpoly, addline, innerto

"""
DELAWARE-60-30-5-2-w-boundary: like DELAWARE-60-30-5-2, but
with boundary constraints (watersheds, rivers, etc) imposed:
    - ECC-60-to-30 (global ocean)
    - 45km (global land)
    - 30km (North Atlantic ocean)
    - 5.km (Delaware, etc watersheds)
    - 2.km (Delaware coastline + PFZ)

"""
# Authors: Darren Engwirda

HERE = os.path.abspath(os.path.dirname(__file__))
TDIR = os.path.join(HERE, "tmp")
ODIR = os.path.join(HERE, "out")

GEOM = [jigsawpy.jigsaw_msh_t()]
SPAC = [jigsawpy.jigsaw_msh_t()]

FULL_SPHERE_RADIUS = constants["SHR_CONST_REARTH"] / 1.E+003


def filt_narivs(feat):

    return feat["properties"]["UP_CELLS"] > 12500


def setgeom():

    geom = jigsawpy.jigsaw_msh_t()
    poly = jigsawpy.jigsaw_msh_t()

#------------------------------------ define JIGSAW geometry

    geom.mshID = "ellipsoid-mesh"
    geom.radii = np.full(
        +3, FULL_SPHERE_RADIUS, dtype=geom.REALS_t)

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

    filename = os.path.join(
        "data",
        "NHD_H_0206_HU4_Shape", "Shape", "WBDHU4.shp")

    loadshp(filename, poly)

    poly.point["coord"] *= np.pi / +180.

    addpoly(geom, poly, +3)

    filename = os.path.join(
        "data",
        "NHD_H_0207_HU4_Shape", "Shape", "WBDHU4.shp")

    loadshp(filename, poly)

    poly.point["coord"] *= np.pi / +180.

    addpoly(geom, poly, +4)

    filename = os.path.join(
        "data",
        "NHD_H_0208_HU4_Shape", "Shape", "WBDHU4.shp")

    loadshp(filename, poly)

    poly.point["coord"] *= np.pi / +180.

    addpoly(geom, poly, +5)

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

    addline(geom, poly, +1)

    GEOM[0] = geom                      # save a "pointer"

    return geom


def setopts():

    opts = jigsawpy.jigsaw_jig_t()

#------------------------------------ define JIGSAW useropts

    opts.hfun_scal = "absolute"
    opts.hfun_hmax = float("inf")       # null spacing lim
    opts.hfun_hmin = float(+0.00)

    opts.mesh_dims = +2                 # 2-dim. simplexes

    opts.bisection = +0                 # call heutristic!

    opts.mesh_rad2 = +1.20
    opts.mesh_eps1 = +0.67

    opts.optm_kern = "cvt+dqdx"

    opts.optm_qlim = +9.5E-01           # tighter opt. tol
    opts.optm_iter = +64
    opts.optm_qtol = +1.0E-05

    return opts


def setspac():

    spac = jigsawpy.jigsaw_msh_t()

#------------------------------------ define spacing pattern

    print("BUILDING MESH SPAC.")

    filename = os.path.join(
        HERE, "..",
        "delaware_60_30_5_2", "tmp", "spac.msh")

    if os.path.isfile(filename):        # load from init.

        jigsawpy.loadmsh(filename, spac)

    else:

        raise Exception(
            "File not found: first run delaware_60_30_5_2"
            " to form initial conditions!")

    SPAC[0] = spac                      # save a "pointer"

    return spac


def setinit():

    init = jigsawpy.jigsaw_msh_t()

#------------------------------------ set initial conditions

    print("BUILDING MESH INIT.")

    filename = os.path.join(
        HERE, "..",
        "delaware_60_30_5_2", "tmp", "mesh.msh")

    if os.path.isfile(filename):        # load from init.

        jigsawpy.loadmsh(filename, init)

    else:

        raise Exception(
            "File not found: first run delaware_60_30_5_2"
            " to form initial conditions!")

    geom = copy.deepcopy(GEOM[0])
    spac = copy.deepcopy(SPAC[0])

    zipnear(init, geom, spac, near=+1.E+01)

    return init
