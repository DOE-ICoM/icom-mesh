
import numpy as np
import jigsawpy
import mpas_tools.mesh.creation.mesh_definition_tools as mdt

"""
GLOBAL-30-10: the "Rossby-radius" configuration - 30-to-10km.
    - RRS-30-to-10km (global ocean)

Authors: Darren Engwirda

"""

FULL_SPHERE_RADIUS = +6.371E+003


def setgeom():

    geom = jigsawpy.jigsaw_msh_t()

#------------------------------------ define JIGSAW geometry

    geom.mshID = "ellipsoid-mesh"
    geom.radii = np.full(
        3, FULL_SPHERE_RADIUS, dtype=geom.REALS_t)

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

#   opts.optm_kern = "cvt+dqdx"

    opts.optm_qlim = +9.5E-01           # tighter opt. tol
    opts.optm_iter = 64
    opts.optm_qtol = +1.0E-05

    return opts


def setspac():

    spac = jigsawpy.jigsaw_msh_t()

#------------------------------------ define spacing pattern

    spac.mshID = "ellipsoid-grid"
    spac.radii = np.full(
        3, FULL_SPHERE_RADIUS, dtype=spac.REALS_t)

    spac.xgrid = np.linspace(
        -1. * np.pi, +1. * np.pi, 360)

    spac.ygrid = np.linspace(
        -.5 * np.pi, +.5 * np.pi, 181)

    vals = mdt.RRS_CellWidthVsLat(
        spac.ygrid * 180. / np.pi, 30., 10.)

    spac.value = np.array(np.tile(
        vals, (1, spac.xgrid.size)), dtype=spac.REALS_t)

    return spac
