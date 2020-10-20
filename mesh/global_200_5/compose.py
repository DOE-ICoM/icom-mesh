
import numpy as np
import jigsawpy
from mpas_tools.cime.constants import constants

"""
GLOBAL-200-5: Quasi-uniform global config. with mid-Atlantic
focused regional refinement:
    - 200km (global ocean)
    - 5km halo in the mid-Atlantic

"""
# Authors: Darren Engwirda

FULL_SPHERE_RADIUS = constants["SHR_CONST_REARTH"] / 1.E+003


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

    opts.bisection = +5                 # call heutristic!

#   opts.optm_kern = "cvt+dqdx"

    opts.optm_qlim = +9.5E-01           # tighter opt. tol
    opts.optm_iter = +64
    opts.optm_qtol = +1.0E-05

    return opts


def setspac():

    spac = jigsawpy.jigsaw_msh_t()

#------------------------------------ define spacing pattern

    xmid = -76.5 * np.pi / 180.0
    ymid = +39.5 * np.pi / 180.0

    spac.mshID = "ellipsoid-grid"
    spac.radii = np.full(
        3, FULL_SPHERE_RADIUS, dtype=spac.REALS_t)

    spac.xgrid = np.linspace(
        -1. * np.pi, +1. * np.pi, 720)

    spac.ygrid = np.linspace(
        -.5 * np.pi, +.5 * np.pi, 361)

    xmat, ymat = \
        np.meshgrid(spac.xgrid[:] , spac.ygrid[:])

    spac.value = +200. - 195. * np.exp(-(
        +12.5 * (xmat - xmid) ** 2 +
        +20.0 * (ymat - ymid) ** 2) ** 2)

    return spac
