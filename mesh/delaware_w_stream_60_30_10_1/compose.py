
import os
import numpy as np

import jigsawpy


def setgeom(args):

    geom = jigsawpy.jigsaw_msh_t()

#------------------------------------ define JIGSAW geometry

    geom.mshID = "ellipsoid-mesh"
    geom.radii = np.full(
        3, args.sphere_radius, dtype=geom.REALS_t)

    return geom


def setinit(args):

    init = jigsawpy.jigsaw_msh_t()

#------------------------------------ set initial conditions

    name = os.path.join(
        "..", "delaware_60_30_10_1", "tmp", "mesh.msh")

    jigsawpy.loadmsh(name, init)

    return init


def setopts(args):

    opts = jigsawpy.jigsaw_jig_t()

#------------------------------------ define JIGSAW useropts

    opts.hfun_scal = "absolute"
    opts.hfun_hmax = float("inf")       # null spacing lim
    opts.hfun_hmin = float(+0.00)

    opts.mesh_dims = +2                 # 2-dim. simplexes

    opts.bisection = +0                 # use single level

    opts.optm_qlim = +9.5E-01           # tighter opt. tol
    opts.optm_iter = +32
    opts.optm_qtol = +1.0E-05

    return opts


def setspac(args):

    spac = jigsawpy.jigsaw_msh_t()

#------------------------------------ define spacing pattern

    spac.mshID = "ellipsoid-grid"
    spac.radii = np.full(
        3, args.sphere_radius, dtype=spac.REALS_t)

    spac.xgrid = np.linspace(
        -1. * np.pi, +1. * np.pi, 360)

    spac.ygrid = np.linspace(
        -.5 * np.pi, +.5 * np.pi, 180) 
    
    spac.value = np.full(
        (180, 360), +1.0E+002, dtype=spac.REALS_t)

    return spac
