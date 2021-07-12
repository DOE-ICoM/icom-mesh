
import os
import time
import math
import numpy as np
import jigsawpy
import netCDF4 as nc
import argparse

from skimage.filters.rank import median, percentile
from skimage.filters import scharr, farid, gaussian
from skimage.morphology import disk

import sys
sys.path.append(os.path.dirname(
                os.path.dirname(
                os.path.dirname(os.path.abspath(__file__)))))
from util.saveesm import saveesm

from mpas_tools.cime.constants import constants
import mpas_tools.mesh.creation.mesh_definition_tools as mdt

"""
BAROTROPIC-TIDES: (various) meshes for barotropic tide runs.
    - Uniform background spacing +
    - Wavelength heuristic.
    - GRAD(elev) heuristic.

"""
# Authors: Darren Engwirda

HERE = os.path.abspath(os.path.dirname(__file__))
TDIR = os.path.join(HERE, "tmp")
ODIR = os.path.join(HERE, "out")

FULL_SPHERE_RADIUS = constants["SHR_CONST_REARTH"] / 1.E+003


def coarsen_spacing_pixels(hmat, down):

    print("Coarsening mesh-spacing pixels...")

    rows = hmat.shape[0] // down
    cols = hmat.shape[1] // down

    htmp = np.full(
        (rows, cols), (np.amax(hmat)), dtype=hmat.dtype)

    for jpos in range(down):
        for ipos in range(down):

            iend = hmat.shape[0] - down + ipos + 1
            jend = hmat.shape[1] - down + jpos + 1

            htmp = np.minimum(
                htmp,
            hmat[ipos:iend:down, jpos:jend:down])

    return htmp


def limit_spacing_gradient(spac, dhdx):

    print("Smoothing h(x) via |dh/dx| limits...")

    opts = jigsawpy.jigsaw_jig_t()

    spac.slope = np.full((+1), dhdx, dtype=spac.FLT32_t)

    opts.verbosity = +1

    jigsawpy.lib.marche(opts, spac)

    return spac


def swe_wavelength_spacing(
        ocnh, nwav, hmin, hmax, halo, plev,
        T_M2=12.*3600., grav=9.806):

    print("Computing wavelength heuristic...")

    vals = T_M2 * np.sqrt(
        grav * np.maximum(5, ocnh)) / nwav / 1000.

    vals[ocnh <= 0.] = hmax
    vals = np.maximum(vals, hmin)
    vals = np.minimum(vals, hmax)

    vals = np.asarray(vals, dtype=np.uint16)
    vals = percentile(
        vals, selem=disk(halo), mask=(ocnh>0.), p0=plev)

    return vals


def elev_sharpness_spacing(
        ocnh, nslp, hmin, hmax, halo, plev, sdev):

    print("Computing GRAD(elev) heuristic...")

    dzdx = scharr(gaussian(np.asarray(
        ocnh, dtype=np.float32), sigma=sdev, mode="wrap"))
    
    dzdx = np.maximum(1.E-08, dzdx) # no divide-by-zero

    vals = np.maximum(
        5., np.abs(ocnh)) / dzdx * 2. * np.pi / nslp

    vals = np.maximum(vals, hmin)
    vals = np.minimum(vals, hmax)

    vals = np.asarray(vals, dtype=np.uint16)
    vals = percentile(
        vals, selem=disk(halo), mask=(ocnh>0.), p0=plev)

    return vals


def compose(args):

    saveesm(msh_file=opts.mesh_file, out_path=ODIR, 
            on_a_sphere=True, 
            sphere_radius=constants["SHR_CONST_REARTH"],
            with_cavities=True)

    raise Exception()



    opts = jigsawpy.jigsaw_jig_t()
    mesh = jigsawpy.jigsaw_msh_t()

    opts.geom_file = os.path.join(TDIR, "geom.msh")
    opts.jcfg_file = os.path.join(TDIR, "opts.jig")
    opts.hfun_file = os.path.join(TDIR, "spac.msh")
    opts.mesh_file = os.path.join(TDIR, "mesh.msh")

    geom = setgeom(args)
    spac = setspac(args)

    jigsawpy.savemsh(opts.geom_file, geom)
    jigsawpy.savemsh(opts.hfun_file, spac)

    jigsawpy.savevtk(os.path.join(TDIR, "spac.vtk"), spac)
    if (args.stop_here == "build-spacing"): return

#------------------------------------ define JIGSAW useropts

    opts.hfun_scal = "absolute"
    opts.hfun_hmax = float("inf")       # null spacing lim
    opts.hfun_hmin = float(+0.00)

    opts.verbosity = +1
    opts.mesh_dims = +2                 # 2-dim. simplexes
   
    opts.optm_kern = "cvt+dqdx"

    opts.optm_iter = 32                 # tighter opt. tol
    opts.optm_qtol = +1.00E-05

    rbar = np.mean(geom.radii)          # bisect heuristic
    hbar = np.mean(spac.value)
    nlev = round(math.log2(
        rbar / math.sin(.4 * math.pi) / hbar)
    )

   #jigsawpy.cmd.jigsaw(opts, mesh)
    jigsawpy.cmd.tetris(opts, nlev - 1, mesh)

    jigsawpy.savevtk(os.path.join(TDIR, "mesh.vtk"), mesh)
    if (args.stop_here == "generate-tria"): return

#------------------------------------ export to MPAS unified

    saveesm(msh_file=opts.mesh_file, out_path=ODIR, 
            on_a_sphere=True, 
            sphere_radius=constants["SHR_CONST_REARTH"],
            with_cavities=True)

def setgeom(args):

    geom = jigsawpy.jigsaw_msh_t()

#------------------------------------ define JIGSAW geometry

    geom.mshID = "ellipsoid-mesh"
    geom.radii = np.full(
        3, FULL_SPHERE_RADIUS, dtype=geom.REALS_t)

    return geom


def setspac(args):

    spac = jigsawpy.jigsaw_msh_t()

#------------------------------------ define spacing pattern

    dhdx = args.spac_dhdx   # max allowable slope in spacing
    hmin = args.spac_hmin   # min mesh-spacing value (km)
    hmax = args.spac_hmax   # max mesh-spacing value (km)
    hbar = args.spac_hbar   # constant spacing value (km)

    halo = args.filt_halo   # DEM pixels per filtering radii
    sdev = args.filt_sdev   # std-dev for gaussian filter   
    plev = args.filt_plev   # median-style filter percentile
    
    nwav = args.ncell_wav   # number of cells per wavelength
    nslp = args.ncell_slp   # number of cells per grad(elev)

    print("Loading elevation assets...")
   
    data = nc.Dataset(args.elev_file, "r")

    ocnh = np.asarray(
        data["ocn_thickness"][:], dtype=spac.FLT32_t)


    print("Computing background h(x)...")

    hmat = np.full(
        (ocnh.shape[:]), (hbar), dtype=spac.FLT32_t)

    if (nwav >= 1):
        hmat = np.minimum(
            hmat, swe_wavelength_spacing(
                ocnh, nwav, hmin, hmax, halo, plev))

    if (nslp >= 1):
        hmat = np.minimum(
            hmat, elev_sharpness_spacing(
                ocnh, nslp, hmin, hmax, halo, plev, 
                sdev))
    
    hmat[ocnh <= 0.] = hmax
    hmat = np.maximum(hmat, hmin)
    hmat = np.minimum(hmat, hmax)

#-- pack h(x) data into jigsaw data-type: average pixel-to-
#-- node, careful with periodic BC's.

    hmat = coarsen_spacing_pixels(hmat, down=4)
    
    spac.mshID = "ellipsoid-grid"       # use the elv. grid
    spac.radii = np.full(
        3, FULL_SPHERE_RADIUS, dtype=spac.REALS_t)

    spac.xgrid = np.linspace(
        -1. * np.pi, +1. * np.pi, hmat.shape[1] + 1)

    spac.ygrid = np.linspace(
        -.5 * np.pi, +.5 * np.pi, hmat.shape[0] + 1)

    R = hmat.shape[0]; C = hmat.shape[1]

    spac.value = np.zeros(
        (R + 1, C + 1), dtype=spac.FLT32_t)

    npos = np.arange(+0, hmat.shape[0] + 1)
    epos = np.arange(-1, hmat.shape[1] - 0)
    spos = np.arange(-1, hmat.shape[0] - 0)
    wpos = np.arange(+0, hmat.shape[1] + 1)
    
    npos[npos >= +R] = R - 1; spos[spos <= -1] = +0
    epos[epos <= -1] = C - 1; wpos[wpos >= +C] = +0

    npos, epos = np.meshgrid(
        npos, epos, sparse=True, indexing="ij")
    spos, wpos = np.meshgrid(
        spos, wpos, sparse=True, indexing="ij")

    spac.value += hmat[npos, epos] * (+1. / 4.)
    spac.value += hmat[npos, wpos] * (+1. / 4.)
    spac.value += hmat[spos, epos] * (+1. / 4.)
    spac.value += hmat[spos, wpos] * (+1. / 4.)

    spac = limit_spacing_gradient(spac, dhdx=dhdx)
    
    return spac


if (__name__ == "__main__"):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "--elev-file", dest="elev_file", type=str,
        required=True, help="Path to DEM pixel file.")

    parser.add_argument(
        "--stop-here", dest="stop_here", type=str,
        required=False, default="build-spacing", 
        help="Stop part way through meshing process:\n"
             "{build-spacing}, generate-tria, generate-mpas")

    parser.add_argument(
        "--spac-dhdx", dest="spac_dhdx", type=float,
        default=0.1000,
        required=False, help="Max spacing gradients: {0.100}")

    parser.add_argument(
        "--spac-hmin", dest="spac_hmin", type=float,
        default=10.,
        required=False, help="Min spacing lim. (km): {10.00}")

    parser.add_argument(
        "--spac-hmax", dest="spac_hmax", type=float,
        default=75.,
        required=False, help="Max spacing lim. (km): {75.00}")

    parser.add_argument(
        "--spac-hbar", dest="spac_hbar", type=float,
        default=60.,
        required=False, help="Constant spacing (km): {60.00}")

    parser.add_argument(
        "--filt-halo", dest="filt_halo", type=int,
        default=50,
        required=False, help="Num. pixels in filter: {50}")

    parser.add_argument(
        "--filt-plev", dest="filt_plev", type=float,
        default=0.3250,
        required=False, help="Filter low-percentile: {0.325}")

    parser.add_argument(
        "--filt-sdev", dest="filt_sdev", type=float,
        default=3.0000,
        required=False, help="Filter std.-deviation: {3.000}")

    parser.add_argument(
        "--ncell_slp", dest="ncell_slp", type=int,
        default=0,
        required=False, help="nCells per grad(elev): {0}")

    parser.add_argument(
        "--ncell_wav", dest="ncell_wav", type=int,
        default=60,
        required=False, help="nCells per wavelength: (60)")

    compose(parser.parse_args())
