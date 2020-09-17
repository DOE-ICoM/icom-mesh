
import numpy as np
from scipy import interpolate
import netCDF4 as nc
import jigsawpy as jp


def tileDEM(dem1, dem2):
    """
    TILEDEM: insert a high-resolution elev. dataset (DEM-2),
    within an enclosing dataset (DEM-1).

    Tiling is blocked such that output is still rectangular.

    Authors: Darren Engwirda

    """

#--------------------- insert DEM-2 as a "tile" within DEM-1

    xx11 = np.min(dem2.xgrid)
    xx22 = np.max(dem2.xgrid)
    yy11 = np.min(dem2.ygrid)
    yy22 = np.max(dem2.ygrid)

    xdel = np.mean(np.diff(dem1.xgrid))
    xx11 = xx11 - 0.5 * xdel
    xx22 = xx22 + 0.5 * xdel

    x1st = dem1.xgrid < xx11
    xend = dem1.xgrid > xx22

    ydel = np.mean(np.diff(dem1.ygrid))
    yy11 = yy11 - 0.5 * ydel
    yy22 = yy22 + 0.5 * ydel

    y1st = dem1.ygrid < yy11
    yend = dem1.ygrid > yy22

#--------------------- interpolate to build blockwise values

    zfun = interpolate.RectBivariateSpline(
        dem1.ygrid,
        dem1.xgrid, dem1.value)

    zz11 = dem1.value[y1st, :]
    zz11 = zz11[:, x1st]

    zz12 = zfun(
        dem1.ygrid[y1st], dem2.xgrid[:])

    zz13 = dem1.value[y1st, :]
    zz13 = zz13[:, xend]

    zz21 = zfun(
        dem2.ygrid[:], dem1.xgrid[x1st])

    zz22 = dem2.value[:, :]

    zz23 = zfun(
        dem2.ygrid[:], dem1.xgrid[xend])

    zz31 = dem1.value[yend, :]
    zz31 = zz31[:, x1st]

    zz32 = zfun(
        dem1.ygrid[yend], dem2.xgrid[:])

    zz33 = dem1.value[yend, :]
    zz33 = zz33[:, xend]

#--------------------- concatenate blockwise values together

    tile = jp.jigsaw_msh_t()

    tile.mshID = "ellipsoid-grid"
    tile.radii = dem1.radii

    tile.xgrid = np.concatenate((
        dem1.xgrid[x1st],
        dem2.xgrid[:],
        dem1.xgrid[xend]))

    tile.ygrid = np.concatenate((
        dem1.ygrid[y1st],
        dem2.ygrid[:],
        dem1.ygrid[yend]))

    tile.value = np.block([
        [zz11, zz12, zz13],
        [zz21, zz22, zz23],
        [zz31, zz32, zz33]])

    return tile


def makeDEM():
    """
    MAKEDEM: assemble a "tiled" elevation dataset centred on
    the Delaware Bay region.

    A high-resolution GEBCO dataset is inserted within an
    ETOPO1 block, which is itself inserted in a "downscaled"
    ETOPO1-based global distribution.   

    Tiling is blocked such that output is still rectangular.

    Output is written to a new ETOPO1-style netCDF file.

    Authors: Darren Engwirda

    """

    dat1 = "../data/ETOPO1_Ice_g_gmt4.nc"
    dat2 = \
        "../data/gebco_2019_N43p5_S34p0_W-80p0_E-69p5.nc"

    dat3 = "../data/etopo_gebco_tiled.nc"

#--------------------- load ETOPO + GEBCO elevation datasets

    etop = jp.jigsaw_msh_t()
    gbco = jp.jigsaw_msh_t()

    print("1. Loading ETOPO dataset...")

    data = nc.Dataset(dat1, "r")

    etop.mshID = "ellipsoid-grid"
    etop.radii = np.full(
        3, +6.371E+003, dtype=etop.REALS_t)
    etop.xgrid = np.array(
        data.variables["x"][:])
    etop.ygrid = np.array(
        data.variables["y"][:])
    etop.value = np.array(
        data.variables["z"][:])

    print("2. Loading GEBCO dataset...")

    data = nc.Dataset(dat2, "r")

    gbco.mshID = "ellipsoid-grid"
    gbco.radii = np.full(
        3, +6.371E+003, dtype=gbco.REALS_t)
    gbco.xgrid = np.array(
        data.variables["lon"][:])
    gbco.ygrid = np.array(
        data.variables["lat"][:])
    gbco.value = np.array(
        data.variables["elevation"][:])

#--------------------- downscale global data outside of halo

    print("3. Tile combined dataset...")

    halo = +10.      # expand outside of GEBCO block by halo

    skip = + 8       # downscale ETOPO by keeping n-th entry

    xx11 = np.min(gbco.xgrid) - halo
    xx22 = np.max(gbco.xgrid) + halo
    yy11 = np.min(gbco.ygrid) - halo
    yy22 = np.max(gbco.ygrid) + halo

    xset = np.logical_and(
        etop.xgrid > xx11,
        etop.xgrid < xx22)
    yset = np.logical_and(
        etop.ygrid > yy11,
        etop.ygrid < yy22)

    xset[+0::skip] = True; xset[-1] = True
    yset[+0::skip] = True; yset[-1] = True

    etop.xgrid = etop.xgrid[xset]
    etop.ygrid = etop.ygrid[yset]
    etop.value = etop.value[yset, :]
    etop.value = etop.value[:, xset]

#--------------------- insert GEBCO block inside ETOPO array

    tile = tileDEM(etop, gbco)

#--------------------- write data as a new netCDF elev. file

    data = nc.Dataset(dat3, "w", format="NETCDF4")

    data.createDimension("x", tile.xgrid.size)
    data.createDimension("y", tile.ygrid.size)

    xdat = data.createVariable("x", "f4", ("x"))
    ydat = data.createVariable("y", "f4", ("y"))
    zdat = data.createVariable("z", "f4", ("y", "x"))

    xdat[:] = tile.xgrid[:]
    ydat[:] = tile.ygrid[:]
    zdat[:, :] = tile.value[:, :]

    data.description = \
        "Elev. data tiled from ETOPO1+GEBCO datasets"

    data.close()

    return


if (__name__ == "__main__"): makeDEM()
