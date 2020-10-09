
import numpy as np


def sphdist(rsph, xlon, ylat, xmid, ymid):
    """
    SPHDIST: return the distance from the points [XLON,YLAT]
    to the point [XMID,YMID] on a sphere of radii RSPH.

    Authors: Darren Engwirda

    """

    dlon = .5 * (xlon - xmid)
    dlat = .5 * (ylat - ymid)

    dist = 2. * rsph * np.arcsin(np.sqrt(
        np.sin(dlat) ** 2 +
        np.sin(dlon) ** 2 * np.cos(ylat) * np.cos(ymid)
    ))

    return dist


def blender(val1, val2, dist, blen, bgap):
    """
    BLENDER: return a blend of VAL1 and VAL2, as a nonlinear
    weighted combination.

    TANH((DIST-LEN) ./ GAP) weighting is used -- blending is
    centred about DIST=LEN, and spans an approx. DIST of GAP.

    Authors: Darren Engwirda

    """

    beta = .5 + .5 * np.tanh(
        np.pi * (dist - blen) / bgap)

    return (+0.0 + beta) * val2 + (+1.0 - beta) * val1
