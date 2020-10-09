
import warnings
from pathlib import Path

from jigsawpy import jigsaw_msh_t, certify


def save_mesh_file(mesh, fptr):
    """
    SAVE-MESH-FILE: save a JIGSAW mesh object to *.DAT file.

    """

    if (mesh.vert2 is not None):
        xpts = mesh.vert2["coord"]
        for ipos in range(mesh.vert2.size):
            fptr.write(
                f"{ipos}\t"
                f"{xpts[ipos, 0]:.17g}\t"
                f"{xpts[ipos, 1]:.17g}\n")

    if (mesh.vert3 is not None and
            mesh.vert3.size != +0):
        warnings.warn(
            "VERT3 elements not supported", Warning)

    if (mesh.edge2 is not None and
            mesh.edge2.size != +0):
        warnings.warn(
            "EDGE2 elements not supported", Warning)

    if (mesh.tria3 is not None):
        cell = mesh.tria3["index"]
        for ipos in range(mesh.tria3.size):
            fptr.write(
                f"{ipos}\t"
                f"{cell[ipos, 0] + 1}\t"
                f"{cell[ipos, 1] + 1}\t"
                f"{cell[ipos, 2] + 1}\n")

    if (mesh.quad4 is not None and
            mesh.quad4.size != +0):
        warnings.warn(
            "QUAD4 elements not supported", Warning)

    if (mesh.tria4 is not None and
            mesh.tria4.size != +0):
        warnings.warn(
            "TRIA4 elements not supported", Warning)

    if (mesh.hexa8 is not None and
            mesh.hexa8.size != +0):
        warnings.warn(
            "HEXA8 elements not supported", Warning)

    if (mesh.wedg6 is not None and
            mesh.wedg6.size != +0):
        warnings.warn(
            "WEDG6 elements not supported", Warning)

    if (mesh.pyra5 is not None and
            mesh.pyra5.size != +0):
        warnings.warn(
            "PYRA5 elements not supported", Warning)

    return


def save_grid_file(mesh, fptr):
    """
    SAVE-GRID-FILE: save a JIGSAW mesh object to *.DAT file.

    """

    if (mesh is not None):
        warnings.warn(
            "--GRID objects not supported", Warning)

    return


def savefvc(name, mesh):
    """
    SAVEFVC: save a JIGSAW MSH object to an FVCOM .DAT file.

    SAVEFVC(NAME, MESH)

    MESH is JIGSAW's primary mesh/grid/geom class. See MSH_t
    for details.

    Data in MESH is written as-needed -- any objects defined
    will be saved to file (if supported).

    Authors: Darren Engwirda

    """

    if (not isinstance(name, str)):
        raise Exception("Incorrect type: NAME.")

    if (not isinstance(mesh, jigsaw_msh_t)):
        raise Exception("Incorrect type: MESH.")

    certify(mesh)

    fext = Path(name).suffix

    if (fext.strip() != ".dat"): name += ".dat"

    kind = mesh.mshID.lower()

    with Path(name).open("w") as fptr:
    #----------------------------------- write JIGSAW object
        if   (kind == "euclidean-mesh"):

            save_mesh_file(mesh, fptr)

        elif (kind == "euclidean-grid"):

            save_grid_file(mesh, fptr)

        elif (kind == "ellipsoid-mesh"):

            save_mesh_file(mesh, fptr)

        elif (kind == "ellipsoid-grid"):

            save_grid_file(mesh, fptr)

        else:
            raise Exception(
                "MESH.mshID is not supported!!")

    return
