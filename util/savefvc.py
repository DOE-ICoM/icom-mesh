
import warnings
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from pathlib import Path

from jigsawpy import jigsaw_msh_t, certify, orient1


def in_tri2(ppos, tri2, test, rtol):
    """
    IN-TRI2: return a T-by-1 array STAT, with STAT[I] = TRUE
    if TEST lies "inside" the I-TH triangle.

    """

    TEST = np.tile(test, (tri2.shape[0], 1))

    sgn1 = orient1(
        ppos[tri2[:, +0], :],
        ppos[tri2[:, +1], :], TEST
    )

    sgn2 = orient1(
        ppos[tri2[:, +1], :],
        ppos[tri2[:, +2], :], TEST
    )

    sgn3 = orient1(
        ppos[tri2[:, +2], :],
        ppos[tri2[:, +0], :], TEST
    )

    return np.logical_and.reduce((
        sgn1 * sgn2 >= -rtol,
        sgn2 * sgn3 >= -rtol,
        sgn3 * sgn1 >= -rtol))


def cullfvc(mesh, seed):
    """
    CULLFVC: clean-up a JIGSAW mesh obj. for export to FVCOM

    CULLFVC(MESH, SEED)

    Keep only those triangles in MESH that can be "reached"
    from the interior point(s) defined in SEED. Here, MESH
    is a standard msh_t object returned by JIGSAW and SEED
    is an S-by-2 array of points to cull from.

    Any edge elements defined in MESH are also culled-away.

    Authors: Darren Engwirda

    """

    certify(mesh)

    edge = np.empty((0, 3), dtype=np.int32)

    if (mesh.edge2 is not None and
            mesh.edge2.size != +0):
    #----------------------------------- destroy EDGE2 cells
        mesh.edge2 = np.empty(
            (+0), dtype=mesh.EDGE2_t)

    if (mesh.tria3 is not None and
            mesh.tria3.size != +0):
    #----------------------------------- connect TRIA3 cells
        cell = mesh.tria3["index"][:]

        indx = np.arange(0, cell.shape[0])
        indx = np.reshape(indx, (indx.size, 1))

        evec = np.sort(cell[:, (0, 1)], axis=1)
        evec = np.concatenate(
            (evec, indx), axis=1)
        edge = np.concatenate(
            (edge, evec), axis=0)

        evec = np.sort(cell[:, (1, 2)], axis=1)
        evec = np.concatenate(
            (evec, indx), axis=1)
        edge = np.concatenate(
            (edge, evec), axis=0)

        evec = np.sort(cell[:, (2, 0)], axis=1)
        evec = np.concatenate(
            (evec, indx), axis=1)
        edge = np.concatenate(
            (edge, evec), axis=0)

    if (edge.size != 0):
    #----------------------------------- get connected cells
        indx = np.lexsort((edge[:, 0], edge[:, 1]))
        edge = edge[indx + 0, :]

        diff = np.diff(edge[:, 0:2], axis=0)

        indx = np.argwhere(np.all(diff==0, axis=1))

        itri = edge[indx + 0, 2]
        jtri = edge[indx + 1, 2]

        rows = np.concatenate(
            (itri, jtri), axis=0).flatten()
        cols = np.concatenate(
            (jtri, itri), axis=0).flatten()
        data = np.ones((rows.size), dtype=np.int32)

        conn = csr_matrix((data, (rows, cols)))

        ncon, cidx = connected_components(
            conn, 
            return_labels=True, directed=False)

    #----------------------------------- reachable via seed?
        keep = np.full(
            mesh.tria3.size, False, dtype=bool)
        used = np.full(
            mesh.point.size, False, dtype=bool)

        pmax = np.max(mesh.point["coord"], axis=0)
        pmin = np.min(mesh.point["coord"], axis=0)
        rtol = +1.0E-12 * np.prod(pmax - pmin)

        for spos in seed:

            indx = np.argwhere(in_tri2(
                mesh.point["coord"],
                mesh.tria3["index"], spos, rtol)
            )

            if (indx.size != 0):

                keep[cidx == cidx[indx[0]]] = True

    #----------------------------------- set reachable parts
        mesh.tria3 = mesh.tria3[keep]

        used[mesh.tria3["index"].flatten()] = True

        redo = np.zeros(
            (mesh.point.size), dtype=np.int32)
        redo[used] = np.arange(
            0, np.count_nonzero(used)
        )

        mesh.tria3["index"] = \
            redo[mesh.tria3["index"]]

        mesh.point = mesh.point[used]

    return


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
