
import numpy as np
from scipy import sparse

import jigsawpy


def in_tri2(ppos, tri2, test, rtol):
    """
    IN-TRI2: return a T-by-1 array STAT, with STAT[I] = TRUE
    if TEST lies "inside" the I-TH triangle.

    """

    TEST = np.tile(test, (tri2.shape[0], 1))

    sgn1 = jigsawpy.orient1(
        ppos[tri2[:, +0], :],
        ppos[tri2[:, +1], :], TEST
    )

    sgn2 = jigsawpy.orient1(
        ppos[tri2[:, +1], :],
        ppos[tri2[:, +2], :], TEST
    )

    sgn3 = jigsawpy.orient1(
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

    """

    jigsawpy.certify(mesh)

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

        conn = sparse.csr_matrix(
            (data, (rows, cols)))

        ncon, cidx = \
            sparse.csgraph.connected_components(
                conn, directed=False, return_labels=True)

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
