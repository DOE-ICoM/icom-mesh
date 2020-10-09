
from jigsawpy import jigsaw_msh_t, certify, savevtk


def saveats(name, mesh):

    savevtk(name + ".vtk", mesh)

    return
