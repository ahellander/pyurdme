#!/usr/bin/env python
import pyurdme
from pyurdme.nsmsolver import NSMSolver
import numpy

from dolfin import *
import dolfin

class CoralReef(pyurdme.URDMEModel):
    """ Model of a coral reef, published in Sandin & McNamara (2012) """

    def __init__(self, name="coral_reef"):
        pyurdme.URDMEModel.__init__(self, name)

        Coral = pyurdme.Species(name="Coral",diffusion_constant=1.0,dimension=2)
        #Coral = pyurdme.Species(name="Coral",diffusion_constant=0.0,dimension=2)
        self.addSpecies([Coral])

        P_cr  = pyurdme.Parameter(name="P_cr", expression=0.01)
        G_c   = pyurdme.Parameter(name="G_c", expression=0.01)
        P_cd  = pyurdme.Parameter(name="P_cd", expression=0.15)
        C_th  = pyurdme.Parameter(name="C_th", expression=900)
        C_max = pyurdme.Parameter(name="C_max", expression=7000)
        

        # A unit square
        # each grid point is 10cm x 10cm, domain is 5m x 5m
        self.mesh = pyurdme.Mesh.SquareMesh(L=500, nx=50, ny=50, periodic=True)

        # Place the A molecules in the voxel nearest the center of the square
        self.placeNear({Coral:1000}, point=[0,0])

        self.timespan(numpy.linspace(0,500,500))

if __name__ == "__main__":
    reef = CoralReef()
    sol = NSMSolver(reef)
    result = sol.run()
    result.toVTK(species='Coral',folder_name="output_coral")
