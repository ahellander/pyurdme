#!/usr/bin/env python
import pyurdme
from pyurdme.nsmsolver import NSMSolver
import numpy
from dolfin import *
import dolfin

#class unitSquareMeshPeriodicBoundary(dolfin.SubDomain):
class unitSquareMeshPeriodicBoundary(SubDomain):
    """ Sub domain for Periodic boundary condition """
    def __init__(self, Lx=1.0, Ly=1.0):
        self.Lx = Lx
        self.Ly = Ly

    def inside(self, x, on_boundary):
        """ Left boundary is "target domain" G """
        # return True if on left or bottom boundary AND NOT on one of the two corners (0, 1) and (1, 0)
        return bool((dolfin.near(x[0], 0) or dolfin.near(x[1], 0)) and
                (not ((dolfin.near(x[0], 0) and dolfin.near(x[1], 1)) or
                        (dolfin.near(x[0], 1) and dolfin.near(x[1], 0)))) and on_boundary)

    def map(self, x, y):
        ''' # Map right boundary G (x) to left boundary H (y) '''
        if dolfin.near(x[0], self.Lx) and dolfin.near(x[1], self.Ly):
            y[0] = x[0] - self.Lx
            y[1] = x[1] - self.Ly
        elif dolfin.near(x[0], self.Lx):
            y[0] = x[0] - self.Lx
            y[1] = x[1]
        else:   # near(x[1], 1)
            y[0] = x[0]
            y[1] = x[1] - self.Ly


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
        self.mesh = pyurdme.Mesh.SquareMesh(L=500, nx=50, ny=50)
        self.mesh.constrained_domain =  unitSquareMeshPeriodicBoundary()

        # Place the A molecules in the voxel nearest the center of the square
        self.placeNear({Coral:1000}, point=[0.,0.])

        self.timespan(numpy.linspace(0,1,50))        

if __name__ == "__main__":
    reef = CoralReef()
    sol = NSMSolver(reef)
    result = sol.run()