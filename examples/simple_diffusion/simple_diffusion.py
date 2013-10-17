""" pyurdme model for a single species diffusing on the unit square. """

import dolfin
from pyurdme.pyurdme import *

def simple_diffusion():
    """ Initial condition is a delta function at the center voxel. 
        The solution should be a Gaussian, up to the point where
        the BC becomes important. """

    model = URDMEModel(name="simple_diffusion")

    D = 0.01
    A = Species(name="A",initial_value=1000000,diffusion_constant=D,dimension=2)

    model.addSpecies([A])

    # A unit square with nx points on the x-axis and ny points on the y-axis
    nx = 40
    ny = 40
    L  = 0.5e-6
    
    #c1 = dolfin.Circle(0,0,1)
    #c2 = dolfin.Circle(0,0,0.3)
    #mesh = dolfin.Mesh(c1-c2,20)
    #dolfin.plot(mesh)
    #dolfin.interactive()
    model.mesh = unitSquareMesh(nx,ny)

    model.timespan(numpy.linspace(0,1,50))
    
    # Place the A molecules in the voxel nearest to the center of the square
    model.placeNear(species=A,point=[0.5,0.5])
    
    return model

if __name__ == '__main__':
    model = simple_diffusion()
    result = urdme(model)

    # Dump a snapshot of the state in paraview format
    file = dolfin.File("A.pvd")
    file << model.sol['A']

