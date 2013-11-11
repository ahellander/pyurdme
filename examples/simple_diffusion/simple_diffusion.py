""" pyurdme model for a single species diffusing on the unit square. """

import dolfin
from pyurdme.pyurdme import *
import pickle
import json

class simple_diffusion(URDMEModel):
    """ Initial condition is a delta function at the center voxel. 
        The solution should be a Gaussian, up to the point where
        the BC becomes important. """

    def __init__(self):
    
        URDMEModel.__init__(self,name="simple_diffusion")

        D = 0.01
        A = Species(name="A",initial_value=1000000,diffusion_constant=D,dimension=2)

        self.addSpecies([A])

        # A unit square
        self.mesh = unitSquareMesh(40,40)
        
        # Place the A molecules in the voxel nearest the center of the square
        self.placeNear(species=A,point=[0.5,0.5])


if __name__ == '__main__':

    model = simple_diffusion()
    
    solver_options = {'tspan':numpy.linspace(0,1,50),
                       'seed':1432423}
    
    #solver_data = model.solverData()

    result = urdme(model)

    # Dump a snapshot of the state in paraview format
    file = dolfin.File("A.pvd")
    file << model.sol['A']

