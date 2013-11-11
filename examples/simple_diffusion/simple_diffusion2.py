""" pyurdme model for a single species diffusing on the unit square. """

import dolfin
from pyurdme.pyurdme import *
import json
import pickle

def simple_diffusion2():
    """ One species diffusing on the boundary of a sphere and one species
        diffusing inside the sphere. """

    model = URDMEModel(name="simple_diffusion2")

    D = 0.01
    A = Species(name="A",initial_value=1000,diffusion_constant=D,dimension=2)
    B = Species(name="B",initial_value=1000,diffusion_constant=0.1*D,dimension=2)


    model.addSpecies([A,B])

    # A unit square with nx points on the x-axis and ny points on the y-axis
    
    c1 = dolfin.Circle(0,0,1)
    #c2 = dolfin.Circle(0,0,0.3)
    mesh = dolfin.Mesh(c1,20)
    model.mesh = Mesh(mesh)
    #dolfin.plot(mesh)
    #dolfin.interactive()
    
    meshextend(model)
    
    model.timespan(numpy.linspace(0,1,50))
    meshextend(model)
    
    # Place the A molecules in the voxel nearest to the center of the square
    model.placeNear(species=A,point=[0.5,0.5])
    model.scatter(species=A,subdomain=1)
    
    return model

if __name__ == '__main__':
    model = simple_diffusion2()
    model.initialize()
    
    params = {'model':model,
              'solver_data':{'tspan':numpy.linspace(0,1,50)}}
    
    solver_data = {'tspan':numpy.linspace(0,1,50)}
    
    #pickle.dumps(model)
    
    result = urdme(model,solver_data)

    # Dump a snapshot of the state in paraview format
    file1 = dolfin.File("A.pvd")
    file1 << model.sol['A']
    file2 = dolfin.File("B.pvd")
    file2 << model.sol['A']

