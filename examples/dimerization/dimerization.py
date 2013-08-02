""" pyurdme model file for reversible dimerization on a surface. """
import os
from pyurdme.urdme import *
import dolfin
import numpy
import scipy.io as spio

def dimerization(model_name=""):
    """ Dimerization. The reversible reaction A+B<->C on a surface. """
    
    if model_name == "":
        model_name = "dimerization"
    model = URDMEModel(name=model_name);

    # Species
    A = Species(name="A",initial_value=50,reaction_radius=1e-9,diffusion_constant=1,dimension=2);
    B = Species(name="B",initial_value=50,reaction_radius=1e-9,diffusion_constant=0.1,dimension=2);
    C = Species(name="C",initial_value=100,reaction_radius=1e-9,diffusion_constant=2,dimension=2);

    model.addSpecies([A,B,C])

    # Parameters
    k1  = Parameter(name="k1",expression=10.0)
    k2  = Parameter(name="k2",expression=0.0)

    model.addParameter([k1,k2])

    # Reactions
    R1 = Reaction(name="R1",reactants={A:1,B:1},products={C:1},massaction=True,rate=k1)
    R2 = Reaction(name="R2",reactants={C:1},products={A:1,B:1},massaction=True,rate=k2)
    model.addReaction([R1,R2])

    model.mesh = read_dolfin_mesh('meshes/surface.xml')
            
    # Distribute the molecues over the mesh
    model.scatter(A,subdomain=1)
    model.scatter(B,subdomain=1)
    model.scatter(C,subdomain=1)

    # Time span of the simulation
    model.tspan = range(100)

    return model


if __name__ == '__main__':
    """ Create a model and assemble the URDME input file. """
    
    model = dimerization()
    model.initialize()
    
    result = urdme(model,solver='nem',solver_path="/Users/andreash/bitbucket/nllattice/", seed=10)
    #result = urdme(model,solver='nsm',seed=10)
    
    U = result["U"]
    
    # Plot using VIPER
    #dolfin.plot(model.sol['C'],wireframe=True)
    #dolfin.plot(model.mesh.mesh,wireframe=True,axes=True)
    #dolfin.interactive()
    #spio.savemat("debugoutput.mat",result)
    
    # Dump solution to file in VTK format for ParaView
    file = dolfin.File("testsolution.pvd")
    file << model.sol['C']
    toXYZ(model,"testsolution.xyz")
