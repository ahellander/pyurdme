""" pyurdme model file for reversible dimerization on a surface. """
import os
from pyurdme.urdme import *
import dolfin
import numpy
import scipy.io as spio

#from pyurdme.gmsh import *

def dimerization(model_name=""):
    """ Dimerization. The reversible reaction A+B<->C on a surface. """
    
    if model_name == "":
        model_name = "dimerization"
    model = URDMEModel(name=model_name);

    # Species
    A = Species(name="A",initial_value=0,reaction_radius=1e-9,diffusion_constant=0.01,dimension=2);
    B = Species(name="B",initial_value=0,reaction_radius=1e-9,diffusion_constant=0.01,dimension=2);
    C = Species(name="C",initial_value=100,reaction_radius=1e-9,diffusion_constant=0.01,dimension=2);

    model.addSpecies([A,B,C])

    # Parameters
    NA  = Parameter(name="NA",expression=6.022e23)
    #k1  = Parameter(name="k1",expression="1.0e6/(1000*NA)")
    #k2  = Parameter(name="k2",expression=1.0)
    k1  = Parameter(name="k1",expression=0.0)
    k2  = Parameter(name="k2",expression=0.0)

    model.addParameter([NA,k1,k2])

    # Reactions
    R1 = Reaction(name="R1",reactants={A:1,B:1},products={C:1},massaction=True,rate=k1)
    R2 = Reaction(name="R2",reactants={C:1},products={A:1,B:1},massaction=True,rate=k2)
    model.addReaction([R1,R2])

    # A square domain with Cartesian discretization
    #model.mesh  = CartesianMesh(geometry="line",side_length=1e-6,hmax=1e-7)

    # We could wrap around Gmsh like this, if we wanted to:
    #model.geometry = gmshGeometry(file='meshes/surface.geo')
    #meshinit(model.geometry)
    model.mesh = read_dolfin_mesh('meshes/surface.xml')

    #mesh,physical_ids=meshInit(geom) -> Shells out and calls Gmsh
    # Or simply load a Gmsh mesh
    #model.mesh = read_gmsh('meshes/surface.msh')
    #meshextend(model)

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
    #result = urdme(model,solver='nem',solver_path="/Users/andreash/bitbucket/nllattice/",seed=10)
    result = urdme(model,solver='nsm',seed=10)
    U = result["U"]
    print numpy.shape(U)
    print numpy.sum(U[:,3])
    #model.serialize(filename="debug.mat")
    # Plot using VIPER
    #dolfin.plot(model.sol['C'],wireframe=True)
    #dolfin.plot(model.mesh.mesh,wireframe=True,axes=True)
    #dolfin.interactive()
    
    spio.savemat("debugoutput.mat",result)
    # Dump solution to file in VTK format for ParaView
    file = dolfin.File("testsolution.pvd")
    file << model.sol['C']
    
    #print result
