""" pyurdme model file for reversible dimerization on a surface. """
import os
from pyurdme.urdme import *
import dolfin
import numpy
import scipy.io as spio

def dimerization(model_name=""):
    """ Dimerization. The reversible reaction A+B<->C on a surface. """
    # Units are scaled to be in \mu m
    
    if model_name == "":
        model_name = "dimerization"
    model = URDMEModel(name=model_name);

    # Species
    A = Species(name="A",initial_value=50,reaction_radius=1e-3,diffusion_constant=1e-2,dimension=2);
    B = Species(name="B",initial_value=50,reaction_radius=0.5e-3,diffusion_constant=1e-2,dimension=2);
    C = Species(name="C",initial_value=100,reaction_radius=0.2e-3,diffusion_constant=1e-2,dimension=2);

    model.addSpecies([A,B,C])

    # Parameters
    k1  = Parameter(name="k1",expression=100.0)
    k2  = Parameter(name="k2",expression=0.0)

    model.addParameter([k1,k2])

    # Reactions
    R1 = Reaction(name="R1",reactants={A:1,B:1},products={C:1},massaction=True,rate=k1)
    R2 = Reaction(name="R2",reactants={C:1},products={A:1,B:1},massaction=True,rate=k2)
    model.addReaction([R1,R2])

    model.mesh = read_dolfin_mesh('meshes/plane.xml')
      
    # Assemble the subdomain vector.

    # Read the facet and interior cell physical domain markers into a Dolfin MeshFunction
    file_in = dolfin.File("meshes/plane_facet_region.xml")
    facet_function = dolfin.MeshFunction("size_t",model.mesh.mesh)
    file_in >> facet_function

    file_in = dolfin.File("meshes/plane_physical_region.xml")
    physical_region = dolfin.MeshFunction("size_t",model.mesh.mesh)
    file_in >> physical_region

    #print model.mesh.mesh.num_facets()
    #topology = model.mesh.mesh.topology()
    #print topology.dim
    #print topology.global_indices(0)
    #exit(-1)

    # These are the physical id lables from the Gmsh file.
    interior = 7
    boundary = 8

    
    # Distribute the molecues over the mesh
    model.scatter(A,subdomain=1)
    model.scatter(B,subdomain=1)
    model.scatter(C,subdomain=1)

    # Time span of the simulation
    model.timespan(numpy.arange(0,25,0.1))

    return model


if __name__ == '__main__':
    """ Create a model and assemble the URDME input file. """
    
    model = dimerization()
    model.initialize()
    

    #result = urdme(model,solver='nem',solver_path="/Users/andreash/bitbucket/nllattice/", seed=10)
    result = urdme(model,solver='nsm',seed=10)
    
    U = result["U"]
    
    print "A" + str(numpy.sum(U[::3,:],axis=0))
    
    print "B" + str(numpy.sum(U[1::3,:],axis=0))
    print "C" + str(numpy.sum(U[2::3,:],axis=0))
    
    # Plot using VIPER
    #dolfin.plot(model.sol['C'],wireframe=True)
    #dolfin.plot(model.mesh.mesh,wireframe=True,axes=True)
    #dolfin.interactive()
    #spio.savemat("debugoutput.mat",result)
    
    # TODO: Wrap and make it be something like this instead
    #urdme.save_solution(model,filename="",format="xyz")
    #urdme.plot_solution(model, species="C",time="0")
    
    # Dump solution to file in VTK format for ParaView
    file = dolfin.File("testsolution.pvd")
    file << model.sol['C']
    #toCSV(model,"solution1")
