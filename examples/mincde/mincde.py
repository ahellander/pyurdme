""" pyURDME model file for the MinCDE example. """

import os
from pyurdme.urdme import *

def mincde(model_name=""):

    if model_name == "":
        model_name = "mincde"
    
    model = URDMEModel(name=model_name)


    # Load mesh
    model.mesh = read_dolfin_mesh("mesh/coli.xml")

    # Some kind of average mesh size
    hmax = model.mesh.mesh.hmax()
    hmin = model.mesh.mesh.hmin()
    h = (hmax+hmin)/2

    # Read the facet and interior cell physical domain markers into a Dolfin MeshFunction 
    file_in = dolfin.File("mesh/coli_facet_region.xml")
    facet_function = dolfin.MeshFunction("size_t",model.mesh.mesh)
    file_in >> facet_function

    file_in = dolfin.File("mesh/coli_physical_region.xml")
    physical_region = dolfin.MeshFunction("size_t",model.mesh.mesh)
    file_in >> physical_region

    #dx = dolfin.Measure("dx")[physical_region]
    #ds = dolfin.Measure("ds")[facet_function]

    # List of Physical domain markers that match those in the  Gmsh .geo file.
    boundary = [73,74,75]
    interior = [76]

    # Species
    # TODO: We need a way to localize species to subdomains/boundaries
    MinD_m     = Species(name="MinD_m",initial_value=0,diffusion_constant=1e-14,dimension=2,active_on=boundary)
    MinD_c_atp = Species(name="MinD_c_atp",initial_value=0,diffusion_constant=2.5e-12)
    MinD_c_adp = Species(name="MinD_c_adp",initial_value=4500,diffusion_constant=2.5e-12)
    MinD_e     = Species(name="MinD_e",initial_value=1575,diffusion_constant=2.5e-12)
    MinDE      = Species(name="MinDE",initial_value=0,diffusion_constant=1e-14,dimension=2,active_on=boundary)
        
    model.addSpecies([MinD_m,MinD_c_atp,MinD_c_adp,MinD_e,MinDE])

    # Parameters
    NA = Parameter(name="NA",expression="6.022e23")
    sigma_d  = Parameter(name="sigma_d",expression=2.5e-8/h)
    sigma_dD = Parameter(name="sigma_dD",expression=0.0016e-18)
    sigma_e  = Parameter(name="sigma_e",expression=0.093e-18)
    sigma_de = Parameter(name="sigma_de",expression=0.7)
    sigma_dt = Parameter(name="sigma_dt",expression=1.0)
        
    model.addParameter([NA,sigma_d,sigma_dD,sigma_e,sigma_de,sigma_dt])

    # Reactions
    R1 = Reaction(name="R1",reactants={MinD_c_atp:1},products={MinD_m:1},massaction=True,rate=sigma_d, restrict_to=boundary)
    R2 = Reaction(name="R2",reactants={MinD_c_atp:1,MinD_m:1},products={MinD_m:2},massaction=True,rate=sigma_dD)
    R3 = Reaction(name="R3",reactants={MinD_m:1,MinD_e:1},products={MinDE:1},massaction=True,rate=sigma_e)
    R4 = Reaction(name="R4",reactants={MinDE:1,MinD_e:1},products={MinD_c_adp:1,MinD_e:1},massaction=True,rate=sigma_de)
    R5 = Reaction(name="R5",reactants={MinD_c_adp:1},products={MinD_c_atp:1},massaction=True,rate=sigma_dt)
    R6 = Reaction(name="R6",reactants={MinDE:1,MinD_c_atp:1},products={MinD_m:1,MinDE:1},massaction=True,rate=sigma_dD)
    
    model.addReaction([R1,R2,R3,R4,R5,R6])
    
    # Distribute molecules randomly over the mesh according to their initial values
    model.scatter(MinD_c_adp)
    model.scatter(MinD_e)

    model.tspan = range(100)
    return model

if __name__=="__main__":
    """ Dump model to a file. """
                     
    model = mincde(model_name="mincde")
    model.initialize()

    result = urdme(model)
    
    file = dolfin.File("testsolution.pvd")
    file << model.sol['MinD_m']
    toXYZ(model,"testsolution.xyz")

    print result
