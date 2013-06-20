""" pyURDME model file for the MinCDE example. """

import os
from pyurdme.urdme import *

def mincde(model_name=""):

	if model_name == "":
		model_name = "mincde"
	model = URDMEModel(name=model_name)

	# Species
	# TODO: We need a way to localize species to subdomains/boundaries
	MinD_m = Species(name="MinD_m",initial_value=0,diffusion_constant=1e-14);
	MinD_c_atp = Species(name="MinD_c_atp",initial_value=0,diffusion_constant=2.5e-12);
	MinD_c_adp = Species(name="MinD_c_adp",initial_value=4500,diffusion_constant=2.5e-12);
	MinD_e = Species(name="MinD_e",initial_value=1575,diffusion_constant=2.5e-12);
	MinDE = Species(name="MinDE",initial_value=0,diffusion_constant=1e-14);
	
	model.addSpecies([MinD_m,MinD_c_atp,MinD_c_adp,MinD_e,MinDE])
	
	# Parameters 
	sigma_d = Parameter(name="sigma_d",expression=2.5e-8)
	sigma_dD = Parameter(name="sigma_dD",expression=0.0016e-18)
	sigma_e = Parameter(name="sigma_e",expression=0.093e-18)
	sigma_de = Parameter(name="sigma_de",expression=0.7)
	sigma_dt = Parameter(name="sigma_dt",expression=1.0)
	
	model.addParameter([sigma_d,sigma_dD,sigma_e,sigma_de,sigma_dt])

	# Reactions
	# TODO: These reactions will serialize without volume dependency. This needs to be adressed in the model.py module.
    # How do we get the voxel volumes into the propensity functions in the .c model file ?
	R1 = Reaction(name="R1",reactants={MinD_c_atp:1},products={MinD_m:1},massaction=True,rate=sigma_d)
    R2 = Reaction(name="R2",reactants={MinD_c_atp:1,MinD_m:1},products={MinD_m:2},massaction=True,rate=sigma_dD)
	R3 = Reaction(name="R3",reactants={MinD_m:1,MinD_e:1},products={MinDE:1},massaction=True,rate=sigma_e)
    R4 = Reaction(name="R4",reactants={MinDE:1,MinD_e:1},products={MinD_c_adp:1,MinD_e:1},massaction=True,rate=sigma_de)
	R5 = Reaction(name="R5",reactants={MinD_c_adp:1},products={MinD_c_atp:1},massaction=True,rate=sigma_dt)
	R6 = Reaction(name="R6",reactants={MinDE:1,MinD_c_atp:1},products={MinD_m:1,MinDE:1},massaction=True,rate=sigma_dD)
	
	model.addReaction([R1,R2,R3,R4,R5,R6])
	
	# Load mesh
    model.mesh = read_gmsh_mesh("coli.xml")
    
	return model

if __name__=="__main__":
	""" Dump model to a file. """
	model = mincde(model_name="mincde")
#model.serialize("testmincde.mat")
	# TODO: This does not work yet. 
	#urdme(model)
	