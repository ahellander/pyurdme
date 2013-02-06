""" pyURMDE model file for the MinCDE example. """

import os
from pyurdme.urdme import *

def mincde(model_name=""):

	if model_name == "":
		model_name = "dimerization"
	model = URDMEModel(name=model_name)

	# Species
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
    sigma_de = Paramter(name="sigma_de",expression=0.7)
    sigma_dt = Parameter(name="sigma_dt",expression=1.0)
    model.addParameter([sigma_d,sigma_dD,sigma_e,sigma_de,sigma_dt])

	# Reactions

	return model

if __name__=="__main__":
	""" Dump model to a file. """
	model = mincde(model_name="mincde")	
	#model.serialize(filename="mincdeinput.mat")
