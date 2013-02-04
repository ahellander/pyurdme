""" pyURMDE model file for the MinCDE example. """

import os
from pyurdme.urdme import *

def mincde(model_name=""):
	
	if model_name == "":
        model_name = "dimerization"
    model = URDMEModel(name=model_name);

	# Species
    MinD_m = Species(name="MinD_m",initial_value=0,diffusion_constant=1e-14);
    MinD_c_atp = Species(name="MinD_c_atp",initial_value=0,diffusion_constant=2.5e-12);
    MinD_c_adp = Species(name="MinD_c_adp",initial_value=4500,diffusion_constant=2.5e-12);
    MinD_e = Species(name="MinD_e",initial_value=1575,diffusion_constant=2.5e-12);
    MinDE = Species(name="MinDE",initial_value=0,diffusion_constant=1e-14);
	
	# Parameters 
	
	# Reactions
	
