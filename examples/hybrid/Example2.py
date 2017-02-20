import os
import pyurdme
import pyurdme.nsmsolver
import dolfin
import numpy
import scipy.io as spio
import matplotlib.pyplot as plt
import math
from scipy import integrate
import datetime as dt

class Example2(pyurdme.URDMEModel):
    """ The reversible reaction A+B <->C in 3D.  """
    
    def __init__(self,voxel_size=0.1e-6, val=1e-20, sigma=1e-9):
        
        pyurdme.URDMEModel.__init__(self,name="Example2")

        gamma = 1e-12

        # Substrate and enzymes
        S1  = pyurdme.Species(name="S1",reaction_radius=sigma,diffusion_constant=gamma)
        S11 = pyurdme.Species(name="S11",reaction_radius=sigma,diffusion_constant=gamma)
        S12 = pyurdme.Species(name="S12",reaction_radius=sigma,diffusion_constant=gamma)
        S2  = pyurdme.Species(name="S2",reaction_radius=sigma,diffusion_constant=gamma)
        self.add_species([S1,S11,S12,S2])

        pi = math.pi
        self.voxel_size = voxel_size
        L = 1e-6
        h = voxel_size
        nx = int(L/h)

        N = nx*nx*nx
        self.mesh = pyurdme.URDMEMesh.generate_cube_mesh(L,nx,nx,nx)
       
        # Microscopic association and disassociation rate
        kr  = pyurdme.Parameter(name="kr",expression=val)
        kd  = pyurdme.Parameter(name="kd",expression=10.0)

        self.add_parameter([kr,kd])
    
        # Reactions
        R1 = pyurdme.Reaction(name="R1",reactants={S1:1},products={S11:1,S12:1},massaction=True, rate=kd)
        R2 = pyurdme.Reaction(name="R2",reactants={S11:1,S12:1},products={S2:1},massaction=True, rate=kr)
        self.add_reaction([R1,R2])
        
        # Scatter the molecules over the mesh
        self.set_initial_condition_scatter({S1:100})

        # Time span of the simulation
        self.timespan(numpy.arange(0,2.0,0.05))
 

if __name__ == '__main__':

    from pyurdme.microsolver import MMMSSolver 

    model = Example2(voxel_size=0.3e-6)
    solver = MMMSSolver(model, min_micro_timestep=1e-4)
    solver.create_input_file("Example2.json")

    #res = solver.propose_mesh_resolution_per_reaction(rel_tol=0.05)
    #solver.partition_system(rel_tol=0.05)

    # To use the meso-micro hybrid solver, simply specify the initial species partitioning.   
    solver.set_modeling_level({"S1":"micro", "S11":"micro", "S12":"meso", "S2":"meso"})
    #microres = solver.run()
    #print microres.get_particles(0,0)
