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


def F(x):
    f = (4*math.log(1/x)-(1-x*x)*(3-x*x))/(4*(1-x*x)*(1-x*x))
    return f

def mesoreac2D(rho,vol,gamma,kr):
    R = math.sqrt(vol/math.pi)
    lam = rho/R
    alpha = kr/(2*math.pi*gamma)
    return math.pi*R*R/kr*(1+alpha*F(lam))

class Example2(pyurdme.URDMEModel):
    """ The reversible reaction A+B<->C in 2D.  """
    
    def __init__(self,voxel_size=0.2e-6, val=1e-18, sigma=1e-9):
        
        pyurdme.URDMEModel.__init__(self,name="Example2")

        gamma = 1e-12

        # Substrate and enzymes
        S1 = pyurdme.Species(name="S1",reaction_radius=sigma,diffusion_constant=gamma)
        S11 = pyurdme.Species(name="S11",reaction_radius=sigma,diffusion_constant=gamma)
        S12 = pyurdme.Species(name="S12",reaction_radius=sigma,diffusion_constant=gamma)
        S2  = pyurdme.Species(name="S2",reaction_radius=sigma,diffusion_constant=gamma)

        self.add_species([S1,S11,S12,S2])

        # Mesoscopic association and dissacosiation rate according to Hellander et. al.
        D = S11.diffusion_constant+S12.diffusion_constant
        rho = S11.reaction_radius+S12.reaction_radius
        beta = 1.5164
        
        # Microscopic association and disassociation rate
        # val = 2*3.1416*A.diffusion_constant
        kr_micro  = pyurdme.Parameter(name="krm",expression=val)

        pi = math.pi

        L = 1e-6
        h = voxel_size
        nx = int(L/h)

        N = nx*nx*nx
        self.mesh = pyurdme.URDMEMesh.generate_cube_mesh(L,nx,nx,nx)
        
        ka_kc = 4.0*pi*rho*D*kr_micro.value/(4.0*pi*rho*D+kr_micro.value)
        ka_meso_val = pyurdme.Parameter(name="ka_meso_val",expression=ka_kc/(1-beta*ka_kc/(6*D*h)))
        ka_meso  = pyurdme.Parameter(name="ka_meso",expression="ka_meso_val")
        kd  = pyurdme.Parameter(name="kd",expression=10.0)
            
        # Reactions
        R1 = pyurdme.Reaction(name="R1",reactants={S1:1},products={S11:1,S12:1},massaction=True, rate=kd)
        R2 = pyurdme.Reaction(name="R2",reactants={S11:1,S12:1},products={S2:1},massaction=True, rate=ka_meso)
    
        self.add_parameter([ka_meso_val,kr_micro,kd,ka_meso])
        self.add_reaction([R1,R2])
        
        # Distribute the molecules over the mesh
        self.set_initial_condition_scatter({S1:100})

        # Time span of the simulation
        self.timespan(numpy.arange(0,2.0,0.05))
 

if __name__ == '__main__':

    from pyurdme.microsolver import MMMSSolver 

    model = Example2(voxel_size=0.2e-6)
    solver = MMMSSolver(model, min_micro_timestep=1e-4)

    # To use the meso-micro hybrid solver, simply specify the species partitioning.   
    solver.set_modeling_level({"S1":"micro", "S11":"micro", "S12":"meso", "S2":"meso"})
    solver.create_input_file("test.txt")
    # TODO1: Automatically determine this partitioning based on a priori error estimate
    # TODO2: Allow the partitioning to be based on subdomain, or based on a URDMEDataFunction
     
    #solver._write_mesh_file("test_mesh.h5")
    #solver.serialize("urdmeinput.mat")
    #solver.create_input_file("test.txt")
    microres = solver.run()
    print microres.get_particles(0,0)
