""" pyURDME model file for the polarization 1D example. """

import os
from pyurdme.pyurdme import *
import dolfin


class PheromoneGradient(pyurdme.URDMEDataFunction):
    L_max = 4
    L_min = 0
    def map(self, x):
        #ligand_c[i] = ( (L_max-L_min)*.5*(1+cos( .5*(i*l - 2*3.14159))) + L_min)*MOLAR ;
        #  x[0] == i*l


class G_protein_cycle_1D(URDMEModel):

    def __init__(self,model_name="polarization3D"):
        URDMEModel.__init__(self,model_name)

        # Species
        # R RL G Ga Gbg Gd
        R   = Species(name="R",  diffusion_constant=5.3e-15,dimension=2)
        RL  = Species(name="RL", diffusion_constant=2.5e-12,dimension=2)
        G   = Species(name="G",  diffusion_constant=2.5e-12,dimension=2)
        Ga  = Species(name="Ga", diffusion_constant=5.3e-15,dimension=2)
        Gbg = Species(name="Gbg",diffusion_constant=5.3e-15,dimension=2)
        Gd  = Species(name="Gd", diffusion_constant=5.3e-15,dimension=2)
        
        self.addSpecies([R,RL,G,Ga,Gbg,Gd])
    

        L = 4*3.14159
        NUM_VOXEL = 200
        self.mesh = pyurdme.Mesh.IntervalMesh(nx=NUM_VOXEL, a=-2*3.14159, b=2*3.14159)
        
        MOLAR = Parameter(name="MOLAR",expression=6.02e-01*((L/NUM_VOXEL)**3))
        SA    = Parameter(name="SA" ,expression=201.056)
        V     = Parameter(name="V" ,expression=33.5)
        k_RL  = Parameter(name="k_RL" ,expression="2e-03/MOLAR")
        k_RLm = Parameter(name="k_RLm" ,expression=1e-02)
        k_Rs  = Parameter(name="k_Rs" ,expression="4.0/SA")
        k_Rd0 = Parameter(name="k_Rd0" ,expression=4e-04)
        k_Rd1 = Parameter(name="k_Rd1" ,expression=4e-04)
        k_G1  = Parameter(name="k_G1" ,expression="1.0*SA")
        k_Ga  = Parameter(name="k_Ga" ,expression="1e-06*SA")
        k_Gd  = Parameter(name="k_Gd" ,expression=0.1)
        
        self.addParameter([MOLAR,SA,V,k_RL,k_RLm,k_Rs,k_Rd0,k_Rd1,k_G1,k_Ga,k_Gd]) 

        # Reactions
        R1 = Reaction(name="R1",reactants={MinD_c_atp:1},products={MinD_m:1},massaction=True,rate=sigma_d, restrict_to=boundary)
        
        self.addReaction([R1,R2,R3,R4,R5,R6])
        
        # Distribute molecules randomly over the mesh according to their initial values
        self.scatter({R:10000})
        self.scatter({G:10000})

        self.timespan(range(1000))


if __name__=="__main__":
    """ Dump model to a file. """
                     
    model = G_protein_cycle_1D()
    result = urdme(model)
    print result
