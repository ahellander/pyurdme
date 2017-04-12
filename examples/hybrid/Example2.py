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
    
    def __init__(self,voxel_size=0.1e-6, val=1e-21, sigma=1e-9):
        
        pyurdme.URDMEModel.__init__(self,name="Example2")

        gamma = 1e-12

        # Substrate and enzymes
        S1  = pyurdme.Species(name="S1",reaction_radius=sigma,diffusion_constant=gamma)
        S2  = pyurdme.Species(name="S2",reaction_radius=sigma,diffusion_constant=gamma)
        S3  = pyurdme.Species(name="S3",reaction_radius=sigma,diffusion_constant=gamma)
        S4  = pyurdme.Species(name="S4",reaction_radius=sigma,diffusion_constant=gamma)
        S11 = pyurdme.Species(name="S11",reaction_radius=sigma,diffusion_constant=gamma)
        S12 = pyurdme.Species(name="S12",reaction_radius=sigma,diffusion_constant=gamma)
        S21 = pyurdme.Species(name="S21",reaction_radius=sigma,diffusion_constant=gamma)
        S22 = pyurdme.Species(name="S22",reaction_radius=sigma,diffusion_constant=gamma)
        S31 = pyurdme.Species(name="S31",reaction_radius=sigma,diffusion_constant=gamma)
        S32 = pyurdme.Species(name="S32",reaction_radius=sigma,diffusion_constant=gamma)

        self.add_species([S1,S2,S3,S4,S11,S12,S21,S22,S31,S32])

        pi = math.pi
        self.voxel_size = voxel_size
        L = 1e-6
        h = voxel_size
        nx = int(L/h)

        N = nx*nx*nx
        self.mesh = pyurdme.URDMEMesh.generate_cube_mesh(L,nx,nx,nx)
       
        # Microscopic association and disassociation rate
        k11  = pyurdme.Parameter(name="k11",expression=1.0)
        k12  = pyurdme.Parameter(name="k12",expression=val)
        k13  = pyurdme.Parameter(name="k13",expression=1.0)
        k21  = pyurdme.Parameter(name="k21",expression=1.0)
        k22  = pyurdme.Parameter(name="k22",expression=val)
        k32  = pyurdme.Parameter(name="k32",expression=val)
        
        
        self.add_parameter([k11,k12,k13,k21,k22,k32])
    
        # Reactions
        R1 = pyurdme.Reaction(name="R1",reactants={S1:1},products={S11:1,S12:1},massaction=True, rate=k11)
        R2 = pyurdme.Reaction(name="R2",reactants={S11:1,S12:1},products={S2:1},massaction=True, rate=k12)
        R3 = pyurdme.Reaction(name="R3",reactants={S2:1},products={S21:1,S22:1},massaction=True, rate=k21)
        R4 = pyurdme.Reaction(name="R4",reactants={S21:1,S22:1},products={S3:1},massaction=True, rate=k22)
        R5 = pyurdme.Reaction(name="R5",reactants={S3:1},products={S31:1,S32:1},massaction=True, rate=k13)
        R6 = pyurdme.Reaction(name="R6",reactants={S31:1,S32:1},products={S4:1},massaction=True, rate=k32)


        self.add_reaction([R1,R2,R3,R4,R5,R6])
        
        # Scatter the molecules over the mesh
        self.set_initial_condition_scatter({S1:100})

        # Time span of the simulation
        self.timespan(numpy.arange(0,2.0,0.05))
 

if __name__ == '__main__':

    from pyurdme.rdsimsolver import RDSIMSolver
    from pyurdme.nsmsolver import NSMSolver
    import numpy
    import time

    model = Example2(voxel_size=0.1e-6)
    solver = RDSIMSolver(model)
    meso_rates = solver.meso_rates(0.1e-18)
        

    print meso_rates
    model_nsm = Example2(voxel_size=0.1e-6,val=meso_rates["k12"])
    
    solver = NSMSolver(model_nsm)
    N = 10
    time_pyurdme = []
    for i in range(N):
        tic = time.time()
        res = solver.run()
        time_pyurdme.append(time.time()-tic)

        try:
            meanS1 = meanS1 + numpy.sum(res.get_species("S1"),axis=1)
        except:
            meanS1 = numpy.sum(res.get_species("S1"),axis=1)
    print meanS1/N
    print "PyURDME: {0}".format(numpy.sum(time_pyurdme[1:-1]))

    solver = RDSIMSolver(model)
    time_rdsim = []
    for i in range(N):
        tic = time.time()
        res = solver.run()
        time_rdsim.append(time.time()-tic)
        try:
            meanS1_rdsim = meanS1_rdsim + res.get_summary_statistic("S1")
        except:
            meanS1_rdsim = res.get_summary_statistic("S1")
    print meanS1_rdsim/float(N)
    print "RDSIM: {0}".format(numpy.sum(time_rdsim[1:-1]))


# import matplotlib.pyplot as plt
#    plt.plot(model.tspan, meanS1,model.tspan,meanS1_rdsim)

#plt.show()
    #res = solver.propose_mesh_resolution_per_reaction(rel_tol=0.05)
    #solver.partition_system(rel_tol=0.05)

    # To use the meso-micro hybrid solver, simply specify the initial species partitioning.   
    #solver.set_modeling_level({"S1":"micro", "S11":"micro", "S12":"meso", "S2":"meso"})
    #microres = solver.run()
    #print microres.get_particles(0,0)
