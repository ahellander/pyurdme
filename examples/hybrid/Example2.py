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
    
    def __init__(self,voxel_size=0.1e-6, vala=1e-21, vald=1.0,  sigma=1e-9,nmol=100):
        
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
        species = [S1,S2,S3,S4,S11,S12,S21,S22,S31,S32]

        self.add_species(species)

        pi = math.pi
        self.voxel_size = voxel_size
        L = 1e-6
        h = voxel_size
        nx = int(L/h)

        N = nx*nx*nx
        self.mesh = pyurdme.URDMEMesh.generate_cube_mesh(L,nx,nx,nx)
       
        # Microscopic association and disassociation rate
        k11  = pyurdme.Parameter(name="k11",expression=vald)
        k12  = pyurdme.Parameter(name="k12",expression=vala)
        k13  = pyurdme.Parameter(name="k13",expression=vald)
        k21  = pyurdme.Parameter(name="k21",expression=vald)
        k22  = pyurdme.Parameter(name="k22",expression=vala)
        k32  = pyurdme.Parameter(name="k32",expression=vala)
        
        
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
        for S in species:
            self.set_initial_condition_scatter({S:nmol})

        # Time span of the simulation
        self.timespan(numpy.arange(0,2.0,0.05))
 

if __name__ == '__main__':

    from pyurdme.rdsimsolver import RDSIMSolver
    from pyurdme.nsmsolver import NSMSolver
    import numpy
    import time

    model = Example2(voxel_size=0.1e-6)
    solver = RDSIMSolver(model, solver_type="RDME")
    
    # res = solver.run()
    #res.get_summary_statistic(["S1"])
    solver._write_mesh_file("Example2_mesh.h5")
    solver.create_input_file("Example2.json")
    solver.serialize("Example2_model.mat")
    
    meso_rates = solver.meso_rates(0.1e-18)
    
    #print meso_rates
    val = meso_rates["k12"]
    oldq= model.listOfParameters["k11"].value/model.listOfParameters["k12"].value
    vala=meso_rates["k12"]
    vald=oldq * vala
    # print vald
    #model_nsm = Example2(voxel_size=0.1e-6,vala=vala, vald=vald)
    
    
    voxel_sizes=[0.4e-6, 0.2e-6, 0.1e-6, 0.05e-6]
    # solver = NSMSolver(model_nsm)
    N = 10
    """time_nsm = []
    for vxs in voxel_sizes:
        mst = 0.0
        model_nsm = Example2(voxel_size=vxs,vala=vala, vald=vald)
        solver = NSMSolver(model_nsm)
        for i in range(N):
            tic = time.time()
            res = solver.run()
            mst += (time.time()-tic)
        time_nsm.append(mst/N)
        print time_nsm


    print "Time NSM: ", time_nsm
    """


    time_rdsim = []
    for vxs in voxel_sizes:
        mst = 0.0
        model = Example2(voxel_size=vxs,vala=vala, vald=vald)
        solver = RDSIMSolver(model)
        for i in range(N):
            tic = time.time()
            res = solver.run()
            res.get_summary_statistic(["S1"])
            mst += (time.time()-tic)
        time_rdsim.append(mst/N)
        print time_rdsim


    print "Time RDSIM: ", time_rdsim
#try:
#           meanS1 = meanS1 + numpy.sum(res.get_species("S1"),axis=1)
#       except:
#           meanS1 = numpy.sum(res.get_species("S1"),axis=1)
#print meanS1/N
#    print "PyURDME: {0}".format(numpy.sum(time_pyurdme[1:-1]))

#  solver = RDSIMSolver(model)
#    time_rdsim = []
#    for i in range(N):
#        tic = time.time()
#        res = solver.run()
#        time_rdsim.append(time.time()-tic)
#        try:
#            meanS1_rdsim = meanS1_rdsim + res.get_summary_statistic("S1")
#        except:
#            meanS1_rdsim = res.get_summary_statistic("S1")
#    print meanS1_rdsim/float(N)
#    print "RDSIM: {0}".format(numpy.sum(time_rdsim[1:-1]))


# import matplotlib.pyplot as plt
#    plt.plot(model.tspan, meanS1,model.tspan,meanS1_rdsim)

#plt.show()
    #res = solver.propose_mesh_resolution_per_reaction(rel_tol=0.05)
    #solver.partition_system(rel_tol=0.05)

    # To use the meso-micro hybrid solver, simply specify the initial species partitioning.   
    #solver.set_modeling_level({"S1":"micro", "S11":"micro", "S12":"meso", "S2":"meso"})
    #microres = solver.run()
    #print microres.get_particles(0,0)
