#!/usr/bin/env python
import pyurdme
from pyurdme.nsmsolver import NSMSolver
import numpy

class CoralReef(pyurdme.URDMEModel):
    """ Model developed by Briggs and Drawert 3/31/2014, based on a 
        non-spatial model by Briggs and Adam.
    """

    def __init__(self, name="coral_reef"):
        pyurdme.URDMEModel.__init__(self, name)

        # Species
        Coral = pyurdme.Species(name="Coral",diffusion_constant=0.0)
        Coral_m = pyurdme.Species(name="Coral_m",diffusion_constant=1.0)
        MA = pyurdme.Species(name="MA", diffusion_constant=0.0)
        MA_m = pyurdme.Species(name="MA_m", diffusion_constant=1.0)
        Turf = pyurdme.Species(name="Turf", diffusion_constant=0.0)
        self.addSpecies([Coral, MA, Coral_m, MA_m, Turf])

        # Parameters
        phi_c  = pyurdme.Parameter(name="phi_c", expression=0.001) #1/year
        phi_m  = pyurdme.Parameter(name="phi_m", expression=0.0001) #1/year
        g_tc  = pyurdme.Parameter(name="g_tc", expression=0.1) #1/year
        g_tm = pyurdme.Parameter(name="g_tm", expression=0.2) #1/year
        Gamma  = pyurdme.Parameter(name="Gamma", expression=0.05)
        dc = pyurdme.Parameter(name="dc", expression=0.05) #1/year
        dm = pyurdme.Parameter(name="dm", expression=0.2) #1/year
        psi_g = pyurdme.Parameter(name="psi_g", expression=0.0)
        #
        mu_c = pyurdme.Parameter(name="mu_c", expression=1.0) #1/year
        mu_m = pyurdme.Parameter(name="mu_m", expression=1.0) #1/year
        self.addParameter([phi_c, phi_m, g_tc , g_tm , Gamma, dc, dm, psi_g, mu_c , mu_m])

        # Reactions:

        # Propagule Recruitment
        # T -> C  : phi_c
        R0 = pyurdme.Reaction(name="R0", reactants={Turf:1}, products={Coral:1}, rate=phi_c)
        # T -> MA : phi_m
        R1 = pyurdme.Reaction(name="R1", reactants={Turf:1}, products={MA:1}, rate=phi_m)
        
        # Death rates
        # C -> T : dc
        R3 = pyurdme.Reaction(name="R3", reactants={Coral:1}, products={Turf:1}, rate=dc)
        # MA -> T : dm
        R4 = pyurdme.Reaction(name="R4", reactants={MA:1}, products={Turf:1}, rate=dm)

        # Growth
        # T + C_m -> C : g_tc * exp(-1*psi_g*MA*100)
        R5 = pyurdme.Reaction(name="R5", reactants={Turf:1, Coral_m:1}, products={Coral:1}, propensity_function="Turf * Coral_m * g_tc * exp(-1.0 * psi_g * MA * 100)/vol")
        # T + MA_m -> MA : g_tm
        R6  = pyurdme.Reaction(name="R6", reactants={Turf:1, MA_m:1}, products={MA:1}, propensity_function="Turf * MA_m * g_tm / vol")
        # MA + C -> 2MA : Gamma*g_tm
        R7 = pyurdme.Reaction(name="R7", reactants={Coral:1}, products={MA:1}, propensity_function="MA * Coral * Gamma * g_tm / vol")

        # Movement, working in meters
        # C -> C + C_m : mu_c + g_tc * exp(-1*psi_g*MA*100)
        R8 = pyurdme.Reaction(name="R8", reactants={}, products={Coral_m:1}, propensity_function="Coral * (mu_c + g_tc * exp(-1.0 * psi_g * MA * 100))")
        # C_m -> 0 :  mu_c
        R9 = pyurdme.Reaction(name="R9", reactants={Coral_m:1}, products={}, rate=mu_c)
        # MA -> MA + MA_m : g_tc + mu_m 
        R10 = pyurdme.Reaction(name="R10", reactants={}, products={MA_m:1}, propensity_function="MA * (g_tc + mu_m)")
        # MA_m -> 0 :  mu_m
        R11 = pyurdme.Reaction(name="R11", reactants={MA_m:1},  products={}, rate=mu_m)

        self.addReaction([R0, R1, R3, R4, R5, R6, R7, R8, R9, R10, R11])


        # A unit square
        # each grid point is 10cm x 10cm, domain is 5m x 5m
        self.mesh = pyurdme.URDMEMesh.SquareMesh(L=5, nx=50, ny=50, periodic=True)

        # Place the A molecules in the voxel nearest the center of the square
        #self.placeNear({Coral:100}, point=[10,10])
        
        self.distributeUniformly({Turf:100})
        self.placeNear({Coral:100}, point=[1,1])
        self.placeNear({Turf:0}, point=[1,1])
        self.placeNear({MA:100}, point=[4,4])
        self.placeNear({Turf:0}, point=[4,4])


        for vndx in range(self.u0.shape[1]):
            tot = 0
            for sndx, sname in enumerate(self.listOfSpecies):
                tot += self.u0[sndx][vndx]
            if tot > 100:
                for sndx, sname in enumerate(self.listOfSpecies):
                    print "u0[{0}][{1}] = {2}".format(sname, vndx, self.u0[sndx][vndx])

        #self.timespan(numpy.linspace(0,500,501)) #500 years
        #self.timespan(numpy.linspace(0,5,72)) #5 years, by months
        self.timespan(numpy.linspace(0,11,66)) #10 years, by 2 months

if __name__ == "__main__":
    reef = CoralReef()
    sol = NSMSolver(reef, report_level=1)
    result = sol.run()
    result.toVTK(species='Coral',folder_name="output_coral")
