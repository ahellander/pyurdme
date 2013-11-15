""" pyURDME model file for the MinCDE example. """

import os
from pyurdme.pyurdme import *

class mincde(URDMEModel):

    
    def __init__():
        URDMEModel.__init__(name="mincde")

        # Species
        # TODO: We need a way to localize species to subdomains/boundaries
        MinD_m     = Species(name="MinD_m",diffusion_constant=1e-14,dimension=2)
        MinD_c_atp = Species(name="MinD_c_atp",diffusion_constant=2.5e-12)
        MinD_c_adp = Species(name="MinD_c_adp",diffusion_constant=2.5e-12)
        MinD_e     = Species(name="MinD_e",diffusion_constant=2.5e-12)
        MinDE      = Species(name="MinDE",diffusion_constant=1e-14,dimension=2)
        
        self.addSpecies([MinD_m,MinD_c_atp,MinD_c_adp,MinD_e,MinDE])

        
        initial_value = {MinD_m:1000,
                         MinD_c_adp:0,
                         MinD_c_atp:4500,
                         MinD_e:1575,
                         MinDE:0}
        
        # Parameters
        
        # Some kind of average mesh size to feed into the propensity functions
        hmax = model.mesh.hmax()
        hmin = model.mesh.hmin()
        h = (hmax+hmin)/2

        NA = Parameter(name="NA",expression="6.022e23")
        sigma_d  = Parameter(name="sigma_d",expression=2.5e-8/h)
        sigma_dD = Parameter(name="sigma_dD",expression="9.0e6/(1000.0*NA)")
        sigma_e  = Parameter(name="sigma_e",expression="5.56e7/(1000.0*NA)")
        sigma_de = Parameter(name="sigma_de",expression=0.7)
        sigma_dt = Parameter(name="sigma_dt",expression=1.0)
            
        self.addParameter([NA,sigma_d,sigma_dD,sigma_e,sigma_de,sigma_dt])

        # Load mesh in dolfins xml format
        self.mesh = read_dolfin_mesh("mesh/coli.xml")
        
        # Read the facet and interior cell physical domain markers into a Dolfin MeshFunction
        file_in = dolfin.File("mesh/coli_facet_region.xml")
        facet_function = dolfin.MeshFunction("size_t",model.mesh.mesh)
        file_in >> facet_function
        
        file_in = dolfin.File("mesh/coli_physical_region.xml")
        physical_region = dolfin.MeshFunction("size_t",model.mesh.mesh)
        file_in >> physical_region
        
        self.subdomains = [facet_function,physical_region]
        
        # List of Physical domain markers that match those in the  Gmsh .geo file.
        boundary = [73,74,75]
        interior = [76]
        
        # Reactions
        R1 = Reaction(name="R1",reactants={MinD_c_atp:1},products={MinD_m:1},massaction=True,rate=sigma_d, restrict_to=boundary)
        R2 = Reaction(name="R2",reactants={MinD_c_atp:1,MinD_m:1},products={MinD_m:2},massaction=True,rate=sigma_dD)
        R3 = Reaction(name="R3",reactants={MinD_m:1,MinD_e:1},products={MinDE:1},massaction=True,rate=sigma_e)
        R4 = Reaction(name="R4",reactants={MinDE:1,MinD_e:1},products={MinD_c_adp:1,MinD_e:1},massaction=True,rate=sigma_de)
        R5 = Reaction(name="R5",reactants={MinD_c_adp:1},products={MinD_c_atp:1},massaction=True,rate=sigma_dt)
        R6 = Reaction(name="R6",reactants={MinDE:1,MinD_c_atp:1},products={MinD_m:1,MinDE:1},massaction=True,rate=sigma_dD,restrict_to=boundary)
        
        self.addReaction([R1,R2,R3,R4,R5,R6])
        
        # Distribute molecules randomly over the mesh according to their initial values
        self.scatter(MinD_m)
        self.scatter(MinD_c_adp)
        self.scatter(MinD_e)

        self.timespan(range(100))

if __name__=="__main__":
    """ Dump model to a file. """
                     
    model = mincde(model_name="mincde")
    model.initialize()
   
    result = urdme(model)
    
    file = dolfin.File("mindm.pvd")
    file << model.sol['MinD_m']
    toXYZ(model,'mindm.xyz',format="VMD")
    #file = dolfin.File("minde.pvd")
    #file << model.sol['MinD_e']

    
    #file = dolfin.File("mincatp.pvd")
    #file << model.sol['MinD_c_atp']

    
    print result



#model.mesh.init()
#
#edge2vertex = model.mesh.mesh.topology()(1,0)
#facet2vertex   = model.mesh.mesh.topology()(2,0)
#cell2vertex  = model.mesh.mesh.topology()(3,0)

#sd = numpy.zeros((1,model.mesh.getNumVoxels()))
#for i in range(physical_region.size()):
#    for vtx in cell2vertex(i):
#        sd[0,vtx] = physical_region[i]

#for i in range(facet_function.size()):
#    for vtx in facet2vertex(i):
#        if facet_function[i] != 0:
#            sd[0,vtx] = facet_function[i]

#model.sd = sd.flatten()

#dx = dolfin.Measure("dx")[physical_region]
#ds = dolfin.Measure("ds")[facet_function]
