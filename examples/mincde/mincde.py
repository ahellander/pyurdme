""" pyURDME model file for the MinCDE example. """

import os
from pyurdme.pyurdme import *
import dolfin


class Membrane(dolfin.SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary

class mincde(URDMEModel):

    def __init__(self,model_name):
        URDMEModel.__init__(self,name="mincde")

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
        
        # Load mesh in dolfins xml format
        self.mesh = read_dolfin_mesh("mesh/coli.xml")
        
        # Read the facet and interior cell physical domain markers into a Dolfin MeshFunction
        
        # TODO:  There is an issue here in that FeniCS dolfin-convert writes the value 0 for all the faces that
        # are not on the boundary. I think we migth have to write our own Gmsh2Dolfin converter.
        #file_in = dolfin.File("mesh/coli_facet_region.xml")
        #facet_function = dolfin.MeshFunction("size_t",self.mesh)
        #file_in >> facet_function
        
        #file_in = dolfin.File("mesh/coli_physical_region.xml")
        #physical_region = dolfin.MeshFunction("size_t",self.mesh)
        #file_in >> physical_region
        
        subdomains = dolfin.MeshFunction("size_t",self.mesh,self.mesh.topology().dim()-1)
        subdomains.set_all(1)
        
        # Mark the boundary points
        membrane = Membrane()
        membrane.mark(subdomains,2)
        boundary = [2]

        self.subdomains = [subdomains]
        # Average mesh size to feed into the propensity functions
        hmax = self.mesh.hmax()
        hmin = self.mesh.hmin()
        h = (hmax+hmin)/2

        # Parameters
        NA = Parameter(name="NA",expression="6.022e23")
        sigma_d  = Parameter(name="sigma_d",expression=2.5e-8/h)
        sigma_dD = Parameter(name="sigma_dD",expression="9.0e6/(1000.0*NA)")
        sigma_e  = Parameter(name="sigma_e",expression="5.56e7/(1000.0*NA)")
        sigma_de = Parameter(name="sigma_de",expression=0.7)
        sigma_dt = Parameter(name="sigma_dt",expression=1.0)
            
        self.addParameter([NA,sigma_d,sigma_dD,sigma_e,sigma_de,sigma_dt])

        # List of Physical domain markers that match those in the  Gmsh .geo file.
        #boundary = [73,74,75,79]
        #interior = [76]
        
        # Reactions
        R1 = Reaction(name="R1",reactants={MinD_c_atp:1},products={MinD_m:1},massaction=True,rate=sigma_d, restrict_to=boundary)
        R2 = Reaction(name="R2",reactants={MinD_c_atp:1,MinD_m:1},products={MinD_m:2},massaction=True,rate=sigma_dD)
        R3 = Reaction(name="R3",reactants={MinD_m:1,MinD_e:1},products={MinDE:1},massaction=True,rate=sigma_e)
        R4 = Reaction(name="R4",reactants={MinDE:1,MinD_e:1},products={MinD_c_adp:1,MinD_e:1},massaction=True,rate=sigma_de)
        R5 = Reaction(name="R5",reactants={MinD_c_adp:1},products={MinD_c_atp:1},massaction=True,rate=sigma_dt)
        R6 = Reaction(name="R6",reactants={MinDE:1,MinD_c_atp:1},products={MinD_m:1,MinDE:1},massaction=True,rate=sigma_dD,restrict_to=boundary)
        
        self.addReaction([R1,R2,R3,R4,R5,R6])
        
        # Restrict to boundary
        self.restrict(MinD_m,boundary)
        self.restrict(MinDE,boundary)
        
        # Distribute molecules randomly over the mesh according to their initial values
        self.scatter({MinD_m:1000},subdomains=boundary)
        self.scatter({MinD_c_atp:4500})
        self.scatter({MinD_e:1575},subdomains=boundary)

        self.timespan(range(100))

if __name__=="__main__":
    """ Dump model to a file. """
                     
    model = mincde(model_name="mincde")
    result = urdme(model)
    
    file = dolfin.File("mindm.pvd")
    # Dump the solution at time 99 (end time) in pvd format
    file << model.sol['MinD_m'][99]
    toXYZ(model,'mindm.xyz',format="VMD")
    
    print result
