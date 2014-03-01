""" pyURDME model file for the MinCDE example. """

import os.path
from pyurdme import pyurdme
import dolfin
import numpy

class Membrane(dolfin.SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary


class Cytosol(dolfin.SubDomain):
    def inside(self,x,on_boundary):
        return not on_boundary


class mincde(pyurdme.URDMEModel):

    def __init__(self,model_name="mincde"):
        pyurdme.URDMEModel.__init__(self,model_name)

        # Species
        # TODO: We need a way to localize species to subdomains/boundaries
        MinD_m     = pyurdme.Species(name="MinD_m",diffusion_constant=1e-14,dimension=2)
        MinD_c_atp = pyurdme.Species(name="MinD_c_atp",diffusion_constant=2.5e-12,dimension=3)
        MinD_c_adp = pyurdme.Species(name="MinD_c_adp",diffusion_constant=2.5e-12,dimension=3)
        MinD_e     = pyurdme.Species(name="MinD_e",diffusion_constant=2.5e-12,dimension=3)
        MinDE      = pyurdme.Species(name="MinDE",diffusion_constant=1e-14,dimension=2)
        
        self.addSpecies([MinD_m,MinD_c_atp,MinD_c_adp,MinD_e,MinDE])
        
        # Make sure that we have the correct path to the mesh file even if we are not executing from the basedir.
        basedir = os.path.dirname(os.path.abspath(__file__))
        self.mesh = pyurdme.Mesh.read_dolfin_mesh(basedir+"/mesh/coli.xml")
        
        interior = dolfin.CellFunction("size_t",self.mesh)
        interior.set_all(1)
        boundary = dolfin.FacetFunction("size_t",self.mesh)
        boundary.set_all(0)
        
        # Mark the boundary points
        membrane = Membrane()
        membrane.mark(boundary,2)
        
        self.addSubDomain(interior)
        self.addSubDomain(boundary)
        
        # Average mesh size to feed into the propensity functions
        hmax = self.mesh.hmax()
        hmin = self.mesh.hmin()
        h = (hmax+hmin)/2

        # Parameters
        NA = pyurdme.Parameter(name="NA",expression="6.022e23")
        sigma_d  = pyurdme.Parameter(name="sigma_d",expression=2.5e-8/hmin)
        sigma_dD = pyurdme.Parameter(name="sigma_dD",expression="9.0e5/(1000.0*NA)")
        sigma_e  = pyurdme.Parameter(name="sigma_e",expression="5.56e7/(1000.0*NA)")
        sigma_de = pyurdme.Parameter(name="sigma_de",expression=0.7)
        sigma_dt = pyurdme.Parameter(name="sigma_dt",expression=1.0)
        
        self.addParameter([NA,sigma_d,sigma_dD,sigma_e,sigma_de,sigma_dt])

        # List of Physical domain markers that match those in the  Gmsh .geo file.
        interior = [1]
        boundary = [2]
        
        # Reactions
        R1 = pyurdme.Reaction(name="R1",reactants={MinD_c_atp:1},products={MinD_m:1},massaction=True,rate=sigma_d, restrict_to=boundary)
        R2 = pyurdme.Reaction(name="R2",reactants={MinD_c_atp:1,MinD_m:1},products={MinD_m:2},massaction=True,rate=sigma_dD)
        R3 = pyurdme.Reaction(name="R3",reactants={MinD_m:1,MinD_e:1},products={MinDE:1},massaction=True,rate=sigma_e)
        R4 = pyurdme.Reaction(name="R4",reactants={MinDE:1},products={MinD_c_adp:1,MinD_e:1},massaction=True,rate=sigma_de)
        R5 = pyurdme.Reaction(name="R5",reactants={MinD_c_adp:1},products={MinD_c_atp:1},massaction=True,rate=sigma_dt)
        R6 = pyurdme.Reaction(name="R6",reactants={MinDE:1,MinD_c_atp:1},products={MinD_m:1,MinDE:1},massaction=True,rate=sigma_dD)
        
        self.addReaction([R1,R2,R3,R4,R5,R6])
        
        # Restrict to boundary
        self.restrict(MinD_m,boundary)
        self.restrict(MinDE,boundary)
        
        # Distribute molecules over the mesh according to their initial values
        self.scatter({MinD_c_adp:4500})
        self.scatter({MinD_e:1575})

        self.timespan(range(50))

if __name__=="__main__":
    """ Dump model to a file. """
                     
    model = mincde(model_name="mincde")
    model.timespan(range(50))
    result = pyurdme.urdme(model)
    jsondoc = result.toTHREEJs("MinD_m",49)
    with open("testmesh.json",'w') as file:
        file.write(jsondoc)

    particlestr = result.printParticlejs(["MinD_m", "MinDE"],49)
    with open("particletest.html",'w') as file:
        file.write(particlestr)



