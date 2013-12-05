""" pyURDME model file for the MinCDE example. """

import os
from pyurdme.pyurdme import *
import dolfin


class Membrane(dolfin.SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary


class Cytosol(dolfin.SubDomain):
    def inside(self,x,on_boundary):
        return not on_boundary


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

        # Load mesh in dolfins xml format

        # This (below) hangs for unknown reason...
        # Build CSG Coli using Dolfin.
        #sphere1 = dolfin.Sphere(dolfin.Point(0,0,2.25),0.5)
        #sphere2 = dolfin.Sphere(dolfin.Point(0,0,-2.25),0.5)
        #cylinder = dolfin.Cylinder(dolfin.Point(0,0,2.25),dolfin.Point(0,0,-2.25),0.5)
        #dolfin.plot(cylinder+sphere1+sphere2)
        #geom = cylinder+sphere1+sphere2
        #self.mesh = dolfin.Mesh(geom,1)
        #dolfin.info(self.mesh)
        #dolfin.plot(self.mesh)
        #dolfin.interactive()
        
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
        
        #subdomains = dolfin.MeshFunction("size_t",self.mesh,1)

        physical_region = dolfin.MeshFunction("size_t",self.mesh,self.mesh.topology().dim())
        physical_region.set_all(1)
        
        facet_function = dolfin.MeshFunction("size_t",self.mesh,self.mesh.topology().dim()-1)
        facet_function.set_all(0)
        # Mark the boundary points
        membrane = Membrane()
        # interior = Cytosol()
        membrane.mark(facet_function,74)
        #interior.mark(subdomains,5)
        #boundary = [2]
        
        #self.subdomains = [subdomains]
        self.subdomains = [physical_region, facet_function]
        # Average mesh size to feed into the propensity functions
        hmax = self.mesh.hmax()
        hmin = self.mesh.hmin()
        h = (hmax+hmin)/2

        # Parameters
        NA = Parameter(name="NA",expression="6.022e23")
        #sigma_d  = Parameter(name="sigma_d",expression=2.5e-8/hmin)
        #sigma_dD = Parameter(name="sigma_dD",expression="9.0e6/(1000.0*NA)")
        #sigma_e  = Parameter(name="sigma_e",expression="5.56e7/(1000.0*NA)")
        #sigma_de = Parameter(name="sigma_de",expression=0.7)
        #sigma_dt = Parameter(name="sigma_dt",expression=1.0)
        
        sigma_d  = Parameter(name="sigma_d",expression=0.0)
        sigma_dD = Parameter(name="sigma_dD",expression=0.0)
        sigma_e  = Parameter(name="sigma_e",expression=0.0)
        sigma_de = Parameter(name="sigma_de",expression=0.0)
        sigma_dt = Parameter(name="sigma_dt",expression=0.0)
        
        self.addParameter([NA,sigma_d,sigma_dD,sigma_e,sigma_de,sigma_dt])

        # List of Physical domain markers that match those in the  Gmsh .geo file.
        boundary = [74]
        interior = [1]
        
        # Reactions
        R1 = Reaction(name="R1",reactants={MinD_c_atp:1},products={MinD_m:1},massaction=True,rate=sigma_d, restrict_to=boundary)
        R2 = Reaction(name="R2",reactants={MinD_c_atp:1,MinD_m:1},products={MinD_m:2},massaction=True,rate=sigma_dD,restrict_to=boundary)
        R3 = Reaction(name="R3",reactants={MinD_m:1,MinD_e:1},products={MinDE:1},massaction=True,rate=sigma_e,restrict_to=boundary)
        R4 = Reaction(name="R4",reactants={MinDE:1},products={MinD_c_adp:1,MinD_e:1},massaction=True,rate=sigma_de,restrict_to=boundary)
        R5 = Reaction(name="R5",reactants={MinD_c_adp:1},products={MinD_c_atp:1},massaction=True,rate=sigma_dt)
        R6 = Reaction(name="R6",reactants={MinDE:1,MinD_c_atp:1},products={MinD_m:1,MinDE:1},massaction=True,rate=sigma_dD,restrict_to=boundary)
        
        self.addReaction([R1,R2,R3,R4,R5,R6])
        
        # Restrict to boundary
        self.restrict(MinD_m,boundary)
        self.restrict(MinDE,boundary)
        
        # Distribute molecules randomly over the mesh according to their initial values
        self.scatter({MinD_m:1000},subdomains=boundary)
        #self.scatter({MinD_c_adp:4500})
        #self.scatter({MinD_e:1575},subdomains=boundary)

        self.timespan(range(100))

if __name__=="__main__":
    """ Dump model to a file. """
                     
    model = mincde(model_name="mincde")
    model.serialize("debug_input.mat")
    result = urdme(model)
<<<<<<< HEAD
    
    file = dolfin.File("mindm.pvd")
    file << model.sol['MinD_m']
    toXYZ(model,'mindm.xyz',file_format="VMD")
    #file = dolfin.File("minde.pvd")
    #file << model.sol['MinD_e']
=======
    U = result["U"]
>>>>>>> mincdesubdomains

    print numpy.sum(U[::5,:],axis=0)
    dumps(model,"MinD_m","mindout")
    
    sd = model.subdomainVector(model.subdomains)
    func = dolfin.Function(model.xmesh.function_space["MinD_m"])
    dofmap = model.xmesh.vertex_to_dof_map["MinD_m"]
    
    numvox = model.mesh.getNumVoxels()
    #func = dolfin.MeshFunction("uint",model.mesh,1)
    func_vector = func.vector()
    #func_vector = func.array()

    # This works!!
    for i in range(numvox):
        dof = dofmap[i]
        vox = dof/5
        func_vector[vox] = sd[i]
    file = dolfin.File("sd.pvd")
    file << func

    toXYZ(model,'mindm.xyz',format="VMD")
    
#print result
