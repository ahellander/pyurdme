#!/usr/bin/env python
""" PyURDME model with one species diffusing in the unit circle and one
    species diffusing on the boundary of the circle. Subdomains are 
    here handled by Dolfin's native subdomain model. """

import dolfin
import pyurdme
import numpy

class MembranePatch(dolfin.SubDomain):
    """ This class defines a Dolfin subdomain. Facets on lower left quadrant of 
        the boundary of the domain will be included. """
    def inside(self,x,on_boundary):
        return on_boundary and x[0] < 0.0 and x[1] < 0.0

class Membrane(dolfin.SubDomain):
    """ This class defines a Dolfin subdomain. Facets on lower left quadrant of
        the boundary of the domain will be included. """
    def inside(self,x,on_boundary):
        return on_boundary

class Cytosol(dolfin.SubDomain):
    """ This class defines a Dolfin subdomain. Facets on lower left quadrant of
        the boundary of the domain will be included. """
    def inside(self,x,on_boundary):
        return not on_boundary


class simple_diffusion2(pyurdme.URDMEModel):
    """ One species diffusing on the boundary of a sphere and one species
        diffusing inside the sphere. """
    
    def __init__(self):
        pyurdme.URDMEModel.__init__(self,name="simple_diffusion2")

        A = pyurdme.Species(name="A",diffusion_constant=0.1,dimension=2)
        B = pyurdme.Species(name="B",diffusion_constant=0.1,dimension=1)

        self.addSpecies([A,B])

        # A circle
        c1 = dolfin.Circle(0,0,1)
        mesh = dolfin.Mesh(c1,20)
        self.mesh = pyurdme.URDMEMesh(mesh)
        
        # A mesh function for the cells
        cell_function = dolfin.CellFunction("size_t",self.mesh)
        cell_function.set_all(1)
        
        # Create a mesh function over then edges of the mesh
        facet_function = dolfin.FacetFunction("size_t",self.mesh)
        facet_function.set_all(0)
        
        # Mark the boundary points
        membrane = Membrane()
        membrane.mark(facet_function,2)
        
        membrane_patch = MembranePatch()
        membrane_patch.mark(facet_function,3)
        
        self.addSubDomain(cell_function)
        self.addSubDomain(facet_function)
        
        k1 = pyurdme.Parameter(name="k1",expression=100.0)
        self.addParameter([k1])
        
        R1 = pyurdme.Reaction(name="R1",reactants={A:1},products={B:1},massaction=True,rate=k1,restrict_to=3)
        self.addReaction([R1])
        
        # Restrict species B to the membrane subdomain
        self.restrict(species=B,subdomains=[2,3])
        self.timespan(numpy.linspace(0,1,50))
        
        # Place the A molecules in the voxel nearest to the center of the square
        self.placeNear({A:10000},point=[0,0])

if __name__ == '__main__':
    
    model = simple_diffusion2()
    result = pyurdme.urdme(model)
    # Dump timeseries in Paraview format
    result.toVTK(species="B",folder_name="Bout")
    result.toVTK(species="A",folder_name="Aout")


