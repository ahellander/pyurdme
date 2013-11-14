""" PyURDME model with one species diffusing in the unit circle and one
    species diffusing on the boundary of the circle. Subdomains are 
    here handled by Dolfin's native subdomain model. """

import dolfin
from pyurdme.pyurdme import *


class Membrane(dolfin.SubDomain):
    """ This class defines a Dolfin subdomain. Facets on lower left quadrant of 
        the boundary of the domain will be included. """
    def inside(self,x,on_boundary):
        return on_boundary and x[0] < 0.0 and x[1] < 0.0


class simple_diffusion2(URDMEModel):
    """ One species diffusing on the boundary of a sphere and one species
        diffusing inside the sphere. """
    
    def __init__(self):
        URDMEModel.__init__(self,name="simple_diffusion2")

        D = 0.01
        A = Species(name="A",diffusion_constant=D,dimension=2)
        B = Species(name="B",diffusion_constant=0.1*D,dimension=1)

        self.addSpecies([A,B])

        # A circle
        c1 = dolfin.Circle(0,0,1)
        mesh = dolfin.Mesh(c1,40)
        self.mesh = Mesh(mesh)
        
        # Create a mesh function over vertices of the mesh
        sub_domains = dolfin.MeshFunction("size_t",self.mesh,self.mesh.topology().dim()-1)
        sub_domains.set_all(1)
        # Mark the boundary points
        mem = Membrane()
        mem.mark(sub_domains,2)
        
        self.subdomains = [sub_domains]
        
        file = dolfin.File("subdomains.pvd")
        file << sub_domains
    
        
        # Restrict species A to the membrane subdomain
        #self.restrict(species=A,subdomain=[2])
        self.timespan(numpy.linspace(0,1,50))
        
        # Place the A molecules in the voxel nearest to the center of the square
        self.placeNear({A:10000},point=[0,0])
        self.scatter({B:1000},subdomain=[2])
        

if __name__ == '__main__':
    
    model = simple_diffusion2()
    result = urdme(model)
    print model.sd

    # Dump a snapshot of the state in paraview format
    file1 = dolfin.File("A.pvd")
    file1 << model.sol['A']
    file2 = dolfin.File("B.pvd")
    file2 << model.sol['B']

