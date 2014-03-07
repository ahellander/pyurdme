""" pyurdme model for a single species diffusing on the unit square. """

import dolfin
import pyurdme
import numpy

class simple_diffusion(pyurdme.URDMEModel):
    """ Initial condition is a delta function at the center voxel. 
        The solution should be a Gaussian, up to the point where
        the BC becomes important. """

    def __init__(self):
    
        pyurdme.URDMEModel.__init__(self,name="simple_diffusion")

        D = 0.01
        A = pyurdme.Species(name="A",diffusion_constant=D,dimension=2)

        self.addSpecies([A])

        # A unit square
        self.mesh = pyurdme.Mesh.unitSquareMesh(40,40)
                
        # Place the A molecules in the voxel nearest the center of the square
        self.placeNear({A:100000},point=[0.5,0.5])

        self.timespan(numpy.linspace(0,1,50))


if __name__ == '__main__':

    model = simple_diffusion()
    result = pyurdme.urdme(model)

    # Dump a snapshot of the state in paraview format. To visualize the solution,
    # open output/trajectory.pvd in ParaView.
    result.toVTK("A", "output")
