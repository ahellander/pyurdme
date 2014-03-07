from pyurdme import pyurdme
import numpy

class CoralReef(pyurdme.URDMEModel):
    """ Model of a coral reef, published in Sandin & McNamara (2012) """

    def __init__(self, name="coral_reef"):
        pyurdme.URDMEModel.__init__(self, name)

        Coral = pyurdme.Species(name="Coral",diffusion_constant=0.0,dimension=2)
        self.addSpecies([Coral])

        P_cr  = pyurdme.Parameter(name="P_cr", expression=0.01)
        G_c   = pyurdme.Parameter(name="G_c", expression=0.01)
        P_cd  = pyurdme.Parameter(name="P_cd", expression=0.15)
        C_th  = pyurdme.Parameter(name="C_th", expression=900)
        C_max = pyurdme.Parameter(name="C_max", expression=7000)
        

        # A unit square
        # each grid point is 10cm x 10cm, domain is 5m x 5m
        self.mesh = pyurdme.Mesh.unitSquareMesh(50,50)

        # Place the A molecules in the voxel nearest the center of the square
        self.placeNear({Coral:1}, point=[0.5,0.5])

        self.timespan(numpy.linspace(0,1,50))        
