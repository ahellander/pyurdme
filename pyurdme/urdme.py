from model import *
import numpy as np
#import dolfin

# A URDME Model is an extenstion of Model
class URDMEModel(Model):

    """ 
        An URDME Model will extend Model with spatial information/Users/andreash/github/ahellander/pyurdme/examples/dimerization/surface.xml
        and methods to create URDME solver input data structures.
        
    """

    def createStoichiometricMatrix(self):
        """ Create a sparse stoichiometric matrix
            from the model's listOfReactions. """
    
    def assemble(self):
        """ Assemble the diffusion matrix. """
    
    def meshextend(self):
        """ Extend the raw mesh datastructure to include
            information about degrees of freedom. """
    
        # Initialize the Dolfin mesh object.
        # model.mesh.init()
        # model.dofs = [...]
    
    def validate(self):
        """ Validate the model data structures. """
    
    def serialize(self):
        """ 
            Serialize the model object to binary file compatible
            with core URDME solvers. 
            
        """
        # Stoichimetric matrix
        # N = createStochiometricMatrix(self)

        # Dependency Graph
        # G =


        # Diffusion matrix
        # D =

        # Initial condition
        #u0 =


        # Volume vector
        # vol =

        # Data vector
        # data = []




class URDMEMesh():
    """ Mesh object. A thin wrapper around the Fenics mesh object. """

        def __init__():
            mesh = dolfin.Mesh()

