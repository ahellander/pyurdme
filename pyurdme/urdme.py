from model import *
import numpy as np

# A URDME Model is an extenstion of Model
class URDMEModel(Model):

    """ 
        An URDME Model will extend Model with spatial information
        and methods to create URDME solver input data structures.
        
    """

    def stoichiometric_matrix(reactions,species):
        """ Create a sparse stoichiometric matrix on the 
            format accepeted by the core URDME solvers. """
        
    
    def serialize(self):
        """ 
            This will be very different from StochKit serialization.
            URDME needs several datastructures (sparse matrices) to
            be created, and written as a binary .mat file. 
            
        """
        N = stoichiometric_matrix(self.listOfReactions,self.listOfSpecies)
