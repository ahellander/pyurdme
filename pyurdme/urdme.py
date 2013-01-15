from model import *
import numpy as np
import scipy.sparse as scisp
import scipy.io as spio
#import dolfin

# A URDME Model is an extenstion of Model
class URDMEModel(Model):

    """ 
        An URDME Model extends Model with spatial information and methods 
        to create URDME solver input data structures.
        
    """

    def createStoichiometricMatrix(self):
        """ Create a sparse stoichiometric matrix
            from the model's listOfReactions. """
        N = numpy.zeros((self.num_species,self.num_reactions))
        for i in range(len(self.listOfReactions)):
            R = self.listOfReactions[i]
            N[i,]

    def createDependencyGraph(self):
        """ Dependecy graph. """
        #TODO: Create a better depencey graph
        GF = np.ones((self.num_reactions,self.num_reactions+self.num_species))
        self.G=scisp.csc_matrix(GF)
    
    def assemble(self):
        """ Assemble the diffusion matrix. """
    
    def meshextend(self):
        """ Extend the raw mesh datastructure to include
            information about degrees of freedom. Initialize
            URDME datastructures that depend on the size 
            of the mesh. """
    
        # Initialize the Dolfin mesh object.
        # model.mesh.init()
        # self.num_voxels  = self.mesh.num_vertices()
        # self.num_species = len(self.listOfSpecies)
        self.u0 = np.zeros((self.num_species,self.num_voxels))
        # model.dofs = [...]
    
    def scatter(self,species,subdomain=None):
        """ Scatter an initial number of molecules over the 
            voxels in a subdomain. """
    
        Sname = species.name
        numS = species.initial_value
    
        for i in range(numS):
            vtx=np.random.randint(0,self.num_voxels)
            self.u0[0,vtx]+=1
    
        # Is the initial condition matrix constructed?
        # Check for exisitence of u0 here.
        
    
    def validate(self):
        """ Validate the model data structures. """
    
    def serialize(self,filename=[]):
        """ 
            Serialize the model object to binary file compatible
            with core URDME solvers. 
            
        """
        # Stoichimetric matrix
        # N = createStochiometricMatrix(self)

        # Dependency Graph
        G = self.createDependencyGraph()


        # Diffusion matrix
        # D =

        # Initial condition
        #u0 =


        # Volume vector
        # vol =

        # Data vector
        # data = []
        filename = filename
        spio.savemat(filename,{'num_species':self.num_species,'u0':self.u0,'G':self.G})


class URDMEMesh():
    """ Mesh object. A thin wrapper around the Fenics mesh object. """

    def __init__():
        mesh = dolfin.Mesh()


class URDMEXmesh():
    """ Extended mesh object. """
    

def urdme():
    """ URDME solver interface, similar to the Matlab function interface. """

    # Shell out and compile the solver  

