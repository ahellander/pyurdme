from model import *
import numpy as np
import scipy.sparse as scisp
import scipy.io as spio

import os

# A URDME Model is an extenstion of Model
class URDMEModel(Model):

    """ 
        An URDME Model extends Model with spatial information and methods 
        to create URDME solver input data structures.
        
    """

    def createStoichiometricMatrix(self):
        """ Create a sparse stoichiometric matrix
            from the model's listOfReactions. """
        ND = np.zeros((self.num_species,self.num_reactions))
        i=0
        for r in self.listOfReactions:
            R = self.listOfReactions[r]
            reactants = R.reactants
            products = R.products
            
            for s in reactants:
                ND[self.species_map[s],i]=-reactants[s]
            for s in products:
                ND[self.species_map[s],i]=products[s]
            i = i+1
    
        self.N = scisp.csc_matrix(ND)


    def createDependencyGraph(self):
        """ Dependecy graph. """
        #TODO: Create a better depencey graph
        GF = np.ones((self.num_reactions,self.num_reactions+self.num_species))
        self.G=scisp.csc_matrix(GF)
    
    def createPropensityFile(self):
        """ Generate a C propensity file for use with the 
            URDME Core C solvers. """
        
    
    def assemble(self):
        """ Assemble the diffusion matrix and volume vector. """
        #M =
    
    
    def meshextend(self):
        """ Extend the raw mesh datastructure to include
            information about degrees of freedom. Initialize
            URDME datastructures that depend on the size 
            of the mesh. """
        # Construct a species map (dict mapping species name to integer index)
        i=0
        self.species_map = {}
        for S in self.listOfSpecies:
            self.species_map[S]=i
            i = i+1;
    
        # Initialize the Dolfin mesh object.
        # model.mesh.init()
        # self.num_voxels  = self.mesh.num_vertices()
        # self.num_species = len(self.listOfSpecies)
        #xmesh = URDMEXmesh()
        #xmesh.dofs =
    
    
    def InitInitialValue(self):
        """ Create all-zeros inital condition matrix. """    
        self.u0 = np.zeros((self.num_species,self.num_voxels))
    
    
    def scatter(self,species,subdomain=None):
        """ Scatter an initial number of molecules over the 
            voxels in a subdomain. """
    
        Sname = species.name
        numS = species.initial_value
        specindx= self.species_map[Sname]
        
        # TODO: USE THE SUBDOMAIN INFO
        for i in range(numS):
            vtx=np.random.randint(0,self.num_voxels)
            self.u0[specindx,vtx]+=1
    
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
        N = self.createStoichiometricMatrix()

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
        spio.savemat(filename,{'num_species':self.num_species,'N':self.N,'u0':self.u0,'G':self.G},oned_as='column')


class URDMEMesh():
    """ Mesh object. A thin wrapper around the Fenics mesh object. """

    def __init__():
        mesh = dolfin.Mesh()


class URDMEXmesh():
    """ Extended mesh object. """
    

def urdme(model=None,solver='nsm'):
    """ URDME solver interface, similar to the Matlab function interface. """

    # Shell out and compile the solver  
    #URDME_ROOT = os.popen('urdme_init -r')
    #print URDME_ROOT
