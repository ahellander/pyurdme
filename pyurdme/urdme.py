from model import *
import numpy as np
import scipy.sparse as scisp
import scipy.io as spio
import subprocess
import os

try:
    import dolfin
except:
    ONLY_CARTESIAN=True

class URDMEModel(Model):
    """ 
        An URDME Model extends Model with spatial information and methods 
        to create URDME solver input.
    """
    def __init__(self,name=""):
	Model.__init__(self)

	self.tspan = None
	self.mesh = None
	self.D = None
		

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
        """ Construct the sparse dependecy graph. """
        #TODO: Create a better dependency graph
        GF = np.ones((self.num_reactions,self.num_reactions+self.num_species))
        self.G=scisp.csc_matrix(GF)
    
    def createPropensityFile(self):
        """ Generate a C propensity file for used to compile the URDME solvers. """
    
    def initializeSubdomainVector(self):
        """ Create URDME 'sd' vector. """
        # TODO: Support arbitrary sd-numbers and more than one subdomain
        self.sd = np.ones((1,self.num_voxels))
    
    def assemble(self):
        """ 
            Assemble the diffusion matrix and volume vector.
            Here, we will need to use Fenics to build function spaces
            and assemble.
        """
    
        if self.mesh.mesh_type == "Cartesian":
            vol,D = self.assembleCartesian()
        elif self.mesh.mesh_type == "Dolfin":
            """ Assemble using Dolfin. """
        # Build function space
        #T = []
        #for s in self.num_species:
        #    T.append(dolfin.FunctionSpace(self.mesh,"CG",1))

        #V =
    
    def assembleCartesian(self):
        """ 
            Assemble the jump-coefficient matrix for the case of
            a Cartesian mesh and a simple geomtry. 
        """
        
        # Allocate space
        Ndofs = self.xmesh.Ndofs
        vol = np.zeros((Ndofs,1))
        DF = np.zeros((Ndofs,Ndofs))
    
        
    def meshextend(self):
        """ 
            Extend the primary mesh with information about degrees of freedom. 
            Initialize URDME datastructures that depend on the size of the mesh.
        """
        
        xmesh = Xmesh()
        
        # Construct a species map (dict mapping model species name to an integer index)
        i=0
        self.species_map = {}
        for S in self.listOfSpecies:
            self.species_map[S]=i
            i = i+1;
        
        dims = np.shape(self.mesh.p)
        xmesh.Ndofs = dims[1]*len(self.species_map)
        xmesh.dof_names = []
        #for i in self.listOfSpecies
        #     xmesh.dof_names[i]=
        xmesh.dofs = np.zeros((1,xmesh.Ndofs))
        
        # Initialize the Dolfin mesh object.
        # model.mesh.init()
        # self.num_voxels  = self.mesh.num_vertices()
        # self.num_species = len(self.listOfSpecies)
        #xmesh = URDMEXmesh()
        #xmesh.dofs =
        self.xmesh = xmesh
    
    
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
        vol,D = assemble(self)

        # Initial condition
        #u0 =


        # Volume vector
        # vol =

        # Data vector
        # data = []
        filename = filename
        spio.savemat(filename,{'num_species':self.num_species,'N':self.N,'u0':self.u0,'G':self.G,'sd':self.sd,'D':D,'vol':vol},oned_as='column')


class Mesh():
    """ A thin wrapper around the Dolfin mesh object. """

    def __init__(self,mesh_type="dolfin"):
        
        self.mesh_type = mesh_type
        
        if mesh_type == "dolfin":
            self.mesh = dolfin.Mesh()
        elif mesh_type == "Cartesian":
            return

def CartesianMesh(geometry=None,side_length=None,hmax=None):
    """
        Create a Cartesian mesh for a line, square or cube.
        This is more or less only useful for debugging and development.
        
    """
    valid_geometries = ["line","square","cube"]
    if geometry not in valid_geometries:
        raise
    
    mesh = Mesh(mesh_type="Cartesian")
    
    # Vertices
    N = int(side_length/hmax)+1

    if geometry == "line":
        dim = 1
        mesh.p = np.zeros((dim,N))
        mesh.e = np.zeros(2);
        mesh.t = np.zeros((2,N-1),dtype=np.int)

        # Points
        for i in range(N):
            x = i*hmax
            mesh.p[0,i]=x
            
        for i in range(N-1):
            mesh.t[0,i]=i
            mesh.t[1,i]=i+1

        # Edges
        mesh.e[0]=0
        mesh.e[1]=N-1


    elif geometry == "square":
        dim = 2
        Nvoxels = pow(N,dim)
        mesh.p = np.zeros((dim,Nvoxels))
        for i in range(N):
            x = i*hmax
            for j in range(N):
                y = j*hmax
                mesh.p[0,i*N+j]=x
                mesh.p[1,i*N+j]=y

    # Boundary segments
    #mesh.e =
    return mesh

def assemble(model):
    """
        Assemble the diffusion matrix and volume vector.
        Here, we will need to use Fenics to build function spaces
        and assemble.
    """
    if model.mesh.mesh_type == "Cartesian":
        # Allocate space
        Ndofs = model.xmesh.Ndofs
        vol = np.zeros((Ndofs,1))
        DF = np.zeros((Ndofs,Ndofs))
        vol = np.zeros(Ndofs)
        dim,nt = np.shape(model.mesh.t)
        
        for t in range(nt):
            nodes = model.mesh.t[:,t]
            # Stupid numpy indexing
            coords = model.mesh.p[:,nodes]
            coords = coords[0]
            h = abs(coords[0]-coords[1])
            h2 = pow(h,2)
            
            vol[nodes[0]]+=h2/2
            vol[nodes[1]]+=h2/2
            
            DF[nodes[0],nodes[1]] = 1.0/h2
            DF[nodes[0],nodes[0]] -= 1.0/h2
            DF[nodes[1],nodes[0]] = 1.0/h2
            DF[nodes[1],nodes[1]] -= 1.0/h2

        D = scisp.csc_matrix(DF)
        return (vol,D)


class Xmesh():
    """ Extended mesh object. """
    #dofs.coords
    #dofs.names
    

def urdme(model=None,solver='nsm'):
    """ URDME solver interface, similar to the Matlab function interface. """

    # Shell out and compile the solver
    
    # Set URDME_ROOT
    URDME_ROOT = subprocess.check_output(['urdme_init','-r'])
    # Trim newline
    URDME_ROOT = URDME_ROOT[:-1]

    #subprocess.call(['export','URDME_ROOT='+URDME_ROOT],stderr=subprocess.STDOUT,shell=True)
    #print subprocess.STDOUT
    #test = subprocess.check_output(['echo','$URDME_ROOT'])
    #print test

    # Compile the solver
    #makefile = 'Makefile.' + solver
    #subprocess.call(['make','-f',URDME_ROOT+'build/'+makefile],stderr=subprocess.STDOUT)
#print subprocess.STDOUT
