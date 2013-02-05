from model import *
import numpy as np
import scipy.sparse as scisp
import scipy.io as spio
import subprocess
import os
import tempfile

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
        ND = np.zeros((self.getNumSpecies(),self.getNumReactions()))
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
        GF = np.ones((self.getNumReactions(),self.getNumReactions()+self.getNumSpecies()))
        self.G=scisp.csc_matrix(GF)
    
    def createPropensityFile(self):
        """ Generate a C propensity file for used to compile the URDME solvers. """
        template = load()
    
    def initializeSubdomainVector(self):
        """ Create URDME 'sd' vector. """
        # TODO: Support arbitrary sd-numbers and more than one subdomain
        self.sd = np.ones((1,self.mesh.getNumVoxels()))
    
    def initializeInitialValue(self):
        """ Create all-zeros inital condition matrix. """
        ns = self.getNumSpecies()
        nv = self.mesh.getNumVoxels()
        self.u0 = np.zeros((ns,nv))
    
    
    def scatter(self,species,subdomain=None):
        """ Scatter an initial number of molecules over the 
            voxels in a subdomain. """
    
        Sname = species.name
        numS = species.initial_value
        specindx= self.species_map[Sname]
        
        if not hasattr(self,"u0"):
            self.initializeInitialValue()
        
        # TODO: USE THE SUBDOMAIN INFO
        for i in range(numS):
            vtx=np.random.randint(0,self.mesh.getNumVoxels())
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
        
        # Subdomain vector
        if not hasattr(self,"sd"):
            self.initializeSubdomainVector()
        
        # Data vector. If not present in model, it defaults to a vector with all
        # elements zero.
        if not hasattr(self,"data"):
            data = np.zeros((1,self.mesh.getNumVoxels()))
        
        # data = []
        filename = filename
        spio.savemat(filename,{'N':self.N,'u0':self.u0,'G':self.G,'sd':self.sd,'D':D,'vol':vol[0::3],'tspan':np.asarray(self.tspan,dtype=np.float),'data':data},oned_as='column')


class Mesh():
    """ A thin wrapper around the Dolfin mesh object. """

    def __init__(self,mesh_type="dolfin"):
        
        self.mesh_type = mesh_type
        
        if mesh_type == "dolfin":
            self.mesh = dolfin.Mesh()
        elif mesh_type == "Cartesian":
            return

    def getNumVoxels(self):
        dims = np.shape(self.p)
        return dims[1]

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

        Ndofs = model.xmesh.Ndofs
        DF  = np.zeros((Ndofs,Ndofs),dtype=np.float)
        vol = np.zeros(Ndofs)
        dim,nt = np.shape(model.mesh.t)
        nodes2dofs = model.xmesh.nodes["dofs"]
        dof_names = model.xmesh.dofs["names"]
        
        for t in range(nt):
            
            nodes = model.mesh.t[:,t]
            
            # Stupid numpy indexing
            coords = model.mesh.p[:,nodes]
            coords = coords[0]
            h = abs(coords[0]-coords[1])
            h2 = pow(h,2)
            
            dofs  = nodes2dofs[:,nodes]
            dims= np.shape(dofs)
            ns = dims[0]

            for j in range(ns):
                dof = dofs[j]
                vol[dof[0]]+=h/2
                vol[dof[1]]+=h/2
                
                # Diffusion constant
                d1 = model.listOfSpecies[dof_names[j]].diffusion_constant

                
                DF[dof[0],dof[1]]  = d1/h2
                DF[dof[0],dof[0]] -= d1/h2
                DF[dof[1],dof[0]]  = d1/h2
                DF[dof[1],dof[1]] -= d1/h2

        D = scisp.csc_matrix(DF)
        return (vol,D)


class Xmesh():
    """ Extended mesh object. """

    def __init__(self):
        self.dofs = {}
        self.nodes = {}
        #dofs.coords = []
        #dofs.names = []
        #nodes.dofs = []

def meshextend(model):
    """
        Extend the primary mesh with information about degrees of freedom.
        Initialize URDME datastructures that depend on the size of the mesh.
        """
    
    xmesh = Xmesh()
    
    # Construct a species map (dict mapping model species name to an integer index)
    i=0
    model.species_map = {}
    for S in model.listOfSpecies:
        model.species_map[S]=i
        i = i+1;
    
    dims = np.shape(model.mesh.p)
    Nvoxels = dims[1]
    xmesh.Ndofs = Nvoxels*len(model.species_map)
    
    #dof_names = np.zeros(len(model.species_map),dtype=np.int)
    dof_names = [""]*len(model.species_map)
    for S in model.species_map:
        dof_names[model.species_map[S]] = S

    dofs = np.zeros((len(model.species_map),Nvoxels))
    nodes = np.zeros(xmesh.Ndofs)
    dof = 0
    for i in range(Nvoxels):
        for S in model.species_map:
             dofs[model.species_map[S],i]=dof
             nodes[dof]=i
             dof+=1
                 
    xmesh.nodes["dofs"]=dofs
    xmesh.dofs["nodes"]=nodes
    xmesh.dofs["names"]=dof_names
   
    model.xmesh = xmesh



def urdme(model=None,solver='nsm'):
    """ URDME solver interface, similar to the Matlab function interface. """

    # Shell out and compile the solver
    
    # Set URDME_ROOT
    URDME_ROOT = subprocess.check_output(['urdme_init','-r'])
    # Trim newline
    URDME_ROOT = URDME_ROOT[:-1]
    URDME_BUILD = URDME_ROOT + '/build/'

    # Compile the solver
    makefile = 'Makefile.' + solver
    subprocess.call(['make','-f',URDME_BUILD+makefile,'URDME_ROOT='+URDME_ROOT,'URDME_MODEL='+'dimerization'],stderr=subprocess.STDOUT)

    # Get temporary input and output files
    infile = tempfile.NamedTemporaryFile(delete=False)
    model.serialize(filename=infile)
    infile.close
    outfile = tempfile.NamedTemporaryFile(delete=False)
    outfile.close()

    # Execute the solver
    subprocess.call(['.urdme/dimerization.nsm',infile.name,'slask.mat'],stderr=subprocess.STDOUT)
    #Load the result.
    #AH: TODO! SciPy fails to read the file (But it loads fine in Matlab)!.
    try:
        result = spio.loadmat('slask.mat')
        i = result["iU"]
        j = result["jU"]
        data = result["sU"]
        M = result["mU"]
        N = result["nU"]
        ij = [i,j]
        model.U = scisp.csc_matrix((data,ij),shape=(M,N))
        result.close()
    except:
        pass

    # Clean up
    subprocess.call(['rm','-rf',infile.name])
    subprocess.call(['rm','-rf',outfile.name])


class URDMEError(Exception):
    pass





