from model import *
import numpy as np
import scipy.sparse as scisp
import scipy.io as spio
import subprocess
import os
import tempfile
import re
import sys
import shutil

import gmsh

# Need a way to read hdf5 files
try:
    import h5py
except:
    pass


try:
    from meshpy.gmsh_reader import *
    #import meshpy
#import meshpy.gmsh_reader
#import meshpy.gmsh_reader.GmshMeshRecieverBase
except Exception,e:
    print "Failed to import meshpy."+e
    raise

try:
    import dolfin
except:
    ONLY_CARTESIAN=True

class MeshImportError(Exception):
    pass

class URDMEModel(Model):
    """ 
        An URDME Model extends Model with spatial information and methods 
        to create URDME solver input.
    """
    def __init__(self,name=""):
	Model.__init__(self,name)

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
    
    def createPropensityFile(self,file_name=None):
        """ Generate a C propensity file for used to compile the URDME solvers. """
        
        template = open(os.path.abspath(os.path.dirname(__file__))+'/data/propensity_file_template.c','r')
        propfile = open(file_name,"w")
        propfilestr = template.read()
        #propfile.write(templatestr)
        
        speciesdef = ""
        i=0
        for S in self.listOfSpecies:
            speciesdef += "#define "+S+" " +"x["+str(i)+"]"+"\n"
            i+=1
        
        propfilestr = propfilestr.replace("__DEFINE_SPECIES__",speciesdef)
        
        propfilestr= propfilestr.replace("__NUMBER_OF_REACTIONS__",str(self.getNumReactions()))
        
        parameters = ""
        for p in self.listOfParameters:
            parameters += "const double "+p+" = " +str(self.listOfParameters[p].value)+";\n"
        propfilestr=propfilestr.replace("__DEFINE_PARAMETERS__",str(parameters))
    
        # Reactions
        funheader = "double __NAME__(const int *x, double t, const double vol, const double *data, int sd)"
    
        funcs = ""
        funcinits = ""
        i = 0
        for R in self.listOfReactions:
            func = ""
            rname=self.listOfReactions[R].name
            func += funheader.replace("__NAME__",rname) + "\n{\n"
            func += "    return " + self.listOfReactions[R].propensity_function + ";"
            func +="\n}"
            funcs += func + "\n\n"
            funcinits += "    ptr["+str(i)+"] = " + rname +";\n"
            i+=1
              
        propfilestr = propfilestr.replace("__DEFINE_REACTIONS__",funcs)
        propfilestr = propfilestr.replace("__DEFINE_PROPFUNS__",funcinits)
                
                
        propfile.write(propfilestr)
        propfile.close()
    
    
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

def read_gmsh_mesh(meshfile):
    """ Read a Gmsh mesh file """
    mr = GmshMeshReceiverBase()
    try:
        mesh = read_gmsh(mr,filename=meshfile)
    except:
        raise MeshImportError("Failed to import mesh: "+filename)

    print mesh
    return mesh

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


def urdme(model=None,solver='nsm',seed=None,report_level=1):
    """ URDME solver interface, analogous to the Matlab URDME function interface. """

    # Set URDME_ROOT
    URDME_ROOT = subprocess.check_output(['urdme_init','-r'])
    # Trim newline
    URDME_ROOT = URDME_ROOT[:-1]
    URDME_BUILD = URDME_ROOT + '/build/'

    # Write the propensity file
    try:
      os.mkdir('.urdme')
    except:
      pass
    
    propfilename= model.name+'_pyurdme_generated_model'
    model.createPropensityFile(file_name='.urdme/'+propfilename+'.c')

    # Build the solver
    makefile = 'Makefile.' + solver
    handle = subprocess.Popen(['make','-f',URDME_BUILD+makefile,'URDME_ROOT='+URDME_ROOT,'URDME_MODEL='+propfilename], stdout = subprocess.PIPE, stderr=subprocess.PIPE)
    if report_level >=1:
      print handle.stdout.read()
      print handle.stderr.read()

    # Get temporary input and output files
    infile = tempfile.NamedTemporaryFile(delete=False)
    model.serialize(filename=infile)
    infile.close
    outfile = tempfile.NamedTemporaryFile(delete=False)
    outfile.close()
    
    # Execute the solver
    
    if seed is not None:
     try: 
      handle = subprocess.Popen(['.urdme/'+propfilename+'.'+solver,infile.name,'slask.mat',str(seed)], stdout = subprocess.PIPE, stderr=subprocess.PIPE)
      if report_level >= 1:
        print handle.stdout.read()
        print handle.stderr.read()
     except:
      return 'Call to URDME failed miserably'
    else:
      handle = subprocess.Popen(['.urdme/'+propfilename+'.'+solver,infile.name,'slask.mat'], stdout = subprocess.PIPE, stderr=subprocess.PIPE)
      if report_level >= 1:
        print handle.stdout.read()
        print handle.stderr.read()

    #Load the result.
    #AH: TODO! SciPy fails to read the file (But it loads fine in Matlab)!.
    try:
        # Pytables (Matlab >= 7.3)
        #file = tables.openFile('slask.mat')
        #U = file.root.U[:]
        #tspan = file.root.tspan[:]
        # SciPy (Matlab <= 7.1)
        #result = spio.loadmat('slask.mat',mat_dtype=True,matlab_compatible=True)
        resultfile = h5py.File('testh5.h5','r')
        U = resultfile['U'].value
        #print U
        #print np.shape(U)
        tspan = resultfile['tspan'].value
        #print tspan
        
        # Clean up
        #shutil.rmtree(infile.name)
        #shutil.rmtree(outfile.name)
        return resultfile
    except Exception,e:
       # Clean up
       subprocess.call(['rm','-rf',infile.name])
       subprocess.call(['rm','-rf',outfile.name])
       raise
       return 'Matfile load failed.'


class URDMEError(Exception):
    pass

if __name__ == '__main__':
    """ Command line interface to URDME. Execute URDME given a model file. """ 





