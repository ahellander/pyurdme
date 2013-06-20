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

import numpy
import scipy.sparse


# Need a way to read hdf5 files
try:
    import h5py
except:
    print "Fatal Error: You are missing the h5py library."
    raise

try:
    import dolfin
    dolfin.parameters["linear_algebra_backend"] = "uBLAS"
except:
    print "Warning: Could not import dolphin. Only simple Cartesain examples will work."
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
        

    def createSystemMatrix(self):
        """ Create one merged system matrix in CCS format for input to the URDME solvers """
                
        stiffness_matrices = self.stiffness_matrices
        mass_matrices = self.mass_matrices
        # Make a dok matrix
        #Ddok = scipy.sparse.dok_matrix((Ndofs,Ndofs))
        i=1;
        Mspecies = len(self.listOfSpecies)
        Ndofs = self.mesh.getNumVoxels()*Mspecies
        S = scipy.sparse.dok_matrix((Ndofs,Ndofs))

        # Create the volume vector from the mass matrices
        vol = numpy.zeros((Ndofs,1))
        spec = 0
        for species,M in mass_matrices.iteritems():
            rows,cols,vals = M.data()
            SM = scipy.sparse.csr_matrix((vals,cols,rows))
            vols = SM.sum(axis=1)
            for j in range(0,len(vols)):
                vol[Mspecies*j+spec,0]=vols[j]
            spec = 0

        vol = vol.flatten()
        
        spec = 0
        for species,K in stiffness_matrices.iteritems():

            rows,cols,vals = K.data()
           
            Kcrs = scipy.sparse.csr_matrix((vals,cols,rows))
            Kdok = Kcrs.todok()
            for entries in Kdok.items():
                ind = entries[0]
                val = entries[1]
                if ind[0] != ind[1] and val > 0.0:
                    val = 0.0
                
                S[Mspecies*ind[0]+spec,Mspecies*ind[1]+spec]=-val/vol[Mspecies*ind[0]+spec]
        
            spec = spec+1

        D = S.tocsc()
        return {'vol':vol,'D':D}

                
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
        
        # Data vector. If not present in model, it defaults to a vector with all elements zero.
        if not hasattr(self,"data"):
            data = np.zeros((1,self.mesh.getNumVoxels()))
        
        # data = []
        filename = filename
        spio.savemat(filename,{'N':self.N,'u0':self.u0,'G':self.G,'sd':self.sd,'D':D,'vol':vol[0::3],'tspan':np.asarray(self.tspan,dtype=np.float),'data':data},oned_as='column')



class Mesh():
    """ A thin wrapper around the Dolfin mesh object. """

    def __init__(self,mesh=None,mesh_type="Dolfin"):
        
        self.mesh_type = mesh_type
        
        if mesh_type == "Dolfin":
            self.mesh = mesh
        elif mesh_type == "Cartesian":
            return

    def getNumVoxels(self):

        return self.mesh.num_vertices()
        #dims = np.shape(self.p)
        #return dims[1]

def read_gmsh_mesh(meshfile):
    """ Read a Gmsh mesh from file. """
    mr = GmshMeshReceiverBase()
    try:
        mesh = read_gmsh(mr,filename=meshfile)
    except:
        raise MeshImportError("Failed to import mesh: "+filename)

    print mesh
    return mesh

def read_dolfin_mesh(filename=None):
    """ Import a mesh in Dolfins native .xml format """
    try:
        dolfin_mesh = dolfin.Mesh(filename)
        mesh = Mesh(mesh=dolfin_mesh,mesh_type="Dolfin")
        return mesh
    except Exception,e:
        raise MeshImportError("Failed to import mesh: "+filename+"\n"+e)


def createCartesianMesh(geometry=None,side_length=None,hmax=None):
    """
    Create a simple Cartesian mesh for a line, square or cube.
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
        Here, we will need to use Dolfin to build function spaces
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
    else:
        # Assemble using Dolfin.
        #
        #
        # Returns: A dict of two dicts (stiffness and mass matrices)
        #          Those dicts in turn has species names as keys and contain
        #          matrices in CSR format (as returned ny Dolfin)
        #
        # TODO: If the mesh is not a Dolfin mesh object, we need to convert to that here...
        
        # Create Function spaces, trial functions and test functions for all the species
        function_space = {}
        trial_functions = {}
        test_functions = {}
        stiffness_matrices = {}
        mass_matrices = {}
        
        for spec in model.listOfSpecies:
            
            species = model.listOfSpecies[spec]
            spec_name = species.name
            
            if species.dimension == 2:
                differential = dolfin.ds
            else:
                differential = dolfin.dx
        
            function_space[spec_name] = dolfin.FunctionSpace(model.mesh.mesh,"Lagrange",1)
            trial_functions[spec_name] = dolfin.TrialFunction(function_space[spec_name])
            test_functions[spec_name] = dolfin.TestFunction(function_space[spec_name])
            a_K = species.diffusion_constant*dolfin.inner(dolfin.nabla_grad(trial_functions[spec_name]), dolfin.nabla_grad(test_functions[spec_name]))*differential
            stiffness_matrices[spec_name] = dolfin.assemble(a_K)
            a_M = trial_functions[spec_name]*test_functions[spec_name]*differential
            mass_matrices[spec_name] = dolfin.assemble(a_M)
        
        return {'K':stiffness_matrices,'M':mass_matrices}


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


def urdme(model=None,solver='nsm',solver_path="",seed=None,report_level=1):
    """ URDME solver interface, analogous to the Matlab URDME function interface. """

    # Set URDME_ROOT
    URDME_ROOT = subprocess.check_output(['urdme_init','-r'])
    # Trim newline
    URDME_ROOT = URDME_ROOT[:-1]
    if solver_path == "":
        URDME_BUILD = URDME_ROOT+'/build/'
    else:
        URDME_BUILD = solver_path+'/build/'
        os.environ['SOLVER_ROOT'] = solver_path

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
      handle = subprocess.Popen(['.urdme/'+propfilename+'.'+solver,infile.name,outfile.name,str(seed)], stdout = subprocess.PIPE, stderr=subprocess.PIPE)
      if report_level >= 1:
        print handle.stdout.read()
        print handle.stderr.read()
     except:
      return 'Call to URDME failed miserably'
    else:
      handle = subprocess.Popen(['.urdme/'+propfilename+'.'+solver,infile.name,outfile.name], stdout = subprocess.PIPE, stderr=subprocess.PIPE)
      if report_level >= 1:
        print handle.stdout.read()
        print handle.stderr.read()

    #Load the result.
    #AH: TODO! SciPy fails to read the file (But it loads fine in Matlab)!.
    try:
        resultfile = h5py.File(outfile.name,'r')
        U = resultfile['U'].value
        tspan = resultfile['tspan'].value
        resultfile.close()
        # Clean up
        os.remove(infile.name)
        os.remove(outfile.name)
        
        return {'U':U,'tspan':tspan}

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





