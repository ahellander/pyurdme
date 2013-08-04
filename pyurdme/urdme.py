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
    print "Fatal Error: pyurdme requires h5py."
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
        
        self.urdme_solver_data = {'initialized':False}
        
        self.tspan = None
        self.mesh = None

    def createStoichiometricMatrix(self):
        """ Create a sparse stoichiometric matrix. """
        ND = np.zeros((self.getNumSpecies(),self.getNumReactions()))
        i=0
        for r in self.listOfReactions:
            R = self.listOfReactions[r]
            reactants = R.reactants
            products  = R.products
            
            for s in reactants:
                ND[self.species_map[s],i]=-reactants[s]
            for s in products:
                ND[self.species_map[s],i]=products[s]
            i = i+1
    
        N = scisp.csc_matrix(ND)
        return N

    def createDependencyGraph(self):
        """ Construct the sparse dependecy graph. """
        #TODO: Create a better dependency graph
        GF = np.ones((self.getNumReactions(),self.getNumReactions()+self.getNumSpecies()))
        G=scisp.csc_matrix(GF)
        return G
    
    def createNewPropensityFile(self,file_name=None):
        """ Generate a C propensity file on the new experimental format. """
        
        template = open(os.path.abspath(os.path.dirname(__file__))+'/data/propensity_file_new_template.c','r')
        propfile = open(file_name,"w")
        propfilestr = template.read()

        propfilestr = propfilestr.replace("__NUMBER_OF_REACTIONS__",str(self.getNumReactions()))
        propfilestr = propfilestr.replace("__NUMBER_OF_SPECIES__",str(len(self.listOfSpecies)))

        
        speciesdef = ""
        for i,sname in enumerate(self.listOfSpecies):
            S = self.listOfSpecies[sname]
            speciesdef += "species *"+sname+";\n\t"
            speciesdef += sname+"= (species *)malloc(sizeof(species));\n\t"
            speciesdef += sname+"->gamma = "+str(S.diffusion_constant)+";\n\t"
            speciesdef += sname+"->sigma = "+str(S.reaction_radius)+";\n\t"
            speciesdef += "ptr["+str(i)+"] = "+sname +";\n\n\t"
            
                                
        propfilestr = propfilestr.replace("__DEFINE_SPECIES__",speciesdef)
        
        
        # Make sure all paramters are evaluated to scalars before we write them to the file.
        self.resolveParameters()
                
        reacstr = ""
                
        for j,sname in enumerate(self.listOfSpecies):
            reacstr += "int "+sname+"="+str(j)+";\n\t"
        
        reacstr += "\n\t"
                
        for i,rname in enumerate(self.listOfReactions):
            R=self.listOfReactions[rname]
            
            reacstr += "reaction *"+rname+";\n\t"
            reacstr += rname+"=(reaction *)malloc(sizeof(reaction));\n\t"
            reacstr += rname+"->order="+str(len(R.reactants))+";\n\t"
            reacstr += rname+"->nr_reactants="+str(len(R.reactants))+";\n\t"
            reacstr += rname+"->nr_products="+str(len(R.products))+";\n\t"
            
            reacstr += rname+"->reactants=(int *)malloc("+rname+"->nr_reactants*sizeof(int));\n\t"
            for j,reactant in enumerate(R.reactants):
                reacstr += rname+"->reactants["+str(j)+"]="+str(reactant)+";\n\t"
            
            reacstr += "\n\t"+rname+"->products=(int *)malloc("+rname+"->nr_products*sizeof(int));\n\t"
            for j,product in enumerate(R.products):
                reacstr += rname+"->products["+str(j)+"]="+str(product)+";\n\t"
    
            reacstr += "\n\t"+rname+"->nr=(int *)calloc("+str(len(self.listOfSpecies))+",sizeof(int));\n\t"
            for j,reactant in enumerate(R.reactants):
                 reacstr += rname+"->nr["+reactant+"]="+str(R.reactants[reactant])+";\n\t"
            for j,product in enumerate(R.products):
                reacstr += rname+"->nr["+product+"]=-"+str(R.reactants[reactant])+";\n\t"

            reacstr += rname+"->k="+str(R.marate.value)+";\n\t"

            reacstr += "\n\tptr["+str(i)+"] = "+rname +";\n\n\t"
                
                
        propfilestr = propfilestr.replace("__DEFINE_REACTIONS__",reacstr)
        
        propfile.write(propfilestr)
        propfile.close()

    
    def createPropensityFile(self,file_name=None):
        """ Generate a C propensity file for used to compile the URDME solvers. """

        template = open(os.path.abspath(os.path.dirname(__file__))+'/data/propensity_file_template.c','r')
        propfile = open(file_name,"w")
        propfilestr = template.read()
        
        speciesdef = ""
        i=0
        for S in self.listOfSpecies:
            speciesdef += "#define "+S+" " +"x["+str(i)+"]"+"\n"
            i+=1
        
        propfilestr = propfilestr.replace("__DEFINE_SPECIES__",speciesdef)
    
        propfilestr = propfilestr.replace("__NUMBER_OF_REACTIONS__",str(self.getNumReactions()))
        
        # Make sure all paramters are evaluated to scalars before we write them to the file.
        self.resolveParameters()
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
        return np.ones((1,self.mesh.getNumVoxels()))
    
    def initializeInitialValue(self):
        """ Create all-zeros inital condition matrix. """
        ns = self.getNumSpecies()
        nv = self.mesh.getNumVoxels()
        self.u0 = np.zeros((ns,nv))
      
    
    def scatter(self,species,subdomain=None):
        """ Scatter an initial number of molecules over the voxels in a subdomain. """
    
        Sname = species.name
        numS = species.initial_value
        if not hasattr(self,'species_map'):
            createSpeciesMap(self)
        specindx= self.species_map[Sname]
        
        if not hasattr(self,"u0"):
            self.initializeInitialValue()
            
        for mol in range(species.initial_value):
            vtx=np.random.randint(0,self.mesh.getNumVoxels())
            self.u0[specindx,vtx]+=1
        

    def createSystemMatrix(self):
        """ Create one merged system matrix in CCS format for input to the URDME solvers """
        
        # Check if the individual stiffness and mass matrices (per species) have been assembled, otherwise assemble
        try:
            stiffness_matrices = self.stiffness_matrices
            mass_matrices = self.mass_matrices
        except:
            matrices = assemble(self)
            self.stiffness_matrices = matrices['K']
            self.mass_matrices = matrices['M']
            stiffness_matrices = self.stiffness_matrices
            mass_matrices = self.mass_matrices
        
        # Make a dok matrix for easy manipulation
        i=1;
        Mspecies = len(self.listOfSpecies)
        Ndofs = self.mesh.getNumVoxels()*Mspecies
        S = scipy.sparse.dok_matrix((Ndofs,Ndofs))

        # Create the volume vector by lumping the mass matrices
        vol = numpy.zeros((Ndofs,1))
        spec = 0
        
        for species,M in mass_matrices.iteritems():
            rows,cols,vals = M.data()
            SM = scipy.sparse.csr_matrix((vals,cols,rows))
            vols = SM.sum(axis=1)
            spec = self.species_map[species]
            for j in range(len(vols)):
                vol[Mspecies*j+spec,0]=vols[j]

        vol = vol.flatten()
        # Assemble one big matrix from the indiviudal stiffness matrices. Multiply by the inverse of
        # the lumped mass matrix, filter out any entries with the wrong sign and renormalize the columns.
        spec = 0
        positive_mass = 0.0
        total_mass = 0.0
                
        for species,K in stiffness_matrices.iteritems():

            rows,cols,vals = K.data()
           
            Kcrs = scipy.sparse.csr_matrix((vals,cols,rows))
            Kdok = Kcrs.todok()

            for entries in Kdok.items():
                ind = entries[0]
                val = entries[1]
                
                if ind[0] != ind[1]:
                    if val > 0.0:
                        positive_mass += val
                        val = 0.0
                    else:
                        total_mass += val
                        
                # The volume can be zero, if the species is not active at the vertex (such as a 2D species at a 3D node)
                if vol[Mspecies*ind[1]+spec]==0:
                    vi = 1
                else:
                    vi = vol[Mspecies*ind[1]+spec]
            
                S[Mspecies*ind[0]+spec,Mspecies*ind[1]+spec]=-val/vi
            
            spec = spec+1

        # Convert to compressed column for compatibility with the URDME solvers.
        D = S.tocsc()
                
        # Renormalize the columns (may not sum to zero since elements may have been filtered out
        sumcol = numpy.zeros((Ndofs,1))
        for i in range(Ndofs):
           col = D.getcol(i)
           for val in col.data:
               if val > 0.0:
                   sumcol[i] += val
                
        D.setdiag(-sumcol.flatten())

    
        print "Fraction of positive off-diagonal entries: " + str(numpy.abs(positive_mass/total_mass))
        return {'vol':vol,'D':D,'relative_positive_mass':positive_mass/total_mass}

                
    def validate(self):
        """ Validate the model data structures. This function should be called 
            prior to writing the model to the input file for the URDME solvers
            since the solvers themselves only do limited error checking. """
    
        # Check that all the columns of the system matrix sums to zero (or close to zero). If not, it does
        # not define a Markov process. 
        maxcolsum = numpy.max(numpy.abs(self.urdme_solver_data['D'].sum(axis=0)))
        if maxcolsum > 1e-10:
            raise InvalidSystemMatrixException("Invalid diffusion matrix. The sum of the columns does not sum to zero. " + str(maxcolsum))
                

    def initialize(self):
        """ Create the datastructures needed by the URDME solvers and put them
            in the urdme_solver_data dict. The content of this dict is what will be dumped
            to the urdme input file.  """
        
        if not self.urdme_solver_data['initialized']:
        
            createSpeciesMap(self)
            
            # Stoichimetric matrix
            N = self.createStoichiometricMatrix()
            self.urdme_solver_data['N'] = N
            # Dependency Graph
            G = self.createDependencyGraph()
            self.urdme_solver_data['G']  = G
            
            # Volume vector
            result =  self.createSystemMatrix()
            vol = result['vol']
            #for v in vol:
            #    print v
            vol = vol[::len(self.listOfSpecies)]
            
            self.urdme_solver_data['vol'] = vol
            D = result['D']
            self.urdme_solver_data['D'] = D
            
            # Subdomain vector
            sd = self.initializeSubdomainVector()
            self.urdme_solver_data['sd'] = sd
            
            # Data vector. If not present in model, it defaults to a vector with all elements zero.
            if not "data" in self.urdme_solver_data:
                data = np.zeros((1,self.mesh.getNumVoxels()))
                self.urdme_solver_data['data'] = data
    
            self.urdme_solver_data['u0'] = self.u0

            tspan= np.asarray(self.tspan,dtype=np.float)
            self.urdme_solver_data['tspan'] = tspan
        
            # Vertex coordinates
            self.urdme_solver_data['p'] = self.mesh.getVoxels()
        
            # Connectivity matrix
            self.urdme_solver_data['K'] = connectivityMatrix(self)

            self.urdme_solver_data['initialized'] = True
    
            
    def serialize(self,filename=None):
        """ 
            Write the datastructures needed by the the core URDME solvers to a .mat input file.
            initialize() must be called prior to calling this function. 
        """
                
        # Validate the data structures before writing them to file. 
        self.validate()
        spio.savemat(filename,self.urdme_solver_data,oned_as='column')


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
    
    def getVoxels(self):
        return self.mesh.coordinates()

def read_gmsh_mesh(meshfile):
    
    """ Read a Gmsh mesh from file. """
    mr = GmshMeshReceiverBase()
    try:
        mesh = read_gmsh(mr,filename=meshfile)
    except:
        raise MeshImportError("Failed to import mesh: "+filename)

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

def connectivityMatrix(model):
    """ Assemble a sparse CCS connectivity matrix for a mesh attached to a model. """


    fs = dolfin.FunctionSpace(model.mesh.mesh,"Lagrange",1)
    trial_function = dolfin.TrialFunction(fs)
    test_function = dolfin.TestFunction(fs)
    a_M = trial_function*test_function*dolfin.dx
    C = dolfin.assemble(a_M)
    rows,cols,vals = C.data()
    C = scipy.sparse.csr_matrix((vals,cols,rows))
    C = C.tocsc()

    return C

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
        
        #try:
        #    xmesh = model.xmesh
        #except:
        #    model.xmesh = meshextend(model)
       
        #function_space = xmesh.function_space
        #trial_functions = xmesh.trial_functions
        #test_functions = xmesh.test_functions
            
        function_space = OrderedDict()
        trial_functions = OrderedDict()
        test_functions = OrderedDict()
        stiffness_matrices = OrderedDict()
        mass_matrices = OrderedDict()
        
        for spec in model.listOfSpecies:
            
            species = model.listOfSpecies[spec]
            spec_name = species.name
            
            if species.dimension == 2:
                # TODO: If the dimension of the mesh is 2 (triangles) and one uses ds,
                # The the mass matrices become strange...
                differential = dolfin.dx
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
    """ Extended mesh object. Contanins dof-mappings, and function spaces etc. """

    def __init__(self):
        self.dofs = {}
        self.nodes = {}
        self.function_space = {}

        #dofs.coords = []
        #dofs.names = []
        #nodes.dofs = []


def createSpeciesMap(model):
    i=0
    model.species_map = {}
    for S in model.listOfSpecies:
        model.species_map[S]=i
        i = i+1;

def toXYZ(model,filename,format="VMD"):
    """ Dump the solution attached to a model as a xyz file. This format can be
        read by e.g. VMD, Jmol and Paraview. """
    
    
    if 'U' not in model.__dict__:
        print "No solution found in the model."
        raise

    #outfile = open(filename,"w")
    dims = numpy.shape(model.U)
    Ndofs = dims[0]
    Mspecies = len(model.listOfSpecies)
    Ncells = Ndofs/Mspecies

    coordinates = model.mesh.getVoxels()
    coordinatestr = coordinates.astype(str)

    if format == "VMD":
        outfile = open(filename,"w")
        filestr = ""
        for i,time in enumerate(model.tspan):
            number_of_atoms = numpy.sum(model.U[:,i])
            filestr += (str(number_of_atoms)+"\n"+"timestep "+str(i) + " time "+str(time)+"\n")
            for j,spec in enumerate(model.listOfSpecies):
                for k in range(Ncells):
                    for mol in range(model.U[k*Mspecies+j,i]):
                        linestr = spec + "\t" + '\t'.join(coordinatestr[k,:]) +"\n"
                        filestr += linestr

        outfile.write(filestr)
        outfile.close()

    elif format == "ParaView":
        foldername = filename
        subprocess.call(["mkdir",foldername])
        for i,time in enumerate(model.tspan):
            outfile = open(foldername+"/"+filename+"."+str(i),"w")
            number_of_atoms = numpy.sum(model.U[:,i])
            filestr = ""
            filestr += (str(number_of_atoms)+"\n"+"timestep "+str(i) + " time "+str(time)+"\n")
            for j,spec in enumerate(model.listOfSpecies):
                for k in range(Ncells):
                    for mol in range(model.U[k*Mspecies+j,i]):
                        linestr = spec + "\t" + '\t'.join(coordinatestr[k,:]) +"\n"
                        filestr += linestr
            outfile.write(filestr)
            outfile.close()


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


def urdme(model=None,solver='nsm',solver_path="", model_file=None, seed=None,report_level=1):
    """ URDME solver interface, analogous to the Matlab URDME interface. """

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
    if model_file == None:
        propfilename= model.name+'_pyurdme_generated_model'
        if solver != "nem":
            model.createPropensityFile(file_name='.urdme/'+propfilename+'.c')
        else:
            model.createNewPropensityFile(file_name='.urdme/'+propfilename+'.c')
    else:
        subprocess.call(['cp',model_file,'.urdme/'+propfilename+'.c'])

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

    subprocess.call(['cp',infile.name,'.urdme/'])

    # Execute the solver
    if seed is not None:
     try: 
      handle = subprocess.Popen(['.urdme/'+propfilename+'.'+solver,infile.name,outfile.name,str(seed)], stdout = subprocess.PIPE, stderr=subprocess.PIPE)
      if report_level >= 1:
        print handle.stdout.read()
        print handle.stderr.read()
     except:
      return 'Call to URDME failed.'
    else:
      handle = subprocess.Popen(['.urdme/'+propfilename+'.'+solver,infile.name,outfile.name], stdout = subprocess.PIPE, stderr=subprocess.PIPE)
      if report_level >= 1:
        print handle.stdout.read()
        print handle.stderr.read()

    subprocess.call(['cp',infile.name,'./debug_input.mat'])
    subprocess.call(['cp',outfile.name,'./debug_output.mat'])

    #Load the result from the hdf5 output file.
    try:
        
        resultfile = h5py.File(outfile.name,'r')
        # print resultfile['U']
        U = resultfile['U']
        U = numpy.array(U)
        # This little hack makes U have the same structure as in the Matlab interface...
        dims = numpy.shape(U)
        U = U.reshape((dims[1],dims[0]))
        U = U.transpose()
        model.U = U
        
        tspan = resultfile['tspan']
        tspan = numpy.array(tspan)
        
        # Create Dolfin Functions for all the species
        model.sol = {}
        for i,spec in enumerate(model.listOfSpecies):
    
            species = model.listOfSpecies[spec]
            spec_name = species.name
            func = dolfin.Function(dolfin.FunctionSpace(model.mesh.mesh,"Lagrange",1))
            func_vector = func.vector()
            dims = U.shape
            
            numvox = model.mesh.getNumVoxels()
            for dof in range(numvox):
                func_vector[dof] = float(U[dof*len(model.listOfSpecies)+i,10])
        
            model.sol[spec_name] = func

        # Clean up
        resultfile.close()
        os.remove(infile.name)
        os.remove(outfile.name)

        return {'U':U,'tspan':tspan}

    except Exception,e:
       # Clean up
       subprocess.call(['rm','-rf',infile.name])
       subprocess.call(['rm','-rf',outfile.name])
       raise


class URDMEError(Exception):
    pass

if __name__ == '__main__':
    """ Command line interface to URDME. Execute URDME given a model file. """ 


class InvalidSystemMatrixException(Exception):
    pass



