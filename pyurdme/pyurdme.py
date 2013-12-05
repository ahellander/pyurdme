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

try:
    import h5py
except:
    print "pyurdme requires h5py."
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
        An URDME Model extends Model with spatial information and methods to create URDME solver input.
        TODO: Documentiation.
    """
    
    def __init__(self,name=""):
        Model.__init__(self,name)
        
        # Currently not used
        self.geometry = None
        
        self.mesh = None
        
        # subdomins is a list of MeshFunctions with subdomain marker information
        self.subdomains = []

        # This dictionary hold information about the subdomains each species is active on
        self.species_to_subdomains = {}
        
        self.tspan = None
    
    def __initializeSpeciesMap(self):
        i=0
        self.species_map = {}
        for S in self.listOfSpecies:
            self.species_map[S]=i
            i = i+1;

    def speciesMap(self):
        """ Get the species map, name to index. """
        if not hasattr(self,'species_map'):
            self.__initializeSpeciesMap()
        
        return self.species_map

    def createStoichiometricMatrix(self):
        """ Generate a stoichiometric matrix in sparse CSC format. """
        
        if self.getNumReactions() > 0:
            ND = np.zeros((self.getNumSpecies(),self.getNumReactions()))
            for i,r in enumerate(self.listOfReactions):
                
                R = self.listOfReactions[r]
                reactants = R.reactants
                products  = R.products
                
                for s in reactants:
                    ND[self.species_map[s],i]-=reactants[s]
                for s in products:
                    ND[self.species_map[s],i]+=products[s]
    
            N = scisp.csc_matrix(ND)
        else:
            N = numpy.zeros((self.getNumSpecies(),self.getNumReactions()))

        return N

    def createDependencyGraph(self):
        """ Construct the sparse dependecy graph. """
        
        #TODO: Automatically create a dependency graph (cannot be optimal, but good enough.)
        GF = np.ones((self.getNumReactions(),self.getNumReactions()+self.getNumSpecies()))
        try:
            G=scisp.csc_matrix(GF)
        except:
            G=GF

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
            
            #print reacstr
            
            reacstr += rname+"->reactants=(int *)malloc("+rname+"->nr_reactants*sizeof(int));\n\t"
            for j,reactant in enumerate(R.reactants):
                reacstr += rname+"->reactants["+str(j)+"]="+str(reactant)+";\n\t"
            
            reacstr += "\n\t"+rname+"->products=(int *)malloc("+rname+"->nr_products*sizeof(int));\n\t"
            for j,product in enumerate(R.products):
                reacstr += rname+"->products["+str(j)+"]="+str(product)+";\n\t"
    
            reacstr += "\n\t"+rname+"->nr=(int *)calloc("+str(len(self.listOfSpecies))+",sizeof(int));\n\t"
                
            for j,reactant in enumerate(R.reactants):
                 reacstr += rname+"->nr["+reactant+"]=-"+str(R.reactants[reactant])+";\n\t"

            for j,product in enumerate(R.products):
                reacstr += rname+"->nr["+product+"]="+str(R.products[product])+";\n\t"
                
            reacstr += rname+"->k="+str(R.marate.value)+";\n\t"

            reacstr += "\n\tptr["+str(i)+"] = "+rname +";\n\n\t"
                
                
        propfilestr = propfilestr.replace("__DEFINE_REACTIONS__",reacstr)
        
        propfile.write(propfilestr)
        propfile.close()

    
    def createPropensityFile(self,file_name=None):
        """ Generate the C propensity file that is used to compile the URDME solvers.
            Only mass action propensities are supported. """
        
        
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
            if self.listOfReactions[R].restrict_to == None:
                func += "    return " + self.listOfReactions[R].propensity_function
                order = len(self.listOfReactions[R].reactants)
                if order == 2:
                    func += "/vol;"
                elif order == 0:
                    func += "*vol;"
                else:
                    func += ";"

            else:
                func += "if("
                for sd in self.listOfReactions[R].restrict_to:
                    func += "sd == "+str(sd)+"||"
                func = func[:-2]
                func += ")\n"
                func += "\treturn " + self.listOfReactions[R].propensity_function
                order = len(self.listOfReactions[R].reactants)
                if order == 2:
                    func += "/vol;"
                elif order == 0:
                    func += "*vol;"
                else:
                    func += ";"

                func += "\nelse"
                func += "\n\treturn 0.0;"
        
            
            func +="\n}"
            funcs += func + "\n\n"
            funcinits += "    ptr["+str(i)+"] = " + rname +";\n"
            i+=1
              
        propfilestr = propfilestr.replace("__DEFINE_REACTIONS__",funcs)
        propfilestr = propfilestr.replace("__DEFINE_PROPFUNS__",funcinits)
                
        propfile.write(propfilestr)
        propfile.close()
    
    def timespan(self, tspan):
        """ Set the time span of simulation. """
        self.tspan = tspan

    def _initialize_species_to_subdomains(self):
        """ Initialize the species mapping to subdomains. The default
            is that a species is active in all the defined subdomains. 
        """
        # The unique elements of the subdomain MeshFunctions
        sds = []
        for subdomain in self.subdomains:
            sds = sds + list(np.unique(subdomain.array()).flatten())
        sds = np.unique(sds)
        sds = list(sds)
        
        sds.remove(0)
        
        # This conversion is necessary for UFL not to choke on the subdomain ids.
        for i,sd in enumerate(sds):
            sds[i]=int(sd)

        # If a species is not present as key in the species_to_subdomain mapping,
        # we label it as active in all subdomains
        for spec_name in self.listOfSpecies:
            species = self.listOfSpecies[spec_name]
            if species not in self.species_to_subdomains.keys():
                self.species_to_subdomains[species] = sds
    
    def restrict(self,species,subdomains):
        self.species_to_subdomains[species] = subdomains


    def subdomainVector(self,subdomains=[]):
        """ Create the 'sd' vector. 'subdomains' is a dolfin FacetFunction,
            and if no subdomain input is specified, they voxels default to
            subdomain 1. """
        
        if not hasattr(self,'xmesh'):
            self.meshextend()
        #fs = dolfin.FunctionSpace(self.mesh,"Lagrange",1)
        # vertex_to_dof_map provides a map between the vertex index and the dof.
        #dofmap=fs.dofmap().vertex_to_dof_map(self.mesh)
        #dofmap = self.xmesh.vertex_to_dof_map["MinD_m"]
        #print dofmap
        self.mesh.init()
        #species_map = self.speciesMap()
        #specindx= species_map["MinD_m"]

        # TODO: Support arbitrary sd-numbers and more than one subdomain
        sd = numpy.zeros((1,self.mesh.getNumVoxels()))
        if subdomains == []:
            self.sd = sd.flatten()
        else:
            for subdomain in subdomains:
                # Map all facet labels to vertex labels
                tovertex = self.mesh.topology()(subdomain.dim(),0)
                for i in range(subdomain.size()):
                    for vtx in tovertex(i):
                        if subdomain[i]!=0: # TODO: Temporary hack to fix issue with Gmesh facet_region files.
                            #dof = dofmap[vtx]
                            #ix = (dof-specindx)/len(species_map)
                            sd[0,vtx] = subdomain[i] #?????????
            
        self.sd = sd.flatten()
        return self.sd


    def initializeInitialValue(self):
        """ Create all-zeros inital condition matrix. """
        ns = self.getNumSpecies()
        nv = self.mesh.getNumVoxels()
        self.u0 = np.zeros((ns,nv))
    
    def meshextend(self):
        """ Extend the primary mesh with information about the degrees of freedom.
            
            TODO: Docs...
            
            """
        
        xmesh = Xmesh()
        
        # Construct a species map (dict mapping model species name to an integer index)
        species_map=self.speciesMap()
        
        # Initialize the function spaces and dof maps.
        for spec in self.listOfSpecies:
            
            species = self.listOfSpecies[spec]
            spec_name = species.name
            spec_index = species_map[spec_name]
            
            xmesh.function_space[spec_name] = dolfin.FunctionSpace(self.mesh,"Lagrange",1)
            # vertex_to_dof_map provides a map between the vertex index and the dof.
            xmesh.vertex_to_dof_map[spec_name]=xmesh.function_space[spec_name].dofmap().dof_to_vertex_map(self.mesh)
            xmesh.vertex_to_dof_map[spec_name]=len(self.listOfSpecies)*xmesh.vertex_to_dof_map[spec_name]+spec_index
            xmesh.dof_to_vertex_map[spec_name] = xmesh.function_space[spec_name].dofmap().vertex_to_dof_map(self.mesh)

        
        
        xmesh.vertex = self.mesh.coordinates()
        self.xmesh = xmesh
        
    
    # Some utility routines to set initial conditions
    def scatter(self,spec_init,subdomains=None):
        """ Scatter an initial number of molecules over the voxels in a subdomain. """
               
        if not hasattr(self,"u0"):
            self.initializeInitialValue()
        
        if not hasattr(self,'xmesh'):
            self.meshextend()

        self._initialize_species_to_subdomains()
        
        if not hasattr(self,'sd'):
            self.subdomainVector(self.subdomains)
                
        for species in spec_init:
            
            if subdomains is None:
                subdomains = self.species_to_subdomains[species]

            spec_name = species.name
            num_spec = spec_init[species]
            species_map = self.speciesMap()
            specindx= species_map[spec_name]

            sd = self.sd
            table = []
            for i,ind in enumerate(sd):
                if ind in subdomains:
                   table.append(i)

            ltab = len(table)
            if ltab == 0:
                raise ModelException("scatter: No voxel in the given subdomains "+str(subdomains)+", check subdomain marking.")
            
            for mol in range(num_spec):
                vtx=np.random.randint(0,ltab)
                ind = table[vtx]
                self.u0[specindx,ind]+=1
                

    def placeNear(self,spec_init, point=None):
        """ Place all molecules of kind species in the voxel nearest a given point. """

        if not hasattr(self,"u0"):
            self.initializeInitialValue()
                
        if not hasattr(self,'xmesh'):
            self.meshextend()

        coords = self.mesh.getVoxels()
        shape = coords.shape

        
        for spec in spec_init:
        
            spec_name = spec.name
            num_spec = spec_init[spec]
        
            
            # Find the voxel with center (vertex) nearest to the point
            p = dolfin.Point(point[0],point[1])
                    
            reppoint = numpy.tile(point,(shape[0],1))
            dist = numpy.sqrt(numpy.sum((coords-reppoint)**2,axis=1))
            ix = numpy.argmin(dist)
           
            species_map = self.speciesMap()
            specindx = species_map[spec_name]
            dofind = self.xmesh.vertex_to_dof_map[spec_name][ix]
            ix = (dofind-specindx)/len(species_map)
            self.u0[specindx,ix]=num_spec



    def createSystemMatrix(self):
        """
            Create the system (diffusion) matrix for input to the URDME solvers. The matrix
            is built by concatenating the individually assembled matrices for each of the species,
            and multiplying with the lumped mass matrix (which define the volume of the voxels).
            
            The dofs in the Dolfin-assembled matrices are reordered so that each column in the
            result matrix corresponds to the vertex numbering in the mesh.
            
            Negative off-diagonal elements in the matrix are set to zero, and the diagonal is renormalized
            in order to assure that the returned matrix is a Markov transition matrix.
            
            Returns a dictionary containing the volumes of the subvolumes, the system diffusion matrix
            and the fraction of the mass of the negative off-diagonal elements that has been filtered out.
            
            """
        
        # Check if the individual stiffness and mass matrices (per species) have been assembled, otherwise assemble them.
        try:
            stiffness_matrices = self.stiffness_matrices
            mass_matrices = self.mass_matrices
        except:
            matrices = assemble(self)
            self.stiffness_matrices = matrices['K']
            self.mass_matrices = matrices['M']
            stiffness_matrices = self.stiffness_matrices
            mass_matrices = self.mass_matrices
        
        # Make a dok matrix of dimension (Ndofs,Ndofs) for easier manipulatio
        
        i=1;
        Mspecies = len(self.listOfSpecies)
        Nvoxels = self.mesh.getNumVoxels()
        Ndofs = Nvoxels*Mspecies
        S = scipy.sparse.dok_matrix((Ndofs,Ndofs))
        
        # Create the volume vector by lumping the mass matrices
        vol = numpy.zeros((Ndofs,1))
        spec = 0
        
        xmesh = self.xmesh
        
        for species,M in mass_matrices.iteritems():
            
            dof2vtx = xmesh.dof_to_vertex_map[species]
            
            rows,cols,vals = M.data()
            SM = scipy.sparse.csr_matrix((vals,cols,rows))
            vols = SM.sum(axis=1)
            
            spec = self.species_map[species]
            for j in range(len(vols)):
                vx = dof2vtx[j]
                dof = Mspecies*vx+spec
                vol[dof,0]=vols[j]

        # This is necessary in order for the array to have the right dimension (Ndofs,1)
        vol = vol.flatten()
        
        # Assemble one big matrix from the indiviudal stiffness matrices. Multiply by the inverse of
        # the lumped mass matrix, filter out any entries with the wrong sign and renormalize the columns.
        spec = 0
        positive_mass = 0.0
        total_mass = 0.0
        
        for species,K in stiffness_matrices.iteritems():
            
            rows,cols,vals = K.data()
            Kcrs = scipy.sparse.csr_matrix((vals,cols,rows))
            #Kcsc = Kscr.tocsc()
            Kdok = Kcrs.todok()
            
            dof2vtx = xmesh.dof_to_vertex_map[species]
            
            
            for entries in Kdok.items():
                ind = entries[0]
                ir = dof2vtx[ind[0]]
                ij = dof2vtx[ind[1]]
                
                val = entries[1]
                
                if ir != ij:
                    if val > 0.0:
                        positive_mass += val
                        val = 0.0
                    else:
                        total_mass += val
                
                
                
                # The volume can be zero, if the species is not active at the vertex (such as a 2D species at a 3D node)
                if vol[Mspecies*ij+spec]==0:
                    vi = 1
                else:
                    vi = vol[Mspecies*ij+spec]
                
                S[Mspecies*ir+spec,Mspecies*ij+spec]=-val/vi
            
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
        
        #print "Fraction of positive off-diagonal entries: " + str(numpy.abs(positive_mass/total_mass))
        return {'vol':vol,'D':D,'relative_positive_mass':positive_mass/total_mass}

    def createSystemMatrix2(self):
        """
            Create the system (diffusion) matrix for input to the URDME solvers. The matrix
            is built by concatenating the individually assembled matrices for each of the species,
            and multiplying with the lumped mass matrix (which define the volume of the voxels).
            
            The dofs are reordered so that each column in the result matrix corresponds to 
            the vertex numbering in the mesh.
            
            Negative off-diagonal elements in the matrix are set to zero, and the diagonal is renormalized
            in order to assure that the returned matrix is a Markov transition matrix. 
            
            Returns a dictionary containing the volumes of the subvolumes, the system diffusion matrix
            and the fraction of the mass of the negative off-diagonal elements that has been filtered out.
            
        """
        
        # Check if the individual stiffness and mass matrices (per species) have been assembled, otherwise assemble them.
        try:
            stiffness_matrices = self.stiffness_matrices
            mass_matrices = self.mass_matrices
        except:
            matrices = assemble(self)
            self.stiffness_matrices = matrices['K']
            self.mass_matrices = matrices['M']
            stiffness_matrices = self.stiffness_matrices
            mass_matrices = self.mass_matrices
        
        # Make a dok matrix for easier manipulation
        i=1;
        Mspecies = len(self.listOfSpecies)
        Ndofs = self.mesh.getNumVoxels()*Mspecies
        S = scipy.sparse.dok_matrix((Ndofs,Ndofs))

        # Create the volume vector by lumping the mass matrices
        vol = numpy.zeros((Ndofs,1))
        spec = 0
        
        xmesh = self.xmesh
    
        for species,M in mass_matrices.iteritems():
            
            dof2vtx = xmesh.dof_to_vertex_map[species]
            
            rows,cols,vals = M.data()
            SM = scipy.sparse.csr_matrix((vals,cols,rows))
            vols = SM.sum(axis=1)
            spec = self.species_map[species]
            for j in range(len(vols)):
                dof = Mspecies*j+spec
                vol[dof,0]=vols[dof2vtx[j]]

        # This is necessary in order for the array to have the right dimension (Ndofs,1) 
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
            
            dof2vtx = xmesh.dof_to_vertex_map[species]


            for entries in Kdok.items():
                ind = entries[0]
                ir = dof2vtx[ind[0]]
                ij = dof2vtx[ind[1]]

                val = entries[1]
                
                if ir != ij:
                    if val > 0.0:
                        positive_mass += val
                        val = 0.0
                    else:
                        total_mass += val
            
            
            
                # The volume can be zero, if the species is not active at the vertex (such as a 2D species at a 3D node)
                if vol[Mspecies*ij+spec]==0:
                    vi = 1
                else:
                    vi = vol[Mspecies*ir+spec]
            
                S[Mspecies*ir+spec,Mspecies*ij+spec]=-val/vi
            
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

        #print "Fraction of positive off-diagonal entries: " + str(numpy.abs(positive_mass/total_mass))
        return {'vol':vol,'D':D,'relative_positive_mass':positive_mass/total_mass}

                
    def validate(self, urdme_solver_data):
        """ Validate the model data structures. 
            
            validate should be called prior to writing the model to the solver input file,
            since the solvers themselves do very limited error checking of the input.
        
        """
    
        # Check that all the columns of the system matrix sums to zero (or close to zero). If not, it does
        # not define a Markov process and the solvers might segfault or produce erraneous results.
        maxcolsum = numpy.max(numpy.abs(urdme_solver_data['D'].sum(axis=0)))
        if maxcolsum > 1e-10:
            raise InvalidSystemMatrixException("Invalid diffusion matrix. The sum of the columns does not sum to zero. " + str(maxcolsum))

    def solverData(self):
        """ Return the datastructures needed by the URDME solvers.
            
           solverData creates and populates a dictionary, urdme_solver_data, 
           containing the mandatory input data structures of the core NSM solver in URDME
           that is derived from the model. The data strucyures are
            
           N    - the stochiometry matrix
           G    - the dependency graph
           vol  - the volume vector
           sd   - the subdomain vector
           data - the data vector
           u0   - the intial condition
           
           This data is also returned, unlike in the Matlab URDME interface
            
           p - the vertex coordinates
           K - a (Nvoxel x Nvoxel) connectivity matrix
            
        """
        
        
        urdme_solver_data = {}
        
        # Stoichimetric matrix
        N = self.createStoichiometricMatrix()
        urdme_solver_data['N'] = N
        # Dependency Graph
        G = self.createDependencyGraph()
        urdme_solver_data['G']  = G
        
        # Volume vector
        result =  self.createSystemMatrix()
        vol = result['vol']
        urdme_solver_data['dofvolumes'] = vol
        
        #TODO: Make use of all dofs values, requires modification of CORE URDME...
        vol = vol[1::len(self.listOfSpecies)]
        
        urdme_solver_data['vol'] = vol
        D = result['D']
        urdme_solver_data['D'] = D
        
        # Subdomain vector
        urdme_solver_data['sd'] = self.subdomainVector(self.subdomains)
        
        # Data vector. If not present in model, it defaults to a vector with all elements zero.
        data = np.zeros((1,self.mesh.getNumVoxels()))
        urdme_solver_data['data'] = data
        
        if not hasattr(self,'u0'):
            self.initializeInitialValue()
        
        urdme_solver_data['u0'] = self.u0
        
        tspan= np.asarray(self.tspan,dtype=np.float)
        urdme_solver_data['tspan'] = tspan
        
        # Vertex coordinates
        urdme_solver_data['p'] = self.mesh.getVoxels()
        
        # Connectivity matrix
        urdme_solver_data['K'] = connectivityMatrix(self)
    
        return urdme_solver_data
    
    def serialize(self,filename=None):
        """ Write the datastructures needed by the the core URDME solvers to a .mat input file. """
        
        urdme_solver_data = self.solverData()
        self.validate(urdme_solver_data)
        spio.savemat(filename,urdme_solver_data,oned_as='column')


class Mesh(dolfin.Mesh):
    """ A URDME mesh extends the Dolfin mesh class. """

    def __init__(self,mesh=None):
        dolfin.Mesh.__init__(self,mesh)
    
    def getNumVoxels(self):
        return self.num_vertices()
    
    def getVoxels(self):
        return self.coordinates()

"""  Wrappers around dolfins built-in simple geometries/meshes.
    
    These following methods will all give regular meshes that will produce discretizations that are
    equivalent to Cartesian grids.

"""

def unitIntervalMesh(nx):
    mesh = dolfin.IntervalMesh(nx,0,1)
    return Mesh(mesh)

def IntervalMesh(nx,a,b):
    mesh = dolfin.IntervalMesh(nx,a,b)
    return Mesh(mesh)

def unitSquareMesh(nx,ny):
    """ Unit Square of with nx,ny points in the respective axes. """
    mesh = dolfin.UnitSquareMesh(nx,ny)
    return Mesh(mesh)

def SquareMesh(L,nx,ny):
    """ Regular mesh of a square with side length L. """
    mesh = dolfin.RectangleMesh(0,0,L,L,nx,ny)
    return Mesh(mesh)
    
def unitCubeMesh(nx,ny,nz):
    """ Unit Square of with nx,ny points in the respective axes. """
    mesh = dolfin.UnitCubeMesh(nx,ny,nz)
    return Mesh(mesh)

#def unitCircle(nx,ny):
#    """ Unit Square of with nx,ny points in the respective axes. """
#    mesh = dolfin.UnitCircleMesh(nx,ny)
#    return Mesh(mesh)

#def unitSphere(nx,ny):
#    """ Unit Square of with nx,ny points in the respective axes. """
#    mesh = dolfin.UnitSquareMesh(nx,ny)
#    return Mesh(mesh)



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
        mesh = Mesh(mesh=dolfin_mesh)
        return mesh
    except Exception,e:
        raise MeshImportError("Failed to import mesh: "+filename+"\n"+str(e))


def connectivityMatrix(model):
    """ Assemble a connectivity matrix in CCS format. """

    fs = dolfin.FunctionSpace(model.mesh,"Lagrange",1)
    trial_function = dolfin.TrialFunction(fs)
    test_function = dolfin.TestFunction(fs)
    a_K = -1*dolfin.inner(dolfin.nabla_grad(trial_function), dolfin.nabla_grad(test_function))*dolfin.dx
    C = dolfin.assemble(a_K)
    rows,cols,vals = C.data()
    C = scipy.sparse.csr_matrix((vals,cols,rows))
    C = C.tocsc()

    return C

def assemble(model):
    """  Assemble the mass and stiffness matrices using Dolfin.
    
         Returns: A dictionary containing two dictionaries, one for the stiffness matrices
                  and one for the mass matrices. Those dictionaries has the species names as keys and
                  the matrices are in CSR format.
    """
    
    if not hasattr(model,'xmesh'):
        model.meshextend()

    model._initialize_species_to_subdomains()
            
    function_space = model.xmesh.function_space
    trial_functions = OrderedDict()
    test_functions = OrderedDict()
    stiffness_matrices = OrderedDict()
    mass_matrices = OrderedDict()

    maxdim = 1
    for spec in model.listOfSpecies:
        dim = model.listOfSpecies[spec].dim()
        if dim > maxdim:
            maxdim = dim

    for spec in model.listOfSpecies:
        trial_functions[spec] = dolfin.TrialFunction(function_space[spec])
        test_functions[spec] = dolfin.TestFunction(function_space[spec])


    weak_form_K = {}
    weak_form_M = {}

    # Set up the weak forms
    for i,subdomain in enumerate(model.subdomains):
        
        # if species.dim() == maxdim:
        # sumbdomain dimension roughly corresponds to rdme_sdlevel
        if subdomain.dim()==maxdim:
            ddx = dolfin.Measure('dx')[subdomain]
        elif subdomain.dim() == maxdim-1:
            ddx = dolfin.Measure('dx')[subdomain]
        else:
            raise ModelException("Three subdomain levels is not supported.")

        for spec_name, species in model.listOfSpecies.items():
            
            spec_dim = species.dim()
            if species.dim() ==  subdomain.dim():
                
                # Find out what subdomains this species is active on
                subdomain_list = model.species_to_subdomains[species]
                
                # Set up the weak forms. We integrate only over those subdomains where the species is active
                for j,sd in enumerate(subdomain_list):
                    if  j==0:
                        weak_form_K[spec_name] = dolfin.inner(dolfin.nabla_grad(trial_functions[spec_name]), dolfin.nabla_grad(test_functions[spec_name]))*ddx(sd)
                        weak_form_M[spec_name] = trial_functions[spec_name]*test_functions[spec_name]*ddx(sd)
                    else:
                        weak_form_K[spec_name] = weak_form_K[spec_name]+dolfin.inner(dolfin.nabla_grad(trial_functions[spec_name]), dolfin.nabla_grad(test_functions[spec_name]))*ddx(sd)
                        weak_form_M[spec_name] = weak_form_M[spec_name]+trial_functions[spec_name]*test_functions[spec_name]*ddx(sd)

    # Assemble the matrices
    for spec,species in model.listOfSpecies.items():
        stiffness_matrices[spec] = dolfin.assemble(weak_form_K[spec])
        # We cannot include the diffusion constant in the assembly, dolfin does not seem to deal well
        # with small diffusion constants (drops small elements)
        stiffness_matrices[spec] = species.diffusion_constant*stiffness_matrices[spec]
        mass_matrices[spec] = dolfin.assemble(weak_form_M[spec])

    return {'K':stiffness_matrices,'M':mass_matrices}



class Xmesh():
    """ Extended mesh object.
        
        Contains function spaces and dof mappings.
    """

    def __init__(self):
        self.coordinates = None
        self.function_space = {}
        self.vertex_to_dof_map = {}
        self.dof_to_vertex_map = {}
       

def dumps(model,species,foldername):
    """ Dump the trajectory of species to a collection of vtk files """
    subprocess.call(["mkdir",foldername])
    func = dolfin.Function(dolfin.FunctionSpace(model.mesh,"Lagrange",1))
    func_vector = func.vector()
    file=dolfin.File(foldername+"/trajectory.pvd")
    numvox = model.mesh.getNumVoxels()

    for i,time in enumerate(model.tspan):
        #outfile = open(foldername+"/"+filename+"."+str(i),"w")
        solvector = (model.sol[species][time]).vector()
        for dof in range(numvox):
            func_vector[dof] = solvector[dof]
        file << func


def toXYZ(model,filename,format="ParaView"):
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

def toCSV(model,filename):
    """ Dump the solution attached to a model as a .csv file. """
    #TODO: Make this work for 2D meshes with only two coordinates.
    
    if 'U' not in model.__dict__:
        print "No solution found in the model."
        raise
    
    dims = numpy.shape(model.U)
    Ndofs = dims[0]
    Mspecies = len(model.listOfSpecies)
    Ncells = Ndofs/Mspecies
    
    coordinates = model.mesh.getVoxels()
    coordinatestr = coordinates.astype(str)
    subprocess.call(["mkdir",filename])
    for i,time in enumerate(model.tspan):
        outfile = open(filename+'/'+filename+str(i)+".csv","w")
        number_of_atoms = numpy.sum(model.U[:,i])
        filestr = "xcoord,ycoord,zcoord,radius,type\n"
        for j,spec in enumerate(model.listOfSpecies):
            for k in range(Ncells):
                for mol in range(model.U[k*Mspecies+j,i]):
                    obj = model.listOfSpecies[spec]
                    reaction_radius = obj.reaction_radius
                    linestr = coordinatestr[k,0]+","+coordinatestr[k,1]+","+coordinatestr[k,2]+","+str(reaction_radius)+","+str(j)+"\n";
                    filestr += linestr
        outfile.write(filestr)
        outfile.close()


def read_solution(filename):

    resultfile = h5py.File(filename,'r')

    U = resultfile['U']
    U = numpy.array(U)
    # This little hack makes U have the same structure as in the Matlab interface...
    dims = numpy.shape(U)
    U = U.reshape((dims[1],dims[0]))
    U = U.transpose()
    
    tspan = resultfile['tspan']
    tspan = numpy.array(tspan).flatten()
    resultfile.close()

    return {'U':U, 'tspan':tspan}


def urdme(model=None,solver='nsm',solver_path="", model_file=None, input_file=None, seed=None,report_level=0):
    """ URDME solver interface.
            
        TODO: Docs...
        
        After sucessful execution, urdme returns a dictionary, result, with the following members
            U:         the raw copy number output in a matrix with dimension (Ndofs, num_time_points)
            tspan:     the time span vector containing the time points that corresponds to the columns in U
            status:    Sucess if the solver executed without error
            stdout:    the standard ouput stream from the call to the core solver
            stderr:    the standard error stream from the call to the core solver
            
    """

    # Set URDME_ROOT. This will fail if URDME is not installed on the system.
    try:
        URDME_ROOT = subprocess.check_output(['urdme_init','-r'])
    except Exception,e:
        print "Could not determine the location of URDME."
        raise
    
    # Trim newline
    URDME_ROOT = URDME_ROOT[:-1]
    if solver_path == "":
        URDME_BUILD = URDME_ROOT+'/build/'
    else:
        URDME_BUILD = solver_path+'/build/'
        os.environ['SOLVER_ROOT'] = solver_path

    # Write the propensity file
    if os.path.isdir('.urdme'):
        shutil.rmtree('.urdme')

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
    handle.wait()

    if report_level >=1:
      print handle.stdout.read()
      print handle.stderr.read()

    if input_file is None:
        # Get temporary input and output files
        infile = tempfile.NamedTemporaryFile(delete=False)

        # Write the model to an input file in .mat format
        model.serialize(filename=infile)
        infile.close()
        infile_name = infile.name
    else:
        infile_name = input_file


    outfile = tempfile.NamedTemporaryFile(delete=False)
    outfile.close()

    # Execute the solver
    if seed is not None:
     try: 
      handle = subprocess.Popen(['.urdme/'+propfilename+'.'+solver,infile_name,outfile.name,str(seed)], stdout = subprocess.PIPE, stderr=subprocess.PIPE)
      handle.wait()
      if report_level >= 1:
        print handle.stdout.read()
        print handle.stderr.read()
     except:
        return {"status":"Failed","stderr":handle.stderr.read(),"stdout":handle.stdout.read()}
    else:
      try:
        handle = subprocess.Popen(['.urdme/'+propfilename+'.'+solver,infile_name,outfile.name], stdout = subprocess.PIPE, stderr=subprocess.PIPE)
        handle.wait()
        if report_level >= 1:
          print handle.stdout.read()
          print handle.stderr.read()
      except:
        return {"status":"Failed","stderr":handle.stderr.read(),"stdout":handle.stdout.read()}

    #Load the result from the hdf5 output file.
    try:
        
        result = read_solution(outfile.name)
        U = result['U']
        tspan = result['tspan']
        model.U = U
        
        # Create Dolfin Functions for all the species
        model.sol = {}
        
        dims = U.shape
        numvox = model.mesh.getNumVoxels()
        # The result is loaded in dolfin Functions, one for each species and time point
        for i,spec in enumerate(model.listOfSpecies):
    
            species = model.listOfSpecies[spec]
            spec_name = species.name
            dof_to_vertex_map = model.xmesh.dof_to_vertex_map[spec]
            vertex_to_dof_map = model.xmesh.vertex_to_dof_map[spec]

            spec_sol = {}
            for j,time in enumerate(tspan):
                func = dolfin.Function(dolfin.FunctionSpace(model.mesh,"Lagrange",1))
                func_vector = func.vector()
            
                for voxel in range(numvox):
                    dof = voxel*len(model.listOfSpecies)+i
                    ix  = vertex_to_dof_map[voxel]
                    dolfvox = (ix-i)/len(model.listOfSpecies)
                    func_vector[dolfvox] = float(U[dof,j])
                
                spec_sol[time] = func
            
            model.sol[spec] = spec_sol


        # Clean up
        if input_file is None:
            os.remove(infile.name)
        os.remove(outfile.name)

        return dict({"Status":"Sucess","stdout":handle.stdout.read(),"stderr":handle.stderr.read()},**result)

    except Exception,e:
       # Clean up
       if input_file is None:
           subprocess.call(['rm','-rf',infile.name])
       subprocess.call(['rm','-rf',outfile.name])
       raise


class URDMEError(Exception):
    pass

if __name__ == '__main__':
    """ Command line interface to URDME. Execute URDME given a model file. """ 

class ModelException(Exception):
    pass


class InvalidSystemMatrixException(Exception):
    pass



