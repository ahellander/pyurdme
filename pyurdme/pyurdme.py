# pylint: disable-msg=C0301
# pylint: disable-msg=C0103

import os
import re
import shutil
import subprocess
import sys
import tempfile
import types
import warnings
import copy

import numpy
import scipy.io
import scipy.sparse

import gmsh
from model import *

import inspect

import IPython.display

try:
    import h5py
except:
    raise Exception("PyURDME requires h5py.")

try:
    import dolfin
    dolfin.parameters["linear_algebra_backend"] = "uBLAS"
except:
    raise Exception("PyURDME requires FeniCS/dolfin.")

import pickle
import json

# Set log level to report only errors or worse
dolfin.set_log_level(dolfin.ERROR)

class URDMEModel(Model):
    """
        An URDME Model extends Model with spatial information and methods to create URDME solver input.
        TODO: Documentiation.
    """

    def __init__(self, name=""):
        Model.__init__(self, name)

        # Currently not used
        self.geometry = None

        self.mesh = None
        self.xmesh = None
        self.stiffness_matrices = None
        self.mass_matrices = None
        
        # subdomins is a list of MeshFunctions with subdomain marker information
        self.subdomains = OrderedDict()

        # This dictionary hold information about the subdomains each species is active on
        self.species_to_subdomains = {}
        self.tspan = None

        # URDMEDataFunction objects to construct the data vector.
        self.listOfDataFunctions = []

        # Volume of each voxel in the dolfin dof ordering (not vetex ordering).
        self.dofvol = None

    def __getstate__(self):
        """ Used by pickle to get state when pickling. Because we
            have Swig wrappers to extension modules, we need to remove some instance variables
            for the object to pickle. """

        #  Filter out any instance variable that is not picklable...
        state = {}
        for key, item in self.__dict__.items():
            try:
                pickle.dumps(item)
                state[key] = item
            except Exception as e:
                if key == "mesh":
                    tmpfile = tempfile.NamedTemporaryFile(suffix=".xml")
                    dolfin.File(tmpfile.name) << item
                    tmpfile.seek(0)
                    state['mesh'] = {}
                    state['mesh']['data'] = tmpfile.read()
                    tmpfile.close()
                    if item.constrained_domain is not None:
                        # Warning: This is black magic.
                        try:
                            cdd = {}
                            cdd['source'] = inspect.getsource(item.constrained_domain.__class__)
                            cdd['name'] = item.constrained_domain.__class__.__name__
                            cdd['dict'] = {}
                            for k,v in item.constrained_domain.__dict__.iteritems():
                                if type(v).__name__ != 'SwigPyObject':
                                    cdd['dict'][k] = v
                            state['mesh']['constrained_domain'] = cdd
                        except Exception as e:
                            sys.stderr.write("error pickling mesh.constrained_domain: {0}\n".format(e))
                            raise e
                    if item.num_dof_voxels is not None:
                        state['mesh']['num_dof_voxels'] = item.num_dof_voxels
                elif key == "subdomains":
                    sddict = OrderedDict()
                    for sdkey, sd_func in item.items():
                        tmpfile = tempfile.NamedTemporaryFile(suffix=".xml")
                        dolfin.File(tmpfile.name) << sd_func
                        tmpfile.seek(0)
                        sddict[sdkey] = tmpfile.read()
                        tmpfile.close()
                    state[key] = sddict
                else:
                    state[key] = None


        return state

    def __setstate__(self, state):
        """ Used by pickle to set state when unpickling. """

        self.__dict__ = state

        if 'mesh' in state:
            # Recreate the mesh
            try:
                fd = tempfile.NamedTemporaryFile(suffix=".xml")
                fdname = fd.name
                fd.write(state['mesh']['data'])
                fd.seek(0)
                mesh = URDMEMesh.read_dolfin_mesh(fdname)
                fd.close()
                if 'constrained_domain' in state['mesh']:
                    # Black magic to match that in __getstate__
                    cdd = state['mesh']['constrained_domain']
                    compiled_class = compile(cdd['source'], 'pyurdme.mesh.constrained_domain', 'exec')
                    eval(compiled_class)
                    compiled_object = eval("{0}()".format(cdd['name']))
                    for k,v in cdd['dict'].iteritems():
                        compiled_object.__dict__[k] = v
                    mesh.constrained_domain = compiled_object
                if 'num_dof_voxels' in state['mesh']:
                    mesh.num_dof_voxels = state['mesh']['num_dof_voxels']
                self.__dict__['mesh'] = mesh
            except Exception as e:
                print "Error unpickling model, could not recreate the mesh."
                raise e

        if 'subdomains' in state:
            # Recreate the subdomain functions
            try:
                sddict = OrderedDict()
                for sdkey, sd_func_str in state["subdomains"].items():
                    fd = tempfile.NamedTemporaryFile(suffix=".xml")
                    fdname = fd.name
                    fd.write(sd_func_str)
                    fd.seek(0)
                    fd_in = dolfin.File(fdname)
                    func = dolfin.MeshFunction("size_t", self.__dict__["mesh"])
                    fd_in >> func
                    sddict[sdkey] = func
                    fd.close()
                self.__dict__["subdomains"] = sddict
            except Exception as e:
                raise Exception("Error unpickling model, could not recreate the subdomain functions"+str(e))
    
        self.meshextend()


    def addDataFunction(self, data_function):
        """ Add a URDMEDataFunction object to this object. """
        if isinstance(data_function, URDMEDataFunction):
            self.listOfDataFunctions.append(data_function)
        else:
            raise Exception("data_function not of type URDMEDataFunction")

    def __initializeSpeciesMap(self):
        i = 0
        self.species_map = {}
        for S in self.listOfSpecies:
            self.species_map[S] = i
            i = i + 1

    def speciesMap(self):
        """ Get the species map, name to index. """
        if not hasattr(self, 'species_map'):
            self.__initializeSpeciesMap()

        return self.species_map

    def addSubDomain(self, subdomain):
        if not subdomain.dim() in self.subdomains.keys():
            self.subdomains[subdomain.dim()] = subdomain
        else:
            raise ModelException("Failed to add subdomain function of dim "+str(subdomain.dim())+". Only one subdomain function of a given dimension is allowed.")

    def createStoichiometricMatrix(self):
        """ Generate a stoichiometric matrix in sparse CSC format. """

        if not hasattr(self, 'species_map'):
            self.__initializeSpeciesMap()
        if self.getNumReactions() > 0:
            ND = numpy.zeros((self.getNumSpecies(), self.getNumReactions()))
            for i, r in enumerate(self.listOfReactions):
                R = self.listOfReactions[r]
                reactants = R.reactants
                products  = R.products

                for s in reactants:
                    ND[self.species_map[s], i] -= reactants[s]
                for s in products:
                    ND[self.species_map[s], i] += products[s]

            N = scipy.sparse.csc_matrix(ND)
        else:
            N = numpy.zeros((self.getNumSpecies(), self.getNumReactions()))

        return N

    def createDependencyGraph(self):
        """ Construct the sparse dependency graph. """
        # We cannot safely generate a dependency graph (without attempting to analyze the propensity string itself)
        # if the model contains custom propensities.
        mass_action_model = True
        for name,reaction in self.listOfReactions.items():
            if not reaction.massaction:
                GF = numpy.ones((self.getNumReactions(), self.getNumReactions() + self.getNumSpecies()))
                mass_action_model = False
    
        if mass_action_model:
            GF = numpy.zeros((self.getNumReactions(), self.getNumReactions() + self.getNumSpecies()))
            species_map = self.speciesMap()
            
            involved_species = []
            reactants = []
            for name, reaction in self.listOfReactions.items():
                temp = []
                temp2 = []
                for s in reaction.reactants:
                    temp.append(species_map[s])
                    temp2.append(species_map[s])
                for s in reaction.products:
                    temp.append(species_map[s])
                involved_species.append(temp)
                reactants.append(temp2)
                    
            species_to_reactions = []
            for species in self.listOfSpecies:
                temp = []
                for j,x in enumerate(reactants):
                    if species_map[species] in x:
                        temp.append(j)
                species_to_reactions.append(temp)
            

            reaction_to_reaction = []
            for name, reaction in self.listOfReactions.items():
                temp = []
                for s in reaction.reactants:
                    if species_to_reactions[species_map[s]] not in temp:
                        temp = temp+species_to_reactions[species_map[s]]
                
                for s in reaction.products:
                    if species_to_reactions[species_map[s]] not in temp:
                        temp = temp+ species_to_reactions[species_map[s]]
                
                temp = list(set(temp))
                reaction_to_reaction.append(temp)
            
            # Populate G
            for j, spec in enumerate(species_to_reactions):
                for s in spec:
                    GF[s,j] = 1
            
            for i,reac in enumerate(reaction_to_reaction):
                for r in reac:
                    GF[r,self.getNumSpecies()+i] = 1

                
        try:
            G = scipy.sparse.csc_matrix(GF)
        except Exception as e:
            G = GF

        return G


    def timespan(self, tspan):
        """ Set the time span of simulation. """
        self.tspan = tspan


    def _initialize_default_subdomain(self):
        """" Create a default subdomain function. The default is all voxels belong
             to subdomain 1.
        """

        subdomain = dolfin.MeshFunction("size_t", self.mesh, self.mesh.topology().dim()-1)
        subdomain.set_all(1)
        self.addSubDomain(subdomain)

    def _initialize_species_to_subdomains(self):
        """ Initialize the species mapping to subdomains. The default
            is that a species is active in all the defined subdomains.
        """

        # If no subdomain function has been set by the user,
        # we need to create a default subdomain here.
        if not self.subdomains:
            self._initialize_default_subdomain()

        # The unique elements of the subdomain MeshFunctions
        sds = []
        for dim, subdomain in self.subdomains.items():
            sds = sds + list(numpy.unique(subdomain.array()).flatten())
        sds = numpy.unique(sds)
        sds = list(sds)

        # This explicit typecast is necessary for UFL not to choke on the subdomain ids.
        for i, sd in enumerate(sds):
            sds[i] = int(sd)

        # This conversion is necessary for UFL not to choke on the subdomain ids.
        for i, sd in enumerate(sds):
            sds[i] = int(sd)
        try:
            sds.remove(0)
        except Exception:
            pass

        # If a species is not present as key in the species_to_subdomain mapping,
        # we label it as active in all subdomains
        for spec_name in self.listOfSpecies:
            species = self.listOfSpecies[spec_name]
            if species not in self.species_to_subdomains.keys():
                self.species_to_subdomains[species] = sds


    def restrict(self, species, subdomains):
        self.species_to_subdomains[species] = subdomains


    def subdomainVector(self, subdomains={}):
        """ Create the 'sd' vector. 'subdomains' is a dolfin FacetFunction,
            and if no subdomain input is specified, they voxels default to
            subdomain 1. """


        # TODO: We need to make sure that the highest dimension is applied
        #       first, otherwise the cell level will overwrite all markings
        #       applied on boundaries.

        if not hasattr(self,'xmesh'):
            self.meshextend()

        self.mesh.init()

        # TODO: Support arbitrary sd-numbers and more than one subdomain
        sd = numpy.zeros((1, self.mesh.getNumVoxels()))
        if subdomains == {}:
            self.sd = sd.flatten()
            print subdomains
        else:
            for dim, subdomain in subdomains.items():
                # Map all facet labels to vertex labels
                tovertex = self.mesh.topology()(dim, 0)
                for i in range(subdomain.size()):
                    for vtx in tovertex(i):
                        if subdomain[i] != 0: # TODO: Temporary hack to fix issue with Gmesh facet_region files.
                            sd[0, vtx] = subdomain[i]

        self.sd = sd.flatten()
        return self.sd


    def initializeInitialValue(self):
        """ Create all-zeros inital condition matrix. """
        ns = self.getNumSpecies()
        if self.xmesh == None:
            self.meshextend()
        nv = self.mesh.getNumVoxels()
        #dims = numpy.shape(self.mesh.FunctionSpace().dofmap().dof_to_vertex_map(self.mesh))
        #nv = dims[0]
        #print nv

        self.u0 = numpy.zeros((ns, nv))

    def meshextend(self):
        """ Extend the primary mesh with information about the degrees of freedom.

            TODO: Docs...

            """

        xmesh = URDMEXmesh()

        # Construct a species map (dict mapping model species name to an integer index)
        species_map = self.speciesMap()

        # Initialize the function spaces and dof maps.
        for spec in self.listOfSpecies:
            
            species = self.listOfSpecies[spec]
            spec_name = species.name
            spec_index = species_map[spec_name]

            xmesh.function_space[spec_name] = self.mesh.FunctionSpace()
            
            xmesh.vertex_to_dof_map[spec_name] = dolfin.vertex_to_dof_map(xmesh.function_space[spec_name])
            xmesh.vertex_to_dof_map[spec_name] = len(self.listOfSpecies) * xmesh.vertex_to_dof_map[spec_name] + spec_index
            xmesh.dof_to_vertex_map[spec_name] = dolfin.dof_to_vertex_map(xmesh.function_space[spec_name])
    

        xmesh.vertex = self.mesh.coordinates()
        self.xmesh = xmesh


    # Some utility routines to set initial conditions
    def scatter(self, spec_init, subdomains=None):
        """ Scatter an initial number of molecules over the voxels in a subdomain. """

        if not hasattr(self,"u0"):
            self.initializeInitialValue()

        if not hasattr(self, 'xmesh'):
            self.meshextend()

        self._initialize_species_to_subdomains()

        if not hasattr(self, 'sd'):
            self.subdomainVector(self.subdomains)

        for species in spec_init:

            if subdomains is None:
                subdomains = self.species_to_subdomains[species]

            spec_name = species.name
            num_spec = spec_init[species]
            species_map = self.speciesMap()
            specindx = species_map[spec_name]

            sd = self.sd
            table = []
            for i, ind in enumerate(sd):
                if ind in subdomains:
                    table.append(i)

            ltab = len(table)
            if ltab == 0:
                raise ModelException("scatter: No voxel in the given subdomains "+str(subdomains)+", check subdomain marking.")

            for mol in range(num_spec):
                vtx = numpy.random.randint(0, ltab)
                ind = table[vtx]
                self.u0[specindx, ind] += 1

    def distributeUniformly(self, spec_init):
        """ Place the same number of molecules of the species in each voxel. """
        if not hasattr(self, "u0"):
            self.initializeInitialValue()

        if not hasattr(self, 'xmesh'):
            self.meshextend()

        species_map = self.speciesMap()
        num_voxels = self.mesh.getNumVoxels()
        for spec in spec_init:
            spec_name = spec.name
            num_spec = spec_init[spec]
            specindx = species_map[spec_name]
            for ndx in range(num_voxels):
                self.u0[specindx, ndx] = num_spec
    
    
    def placeNear(self, spec_init, point=None):
        """ Place all molecules of kind species in the voxel nearest a given point. """

        if not hasattr(self, "u0"):
            self.initializeInitialValue()

        if not hasattr(self, 'xmesh'):
            self.meshextend()

        coords = self.mesh.getVoxels()
        shape = coords.shape


        for spec in spec_init:
            spec_name = spec.name
            num_spec = spec_init[spec]

            # Find the voxel with center (vertex) nearest to the point
            reppoint = numpy.tile(point, (shape[0], 1))
            dist = numpy.sqrt(numpy.sum((coords-reppoint)**2, axis=1))
            ix = numpy.argmin(dist)

            species_map = self.speciesMap()
            specindx = species_map[spec_name]
            #dofind = self.xmesh.vertex_to_dof_map[spec_name][ix]
            #ix = (dofind - specindx) / len(species_map)
            self.u0[specindx, ix] = num_spec



    def createSystemMatrix(self):
        """
            Create the system (diffusion) matrix for input to the URDME solvers. The matrix
            is built by concatenating the individually assembled matrices for each of the species,
            and multiplying with the lumped mass matrix (which define the volume of the voxels).

            Negative off-diagonal elements in the matrix are set to zero, and the diagonal is renormalized
            in order to assure that the returned matrix is a Markov transition matrix.

            Returns a dictionary containing the volumes of the subvolumes, the system diffusion matrix
            and the fraction of the mass of the negative off-diagonal elements that has been filtered out.

            """

        # Check if the individual stiffness and mass matrices (per species) have been assembled, otherwise assemble them.
        if self.stiffness_matrices is not None and self.mass_matrices is not None:
            stiffness_matrices = self.stiffness_matrices
            mass_matrices = self.mass_matrices
        else:
            if self.mesh is None:
                raise ModelException("This model has no mesh, can not create system matrix.")
            matrices = self.assemble()
            self.stiffness_matrices = matrices['K']
            self.mass_matrices = matrices['M']
            stiffness_matrices = self.stiffness_matrices
            mass_matrices = self.mass_matrices

        # Make a dok matrix of dimension (Ndofs,Ndofs) for easier manipulatio

        i = 1
        Mspecies = len(self.listOfSpecies)
        if Mspecies == 0:
            raise ModelException("The model has no species, can not create system matrix.")
        # Use dolfin 'dof' number of voxels, not the number of verticies
        Nvoxels = self.mesh.getNumDofVoxels()
        Ndofs = Nvoxels*Mspecies
        S = scipy.sparse.dok_matrix((Ndofs, Ndofs))

        # Create the volume vector by lumping the mass matrices
        vol = numpy.zeros((Ndofs, 1))
        spec = 0

        xmesh = self.xmesh

        for species, M in mass_matrices.iteritems():

            #dof2vtx = xmesh.dof_to_vertex_map[species]
            rows, cols, vals = M.data()
            SM = scipy.sparse.csr_matrix((vals, cols, rows))
            vols = SM.sum(axis=1)
            
            spec = self.species_map[species]
            for j in range(len(vols)):
                #vx = dof2vtx[j]  # need to use dof ordering
                vx = j
                dof = Mspecies*vx+spec
                vol[dof, 0] = vols[j]
            
        # This is necessary in order for the array to have the right dimension (Ndofs,1)
        vol = vol.flatten()

        # Assemble one big matrix from the indiviudal stiffness matrices. Multiply by the inverse of
        # the lumped mass matrix, filter out any entries with the wrong sign and renormalize the columns.
        spec = 0
        positive_mass = 0.0
        total_mass = 0.0

#try:
#           sd = self.sd
#        except:
        sd = self.subdomainVector(self.subdomains)
        sd_vec_dof = numpy.zeros(self.mesh.getNumDofVoxels())
        vertex_to_dof = dolfin.vertex_to_dof_map(self.mesh.FunctionSpace())
        for ndx, sd_val in enumerate(sd):
            sd_vec_dof[vertex_to_dof[ndx]] = sd_val
        sd = sd_vec_dof

        for species, K in stiffness_matrices.iteritems():

            rows, cols, vals = K.data()
            Kcrs = scipy.sparse.csr_matrix((vals, cols, rows))
            Kdok = Kcrs.todok()

            #dof2vtx = xmesh.dof_to_vertex_map[species]

            for entries in Kdok.items():

                ind = entries[0]
                ir = ind[0]
                ij = ind[1]

                # Use Dolfin dof ordering
                # Depricated: Permutation to make the matrix ordering match that of sd, u0. (Dolfin dof -> URDME dof)
                #ir = dof2vtx[ind[0]]
                #ij = dof2vtx[ind[1]]

                val = entries[1]

                if ir != ij:

                    # Check if this is an edge that the species should diffuse along,
                    # if not, set the diffusion coefficient along this edge to zero. This is
                    # equivalent to how boundary species are handled in the current Matlab interface.
                    if sd[ir] not in self.species_to_subdomains[self.listOfSpecies[species]]:
                        val = 0.0

                    if val > 0.0:
                        positive_mass += val
                        val = 0.0
                    else:
                        total_mass += val

                # The volume can be zero, if the species is not active at the vertex (such as a 2D species at a 3D node)
                if vol[Mspecies*ij+spec] == 0:
                    vi = 1
                else:
                    vi = vol[Mspecies*ij+spec]
                
                S[Mspecies*ir+spec, Mspecies*ij+spec] = -val/vi

            spec = spec + 1

        # Convert to compressed column for compatibility with the URDME solvers.
        D = S.tocsc()

        # Renormalize the columns (may not sum to zero since elements may have been filtered out
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sumcol = numpy.zeros((Ndofs, 1))
            for i in range(Ndofs):
                col = D.getcol(i)
                for val in col.data:
                    if val > 0.0:
                        sumcol[i] += val

            D.setdiag(-sumcol.flatten())        

        return {'vol':vol, 'D':D, 'relative_positive_mass':positive_mass/total_mass}


    def validate(self, urdme_solver_data):
        """ Validate the model data structures.

            validate should be called prior to writing the model to the solver input file,
            since the solvers themselves do very limited error checking of the input.

        """

        for spec_name, species in self.listOfSpecies.items():
            if 0 in self.species_to_subdomains[species]:
                raise ModelException("Subdomain number 0 is reserved. Please check your model.")

        # Check that all the columns of the system matrix sums to zero (or close to zero). If not, it does
        # not define a Markov process and the solvers might segfault or produce erraneous results.
        colsum = numpy.abs(urdme_solver_data['D'].sum(axis=0))
        colsum = colsum.flatten()
        maxcolsum = numpy.argmax(colsum)
        if colsum[0,maxcolsum] > 1e-10:
            D = urdme_solver_data["D"]
            raise InvalidSystemMatrixException("Invalid diffusion matrix. The sum of the columns does not sum to zero. " + str(maxcolsum) + str(colsum[0,maxcolsum]))


    def connectivityMatrix(self):
        """ Assemble a connectivity matrix in CCS format. """
        
        fs = self.mesh.FunctionSpace()
        trial_function = dolfin.TrialFunction(fs)
        test_function = dolfin.TestFunction(fs)
        a_K = -1*dolfin.inner(dolfin.nabla_grad(trial_function), dolfin.nabla_grad(test_function)) * dolfin.dx
        K = dolfin.assemble(a_K)
        rows, cols, vals = K.data()
        Kcrs = scipy.sparse.csc_matrix((vals, cols, rows))
        return Kcrs
        #
        ## Permutation dolfin dof -> URDME dof
        #Kdok = Kcrs.todok()
        #
        #nv  = self.mesh.num_vertices()
        #S = scipy.sparse.dok_matrix((nv,nv))
        #
        #dof2vtx = dolfin.vertex_to_dof_map(fs)
        #
        #for entries in Kdok.items():
        #
        #    ind = entries[0]
        #    ir = ind[0]
        #    ij = ind[1]
        #
        #    # Permutation to make the matrix ordering match that of sd, u0. (Dolfin dof -> URDME dof)
        #    ir = dof2vtx[ind[0]]
        #    ij = dof2vtx[ind[1]]
        #
        #    val = entries[1]
        #
        #
        #    S[ir,ij] = val
        #
        #
        #
        #C = S.tocsc()
        #return C


    def solverData(self):
        """ Return the datastructures needed by the URDME solvers.

           solverData creates and populates a dictionary, urdme_solver_data,
           containing the mandatory input data structures of the core NSM solver in URDME
           that is derived from the model. The data strucyures are

           D    - the Diffusion matrix
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
        
        num_species = self.getNumSpecies()

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
        self.dofvol = vol[::len(self.listOfSpecies)]
        urdme_solver_data['vol'] = self.dofvol
        
        D = result['D']
        urdme_solver_data['D'] = D
        
        #
        num_dofvox = self.dofvol.shape[0]

        # Get vertex to dof ordering
        vertex_to_dof = dolfin.vertex_to_dof_map(self.mesh.FunctionSpace())
        dof_to_vertex = dolfin.dof_to_vertex_map(self.mesh.FunctionSpace())
        
        vertex_to_dof_to_vertex = dof_to_vertex[vertex_to_dof]
        
        # Subdomain vector
        # convert to dof ordering
        sd_vec_dof = numpy.zeros(num_dofvox)
        for ndx, sd_val in enumerate(self.subdomainVector(self.subdomains)):
            sd_vec_dof[vertex_to_dof[ndx]] = sd_val
        urdme_solver_data['sd'] = sd_vec_dof
        
        # Data vector. If not present in model, it defaults to a vector with all elements zero.
        # convert to dof ordering
        data = numpy.zeros((1, num_dofvox))
        if len(self.listOfDataFunctions) > 0:
            data = numpy.zeros((len(self.listOfDataFunctions), num_dofvox))
            coords = self.mesh.coordinates()
            for ndf, df in enumerate(self.listOfDataFunctions):
                for ndx in range(len(coords)):
                    vox_coords = numpy.zeros(3)
                    for cndx in range(len(coords[ndx])):
                        vox_coords[cndx] = coords[ndx][cndx]
                    data[ndf][vertex_to_dof[ndx]] = df.map(vox_coords)
            
        urdme_solver_data['data'] = data

        if not hasattr(self,'u0'):
            self.initializeInitialValue()

        # Initial Conditions, convert to dof ordering
        u0_dof = numpy.zeros((num_species, num_dofvox))
        for vox_ndx in range(self.mesh.getNumVoxels()):
            dof_ndx = vertex_to_dof[vox_ndx]
            # With periodic BCs the same dof_ndx voxel will get written to twice
            # which may overwrite the value.  We need to check for this case.
            if vertex_to_dof_to_vertex[vox_ndx] != vox_ndx:
                vox_ndx2 = vertex_to_dof_to_vertex[vox_ndx]
                for cndx in range(num_species):
                    if self.u0[cndx, vox_ndx] == 0 or self.u0[cndx, vox_ndx] == self.u0[cndx, vox_ndx2]:
                        u0_dof[cndx, dof_ndx] = self.u0[cndx, vox_ndx]
                    elif self.u0[cndx, vox_ndx2] == 0 and vox_ndx < vox_ndx2:
                        self.u0[cndx, vox_ndx2] = self.u0[cndx, vox_ndx]
                    elif self.u0[cndx, vox_ndx2] == 0 and vox_ndx > vox_ndx2:
                        u0_dof[cndx, dof_ndx] = self.u0[cndx, vox_ndx]
                    else:
                        sys.stderr.write("Warning: the initial condition for species {0} in voxel {1} will be discarded due to periodic boundary conditions.\n".format(self.listOfSpecies.keys()[cndx], vox_ndx))
            else:
                for cndx in range(num_species):
                    u0_dof[cndx, dof_ndx] = self.u0[cndx, vox_ndx]
        urdme_solver_data['u0'] = u0_dof

        tspan = numpy.asarray(self.tspan, dtype=numpy.float)
        urdme_solver_data['tspan'] = tspan

        # Vertex coordinates
        # convert to dof ordering
        p_dof = numpy.zeros((num_dofvox, 3))
        for vox_ndx, row in enumerate(self.mesh.getVoxels()):
            p_dof[vertex_to_dof[vox_ndx],:len(row)] = row
        urdme_solver_data['p'] = p_dof

        # Connectivity matrix
        urdme_solver_data['K'] = self.connectivityMatrix()

        urdme_solver_data['report']=0

        return urdme_solver_data


    def serialize(self, filename=None, report_level=0):
        """ Write the datastructures needed by the the core URDME solvers to a .mat input file. """
        urdme_solver_data = self.solverData()
        urdme_solver_data['report'] = report_level
        self.validate(urdme_solver_data)
        scipy.io.savemat(filename, urdme_solver_data, oned_as='column')


    def assemble(self):
        """  Assemble the mass and stiffness matrices using Dolfin.

            Returns: A dictionary containing two dictionaries, one for the stiffness matrices
            and one for the mass matrices. Those dictionaries has the species names as keys and
            the matrices are in CSR format.
            """

        if self.xmesh == None:
            self.meshextend()

        self._initialize_species_to_subdomains()

        function_space = self.xmesh.function_space
        trial_functions = OrderedDict()
        test_functions = OrderedDict()
        stiffness_matrices = OrderedDict()
        mass_matrices = OrderedDict()

        # The maximum dimension that a species is active on (currently not used)
        maxdim = 1
        for spec in self.listOfSpecies:
            dim = self.listOfSpecies[spec].dim()
            if dim > maxdim:
                maxdim = dim

        for spec in self.listOfSpecies:
            trial_functions[spec] = dolfin.TrialFunction(function_space[spec])
            test_functions[spec] = dolfin.TestFunction(function_space[spec])


        weak_form_K = {}
        weak_form_M = {}
        
        ndofs = None

        # Set up the forms
        for spec_name, species in self.listOfSpecies.items():

            # Find out what subdomains this species is active on
            #subdomain_list = self.species_to_subdomains[species]
            weak_form_K[spec_name] = dolfin.inner(dolfin.nabla_grad(trial_functions[spec_name]), dolfin.nabla_grad(test_functions[spec_name]))*dolfin.dx
            weak_form_M[spec_name] = trial_functions[spec_name]*test_functions[spec_name]*dolfin.dx

        # Assemble the matrices
        for spec_name, species in self.listOfSpecies.items():
            stiffness_matrices[spec_name] = dolfin.assemble(weak_form_K[spec_name])
            if ndofs is None:
                ndofs = stiffness_matrices[spec_name].size(0)
                self.mesh.setNumDofVoxels(ndofs)
            
            # We cannot include the diffusion constant in the assembly, dolfin does not seem to deal well
            # with small diffusion constants (drops small elements)
            stiffness_matrices[spec_name] = species.diffusion_constant * stiffness_matrices[spec_name]
            mass_matrices[spec_name] = dolfin.assemble(weak_form_M[spec_name])


        return {'K':stiffness_matrices, 'M':mass_matrices}

    def run(self, solver='nsm', seed=None, report_level=0):
        """ Simulate the model.
        
        Args:
            solver: A str or class type that is a subclass of URDMESolver.  Default: NSM solver.
            seed: An int, the random seed given to the solver.
            report_level: An int, Level of output from the solver: 0, 1, or 2. Default: 0.
        Returns:
            A URDMEResult object with the results of the simulation.
        """
        
        #If solver is a subclass of URDMESolver, use it directly.
        if isinstance(solver, (type, types.ClassType)) and  issubclass(solver, URDMESolver):
            sol = solver(self, report_level=report_level)
        elif type(solver) is str:
            if solver == 'nsm':
                from nsmsolver import NSMSolver
                sol = NSMSolver(self, report_level=report_level)
            else:
                raise URDMEError("Unknown solver: {0}".format(solver_name))
        else:
            raise URDMEError("solver argument to urdme() must be a string or a URDMESolver class object.")

        return sol.run(seed)



class URDMEMesh(dolfin.Mesh):
    """ A URDME mesh extends the Dolfin mesh class. """

    def __init__(self, mesh=None):
        self.constrained_domain = None
        dolfin.Mesh.__init__(self, mesh)
        self.function_space = None
        self.num_dof_voxels = None


    def addPeriodicBoundaryCondition(self, domain):
        self.constrained_domain = domain

    def FunctionSpace(self):
        if self.function_space is not None:
            return self.function_space
        else:
            if self.constrained_domain is not None:
                fs = dolfin.FunctionSpace(self, "Lagrange", 1, constrained_domain=self.constrained_domain)
            else:
                fs = dolfin.FunctionSpace(self, "Lagrange", 1)
            self.function_space = fs
            return fs

    def getNumVoxels(self):
        return self.num_vertices()

    def setNumDofVoxels(self, num):
        self.num_dof_voxels = num

    def getNumDofVoxels(self):
        if self.num_dof_voxels is None:
            raise URDMEError('NumDofVoxels is not set')
        return self.num_dof_voxels

    def getVoxels(self):
        return self.coordinates()

    """  Wrappers around dolfins built-in simple geometries/meshes.

        These following methods will all give regular meshes that will produce discretizations that are
        equivalent to Cartesian grids.

    """
    
    def meshSize(self):
        """ Estimate of mesh size at each vertex. """
        coordinates = self.coordinates()
        
        # Compute the circumradius of the cells
        cr = []
        for i in range(self.num_cells()):
            cell = dolfin.Cell(self, i)
            cr.append(cell.diameter()/2.0)

        # Compute the mean for each vertex based on all incident cells
        vtx2cell = self.topology()(0,self.topology().dim())
        vtxh = []
        for i in range(self.num_vertices()):
            v2c = vtx2cell(i)
            h = 0.0
            for indx in v2c:
                h += cr[indx]
            h = h/len(v2c)
            vtxh.append(h)

        return vtxh


    def normalizedCoordinates(self):
        """ Return vertex coordinates centered at origo. """
        
        # Compute mesh centroid
        vtx = self.coordinates()
        centroid = numpy.mean(vtx,axis=0)
        # Shift so the centroid is now origo
        normalized_vtx = numpy.zeros(numpy.shape(vtx))
        for i,v in enumerate(vtx):
            normalized_vtx[i,:] = v - centroid

        return normalized_vtx

    def scaledNormalizedCoordinates(self):
        """ Return vertex coordinates scaled to the interval (-1,1) and centered at origo. """
        # Scale the verices so the max dimension is in the range (-1,1) to be compatible with the browser display
        vtx = self.coordinates()
        maxvtx = numpy.max(numpy.amax(vtx,axis=0))
        factor = 1/maxvtx
        vtx = factor*vtx
        
        # Compute mesh centroid
        centroid = numpy.mean(vtx,axis=0)
        # Shift so the centroid is now origo
        normalized_vtx = numpy.zeros(numpy.shape(vtx))
        for i,v in enumerate(vtx):
            normalized_vtx[i,:] = v - centroid
        
        
        return factor, normalized_vtx
    
    def scaledCoordinates(self):
        """ Return vertex coordinates scaled to the interval (-1,1). """
        # Scale the verices so the max dimension is in the range (-1,1) to be compatible with the browser display
        vtx = self.coordinates()
        maxvtx = numpy.max(numpy.amax(vtx,axis=0))
        factor = 1/maxvtx
        return factor, factor*vtx


    @classmethod
    def unitIntervalMesh(cls, nx, periodic=False):
        return cls.IntervalMesh(nx=nx, a=0, b=1, periodic=periodic)
    
    @classmethod
    def unitSquareMesh(cls, nx, ny, periodic=False):
        """ Unit Square of with nx,ny points in the respective axes. """
        return cls.SquareMesh(L=1, nx=nx, ny=ny, periodic=periodic)

    @classmethod
    def unitCubeMesh(cls, nx, ny, nz, periodic=False):
        """ Unit Cube of with nx,ny points in the respective axes. """
        return cls.CubeMesh(nx=nx, ny=ny, nz=nz, periodic=periodic)

    @classmethod
    def IntervalMesh(cls, nx, a, b, periodic=False):
        mesh = dolfin.IntervalMesh(nx, a, b)
        ret = URDMEMesh(mesh)
        if isinstance(periodic, bool) and periodic:
            ret.addPeriodicBoundaryCondition(IntervalMeshPeriodicBoundary(a=a, b=b))
        elif isinstance(periodic, dolfin.SubDomain):
            ret.addPeriodicBoundaryCondition(periodic)
        return ret

    @classmethod
    def SquareMesh(cls, L, nx, ny, periodic=False):
        """ Regular mesh of a square with side length L. """
        mesh = dolfin.RectangleMesh(0, 0, L, L, nx, ny)
        ret = URDMEMesh(mesh)
        if isinstance(periodic, bool) and periodic:
            ret.addPeriodicBoundaryCondition(SquareMeshPeriodicBoundary(Lx=L, Ly=L))
        elif isinstance(periodic, dolfin.SubDomain):
            ret.addPeriodicBoundaryCondition(periodic)
        return ret

    @classmethod
    def CubeMesh(cls, L, nx, ny, nz, periodic=False):
        """ Cube with nx,ny points in the respective axes. """
        mesh = dolfin.BoxMesh(0, 0, 0, L, L, L, nx, ny, nz)
        ret = URDMEMesh(mesh)
        if isinstance(periodic, bool) and periodic:
            ret.addPeriodicBoundaryCondition(CubeMeshPeriodicBoundary(Lx=L, Ly=L, Lz=L))
        elif isinstance(periodic, dolfin.SubDomain):
            ret.addPeriodicBoundaryCondition(periodic)
        return ret





    #@classmethod
    #def unitCircle(cls, nx,ny):
    #    """ Unit Square of with nx,ny points in the respective axes. """
    #    mesh = dolfin.UnitCircleMesh(nx,ny)
    #    return Mesh(mesh)

    #@classmethod
    #def unitSphere(cls, nx,ny):
    #    """ Unit Square of with nx,ny points in the respective axes. """
    #    mesh = dolfin.UnitSquareMesh(nx,ny)
    #    return Mesh(mesh)'t



    @classmethod
    def read_gmsh_mesh(cls, meshfile):
        """ Read a Gmsh mesh from file. """
        mr = GmshMeshReceiverBase()
        try:
            mesh = read_gmsh(mr, filename=meshfile)
        except:
            raise MeshImportError("Failed to import mesh: " + filename)

        return mesh

    @classmethod
    def read_dolfin_mesh(cls, filename=None, colors = []):
        """ Import a mesh in Dolfins native .xml format """

        try:
            dolfin_mesh = dolfin.Mesh(filename)
            mesh = URDMEMesh(mesh=dolfin_mesh)
            return mesh
        except Exception as e:
            raise MeshImportError("Failed to import mesh: " + filename+"\n" + str(e))

    def toTHREEJs(self, colors = None):
        """ return a Json string of the mesh in THREE Js format. 
            
            If a colors list is specified, it should have the num_voxels entries
            
        """
        document = {}
        document["metadata"] = {"formatVersion":3}
        gfdg,vtx = self.scaledNormalizedCoordinates()
        
        if self.topology().dim() == 2:
            # 2D
            num_elements = self.num_cells()
            # This is a fix for the built-in 2D meshes that only have x,y-coordinates.
            dims = numpy.shape(vtx)
            if dims[1] == 2:
                vtxx = numpy.zeros((dims[0],3))
                for i, v in enumerate(vtx):
                    vtxx[i,:]=(list(v)+[0])
                vtx = vtxx
        else:
            # 3D
            num_elements = self.num_facets()


        materials = [ {
                     "DbgColor" : 15658734,
                     "DbgIndex" : 0,
                     "DbgName" : "dummy",
                     "colorDiffuse" : [ 1, 0, 0 ],
                     } ]
                     
        document["materials"] = materials
        document["vertices"] = list(vtx.flatten())
        
        if colors == None:
            # Default color is blue
            colors = [255]*self.num_vertices()
        
        document["colors"] = colors
        
        self.init(2,0)
        connectivity = self.topology()(2,0)
        faces = []
        
       
        
        for i in range(num_elements):
            face = connectivity(i)
            f = []
            for ind in face:
                if int(ind) >= self.num_vertices():
                    raise Exception("Out of bounds")

                f.append(int(ind))
            faces += ([128]+f+f)
            document["faces"] = list(faces)
        
        #Test that we can index into vertices
        vertices = document["vertices"]
        
        return json.dumps(document)

    def _ipython_display_(self, filename=None):
        jstr = self.toTHREEJs()
        hstr = None
        with open(os.path.dirname(os.path.abspath(__file__))+"/data/three.js_templates/mesh.html",'r') as fd:
            hstr = fd.read()
        if hstr is None:
            raise Exception("could note open template mesh.html")
        hstr = hstr.replace('###PYURDME_MESH_JSON###',jstr)
        # Create a random id for the display div. This is to avioid multiple plots ending up in the same
        # div in Ipython notebook
        import uuid
        displayareaid=str(uuid.uuid4())
        hstr = hstr.replace('###DISPLAYAREAID###',displayareaid)
        html = '<div id="'+displayareaid+'" class="cell"></div>'

        if filename != None:
            with open(filename, 'w') as fd:
                fd.write("""
<html>
    <head>
        <title>PyURDME Result</title> <style>canvas { width: 100%; height: 100% }</style> </head>

        <body>
""")
                fd.write(html+hstr)
                fd.write("""
        </body>

</html>""")
        else:
            IPython.display.display(IPython.display.HTML(html+hstr))

class URDMEXmesh():
    """ Extended mesh object.

        Contains function spaces and dof mappings.
    """

    def __init__(self):
        self.coordinates = None
        self.function_space = {}
        self.vertex_to_dof_map = {}
        self.dof_to_vertex_map = {}


class URDMEResult(dict):
    """ Result object for a URDME simulation, extends the dict object. """

    def __init__(self, model=None, filename=None, loaddata=False):
        self.model = model
        self.sol = None
        self.U = None
        self.tspan = None
        self.data_is_loaded = False
        self.sol_initialized = False
        self.filename = filename
        if filename is not None and loaddata:
            self.read_solution()


    def get_endtime_model(self):
        """ Return a URDME model object with the initial conditions set to the final time point of the
            result object.
        """
        if self.model is None:
            raise Exception("can not continue a result with no model")
        # create a soft copy
        model2 = copy.copy(self.model)
        # set the initial conditions 
        model2.u0 = numpy.zeros(self.model.u0.shape)
        for s, sname in enumerate(self.model.listOfSpecies):
            model2.u0[s,:] = self.getSpecies(sname, timepoints=-1)
        return model2



    def __getstate__(self):
        """ Used by pickle to get state when pickling. We need to read the contents of the
        output file since we can't pickel file objects. """

        try:
            with open(self.filename,mode='rb') as fh:
                filecontents = fh.read()
        except Exception as e:
            raise Exception(("Error pickling model. Failed to read result file:",str(e)))
        
        state = self.__dict__
        state["filecontents"] = filecontents
        for key, item in state.items():
            try:
                pickle.dumps(item)
            except Exception as e:
                raise Exception(("Failed to pickle URDMEResult:", str(e)))

        return state


    def __setstate__(self, state):
        """ Used by pickle to set state when unpickling. """
        
        # If the object contains filecontents, write those to a new tmp file.
        try:
            filecontents = state.pop("filecontents",None)
            fd = tempfile.NamedTemporaryFile(delete=False)
            with open(fd.name, mode='wb') as fh:
                fh.write(filecontents)
            state["filename"] = fd.name
        except Exception as e:
            print "Error unpickling model, could not recreate the solution file."
            raise e

        for k,v in state.items():
            self.__dict__[k] = v

    def reorderDofToVoxel(self, M, num_species=None):
        """ Reorder the colums of M from dof ordering to vertex ordering. """
        
        fs = self.model.mesh.FunctionSpace()
        v2d = dolfin.vertex_to_dof_map(fs)
        if len(M.shape) == 1:
            num_timepoints = 1
        else:
            num_timepoints = M.shape[0]
        num_vox = self.model.mesh.getNumVoxels()
        if num_species is None:
            num_species = self.model.getNumSpecies()
        num_dofs = num_vox*num_species
        C = numpy.zeros((num_timepoints, num_dofs), dtype=numpy.float64)

#        reorder_map = numpy.zeros((num_vox*num_species), dtype=numpy.int)
#        for vox_ndx in range(num_vox):
#            for cndx in range(num_species):
#                reorder_map[vox_ndx*num_species+cndx] = v2d[vox_ndx]*num_species+cndx
#        if len(M.shape) == 1:
#            C[:,:] = M[reorder_map]
#        else:
#            C[:,:] = M[:, reorder_map]
        #for t in range(num_timepoints):
        for vox_ndx in range(num_vox):
            for cndx in range(num_species):
                try:
                    if len(M.shape) == 1:
                        C[:, vox_ndx*num_species+cndx] = M[v2d[vox_ndx]*num_species+cndx]
                    else:
                        C[:, vox_ndx*num_species+cndx] = M[:, v2d[vox_ndx]*num_species+cndx]
                except IndexError as e:
                    import traceback
                    #traceback.print_stack()
                    print traceback.format_exc()
                    print "C.shape: ", C.shape
                    print "M.shape: ", M.shape
                    print "num_timepoints: ", num_timepoints
                    print "vox_ndx={0},num_species={1},cndx={2}".format(vox_ndx,num_species,cndx)
                    print "v2d[vox_ndx]={0}".format(v2d[vox_ndx])
                    print "vox_ndx*num_species+cndx={0}".format(vox_ndx*num_species+cndx)
                    print "v2d[vox_ndx]*num_species+cndx={0}".format(v2d[vox_ndx]*num_species+cndx)
                    raise e
        return C

    def read_solution(self):
        """ Read the tspan and U matrix into memory. """
        
        resultfile = h5py.File(self.filename, 'r')
        U = resultfile['U']
        U = numpy.array(U)
        
        tspan = resultfile['tspan']
        tspan = numpy.array(tspan).flatten()
        resultfile.close()
        
        # Reorder the dof from dof ordering to voxel ordering
        U = self.reorderDofToVoxel(U)

        self.U = U
        self.tspan = tspan
        self.data_is_loaded = True
 
    def getSpecies(self, species, timepoints="all", concentration=False):
        """ Returns a slice (view) of the output matrix U that contains one species for the timepoints
            specified by the time index array. The default is to return all timepoints. 
            
            Data is loaded by slicing directly in the hdf5 dataset, i.e. it the entire
            content of the file is not loaded in memory and the U matrix 
            is never added to the object.
            
            if concentration is False (default), the integer, raw, trajectory data is returned,
            if set to True, the concentration (=copy_number/volume) is returned.
            
        """
        
        if isinstance(species, Species):
            spec_name = species.name
        else:
            spec_name = species
        
        species_map = self.model.speciesMap()
        num_species = self.model.getNumSpecies()
        spec_indx = species_map[spec_name]
        
        resultfile = h5py.File(self.filename, 'r')
        #Ncells = self.model.mesh.num_vertices()  # Need dof ordering numVoxels
        U = resultfile['U']
        Ncells = U.shape[1]/num_species
        
        if timepoints  ==  "all":
            Uslice= U[:,(spec_indx*Ncells):(spec_indx*Ncells+Ncells)]
        else:
            Uslice = U[timepoints,(spec_indx*Ncells):(spec_indx*Ncells+Ncells)]
        
        if concentration:
            Uslice = self._copynumber_to_concentration(Uslice)
        
        # Reorder the dof from dof ordering to voxel ordering
        Uslice = self.reorderDofToVoxel(Uslice, num_species=1)
        
        # Make sure we return 1D slices as flat arrays
        dims = numpy.shape(Uslice)
        if dims[0] == 1:
            Uslice = Uslice.flatten()
        
        resultfile.close()
        return Uslice
            
    def __setattr__(self, k, v):
        if k in self.keys():
            self[k] = v
        elif not hasattr(self, k):
            self[k] = v
        else:
            raise AttributeError, "Cannot set '%s', cls attribute already exists" % ( k, )

    def __setupitems__(self, k):
        if k == 'sol' and not self.sol_initialized:
            self._initialize_sol()
        elif (k == 'U' or k == 'tspan') and not self.data_is_loaded:
            if self.filename is None:
                raise AttributeError("This result object has no data file.")
            self.read_solution()

    def __getitem__(self, k):
        self.__setupitems__(k)
        if k in self.keys():
            return self.get(k)
        raise KeyError("Object has no attribute {0}".format(k))

    def __getattr__(self, k):
        self.__setupitems__(k)
        if k in self.keys():
            return self.get(k)
        raise AttributeError("Object has no attribute {0}".format(k))

    def __del__(self):
        """ Deconstructor. """
            #   if not self.data_is_loaded:
        try:
            # Clean up data file
            os.remove(self.filename)
        except OSError as e:
            #print "URDMEResult.__del__: Could not delete result file'{0}': {1}".format(self.filename, e)
            pass


    def _initialize_sol(self):
        """ Initialize the sol variable. This is a helper function to toVTK. """
        
        # Create Dolfin Functions for all the species
        sol = {}

        if self.model is None:
            raise URDMEError("URDMEResult.model must be set before the sol attribute can be accessed.")
        numvox = self.model.mesh.num_vertices()
        fs = self.model.mesh.FunctionSpace()
        vertex_to_dof_map = dolfin.vertex_to_dof_map(fs)
        dof_to_vertex_map = dolfin.dof_to_vertex_map(fs)

        # The result is loaded in dolfin Functions, one for each species and time point
        for i, spec in enumerate(self.model.listOfSpecies):

            species = self.model.listOfSpecies[spec]
            spec_name = species.name

            spec_sol = {}
            for j, time in enumerate(self.tspan):
                
                func = dolfin.Function(fs)
                func_vector = func.vector()

                S = self.getSpecies(spec, [j])

                for voxel in range(numvox):
                    ix  = vertex_to_dof_map[voxel]
                    try:
                        func_vector[ix] = float(S[voxel]/self.model.dofvol[ix])

                    except IndexError as e:
                        print "func_vector.size(): ", func_vector.size()
                        print "dolfvox: ",dolfvox
                        print "S.shape: ",S.shape
                        print "voxel: ",voxel
                        print "vertex_to_dof_map[voxel]", vertex_to_dof_map[voxel]
                        print "self.model.dofvol.shape: ", self.model.dofvol.shape
                        raise e

                spec_sol[time] = func

            sol[spec] = spec_sol
        self.sol = sol
        self.sol_initialized = True
        return sol
            
    
    def toVTK(self, species, folder_name):
        """ Dump the trajectory to a collection of vtk files in the folder folder_name (created if non-existant). """
        
        self._initialize_sol()
        subprocess.call(["mkdir", "-p", folder_name])
        func = dolfin.Function(self.model.mesh.FunctionSpace())
        func_vector = func.vector()
        fd = dolfin.File(folder_name+"/trajectory.pvd")
        numvox = self.model.mesh.getNumDofVoxels()

        for i, time in enumerate(self.tspan):
            solvector = (self.sol[species][time]).vector()
            for dof in range(numvox):
                func_vector[dof] = solvector[dof]
            fd << func

    def toXYZ(self, filename, species=None, file_format="VMD"):
        """ Dump the solution attached to a model as a xyz file. This format can be
            read by e.g. VMD, Jmol and Paraview. """

        if self.U is None:
            raise URDMEError("No solution found in the model.")

        #outfile = open(filename,"w")
        dims = numpy.shape(self.U)
        Ndofs = dims[0]
        Mspecies = len(self.model.listOfSpecies)
        Ncells = Ndofs / Mspecies

        coordinates = self.model.mesh.getVoxels()
        coordinatestr = coordinates.astype(str)

        if species == None:
            species = list(self.model.listOfSpecies.keys())

        if file_format == "VMD":
            outfile = open(filename, "w")
            filestr = ""
            for i, time in enumerate(self.tspan):
                number_of_atoms = numpy.sum(self.U[:, i])
                filestr += (str(number_of_atoms) + "\n" + "timestep " + str(i) + " time " + str(time) + "\n")
                for j, spec in enumerate(species):
                    for k in range(Ncells):
                        for mol in range(self.U[k * Mspecies + j, i]):
                            # Sample a random position in a sphere of radius computed from the voxel volume
                            # TODO: Sample volume
                            linestr = spec + "\t" + '\t'.join(coordinatestr[k, :]) + "\n"
                            filestr += linestr

            outfile.write(filestr)
            outfile.close()

        elif file_format == "ParaView":
            foldername = filename
            os.mkdir(foldername)
            for i, time in enumerate(self.tspan):
                outfile = open(foldername + "/" + filename + "." + str(i), "w")
                number_of_atoms = numpy.sum(self.U[:, i])
                filestr = ""
                filestr += (str(number_of_atoms) + "\n" + "timestep " + str(i) + " time " + str(time) + "\n")
                for j, spec in enumerate(self.model.listOfSpecies):
                    for k in range(Ncells):
                        for mol in range(model.U[k * Mspecies + j, i]):
                            linestr = spec + "\t" + '\t'.join(coordinatestr[k, :]) + "\n"
                            filestr += linestr
                outfile.write(filestr)
                outfile.close()



    def printParticlejs(self,species,time_index):
    
        import random
        
        with open(os.path.dirname(os.path.abspath(__file__))+"/data/three.js_templates/particles.html",'r') as fd:
            template = fd.read()
        
        factor, coordinates = self.model.mesh.scaledNormalizedCoordinates()
        dims = numpy.shape(coordinates)
        if dims[1]==2:
            is3d = 0
            vtxx = numpy.zeros((dims[0],3))
            for i, v in enumerate(coordinates):
                vtxx[i,:]=(list(v)+[0])
            coordinates = vtxx
        else:
            is3d = 1
        
        h = self.model.mesh.meshSize()
        
        x=[];
        y=[];
        z=[];
        c=[];
        radius = []

        total_num_particles = 0
        colors = ["blue","red","yellow", "green"]
        
        for j,spec in enumerate(species):
            
            timeslice = self.getSpecies(spec, 0)
            #timeslice = US[time_index,:]
            ns = numpy.sum(timeslice)
            total_num_particles += ns

            for i, particles in enumerate(timeslice):
                # "Radius" of voxel
                hix = h[i]*factor
                hiy = hix;
                hiz = hix*is3d
                
                for particle in range(particles):
                    x.append((coordinates[i,0]+random.uniform(-1,1)*hix))
                    y.append((coordinates[i,1]+random.uniform(-1,1)*hiy))
                    z.append((coordinates[i,2]+random.uniform(-1,1)*hiz))
                    try:
                        radius.append(self.model.listOfSpecies[spec].reaction_radius)
                    except:
                        radius.append(0.01)
                    
                    c.append(colors[j])
    
        docstr = template.replace("__NUM__MOLECULES__", str(total_num_particles))
        docstr = docstr.replace("__X__",str(x))
        docstr = docstr.replace("__Y__",str(y))
        docstr = docstr.replace("__Z__",str(z))
        docstr = docstr.replace("__COLOR__",str(c))
        docstr = docstr.replace("__RADIUS__",str(radius))


        return docstr


    def toTHREEJs(self, species, time_index):
        """ Return a json serialized document that can 
            be read and visualized by three.js.
        """
        
        colors = self._compute_solution_colors(species,time_index)
        return self.model.mesh.toTHREEJs(colors=colors)

    def _copynumber_to_concentration(self,copy_number_data):
        """ Scale compy numbers to concentrations (in unit mol/volume),
            where the volume unit is defined by the user input.
        """
        
        fs = self.model.mesh.FunctionSpace()
        v2d = dolfin.vertex_to_dof_map(fs)
        shape = numpy.shape(copy_number_data)
        if len(shape) == 1:
            shape = (1,shape[0])

        scaled_sol = numpy.zeros(shape)
        scaled_sol[:,:] = copy_number_data
        dims = numpy.shape(scaled_sol)
        
        for t in range(dims[0]):
            timeslice = scaled_sol[t,:]
            for i,cn in enumerate(timeslice):
                scaled_sol[t, i] = float(cn)/(6.022e23*self.model.dofvol[v2d[i]])

        return scaled_sol


    def _compute_solution_colors(self,species, time_index):
        """ Create a color list for species at time. """
        
        timeslice = self.getSpecies(species,time_index, concentration = True)
        import matplotlib.cm
        
        # Get RGB color map proportinal to the concentration.
        cm = matplotlib.cm.ScalarMappable()
        crgba= cm.to_rgba(timeslice, bytes = True)
                                     
        # Convert RGB to HEX
        colors= []
        for row in crgba:
            colors.append(self._rgb_to_hex(tuple(list(row[1:]))))

        # Convert Hex to Decimal
        for i,c in enumerate(colors):
            colors[i] = int(c,0)

        return colors

    def _rgb_to_hex(self, rgb):
        return '0x%02x%02x%02x' % rgb


    def display(self,species,time_index):

        jstr = self.toTHREEJs(species,time_index)
        hstr = None
        with open(os.path.dirname(os.path.abspath(__file__))+"/data/three.js_templates/solution.html",'r') as fd:
            hstr = fd.read()
        if hstr is None:
            raise Exception("could note open template mesh.html")
        hstr = hstr.replace('###PYURDME_MESH_JSON###',jstr)
        
        # Create a random id for the display div. This is to avioid multiple plots ending up in the same
        # div in Ipython notebook
        import uuid
        displayareaid=str(uuid.uuid4())
        hstr = hstr.replace('###DISPLAYAREAID###',displayareaid)
        
        html = '<div id="'+displayareaid+'" class="cell"></div>'
        IPython.display.display(IPython.display.HTML(html+hstr))


class URDMESolver:
    """ Abstract class for URDME solvers. """

    def __init__(self, model, solver_path=None, report_level=0, model_file=None):
        """ Constructor. """
        if not isinstance(model, URDMEModel):
            raise URDMEError("URDMEsolver constructors must take a URDMEModel as an argument.")
        if not issubclass(self.__class__, URDMESolver):
            raise URDMEError("Solver classes must be a subclass of URDMESolver.")
        if not hasattr(self, 'NAME'):
            raise URDMEError("Solver classes must implement a NAME attribute.")

        self.model = model
        self.is_compiled = False
        self.report_level = report_level
        self.model_file = model_file
        self.infile_name = None
        self.delete_infile = False
        self.model_name = self.model.name
        self.solver_base_dir = None

        # For the remote execution
        self.temp_urdme_root = None
        
        self.URDME_ROOT =  os.path.dirname(os.path.abspath(__file__))+"/urdme"

        #print "solver_path={0}".format(solver_path)
        if solver_path is None or solver_path == "":
            self.URDME_BUILD = self.URDME_ROOT + '/build/'
        else:
            self.URDME_BUILD = solver_path + '/build/'
            os.environ['SOLVER_ROOT'] = solver_path

    def __getstate__(self):
        """ Save the state of the solver, saves all instance variables
            and reads all the files necessary to compile the solver off
            of the file system and stores it in a separate state variable.
            If the solver model files is specified, it saves that too.
            This is used by Pickle.
        """
        ret = {}
        # Save the instance variables
        ret['vars'] = self.__dict__.copy()
        # The model object is not picklabe due to the Swig-objects from Dolfin
        #ret['vars']['model'] = None
        ret['vars']['is_compiled'] = False
        # Create temp root
        tmproot = tempfile.mkdtemp()
        # Get the propensity file
        model_file = tmproot+'/'+self.model_name + '_pyurdme_generated_model'+ '.c'
        ret['model_file'] = os.path.basename(model_file)
        if self.model_file == None:
            self.createPropensityFile(file_name=model_file)
        else:
            subprocess.call('cp '+self.model_file+' '+model_file, shell=True)
        # Get the solver source files
        os.mkdir(tmproot+'/include')
        os.mkdir(tmproot+'/src')
        os.mkdir(tmproot+'/src/'+self.NAME)
        #TODO: what if solverdir is not the same as URDME_ROOT ?
        subprocess.call('cp '+self.URDME_ROOT+'/src/*.c '+tmproot+'/src/', shell=True)
        subprocess.call('cp '+self.URDME_ROOT+'/src/'+self.NAME+'/*.* '+tmproot+'/src/'+self.NAME+'/', shell=True)
        subprocess.call('cp '+self.URDME_ROOT+'/include/*.h '+tmproot+'/include/', shell=True)
        #TODO: get the include files from solvers not in the default path (none currently implement this).
        # Get the Makefile
        os.mkdir(tmproot+'/build')
        subprocess.call('cp '+self.URDME_BUILD+'Makefile.'+self.NAME+' '+tmproot+'/build/Makefile.'+self.NAME, shell=True)
        # Get the input file
        input_file = tmproot+'/model_input.mat'
        ret['input_file'] = os.path.basename(input_file)
        self.model.serialize(filename=input_file, report_level=self.report_level)
        ##
        origwd = os.getcwd()
        os.chdir(tmproot)
        tarname = tmproot+'/'+self.NAME+'.tar.gz'
        subprocess.call('tar -czf '+tarname+' src include build '+os.path.basename(input_file)+' '+os.path.basename(model_file), shell=True)
        with open(tarname, 'r') as f:
            ret['SolverFiles'] = f.read()
        os.chdir(origwd)
        shutil.rmtree(tmproot)
        # return the state
        return ret

    def __setstate__(self, state):
        """ Set all instance variables for the object, and create a unique temporary
            directory to store all the solver files.  URDME_BUILD is set to this dir,
            and is_compiled is always set to false.  This is used by Pickle.
        """
        # 0. restore the instance variables
        for key, val in state['vars'].iteritems():
            self.__dict__[key] = val
        # 1. create temporary directory = URDME_ROOT
        self.temp_urdme_root = tempfile.mkdtemp()
        self.URDME_ROOT = self.temp_urdme_root
        self.URDME_BUILD = self.temp_urdme_root+'/build/'
        origwd = os.getcwd()
        os.chdir(self.temp_urdme_root)
        tarname = self.temp_urdme_root+'/'+self.NAME+'.tar.gz'
        with open(tarname, 'wd') as f:
            f.write(state['SolverFiles'])
        subprocess.call('tar -zxf '+tarname, shell=True)
        os.chdir(origwd)
        # Model File
        self.model_file = self.temp_urdme_root+'/'+state['model_file']
        # Input File
        self.infile_name = self.temp_urdme_root+'/'+state['input_file']



    def __del__(self):
        """ Deconstructor.  Removes the compiled solver."""
        if self.delete_infile:
            try:
                os.remove(self.infile_name)
            except OSError as e:
                print "Could not delete '{0}'".format(self.infile_name)
        if self.solver_base_dir is not None:
            try:
                shutil.rmtree(self.solver_base_dir)
            except OSError as e:
                print "Could not delete '{0}'".format(self.solver_base_dir)
        if self.temp_urdme_root is not None:
            try:
                shutil.rmtree(self.temp_urdme_root)
            except OSError as e:
                print "Could not delete '{0}'".format(self.temp_urdme_root)


    def compile(self):
        """ Compile the model."""

        # Create a unique directory each time call to compile.
        self.solver_base_dir = tempfile.mkdtemp()
        self.solver_dir = self.solver_base_dir + '/.urdme/'
        #print "URDMESolver.compile()  self.solver_dir={0}".format(self.solver_dir)

        if self.report_level >= 1:
            print "Compiling Solver"

        if os.path.isdir(self.solver_dir):
            try:
                shutil.rmtree(self.solver_dir)
            except OSError as e:
                pass
        try:
            os.mkdir(self.solver_dir)
        except Exception as e:
            pass

        # Write the propensity file
        self.propfilename = self.model_name + '_pyurdme_generated_model'
        if self.model_file == None:
            if self.report_level > 1:
                prop_file_name=self.solver_dir + self.propfilename + '.c'
                print "Creating propensity file {0}".format(prop_file_name)
            self.createPropensityFile(file_name=prop_file_name)
        else:
            cmd = " ".join(['cp', self.model_file, self.solver_dir + self.propfilename + '.c'])
            if self.report_level > 1:
                print cmd
            subprocess.call(cmd)

        # Build the solver
        makefile = 'Makefile.' + self.NAME
        cmd = " ".join([ 'cd', self.solver_base_dir , ';', 'make', '-f', self.URDME_BUILD + makefile, 'URDME_ROOT=' + self.URDME_ROOT, 'URDME_MODEL=' + self.propfilename])
        if self.report_level > 1:
            print "cmd: {0}\n".format(cmd)
        try:
            handle = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            return_code = handle.wait()
        except OSError as e:
            print "Error, execution of compilation raised an exception: {0}".format(e)
            print "cmd = {0}".format(cmd)
            raise URDMEError("Compilation of solver failed")

        if return_code != 0:
            try:
                print handle.stdout.read()
                print handle.stderr.read()
            except Exception as e:
                pass
            raise URDMEError("Compilation of solver failed, return_code={0}".format(return_code))

        if self.report_level > 1:
            print handle.stdout.read()
            print handle.stderr.read()


        self.is_compiled = True


    def run_ensemble(self, number_of_trajectories, seed=None, input_file=None, loaddata=False):
        """ Run multiple simulations of the model.

        number_of_trajectories: How many trajectories should be run.
        seed: the random number seed.
        input_file: the filename of the solver input data file .
        loaddata: boolean, should the result object load the data into memory on creation.

        Returns:
            A list of URDMEResult objects.
        """
        result = []
        for ndx in range(number_of_trajectories):
            if seed is None:
                result.append(self.run(input_file=input_file, loaddata=loaddata))
            else:
                result.append(self.run(seed=seed+ndx, input_file=input_file, loaddata=loaddata))
        return result

    def run(self, seed=None, input_file=None, loaddata=False):
        """ Run one simulation of the model.

        seed: the random number seed.
        input_file: the filename of the solver input data file .
        loaddata: boolean, should the result object load the data into memory on creation.

        Returns:
            URDMEResult object.
        """
        # Check if compiled, call compile() if not.
        if not self.is_compiled:
            self.compile()

        if input_file is None:
            if self.infile_name is None or not os.path.exists(self.infile_name):
                # Get temporary input and output files
                infile = tempfile.NamedTemporaryFile(delete=False)

                # Write the model to an input file in .mat format
                self.model.serialize(filename=infile, report_level=self.report_level)
                infile.close()
                self.infile_name = infile.name
                #self.delete_infile = True
        else:
            self.infile_name = input_file
            self.delete_infile = False

        outfile = tempfile.NamedTemporaryFile(delete=False)
        outfile.close()

        if not os.path.exists(self.infile_name):
            raise URDMEError("input file not found.")

        # Execute the solver
        urdme_solver_cmd = [self.solver_dir + self.propfilename + '.' + self.NAME, self.infile_name, outfile.name]
        if seed is not None:
            urdme_solver_cmd.append(str(seed))
        if self.report_level >= 1:
            print 'cmd: {0}\n'.format(urdme_solver_cmd)
        try:
            if self.report_level >= 1:  #stderr & stdout to the terminal
                handle = subprocess.Popen(urdme_solver_cmd)
            else:
                handle = subprocess.Popen(urdme_solver_cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            return_code = handle.wait()
        except OSError as e:
            print "Error, execution of solver raised an exception: {0}".format(e)
            print "urdme_solver_cmd = {0}".format(urdme_solver_cmd)
            raise URDMEError("Solver execution failed")

        if return_code != 0:
            print outfile.name
            print return_code
            if self.report_level >= 1:
                print handle.stderr.read(), handle.stdout.read()
            print "urdme_solver_cmd = {0}".format(urdme_solver_cmd)
            raise URDMEError("Solver execution failed")

        #if self.report_level >= 1:
        #    print handle.stdout.read()
        #   print handle.stderr.read()

        #Load the result from the hdf5 output file.
        try:
            result = URDMEResult(self.model, outfile.name, loaddata=loaddata)
            result["Status"] = "Sucess"
            return result
        except Exception as e:
            exc_info = sys.exc_info()
            # Clean up
            #if self.delete_infile:
            #    os.remove(self.infile_name)
            os.remove(outfile.name)
            raise exc_info[1], None, exc_info[2]



    def createPropensityFile(self, file_name=None):
        """ Generate the C propensity file that is used to compile the URDME solvers.
            Only mass action propensities are supported.

        """

        template = open(os.path.abspath(os.path.dirname(__file__)) + '/data/propensity_file_template.c', 'r')
        propfile = open(file_name, "w")
        propfilestr = template.read()

        speciesdef = ""
        i = 0
        for S in self.model.listOfSpecies:
            speciesdef += "#define " + S + " " + "x[" + str(i) + "]" + "\n"
            speciesdef += "#define " + S + "_INDEX " +  str(i) + "\n"
            i += 1

        propfilestr = propfilestr.replace("__DEFINE_SPECIES__", speciesdef)

        propfilestr = propfilestr.replace("__NUMBER_OF_REACTIONS__", str(self.model.getNumReactions()))
        propfilestr = propfilestr.replace("__NUMBER_OF_SPECIES__", str(self.model.getNumSpecies()))
        propfilestr = propfilestr.replace("__NUMBER_OF_VOXELS__", str(self.model.mesh.getNumVoxels()))

        # Create defines for the DataFunctions.
        data_fn_str = ""
        i = 0
        for d in self.model.listOfDataFunctions:
            if d.name is None:
                raise URDMEError("DataFunction {0} does not have a name attributed defined.".format(i))
            data_fn_str += "#define " + d.name + " data[" + str(i) + "]\n"
            i += 1
        propfilestr = propfilestr.replace("__DEFINE_DATA_FUNCTIONS__", str(data_fn_str))

        # Make sure all paramters are evaluated to scalars before we write them to the file.
        self.model.resolveParameters()
        parameters = ""
        for p in self.model.listOfParameters:
            parameters += "const double " + p + " = " + str(self.model.listOfParameters[p].value) + ";\n"
        propfilestr = propfilestr.replace("__DEFINE_PARAMETERS__", str(parameters))


        # Reactions
        funheader = "double __NAME__(const int *x, double t, const double vol, const double *data, int sd)"
        #funheader = "double __NAME__(const int *x, double t, const double vol, const double *data, int sd, int voxel, int *xx, const size_t *irK, const size_t *jcK, const double *prK)"

        funcs = ""
        funcinits = ""
        i = 0
        for R in self.model.listOfReactions:
            func = ""
            rname = self.model.listOfReactions[R].name
            func += funheader.replace("__NAME__", rname) + "\n{\n"
            if self.model.listOfReactions[R].restrict_to == None:
                func += self.model.listOfReactions[R].propensity_function

            else:
                func += "if("
                if isinstance(self.model.listOfReactions[R].restrict_to, list):
                    for sd in self.model.listOfReactions[R].restrict_to:
                        func += "sd == " + str(sd) + "||"
                    func = func[:-2]
                elif isinstance(self.model.listOfReactions[R].restrict_to, int):
                    func += "sd == " +  str(self.model.listOfReactions[R].restrict_to)
                else:
                    raise URDMEError("When restricting reaction to subdomains, you must specify either a list or an int")
                func += "){\n"
                func += self.model.listOfReactions[R].propensity_function

                func += "\n}else{"
                func += "\n\treturn 0.0;}"


            func += "\n}"
            funcs += func + "\n\n"
            funcinits += "    ptr[" + str(i) + "] = " + rname + ";\n"
            i += 1

        propfilestr = propfilestr.replace("__DEFINE_REACTIONS__", funcs)
        propfilestr = propfilestr.replace("__DEFINE_PROPFUNS__", funcinits)

        propfile.write(propfilestr)
        propfile.close()



def urdme(model=None, solver='nsm', solver_path="", model_file=None, input_file=None, seed=None, report_level=0):
    """ URDME solver interface.

        Similar to model.run() the urdme() function provides an interface that is backwards compatiable with the
        previous URDME implementation.

        After sucessful execution, urdme returns a URDMEResults object with the following members:
        U:         the raw copy number output in a matrix with dimension (Ndofs, num_time_points)
        tspan:     the time span vector containing the time points that corresponds to the columns in U
        status:    Sucess if the solver executed without error
        stdout:    the standard ouput stream from the call to the core solver
        stderr:    the standard error stream from the call to the core solver

    """
    
    
    #If solver is a subclass of URDMESolver, use it directly.
    if isinstance(solver, (type, types.ClassType)) and  issubclass(solver, URDMESolver):
        sol = solver(model, solver_path, report_level, model_file=model_file)
    elif type(solver) is str:
        if solver == 'nsm':
            from nsmsolver import NSMSolver
            sol = NSMSolver(model, solver_path, report_level, model_file=model_file)
        elif solver == 'nem':
            from nemsolver import NEMSolver
            sol = NEMSolver(model, solver_path, report_level, model_file=model_file)
        else:
            raise URDMEError("Unknown solver: {0}".format(solver_name))
    else:
        raise URDMEError("solver argument to urdme() must be a string or a URDMESolver class object.")

    sol.compile()
    return sol.run(seed, input_file=input_file)


class URDMEDataFunction():
    """ Abstract class used to constuct the URDME data vector. """
    name = None
    def __init__(self, name=None):
        if name is not None:
            self.name = name
        if self.name is None:
            raise Exception("URDMEDataFunction must have a 'name'")
        
    def map(self, x):
        """ map() takes the coordinate 'x' and returns a list of doubles.
        Args:
            x: a list of 3 ints.
        Returns:
            a doubles.
        """
        raise Exception("URDMEDataFunction.map() not implemented.")


class MeshImportError(Exception):
    """ Exception to raise when encourntering and error importing a mesh. """
    pass

class URDMEError(Exception):
    pass

class ModelException(Exception):
    pass

class InvalidSystemMatrixException(Exception):
    pass


class IntervalMeshPeriodicBoundary(dolfin.SubDomain):
    def __init__(self, a=0.0, b=1.0):
        """ 1D domain from a to b. """
        dolfin.SubDomain.__init__(self)
        self.a = a
        self.b = b

    def inside(self, x, on_boundary):
        return not bool((dolfin.near(x[0], self.b)) and on_boundary)

    def map(self, x, y):
        if dolfin.near(x[0], self.b):
            y[0] = self.a + (x[0] - self.b)

class SquareMeshPeriodicBoundary(dolfin.SubDomain):
    """ Sub domain for Periodic boundary condition """
    def __init__(self, Lx=1.0, Ly=1.0):
        dolfin.SubDomain.__init__(self)
        self.Lx = Lx
        self.Ly = Ly

    def inside(self, x, on_boundary):
        """ Left boundary is "target domain" G """
        # return True if on left or bottom boundary AND NOT on one of the two corners (0, 1) and (1, 0)
        return bool((dolfin.near(x[0], 0) or dolfin.near(x[1], 0)) and
                (not ((dolfin.near(x[0], 0) and dolfin.near(x[1], 1)) or
                        (dolfin.near(x[0], 1) and dolfin.near(x[1], 0)))) and on_boundary)

    def map(self, x, y):
        ''' # Map right boundary G (x) to left boundary H (y) '''
        if dolfin.near(x[0], self.Lx) and dolfin.near(x[1], self.Ly):
            y[0] = x[0] - self.Lx
            y[1] = x[1] - self.Ly
        elif dolfin.near(x[0], self.Lx):
            y[0] = x[0] - self.Lx
            y[1] = x[1]
        else:   # near(x[1], 1)
            y[0] = x[0]
            y[1] = x[1] - self.Ly

class CubeMeshPeriodicBoundary(dolfin.SubDomain):
    """ Sub domain for Periodic boundary condition """
    def __init__(self, Lx=1.0, Ly=1.0, Lz=1.0):
        dolfin.SubDomain.__init__(self)
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz

    def inside(self, x, on_boundary):
        """ Left boundary is "target domain" G """
        # return True if on left or bottom boundary AND NOT on one of the two corners (0, 1) and (1, 0)
        return bool(
                (dolfin.near(x[0], 0) or dolfin.near(x[1], 0) or dolfin.near(x[3], 0))
                and (not (
                        (dolfin.near(x[0], 1) and dolfin.near(x[1], 0) and dolfin.near(x[1], 0)) or
                        (dolfin.near(x[0], 0) and dolfin.near(x[1], 1) and dolfin.near(x[1], 0)) or
                        (dolfin.near(x[0], 0) and dolfin.near(x[1], 0) and dolfin.near(x[1], 1))
                    ))
                and on_boundary
               )

    def map(self, x, y):
        ''' # Map right boundary G (x) to left boundary H (y) '''
        if dolfin.near(x[0], self.Lx) and dolfin.near(x[1], self.Ly) and dolfin.near(x[2], self.Lz):
            y[0] = x[0] - self.Lx
            y[1] = x[1] - self.Ly
            y[2] = x[2] - self.Lz
        elif dolfin.near(x[0], self.Lx) and dolfin.near(x[1], self.Ly):
            y[0] = x[0] - self.Lx
            y[1] = x[1] - self.Ly
            y[2] = x[2]
        elif dolfin.near(x[0], self.Lx) and dolfin.near(x[2], self.Lz):
            y[0] = x[0] - self.Lx
            y[1] = x[1]
            y[2] = x[2] - self.Lz
        elif dolfin.near(x[1], self.Ly) and dolfin.near(x[2], self.Lz):
            y[0] = x[0]
            y[1] = x[1] - self.Ly
            y[2] = x[2] - self.Lz
        elif dolfin.near(x[0], self.Lx):
            y[0] = x[0] - self.Lx
            y[1] = x[1]
            y[2] = x[2]
        elif dolfin.near(x[1], self.Ly):
            y[0] = x[0]
            y[1] = x[1] - self.Ly
            y[2] = x[2]
        elif dolfin.near(x[2], self.Lz):
            y[0] = x[0]
            y[1] = x[1]
            y[2] = x[2] - self.Lz
        else:
            y[0] = x[0]
            y[1] = x[1]
            y[2] = x[2]




if __name__ == '__main__':
    """ Command line interface to URDME. Execute URDME given a model file. """

