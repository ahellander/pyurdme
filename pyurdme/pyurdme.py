# pylint: disable-msg=C0301
# pylint: disable-msg=C0103

from model import *

import numpy as np
import os
import re
import scipy.io as spio
import scipy.sparse as scisp
import shutil
import subprocess
import sys
import tempfile

import gmsh
import numpy
import scipy.sparse

from nsmsolver import NSMSolver
from nemsolver import NEMSolver

try:
    import h5py
except:
    raise Exception("pyurdme requires h5py.")

try:
    import dolfin
    dolfin.parameters["linear_algebra_backend"] = "uBLAS"
except Exception:
    print "Warning: Could not import dolphin. Only simple Cartesain examples will work."
    ONLY_CARTESIAN = True



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

        # subdomins is a list of MeshFunctions with subdomain marker information
        self.subdomains = OrderedDict()

        # This dictionary hold information about the subdomains each species is active on
        self.species_to_subdomains = {}

        self.tspan = None
        self.vol = None

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

    def addSubDomain(self,subdomain):
        if not subdomain.dim() in self.subdomains.keys():
            self.subdomains[subdomain.dim()] = subdomain
        else:
            raise ModelException("Failed to add subdomain function of dim "+str(subdomain.dim())+". Only one subdomain function of a given dimension is allowed.")

    def createStoichiometricMatrix(self):
        """ Generate a stoichiometric matrix in sparse CSC format. """

        if not hasattr(self, 'species_map'):
            self.__initializeSpeciesMap()
        if self.getNumReactions() > 0:
            ND = np.zeros((self.getNumSpecies(), self.getNumReactions()))
            for i, r in enumerate(self.listOfReactions):
                R = self.listOfReactions[r]
                reactants = R.reactants
                products  = R.products

                for s in reactants:
                    ND[self.species_map[s], i] -= reactants[s]
                for s in products:
                    ND[self.species_map[s], i] += products[s]

            N = scisp.csc_matrix(ND)
        else:
            N = numpy.zeros((self.getNumSpecies(), self.getNumReactions()))

        return N

    def createDependencyGraph(self):
        """ Construct the sparse dependecy graph. """

        #TODO: Automatically create a dependency graph (cannot be optimal, but good enough.)
        GF = np.ones((self.getNumReactions(), self.getNumReactions() + self.getNumSpecies()))
        try:
            G = scisp.csc_matrix(GF)
        except:
            G = GF

        return G

    



    def timespan(self, tspan):
        """ Set the time span of simulation. """
        self.tspan = tspan

    def _initialize_species_to_subdomains(self):
        """ Initialize the species mapping to subdomains. The default
            is that a species is active in all the defined subdomains.
        """
        # The unique elements of the subdomain MeshFunctions
        sds = []
        for dim, subdomain in self.subdomains.items():
            sds = sds + list(np.unique(subdomain.array()).flatten())
        sds = np.unique(sds)
        sds = list(sds)

        # This explicit typecast is necessary for UFL not to choke on the subdomain ids.
        for i, sd in enumerate(sds):
            sds[i] = int(sd)
        
        # This conversion is necessary for UFL not to choke on the subdomain ids.
        for i,sd in enumerate(sds):
            sds[i]=int(sd)
        try:
            sds.remove(0)
        except:
            pass
        
        # If a species is not present as key in the species_to_subdomain mapping,
        # we label it as active in all subdomains
        for spec_name in self.listOfSpecies:
            species = self.listOfSpecies[spec_name]
            if species not in self.species_to_subdomains.keys():
                self.species_to_subdomains[species] = sds


    def restrict(self, species, subdomains):
        self.species_to_subdomains[species] = subdomains


    def subdomainVector(self,subdomains={}):
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
        sd = numpy.zeros((1,self.mesh.getNumVoxels()))
        if subdomains == {}:
            self.sd = sd.flatten()
            print subdomains
        else:
            for dim,subdomain in subdomains.items():
                # Map all facet labels to vertex labels
                tovertex = self.mesh.topology()(dim,0)
                for i in range(subdomain.size()):
                    for vtx in tovertex(i):
                        if subdomain[i]!=0: # TODO: Temporary hack to fix issue with Gmesh facet_region files.
                            sd[0,vtx] = subdomain[i]
            
        self.sd = sd.flatten()
        return self.sd


    def initializeInitialValue(self):
        """ Create all-zeros inital condition matrix. """
        ns = self.getNumSpecies()
        nv = self.mesh.getNumVoxels()
        self.u0 = np.zeros((ns, nv))

    def meshextend(self):
        """ Extend the primary mesh with information about the degrees of freedom.

            TODO: Docs...

            """

        xmesh = Xmesh()

        # Construct a species map (dict mapping model species name to an integer index)
        species_map = self.speciesMap()

        # Initialize the function spaces and dof maps.
        for spec in self.listOfSpecies:

            species = self.listOfSpecies[spec]
            spec_name = species.name
            spec_index = species_map[spec_name]

            xmesh.function_space[spec_name] = dolfin.FunctionSpace(self.mesh, "Lagrange", 1)
            # vertex_to_dof_map provides a map between the vertex index and the dof.
            xmesh.vertex_to_dof_map[spec_name] = xmesh.function_space[spec_name].dofmap().dof_to_vertex_map(self.mesh)
            xmesh.vertex_to_dof_map[spec_name] = len(self.listOfSpecies) * xmesh.vertex_to_dof_map[spec_name] + spec_index
            xmesh.dof_to_vertex_map[spec_name] = xmesh.function_space[spec_name].dofmap().vertex_to_dof_map(self.mesh)

        xmesh.vertex = self.mesh.coordinates()
        self.xmesh = xmesh


    # Some utility routines to set initial conditions
    def scatter(self,spec_init,subdomains=None):
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
                vtx = np.random.randint(0, ltab)
                ind = table[vtx]
                self.u0[specindx,ind]+=1
                

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
            matrices = self.assemble()
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
        
        try:
            sd = self.sd
        except:
            sd = self.subdomainVector(self.subdomains)
       
        for species,K in stiffness_matrices.iteritems():
            
            rows,cols,vals = K.data()
            Kcrs = scipy.sparse.csr_matrix((vals,cols,rows))
            Kdok = Kcrs.todok()
            
            dof2vtx = xmesh.dof_to_vertex_map[species]
            
            for entries in Kdok.items():
                
                ind = entries[0]
                ir = ind[0]
                ij = ind[1]
                
                # Permutation to make the matrix ordering match that of sd, u0. (Dolfin dof -> URDME dof)
                ir = dof2vtx[ind[0]]
                ij = dof2vtx[ind[1]]
                
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


    def validate(self, urdme_solver_data):
        """ Validate the model data structures.

            validate should be called prior to writing the model to the solver input file,
            since the solvers themselves do very limited error checking of the input.

        """

        for spec_name,species in self.listOfSpecies.items():
            if 0 in self.species_to_subdomains[species]:
                raise ModelException("Subdomain number 0 is reserved. Please check your model.")

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
        self.vol = vol[1::len(self.listOfSpecies)]

        urdme_solver_data['vol'] = self.vol
        D = result['D']
        urdme_solver_data['D'] = D

        # Subdomain vector
        urdme_solver_data['sd'] = self.subdomainVector(self.subdomains)

        # Data vector. If not present in model, it defaults to a vector with all elements zero.
        data = np.zeros((1, self.mesh.getNumVoxels()))
        urdme_solver_data['data'] = data

        if not hasattr(self,'u0'):
            self.initializeInitialValue()

        urdme_solver_data['u0'] = self.u0

        tspan = np.asarray(self.tspan, dtype=np.float)
        urdme_solver_data['tspan'] = tspan

        # Vertex coordinates
        urdme_solver_data['p'] = self.mesh.getVoxels()

        # Connectivity matrix
        urdme_solver_data['K'] = self.connectivityMatrix()

        #rows,cols,vals = self.stiffness_matrices["MinD_m"].data()
        #SM = scipy.sparse.csr_matrix((vals,cols,rows))
        #urdme_solver_data["Kmindm"] = SM.tocsc()

        return urdme_solver_data


    def serialize(self, filename=None):
        """ Write the datastructures needed by the the core URDME solvers to a .mat input file. """

        urdme_solver_data = self.solverData()
        self.validate(urdme_solver_data)
        spio.savemat(filename, urdme_solver_data, oned_as='column')


    def connectivityMatrix(self):
        """ Assemble a connectivity matrix in CCS format. """
        
        fs = dolfin.FunctionSpace(self.mesh, "Lagrange", 1)
        trial_function = dolfin.TrialFunction(fs)
        test_function = dolfin.TestFunction(fs)
        a_K = -1*dolfin.inner(dolfin.nabla_grad(trial_function), dolfin.nabla_grad(test_function)) * dolfin.dx
        C = dolfin.assemble(a_K)
        rows, cols, vals = C.data()
        C = scipy.sparse.csr_matrix((vals, cols, rows))
        C = C.tocsc()
        return C

    def assemble(self):
        """  Assemble the mass and stiffness matrices using Dolfin.
            
            Returns: A dictionary containing two dictionaries, one for the stiffness matrices
            and one for the mass matrices. Those dictionaries has the species names as keys and
            the matrices are in CSR format.
            """
        
        if not hasattr(self, 'xmesh'):
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
        
        # Set up the forms
        for spec_name, species in self.listOfSpecies.items():
            
            # Find out what subdomains this species is active on
            subdomain_list = self.species_to_subdomains[species]
            weak_form_K[spec_name] = dolfin.inner(dolfin.nabla_grad(trial_functions[spec_name]), dolfin.nabla_grad(test_functions[spec_name]))*dolfin.dx
            weak_form_M[spec_name] = trial_functions[spec_name]*test_functions[spec_name]*dolfin.dx

        # Assemble the matrices
        for spec_name,species in self.listOfSpecies.items():
            stiffness_matrices[spec_name] = dolfin.assemble(weak_form_K[spec_name])
            # We cannot include the diffusion constant in the assembly, dolfin does not seem to deal well
            # with small diffusion constants (drops small elements)
            stiffness_matrices[spec_name] = species.diffusion_constant * stiffness_matrices[spec_name]
            mass_matrices[spec_name] = dolfin.assemble(weak_form_M[spec_name])
        
        
        return {'K':stiffness_matrices, 'M':mass_matrices}





class Mesh(dolfin.Mesh):
    """ A URDME mesh extends the Dolfin mesh class. """

    def __init__(self, mesh=None):
        dolfin.Mesh.__init__(self, mesh)

    def getNumVoxels(self):
        return self.num_vertices()

    def getVoxels(self):
        return self.coordinates()

    """  Wrappers around dolfins built-in simple geometries/meshes.

        These following methods will all give regular meshes that will produce discretizations that are
        equivalent to Cartesian grids.

    """

    @classmethod
    def unitIntervalMesh(cls, nx):
        mesh = dolfin.IntervalMesh(nx, 0, 1)
        return Mesh(mesh)

    @classmethod
    def IntervalMesh(cls, nx, a, b):
        mesh = dolfin.IntervalMesh(nx, a, b)
        return Mesh(mesh)

    @classmethod
    def unitSquareMesh(cls, nx, ny):
        """ Unit Square of with nx,ny points in the respective axes. """
        mesh = dolfin.UnitSquareMesh(nx, ny)
        return Mesh(mesh)

    @classmethod
    def SquareMesh(cls, L, nx, ny):
        """ Regular mesh of a square with side length L. """
        mesh = dolfin.RectangleMesh(0, 0, L, L, nx, ny)
        return Mesh(mesh)

    @classmethod
    def unitCubeMesh(cls, nx, ny, nz):
        """ Unit Square of with nx,ny points in the respective axes. """
        mesh = dolfin.UnitCubeMesh(nx, ny, nz)
        return Mesh(mesh)

    #@classmethod
    #def unitCircle(cls, nx,ny):
    #    """ Unit Square of with nx,ny points in the respective axes. """
    #    mesh = dolfin.UnitCircleMesh(nx,ny)
    #    return Mesh(mesh)

    #@classmethod
    #def unitSphere(cls, nx,ny):
    #    """ Unit Square of with nx,ny points in the respective axes. """
    #    mesh = dolfin.UnitSquareMesh(nx,ny)
    #    return Mesh(mesh)



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
    def read_dolfin_mesh(cls, filename=None):
        """ Import a mesh in Dolfins native .xml format """

        try:
            dolfin_mesh = dolfin.Mesh(filename)
            mesh = Mesh(mesh=dolfin_mesh)
            return mesh
        except Exception as e:
            raise MeshImportError("Failed to import mesh: " + filename+"\n" + str(e))


class Xmesh():
    """ Extended mesh object.

        Contains function spaces and dof mappings.
    """

    def __init__(self):
        self.coordinates = None
        self.function_space = {}
        self.vertex_to_dof_map = {}
        self.dof_to_vertex_map = {}

class URDMEResult(dict):
    """ Result object for a URDME simulation, extendes the dict object. """
    
    def __init__(self, model, filename=None):
        self.model = model
        self.sol = None
        self.U = None
        self.tspan = None
        if filename is not None:
            self.read_solution(filename)
    
    def __setattr__(self, k, v):
        if k in self.keys():
            self[k] = v
        elif not hasattr(self, k):
            self[k] = v
        else:
            raise AttributeError, "Cannot set '%s', cls attribute already exists" % ( k, )
    
    def __getattr__(self, k):
        if k == 'sol':
            ret = self.get('sol')
            if ret is None:
                return self.initialize_sol()
            return ret
        if k in self.keys():
            return self[k]
        raise AttributeError
    
    def _initialize_sol(self):
        """ Initialize the sol variable. """
        # Create Dolfin Functions for all the species
        sol = {}
        
        dims = self.U.shape
        numvox = self.model.mesh.getNumVoxels()
        # The result is loaded in dolfin Functions, one for each species and time point
        for i,spec in enumerate(self.model.listOfSpecies):
            
            species = self.model.listOfSpecies[spec]
            spec_name = species.name
            dof_to_vertex_map = self.model.xmesh.dof_to_vertex_map[spec]
            vertex_to_dof_map = self.model.xmesh.vertex_to_dof_map[spec]
            
            spec_sol = {}
            for j,time in enumerate(self.tspan):
                func = dolfin.Function(dolfin.FunctionSpace(self.model.mesh,"Lagrange",1))
                func_vector = func.vector()
                
                for voxel in range(numvox):
                    dof = voxel*len(self.model.listOfSpecies)+i
                    ix  = vertex_to_dof_map[voxel]
                    dolfvox = (ix-i)/len(self.model.listOfSpecies)
                    func_vector[dolfvox] = float(self.U[dof,j]/self.model.vol[voxel])
                
                spec_sol[time] = func
            
            sol[spec] = spec_sol
        self.sol = sol
        return sol

    def dumps(self,species,folder_name):
        """ Dump the trajectory of species to a collection of vtk files """
        self._initialize_sol()
        subprocess.call(["mkdir","-p",folder_name])
        func = dolfin.Function(dolfin.FunctionSpace(self.model.mesh,"Lagrange",1))
        func_vector = func.vector()
        file=dolfin.File(folder_name+"/trajectory.pvd")
        numvox = self.model.mesh.getNumVoxels()

        for i,time in enumerate(self.tspan):
            solvector = (self.sol[species][time]).vector()
            for dof in range(numvox):
                func_vector[dof] = solvector[dof]
            file << func



    def toXYZ(self, filename, file_format="ParaView"):
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

        if file_format == "VMD":
            outfile = open(filename, "w")
            filestr = ""
            for i, time in enumerate(self.tspan):
                number_of_atoms = numpy.sum(self.U[:, i])
                filestr += (str(number_of_atoms) + "\n" + "timestep " + str(i) + " time " + str(time) + "\n")
                for j, spec in enumerate(self.model.listOfSpecies):
                    for k in range(Ncells):
                        for mol in range(self.U[k * Mspecies + j, i]):
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

    def toCSV(self, filename):
        """ Dump the solution attached to a model as a .csv file. """
        #TODO: Make this work for 2D meshes with only two coordinates.

        if self.U is None:
            raise URDMEError("No solution found in the model.")

        dims = numpy.shape(self.U)
        Ndofs = dims[0]
        Mspecies = len(self.model.listOfSpecies)
        Ncells = Ndofs/Mspecies

        coordinates = self.model.mesh.getVoxels()
        coordinatestr = coordinates.astype(str)
        subprocess.call(["mkdir","-p", filename])
        for i, time in enumerate(self.tspan):
            outfile = open(filename + '/' + filename + str(i) + ".csv", "w")
            number_of_atoms = numpy.sum(self.U[:, i])
            filestr = "xcoord,ycoord,zcoord,radius,type\n"
            for j, spec in enumerate(self.model.listOfSpecies):
                for k in range(Ncells):
                    for mol in range(self.U[k * Mspecies + j, i]):
                        obj = self.model.listOfSpecies[spec]
                        reaction_radius = obj.reaction_radius
                        linestr = coordinatestr[k, 0] + "," + coordinatestr[k, 1] + "," + coordinatestr[k, 2] + "," + str(reaction_radius) + "," + str(j) + "\n"
                        filestr += linestr
            outfile.write(filestr)
            outfile.close()


    def read_solution(self, filename):

        resultfile = h5py.File(filename, 'r')

        U = resultfile['U']
        U = numpy.array(U)
        # This little hack makes U have the same structure as in the Matlab interface...
        dims = numpy.shape(U)
        U = U.reshape((dims[1], dims[0]))
        U = U.transpose()

        tspan = resultfile['tspan']
        tspan = numpy.array(tspan).flatten()
        resultfile.close()
        self.U = U
        self.tspan = tspan


class URDMESolver:
    """ Abstract class for URDME solvers. """
    
    @classmethod
    def get_solver_instance(cls, solver_name, model, solver_path=None):
        """ Lookup table for solvers regitered with pyurdme. """
        if solver_name == 'nsm':
            return NSMSolver(model,solver_path)
        elif solver_name == 'nem':
            return NEMSolver(model,solver_path)
        else:
            raise URDMEError("Unknown solver: {0}".format(solver_name))
            
    
    def __init__(self, model, solver_path=None):
        """ Constructor. """
        if not isinstance(model,URDMEModel):
          raise URDMEError("URDMEsolver constructors must take a URDMEModel as an argument.")
        if not hasattr(self, 'NAME'):
            raise URDMEError("")

        self.model = model
        self.is_compiled = False

        # Set URDME_ROOT. This will fail if URDME is not installed on the system.
        try:
            urdme_init = subprocess.check_output(['which','urdme_init']).strip()
            path = os.readlink(urdme_init)
            if not path.endswith('/urdme/bin/urdme_init'):
                raise Exception('path={0}\n'.format(path))
            self.URDME_ROOT = path.replace('bin/urdme_init','')
        except Exception as e:
            raise URDMEError("Could not determine the location of URDME.")
        
        if solver_path is None:
            self.URDME_BUILD = self.URDME_ROOT + '/build/'
        else:
            self.URDME_BUILD = solver_path + '/build/'
            os.environ['SOLVER_ROOT'] = solver_path
    
    
    def compile(self, model_file=None):
        """ Compile the model."""
        # Write the propensity file
        if os.path.isdir('.urdme'):
            shutil.rmtree('.urdme')
        
        try:
            os.mkdir('.urdme')
        except Exception as e:
            pass
        
        propfilename = self.model.name + '_pyurdme_generated_model'
        if model_file == None:
            propfilename = self.model.name + '_pyurdme_generated_model'
            self.createPropensityFile(file_name='.urdme/' + propfilename + '.c')
        else:
            subprocess.call(['cp', model_file, '.urdme/' + propfilename + '.c'])
        
        # Build the solver
        makefile = 'Makefile.' + self.NAME
        cmd = ['make', '-f', self.URDME_BUILD + makefile, 'URDME_ROOT=' + self.URDME_ROOT, 'URDME_MODEL=' + propfilename]
        if report_level >= 1:
            print "cmd: {0}\n".format(cmd)
        handle = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
        return_code = handle.wait()
        
        if return_code != 0:
            print handle.stdout.read()
            print handle.stderr.read()
            raise URDMEError("Compilation of solver failed")
        
        if report_level >= 1:
            print handle.stdout.read()
            print handle.stderr.read()

    
    
    def run(self, input_file=None, seed=None, report_level=0):
        """ Run one simulation of the model.
            
        Returns:
            URDMEResult object
        """
        if input_file is None:
            # Get temporary input and output files
            infile = tempfile.NamedTemporaryFile(delete=False)
            
            # Write the model to an input file in .mat format
            self.model.serialize(filename=infile)
            infile.close()
            infile_name = infile.name
        else:
            infile_name = input_file
        
        
        outfile = tempfile.NamedTemporaryFile(delete=False)
        outfile.close()
        
        # Execute the solver
        urdme_solver_cmd = ['.urdme/' + propfilename + '.' + self.NAME , infile_name , outfile.name]
        if seed is not None:
            urdme_solver_cmd.append(str(seed))
        if report_level >= 1:
            print 'cmd: {0}\n'.format(urdme_solver_cmd)
        handle = subprocess.Popen(urdme_solver_cmd)
        return_code = handle.wait()
        if return_code != 0:
            raise URDMEError("Solver execution failed")
        
        #Load the result from the hdf5 output file.
        try:
            result = URDMEResult(model, outfile.name)
            
            # Clean up
            if input_file is None:
                os.remove(infile.name)
            os.remove(outfile.name)
            
            result["Status"] = "Sucess"
            return result
        
        except Exception as e:
            exc_info = sys.exc_info()
            # Clean up
            if input_file is None:
                os.remove(infile.name)
            os.remove(outfile.name)
            raise exc_info[1], None, exc_info[2]    

    
    
    def createPropensityFile(self, file_name=None):
        """ Generate the C propensity file that is used to compile the URDME solvers.
            Only mass action propensities are supported. """
        
        template = open(os.path.abspath(os.path.dirname(__file__)) + '/data/propensity_file_template.c', 'r')
        propfile = open(file_name, "w")
        propfilestr = template.read()
        
        speciesdef = ""
        i = 0
        for S in self.model.listOfSpecies:
            speciesdef += "#define " + S + " " + "x[" + str(i) + "]" + "\n"
            i += 1
        
        propfilestr = propfilestr.replace("__DEFINE_SPECIES__", speciesdef)
        
        propfilestr = propfilestr.replace("__NUMBER_OF_REACTIONS__", str(self.model.getNumReactions()))
        
        # Make sure all paramters are evaluated to scalars before we write them to the file.
        self.model.resolveParameters()
        parameters = ""
        for p in self.model.listOfParameters:
            parameters += "const double " + p + " = " + str(self.model.listOfParameters[p].value) + ";\n"
        propfilestr = propfilestr.replace("__DEFINE_PARAMETERS__", str(parameters))
        
        # Reactions
        funheader = "double __NAME__(const int *x, double t, const double vol, const double *data, int sd)"
        
        funcs = ""
        funcinits = ""
        i = 0
        for R in self.model.listOfReactions:
            func = ""
            rname = self.model.listOfReactions[R].name
            func += funheader.replace("__NAME__", rname) + "\n{\n"
            if self.model.listOfReactions[R].restrict_to == None:
                func += "    return " + self.model.listOfReactions[R].propensity_function
                order = len(self.model.listOfReactions[R].reactants)
                if order == 2:
                    func += "/vol;"
                elif order == 0:
                    func += "*vol;"
                else:
                    func += ";"
            
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
                func += ")\n"
                func += "\treturn " + self.model.listOfReactions[R].propensity_function
                order = len(self.model.listOfReactions[R].reactants)
                if order == 2:
                    func += "/vol;"
                elif order == 0:
                    func += "*vol;"
                else:
                    func += ";"
                
                func += "\nelse"
                func += "\n\treturn 0.0;"
            
            
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
        
        TODO: Docs...
        
        After sucessful execution, urdme returns a dictionary, result, with the following members
        U:         the raw copy number output in a matrix with dimension (Ndofs, num_time_points)
        tspan:     the time span vector containing the time points that corresponds to the columns in U
        status:    Sucess if the solver executed without error
        stdout:    the standard ouput stream from the call to the core solver
        stderr:    the standard error stream from the call to the core solver
        
        """
    solver = URDMESolver.get_solver_instance(solver,model,solver_path)
    solver.compile(model_file)
    return solver.run(input_file, seed, report_level)



class MeshImportError(Exception):
    """ Exception to raise when encourntering and error importing a mesh. """
    pass

class URDMEError(Exception):
    pass

class ModelException(Exception):
    pass

class InvalidSystemMatrixException(Exception):
    pass



if __name__ == '__main__':
    """ Command line interface to URDME. Execute URDME given a model file. """

