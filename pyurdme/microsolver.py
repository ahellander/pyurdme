import pyurdme
import dolfin
import os
import tempfile
import subprocess
import shutil
import numpy
from pyurdme import get_N_HexCol
import json
import uuid
import h5py

try:
    # This is only needed if we are running in an Ipython Notebook
    import IPython.display
except:
    pass

class MMMSSolver(pyurdme.URDMESolver):
    """ Micro solver class.TODO: Add short description. """
    
    NAME = 'mmms'
    
    def __init__(self, model, solver_path=None, report_level=0, model_file=None, sopts=None):
        pyurdme.URDMESolver.__init__(self,model,solver_path,report_level,model_file,sopts)

        self.solver_name = 'hybrid'
        self.solver_path = ''
        self.urdme_infile_name = ''
    
    def __getstate__(self):
        """ TODO: Implement"""
    
    def __getstate__(self):
        """ TODO: Implement"""
    
    def __del__(self):
        """ Remove the temporary output folder """
        if self.delete_infile:
            try:
                os.remove(self.infile_name)
            except OSError as e:
                print "Could not delete '{0}'".format(self.infile_name)
    #shutil.rmtree(self.outfolder_name)
    

    def _write_mesh_file(self, filename=None):
        """ Write the mesh data to a HDF5 file that the mmmms solver understands. """
        
        meshfile = h5py.File(filename,"w")
        meshfile.create_group("mesh")
        meshfile.create_group("boundarymesh")
        grp = meshfile["mesh"] 
    
        mesh = self.model.mesh
        mesh.init()
        
        cells = mesh.cells()
        grp.create_dataset("t", data = cells)
        vertices = mesh.coordinates()
        grp.create_dataset("p",data=vertices)

        # Create the bounday mesh triangle entities 
        boundary_mesh = dolfin.BoundaryMesh(mesh, "exterior")

        #print numpy.shape(boundary_mesh.coordinates())
       # print numpy.shape(boundary_mesh.cells())
    
        # Vertex map
        vm = boundary_mesh.entity_map(0).array()
        #print vm
        
        bndtri= []
        for tri in boundary_mesh.cells():
            bndtri.append([vm[v] for v in tri])

        grp.create_dataset("boundaryfacets",data=numpy.array(bndtri))

        # Datastructure, vertex->tetrahedron
        mesh.init(0,0)
        f = dolfin.MeshEntity(mesh, 0,0)
        v2t = []
        for v in dolfin.vertices(mesh):
            #for c in dolfin.cells(v):
            v2t.append([ci.index() for ci in dolfin.cells(v)])
                #print v.index(), c.index()
        v2t = numpy.array(v2t)        
        lmax = 0;
        for l in numpy.array(v2t):
            if len(l) > lmax:
                lmax = len(l)
        temp = -1*numpy.ones((len(v2t),lmax))

        for i,l in enumerate(v2t):
            for j,v in enumerate(l):
                temp[i,j] = v;

        grp.create_dataset("vertex2cells",data=temp)       

        #for cell in dolfin.cells(mesh):
        #    for face in dolfin.faces(cell):
        #        print [v for v in face.entities(0)] 
        
        #print cells
        #D = mesh.topology().dim()
        #fct = dolfin.facets(mesh)

        #edges = []
        #for facet in fct:
        #   edges.append([v for v in facet.entities(D)])
        #print numpy.array(edges)

        #grp.create_dataset("edges",data=numpy.array(edges))
        

    def create_input_file(self, filename):
    
        input_file = open(filename,'w')
        input_file.write("NAME "+self.model.name+"\n")
        
        
        # Write the model dimension
        (np,dim) = numpy.shape(self.model.mesh.coordinates())
        input_file.write("DIMENSION {0}\n".format(dim))
        input_file.write("BOUNDARY 0 1e-6 0 1e-6 0 1e-6\n")
        
        self.model.resolve_parameters()
        params = ""
        for i, pname in enumerate(self.model.listOfParameters):
            P = self.model.listOfParameters[pname]
            params += "PARAMETER {0} {1}\n".format(pname, str(P.value))
        input_file.write(params)

        speciesdef = ""
        
        initial_data = self.model.u0
        spec_map = self.model.get_species_map()
        
        for i, sname in enumerate(self.model.listOfSpecies):
            S = self.model.listOfSpecies[sname]
            sum_mol = numpy.sum(self.model.u0[spec_map[sname],:])
            speciesdef += "SPECIES {0} {1} {2} {3} {4}\n".format(sname, str(S.diffusion_constant), str(S.reaction_radius), sum_mol, "MICRO")
            
        input_file.write(speciesdef)
    
        reacstr = ""
        for i, rname in enumerate(self.model.listOfReactions):
            R = self.model.listOfReactions[rname]
            reacstr += "REACTION "

            for j, reactant in enumerate(R.reactants):
                reacstr += str(reactant)+" "
            
            reacstr += "> "
            
            for j, product in enumerate(R.products):
                reacstr += str(product)+" "
            reacstr += " {0}\n".format(str(R.marate.value))
            
        input_file.write(reacstr)

        input_file.write("T {0}\n".format(str(self.model.tspan[-1])))
        nint = len(self.model.tspan)
        input_file.write("NINT {0}\n".format(str(int(nint))))

    def run(self, number_of_trajectories=1, seed=None, input_file=None, loaddata=False):
        """ Run one simulation of the model.
            
            number_of_trajectories: How many trajectories should be run.
            seed: the random number seed (incimented by one for multiple runs).
            input_file: the filename of the solver input data file .
            loaddata: boolean, should the result object load the data into memory on creation.
            
            Returns:
            URDMEResult object.
            or, if number_of_trajectories > 1
            a list of URDMEResult objects
        """
        
        # Create the urdme input file that contains connectivity matrices etc
        if self.urdme_infile_name is None or not os.path.exists(self.urdme_infile_name):
            # Get temporary input and output files
            urdme_infile = tempfile.NamedTemporaryFile(delete=False, dir=os.environ.get('PYURDME_TMPDIR'))
                
            # Write the model to an input file in .mat format
            self.serialize(filename=urdme_infile, report_level=self.report_level)
            urdme_infile.close()
            self.urdme_infile_name = urdme_infile.name
            self.delete_infile = True

        if not os.path.exists(self.urdme_infile_name):
            raise URDMEError("input file not found.")
        
        # Generate the input file containing the microsolver specific 
        infile = tempfile.NamedTemporaryFile(delete=False, dir=os.environ.get('PYURDME_TMPDIR'))
        self.infile_name = infile.name
        self.create_input_file(infile.name)

        #  print infile.read()
        with open('test.txt','w') as fh:
            fh.write(infile.read())

        infile.close()

        # Create the mesh input file
        mesh_infile = tempfile.NamedTemporaryFile(delete=False, dir=os.environ.get('PYURDME_TMPDIR'))
        self._write_mesh_file(mesh_infile.name)
        self.mesh_infile_name=mesh_infile.name
        mesh_infile.close()

        if not os.path.exists(self.mesh_infile_name):
            raise URDMEError("Mesh input file not found.")
        
        # Generate output file
        outfile = tempfile.NamedTemporaryFile(delete=False, dir=os.environ.get('PYURDME_TMPDIR'))
        outfile.close()        
        
        if self.report_level > 2:
            print model_str
        if number_of_trajectories > 1:
            result_list = []

        solver_str=os.path.dirname(__file__)+"/mmms/bin/mmms"

        solver_cmd = [solver_str,self.infile_name, self.urdme_infile_name,self.mesh_infile_name, outfile.name]
	print solver_cmd
        handle = subprocess.Popen(solver_cmd)#, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        handle.wait()
        
        
        try:
            result = MICROResult(self.model,outfile.name)
            return result
        except IOError,e:
            print handle.stderr.read()
            raise

class MICROResult():

    def __init__(self,model, filename=None):
        self.filename=filename
        self.model = model
    
    def get_particles(self,species, time_index):
        """ Return a dict with the unique ids and positions of all particles of type species
            at time point time. """

        file = h5py.File(self.filename,'r')
        return {
                'unique_ids':numpy.array(file.get("Trajectories/0/Type_{0}/unique_ids_{1}".format(species, time_index))),
                'positions':numpy.array(file.get("Trajectories/0/Type_{0}/positions_{1}".format(species,time_index)))
                }

    def get_species(self,spec_name, time_index):
        """ Return copy number of spec_name in each voxel of the mesh. This function mimics the
            functionality of the mesoscopic solvers, so we insert particles in the voxels
            based on their position.  """

        return None
        #with open(self.output_folder_name+"/0_pos.txt", 'r') as fh:
        #    print fh.read()

    def get_summary_statistic(self, species, time_indices=None):
        """ Return the sum of molecules of a species for set of time points.
            If the result file contains multiple trajectories, then the 
            mean is taken over the realizations. 

            TODO: Implement mean value. 

        """
        
        if time_indices == None:
            tind = range(len(self.model.tspan))

        num_mol = []
        for ti in tind:
            r = self.get_particles(species, ti)
            (nm,dim) = numpy.shape(r['unique_ids'])
            num_mol.append(nm)

        return numpy.array(num_mol)

    def _export_to_particle_js(self,species,time_index, colors=None):
        """ Create a html string for displaying the particles as small spheres. """
        import random
        with open(os.path.dirname(pyurdme.__file__)+"/data/three.js_templates/particles.html",'r') as fd:
            template = fd.read()
        
        
        x=[]
        y=[]
        z=[]
        c=[]
        radius = []
        
        if colors == None:
            colors =  get_N_HexCol(len(species))
        
 
        spec_map = self.model.get_species_map()

        for j,spec in enumerate(species):
            spec_ind=spec_map[spec]
            particles = self.get_particles(spec_ind, time_index)
            vtx = particles['positions']
            maxvtx = numpy.max(numpy.amax(vtx,axis=0))
            factor = 1/maxvtx
            vtx = factor*vtx
            centroid = numpy.mean(vtx,axis=0)
            # Shift so the centroid is now origo
            normalized_vtx = numpy.zeros(numpy.shape(vtx))
            for i,v in enumerate(vtx):
                normalized_vtx[i,:] = v - centroid
            vtx=normalized_vtx

            for v in vtx:
                    x.append(v[0])
                    y.append(v[1])
                    z.append(v[2])
                    c.append(colors[j])
                    radius.append(0.01)
        
        template = template.replace("__X__",str(x))
        template = template.replace("__Y__",str(y))
        template = template.replace("__Z__",str(z))
        template = template.replace("__COLOR__",str(c))
        template = template.replace("__RADIUS__",str(radius))
        
        return template


    def display_particles(self,species, time_index):
        hstr = self._export_to_particle_js(species, time_index)
        displayareaid=str(uuid.uuid4())
        hstr = hstr.replace('###DISPLAYAREAID###',displayareaid)
        
        html = '<div id="'+displayareaid+'" class="cell"></div>'
        IPython.display.display(IPython.display.HTML(html+hstr))

    def __del__(self):
        """ Deconstructor. """
            #   if not self.data_is_loaded:
        try:
            # Clean up data file
            os.remove(self.filename)
        except OSError as e:
            #print "URDMEResult.__del__: Could not delete result file'{0}': {1}".format(self.filename, e)
            pass




