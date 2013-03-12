""" Wrappers for selected Gmsh functionality. """
import subprocess
import numpy as np

class gmshGeometry():
    """ Represent a Gmsh geometry. Wrapper around a Gmesh .geo file. """
    def __init__(self,file=None):
        geomfile = open(file)
        self.geometry_file = geomfile.read()
        geomfile.close()

class GmshMesh():

    def __init__(self,meshfile=None):
        importmesh(meshfile)

        def refine(self):
            """ Refine the mesh. """
            
def importMesh(meshfile):
    """ Import a Gmsh .msh file. """
    meshf = open(meshfile,'r')
    lines = meshf.readlines()
    
    # Read number of nodes
    for lineno in range(len(lines)):
        if lines[lineno]=="$Nodes\n":
            break
    lineno=lineno+1
    number_of_nodes=int(lines[lineno])

    # Read the node values
    p = np.zeros((3,number_of_nodes))
    # !!!! IT IS READING ONE NODE TOO LITTLE.
    for i in range(number_of_nodes):
        line = lines[lineno+i+1]
        vtx = np.asmatrix(line[:-1])
        p[0,i]=vtx[0,1]
        p[1,i]=vtx[0,2]
        p[2,i]=vtx[0,3]

    lineno = lineno+number_of_nodes
    if lines[lineno]!="$EndNodes":
        print lines[lineno]
        

        

def meshinit(geometry):
    # Write temporary geometry file
    geometry_file_name='.gmsh_temporary_geometry_file.geo'
    geometry_file=open(geometry_file_name,'w')
    geometry_file.write(geometry.geometry_file)
    geometry_file.close()
    
    mesh_outfile_name = '.gmsh_temporary_mesh_file.msh'
    mesh_outfile = open(mesh_outfile_name,'w')
    mesh_outfile.close()
    
    subprocess.call(['gmsh','-3','-optimize',geometry_file_name,'-o',mesh_outfile_name])

    importMesh(mesh_outfile_name)
    
    # Clean up
    #subprocess.call(['rm','-rf',geometry_file_name])
    #subprocess.call(['rm','-rf',mesh_outfile_name])

    
def toDolfin(mesh):
    """ Serialize a Gmsh mesh to Dolfin format. """
    # TBD.  
