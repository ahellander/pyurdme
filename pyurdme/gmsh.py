""" Wrappers for selected Gmsh functionality. """
import subprocess

class GmshGeo():
	""" Represent a Gmsh geometry. Wrapper arounf a
	    Gmesh .geo script file. """			
	def __init__(self,geofile=None):
		self.geofile = geofile
		
				
class GmshMesh():

	def __init__(self,meshfile=None):
		importmesh(meshfile)	
		
	def refine(self):
		""" Refine the mesh. """	
	
def imoportmesh(meshfile):
	""" Import a Gmsh .msh file. """	
	# TBD. 


def meshinit(geom):
	""" Call Gmsh and initialize a mesh from the .geo file. """
	# TBD	
	subprocess.call(['gmsh -m '])

def toDolfin(mesh):
	""" Serialize a Gmsh mesh to Dolfin format"""
	# TBD.  	
