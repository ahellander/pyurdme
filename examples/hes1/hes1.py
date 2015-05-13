import matplotlib.pyplot as plt
import os.path
import pyurdme
import dolfin
import numpy
import scipy

class Nucleus(dolfin.SubDomain):
    def inside(self,x,on_boundary):
        return dolfin.between(x[0]**2+x[1]**2+x[2]**2,(0,3.**2))

class MeshSize(pyurdme.URDMEDataFunction):
    def __init__(self,mesh):
        pyurdme.URDMEDataFunction.__init__(self,name="MeshSize")
        self.mesh = mesh
        self.h = mesh.get_mesh_size()

    def map(self,x):
        ret = self.h[self.mesh.closest_vertex(x)]
        return ret

class Hes1(pyurdme.URDMEModel):
    def __init__(self,model_name="hes1"):
        pyurdme.URDMEModel.__init__(self, model_name)

        #Species
        Pf = pyurdme.Species(name="Pf",diffusion_constant=0.,dimension=3)
        Po = pyurdme.Species(name="Po",diffusion_constant=0.,dimension=3)
        mRNA = pyurdme.Species(name="mRNA",diffusion_constant=6.e-1,dimension=3)
        protein = pyurdme.Species(name="protein",diffusion_constant=6.e-1,dimension=3)

        self.add_species([Pf,Po,mRNA,protein])

        #Domains
        basedir = os.path.dirname(os.path.abspath(__file__))
        self.mesh = pyurdme.URDMEMesh.read_mesh(basedir+"/mesh/cell_coarse.msh")
        
        volumes = dolfin.MeshFunction("size_t",self.mesh,0)
        volumes.set_all(2)
        nucleus = Nucleus()
        nucleus.mark(volumes,1)
        volumes[self.mesh.closest_vertex([0,0,0])] = 3

        self.add_subdomain(volumes)

        h = self.mesh.get_mesh_size()
        self.add_data_function(MeshSize(self.mesh))

        #Parameters
        k1 = pyurdme.Parameter(name="k1",expression=1.e9)
        k2 = pyurdme.Parameter(name="k2",expression=0.1)
        alpha_m = pyurdme.Parameter(name="alpha_m",expression=3.)
        alpha_m_gamma = pyurdme.Parameter(name="alpha_m_gamma",expression=3./30.)
        alpha_p = pyurdme.Parameter(name="alpha_p",expression=1.)
        mu_m = pyurdme.Parameter(name="mu_m",expression=0.015)
        mu_p = pyurdme.Parameter(name="mu_p",expression=0.043)
        
        self.add_parameter([k1,k2,alpha_m,alpha_m_gamma,alpha_p,mu_m,mu_p])

        #Domains markers
        nucleus = [1]
        cytoplasm = [2]
        promoter_site = [3]

        #Reactions
        R1 = pyurdme.Reaction(name="R1",reactants={Pf:1,protein:1},products={Po:1},massaction=True,rate=k1, restrict_to=promoter_site)
        R2 = pyurdme.Reaction(name="R2",reactants={Po:1},products={Pf:1,protein:1},massaction=True,rate=k2, restrict_to=promoter_site)
        R3 = pyurdme.Reaction(name="R3",reactants={Pf:1},products={Pf:1,mRNA:1},massaction=True,rate=alpha_m, restrict_to=promoter_site)
        R4 = pyurdme.Reaction(name="R4",reactants={Po:1},products={Po:1,mRNA:1},massaction=True,rate=alpha_m_gamma, restrict_to=promoter_site)
        R5 = pyurdme.Reaction(name="R5",reactants={mRNA:1},products={mRNA:1,protein:1},massaction=True,rate=alpha_p,restrict_to=cytoplasm)
        R6 = pyurdme.Reaction(name="R6",reactants={mRNA:1},products={},massaction=True,rate=mu_m)
        R7 = pyurdme.Reaction(name="R7",reactants={protein:1},products={},massaction=True,rate=mu_p)

        self.add_reaction([R1,R2,R3,R4,R5,R6,R7])

        #Restrict to promoter_site
        self.restrict(Po,promoter_site)
        self.restrict(Pf,promoter_site)

        #Distribute molecules over the mesh
        #self.set_initial_condition_place_near({Pf:1},[0,0,0])
        self.set_initial_condition_scatter({Pf:1},promoter_site)
        self.set_initial_condition_scatter({protein:60},cytoplasm)
        self.set_initial_condition_scatter({mRNA:10},nucleus)

        self.timespan(range(1200))

if __name__=="__main__":
    model = Hes1(model_name="hes1")
    result = model.run(report_level=1)

    protein = result.get_species("protein")
    proteinsum = numpy.sum(protein,axis=1)
    plt.plot(model.tspan,proteinsum,'r')
    mRNA = result.get_species("mRNA")
    mRNAsum=numpy.sum(mRNA[:],axis=1)
    plt.plot(model.tspan,mRNAsum,'b')
    plt.show()

    #print 'Writing species "protein" to folder "proteinOut"'
    #result.export_to_vtk(species='protein',folder_name='proteinOut')
