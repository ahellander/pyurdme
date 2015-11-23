""" pyurdme model for a single species diffusing on the unit square. """

import dolfin
import mshr
import pyurdme
import numpy
import math
import sys

class SimpleCrowding(pyurdme.URDMEModel):
    """ Initial condition is a delta function at the center voxel.
        The solution should be a Gaussian, up to the point where
        the BC becomes important. """
    
    def __init__(self,name="diffusion",rho=0.1,nx=40):
        
    
        pyurdme.URDMEModel.__init__(self,name=name)

        D  = 0.5
        Rc = 1
        Rp = 5
        A = pyurdme.Species(name="A",reaction_radius=Rp,diffusion_constant=D)
        
        # Crowding agents
        B = pyurdme.Species(name="B",reaction_radius=Rc,diffusion_constant=0)
        
        L = 100
        V = L**2
        N = int(V*rho/(math.pi*Rc**2))
        self.add_species([A,B])

        # A (1x1) micro meter sidelength square
        domain = mshr.Rectangle(dolfin.Point(0,0),dolfin.Point(L,L))
        self.mesh = pyurdme.URDMEMesh(mesh=mshr.generate_mesh(domain,nx))

        # Place one A molecules in the voxel closest to the center of the square
        self.set_initial_condition_place_near({A:1},point=numpy.array([0.5,0.5])*L)
        
        # Distribute crowders randomly in the domain.
        self.set_initial_condition_scatter({B:N})
        self.timespan(numpy.linspace(0,50,1000))



def msd(result):
    A = result.get_species("A")
    
    factor, p = result.model.mesh.get_scaled_normalized_coordinates()
    
    dims = numpy.shape(A)
    msd = numpy.zeros((1,dims[0]))
    
    dist = numpy.linalg.norm(p,ord=2,axis=1)
    dist = numpy.power(dist,2)
    
    for time in range(dims[0]):
        msd[0,time]=numpy.dot(dist,A[time,:])/numpy.sum(A[time, :])
    
    msd[0,0] = 0

    return msd

if __name__ == '__main__':

    model = SimpleCrowding(rho=0.3,nx=40)
    
    N = 1
    
    from pyurdme import nsmsolver,nemsolver
    
    solver1 = nsmsolver.NSMSolver(model)
    solver2 = nemsolver.NEMSolver(model, report_level=2)
    import time
    tic = time.time()
    res = solver2.run()
    print "NEM time:",time.time() - tic
    A = res.get_species("A",0)
    B = res.get_species("B",0)
    print "Non-zeros in A:",float(len(numpy.flatnonzero(A)))/len(A)    #print B
    print "Non-zeros in B:",float(len(numpy.flatnonzero(B)))/len(B)
   
    exit(0)
    
    tic = time.time()
    for i in xrange(N):
        if i == 0:
            result = solver2.run()
            msd2 = msd(result)
        else:
            msd2 += msd(result)
        print i

    msd2 = msd2.flatten()/N
    print msd2

    print "Time nem", str(time.time()-tic)

    tspan = result.model.tspan

    from scipy import stats

    for i, val in enumerate(msd2):
        if val > 0.0:
            start = i
            break

    alpha2 = stats.linregress(numpy.log10(tspan[start::]),numpy.log10(msd2[start::]))
    print "NEM: alpha: ", alpha2[0], "D:", pow(10,alpha2[1])/4


