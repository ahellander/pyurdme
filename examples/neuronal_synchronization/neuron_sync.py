#!/usr/bin/env python
""" pyURDME model file for the Locke 2008 model. """

import os
import pyurdme
import dolfin
import math

import matplotlib.pyplot as plt
import numpy


# Sub domain for Periodic boundary condition
class PeriodicBoundary1D(dolfin.SubDomain):
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



class neuron_sync_1D(pyurdme.URDMEModel):

    def __init__(self, model_name="neuron_sync", diffusion=0.1, nano_molar=None):
        pyurdme.URDMEModel.__init__(self,model_name)

        # Species
        X  = pyurdme.Species(name="X",  diffusion_constant=0.0)
        Y  = pyurdme.Species(name="Y", diffusion_constant=0.0)
        Z  = pyurdme.Species(name="Z",  diffusion_constant=0.0)
        V  = pyurdme.Species(name="V", diffusion_constant=diffusion)
        
        self.addSpecies([X, Y, Z, V])
   
        # Create Mesh
        self.mesh = pyurdme.Mesh.IntervalMesh(nx=5, a=0, b=5)
        self.mesh.addPeriodicBoundaryCondition(PeriodicBoundary1D(a=0, b=5))

        v1 = pyurdme.Parameter(name="v1" ,expression=6.8355) #nM/h
        K1 = pyurdme.Parameter(name="K1" ,expression=2.7266) #nM
        n  = pyurdme.Parameter(name="n"  ,expression=5.6645)
        v2 = pyurdme.Parameter(name="v2" ,expression=8.4297) #nM/h
        K2 = pyurdme.Parameter(name="K2" ,expression=0.2910) #nM
        k3 = pyurdme.Parameter(name="k3" ,expression=0.1177) #1/h
        v4 = pyurdme.Parameter(name="v4" ,expression=1.0841) #nM/h
        K4 = pyurdme.Parameter(name="K4" ,expression=8.1343) #nM
        k5 = pyurdme.Parameter(name="k5" ,expression=0.3352) #1/h
        v6 = pyurdme.Parameter(name="v6" ,expression=4.6645) #nM/h
        K6 = pyurdme.Parameter(name="K6" ,expression=9.9849) #nM
        k7 = pyurdme.Parameter(name="k7" ,expression=0.2282) #1/h
        v8 = pyurdme.Parameter(name="v8" ,expression=3.5216) #nM/h
        K8 = pyurdme.Parameter(name="K8" ,expression=7.4519) #nM
        vc = pyurdme.Parameter(name="vc" ,expression=6.7924) #nM/h
        Kc = pyurdme.Parameter(name="Kc" ,expression=4.8283) #nM
        K  = pyurdme.Parameter(name="K"  ,expression=1.0)
        L  = pyurdme.Parameter(name="L"  ,expression=0.0)
        
        self.addParameter([ v1, K1, n, v2, K2, k3, v4, K4, k5, v6, K6, k7, v8, K8, vc, Kc, K, L ]) 

        # Convert ODEs into stochastic reactions
        # since we have explicit diffusion F_i = V_i
        # (1) dX/dt = v1*K1^n/(K1^n + Z^n) - v2*X/(K2 + X) + vc*K*V/(Kc + K*V) + L
        # (2) dY/dt = k3*X - v4*Y/(K4 + Y)
        # (3) dZ/dt = k5*Y - v6*Z/(K6 + Z)
        # (4) dV/dt = -D\nambla^2*V + k7*X - v8*V/(K8 + V)

        # Reactions
        R1c = pyurdme.Reaction(name="R1c", reactants={}, products={X:1}, 
            propensity_function="v1*pow(K1,n)/(pow(K1,n) + pow(Z/vol,n))*vol + vc*K*V/(Kc + K*V/vol) + L*vol")
        R1d = pyurdme.Reaction(name="R1d", reactants={X:1}, products={}, 
            propensity_function="v2*X/(K2 + X/vol)")
        R2c = pyurdme.Reaction(name="R2c", reactants={}, products={Y:1}, 
            propensity_function="k3*X")
        R2d = pyurdme.Reaction(name="R2d", reactants={Y:1}, products={}, 
            propensity_function="v4*Y/(K4 + Y/vol)")
        R3c = pyurdme.Reaction(name="R3c", reactants={}, products={Z:1}, 
            propensity_function="k5*Y")
        R3d = pyurdme.Reaction(name="R3d", reactants={Z:1}, products={}, 
            propensity_function="v6*Z/(K6 + Z/vol)")
        R4c = pyurdme.Reaction(name="R4c", reactants={}, products={V:1}, 
            propensity_function="k7*X")
        R4d = pyurdme.Reaction(name="R4d", reactants={V:1}, products={}, 
            propensity_function="v8*V/(K8 + V/vol)")

        self.addReaction([R1c, R1d, R2c, R2d, R3c, R3d, R4c, R4d])
        
        # Set the Stochasticty
        NANO_MOLAR = 100.0
        if nano_molar is not None:
            NANO_MOLAR = nano_molar
        # Initial Conditions 
        #self.distributeUniformly({X:numpy.floor(2.5*NANO_MOLAR)})
        #self.distributeUniformly({Y:numpy.floor(2.5*NANO_MOLAR)})
        #self.distributeUniformly({Z:numpy.floor(2.5*NANO_MOLAR)})
        self.distributeUniformly({V:numpy.floor(1.0*NANO_MOLAR)})

        # Set time range and sample points
        self.timespan(range(101))


if __name__=="__main__":
    """ Dump model to a file. """
                     
    model = neuron_sync_1D()
    result = pyurdme.urdme(model)
    print result

    x_vals = model.mesh.coordinates()[:, 0]
    tspan = result.tspan
    #X = result.getSpecies("V")
    #plt.plot(tspan, X[:,0], tspan, X[:,1], tspan, X[:,2])
    X = result.getSpecies("X")
    Y = result.getSpecies("Y")
    Z = result.getSpecies("Z")
    V = result.getSpecies("V")
    plt.plot(tspan, X[:,0], tspan, Y[:,0], tspan, Z[:, 0], tspan, V[:,0])
    plt.legend(('X','Y','Z', 'V'))
    plt.show()
