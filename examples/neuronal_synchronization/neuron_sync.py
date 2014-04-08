#!/usr/bin/env python
""" pyURDME model file for Giovanni's three state  model.
cite: Schroder JBR (2012)
"""

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
        self.mesh = pyurdme.URDMEMesh.IntervalMesh(nx=5, a=0, b=5)
        #self.mesh.addPeriodicBoundaryCondition(PeriodicBoundary1D(a=0, b=5))


        alpha = pyurdme.Parameter(name="alpha" ,expression=0.1) 
        KI = pyurdme.Parameter(name="KI" ,expression=1.0) 
        KM = pyurdme.Parameter(name="KM" ,expression=0.5) 
        ks = pyurdme.Parameter(name="ks" ,expression=0.417) 
        vd = pyurdme.Parameter(name="vd" ,expression=1.167) 
        KD = pyurdme.Parameter(name="KD" ,expression=0.13) 
        k1 = pyurdme.Parameter(name="k1" ,expression=0.417) 
        k2 = pyurdme.Parameter(name="k2" ,expression=0.5) 
        n  = pyurdme.Parameter(name="n"  ,expression=4)
        #vm = normrnd(0.416,0.04,N,1);
        vm = pyurdme.Parameter(name="vm" ,expression=0.416)
        L  = pyurdme.Parameter(name="L"  ,expression=0.0)

        
        self.addParameter([ v1, K1, n, v2, K2, k3, v4, K4, k5, v6, K6, k7, v8, K8, vc, Kc, K, L ]) 

        # Convert ODEs into stochastic reactions
        #x1(t,osc1) = x1(t-1,osc1) + dt * ( v(t-1,osc1) * KI^n / (KI^n + x3(t-1,osc1)^n) - vm(osc1) * x1(t-1,osc1)/(KM + x1(t-1,osc1)) );
        #x2(t,osc1) = x2(t-1,osc1) + dt * ( ks * x1(t-1,osc1) - vd * x2(t-1,osc1) / (KD + x2(t-1,osc1)) - k1 * x2(t-1,osc1) + k2 * x3(t-1,osc1) );
        #x3(t,osc1) = x3(t-1,osc1) + dt * ( k1 * x2(t-1,osc1) - k2 * x3(t-1,osc1) );
        #v(t,osc1) = -Diff_Cost/(ell*ell)*dt*(right_v + left_v - 2*v(t-1, osc1)) + 0.83 + L(osc1) + alpha(osc1) * x1(t-1,osc1) - v_deg*v(t-1, osc1);
                        
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
        NANO_MOLAR = 10000.0
        if nano_molar is not None:
            NANO_MOLAR = nano_molar
        # Initial Conditions 
        self.distributeUniformly({X:numpy.floor(2.5*NANO_MOLAR)})
        self.distributeUniformly({Y:numpy.floor(2.5*NANO_MOLAR)})
        self.distributeUniformly({Z:numpy.floor(2.5*NANO_MOLAR)})
        self.distributeUniformly({V:numpy.floor(1.0*NANO_MOLAR)})

        # Set time range and sample points
        self.timespan(range(201))


if __name__=="__main__":
    """ Dump model to a file. """
                     
    model = neuron_sync_1D()
    result = pyurdme.urdme(model, report_level=1)
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
