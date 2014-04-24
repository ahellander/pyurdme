#!/usr/bin/env python
import pyurdme
from examples.cylinder_demo.cylinder_demo3D import cylinderDemo3D
from examples.simple_diffusion import simple_diffusion
import pyurdme.nsmsolver

import numpy
import unittest

class SimpleDiffusion(pyurdme.URDMEModel):
    """ Initial condition is a delta function at the center voxel.
        The solution should be a Gaussian, up to the point where
        the BC becomes important. """
    
    def __init__(self):
        
        pyurdme.URDMEModel.__init__(self,name="simple_diffusion")
        
        D = 0.01
        A = pyurdme.Species(name="A",diffusion_constant=D,dimension=2)
        
        self.addSpecies([A])
        
        # A unit square
        self.mesh = pyurdme.URDMEMesh.unitSquareMesh(40,40)
        
        # Place the A molecules in the voxel nearest the center of the square
        self.placeNear({A:10000},point=[0.5,0.5])
        
        self.timespan(numpy.linspace(0,5,200))




class TestSolverFunctions(unittest.TestCase):

    def setUp(self):
        self.model = SimpleDiffusion()

    def test_same_seed(self):
        """ Test that the output is the same if the same seed is used.  """
        solver = pyurdme.nsmsolver.NSMSolver(self.model)
        result1 = solver.run(seed=1)
        result2 = solver.run(seed=1)
        A1 = result1.getSpecies("A")
        A2 = result2.getSpecies("A")
        
        self.assertFalse((A1-A2).any())

    def test_different_seeds(self):
        """ Test that the output is different if different seeds are used. """
        solver = pyurdme.nsmsolver.NSMSolver(self.model)
        result1 = solver.run(seed=1)
        result2 = solver.run(seed=100)
        A1 = result1.getSpecies("A")
        A2 = result2.getSpecies("A")
        self.assertFalse((A1-A2).any())


    def test_default_seeds(self):
        """ Test that the output is different if no seed is given (default set on C level). """
        solver = pyurdme.nsmsolver.NSMSolver(self.model)
        result1 = solver.run()
        result2 = solver.run()
        A1 = result1.getSpecies("A")
        A2 = result2.getSpecies("A")
        self.assertFalse((A1-A2).any())


if __name__ == '__main__':
    unittest.main()


