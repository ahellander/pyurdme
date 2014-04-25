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
        self.mesh = pyurdme.URDMEMesh.unitSquareMesh(10,10)
        
        # Place the A molecules in the voxel nearest the center of the square
        #self.placeNear({A:10000},point=[0.5,0.5])
        self.scatter({A:10000})
        self.timespan(numpy.linspace(0,5,200))



class TestSolverFunctionality(unittest.TestCase):

    def setUp(self):
        self.model = SimpleDiffusion()
    
    
    def test_solver_io(self):
        """ Test that the initial value in the solver output file is the same as the input initial value. """
        model = SimpleDiffusion()
        result = model.run()
        A = result.getSpecies("A",0)
        self.assertFalse((A-model.u0).any())
    
    def test_no_events(self):
        """ Test that nothing happens if the diffusion is set to zero. """
        model = SimpleDiffusion()
        model.listOfSpecies["A"].diffusion_constant = 0.0
        result = model.run()
        A = result.getSpecies("A")
        self.assertFalse((numpy.mean(A,axis=0)-model.u0).any())
    
    
    def test_same_seed(self):
        """ Test that the output is the same if the same seed is used, edxplicit solver creation  """
        solver = pyurdme.nsmsolver.NSMSolver(self.model)
        result1 = solver.run(seed=1)
        result2 = solver.run(seed=1)
        self.assertEqual(result1,result2)
    
    def test_same_seed2(self):
        """ Test that the output is the same if the same seed is used, use model.run()  """
    
        result1 = self.model.run(seed=1)
        result2 = self.model.run(seed=1)
        self.assertEqual(result1,result2)
    
    def test_different_seeds(self):
        """ Test that the output is different if different seeds are used. """
        solver = pyurdme.nsmsolver.NSMSolver(self.model)
        result1 = solver.run(seed=1)
        result2 = solver.run(seed=100)
        self.assertNotEqual(result1,result2)


    def test_default_seed(self):
        """ Test that the output is different if no seed is given (default set on C level). """
        solver = pyurdme.nsmsolver.NSMSolver(self.model)
        result1 = solver.run()
        result2 = solver.run()
        self.assertNotEqual(result1,result2)


    def test_run_aggregate(self):
        """ cccxxxyyy """


if __name__ == '__main__':
    unittest.main()


