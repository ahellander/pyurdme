#!/usr/bin/env python
import pyurdme
from examples.mincde.mincde import mincde
import scipy.fftpack

import numpy
import pickle
import unittest



class TestMinCDE(unittest.TestCase):


    def test_mincde_oscillation_period(self):
        """ Check that the MinCDE model is producing oscillation of the right period. """
        model = mincde()
        result = model.run()
        mindm = result.get_species("MinD_m")
        mindmsum = numpy.sum(mindm[:,idx],axis=1)
        mindfft = scipy.fftpack.fft(mindmsum[200:]-mean(mindmsum[200:]))
        N = len(mindfft)
        T = model.tspan[-1] - model.tspan[200]
        mindpsd = numpy.abs(mindfft[:floor((N-1)/2)]) 
        mindfreq = numpy.arange(len(mindpsd), dtype=float)/T
        mind_max_period = 1/mindfreq[1+numpy.argmax(mindpsd[1:])]
        self.assertTrue(mind_max_period > 60)
        self.assertTrue(mind_max_period < 70)


if __name__ == '__main__':
    unittest.main()


