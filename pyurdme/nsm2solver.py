""" NSM2 solver. """
import pyurdme

class NSM2Solver(pyurdme.URDMESolver):
    """ NSM2 solver class. """
    NAME = 'nsm2'

    def __init__(self, *args, **kwargs):
        pyurdme.URDMESolver.__init__(self, *args, **kwargs)
        self.is_compiled = True
        self.propfilename = "solver"
        self.solver_dir = self.URDME_ROOT + '/bin/'
