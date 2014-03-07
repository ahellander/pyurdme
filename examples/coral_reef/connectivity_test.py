from dolfin import *
import dolfin
import matplotlib.pyplot as plt



# Sub domain for Periodic boundary condition
class PeriodicBoundary2D(SubDomain):

    # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        # return True if on left or bottom boundary AND NOT on one of the two corners (0, 1) and (1, 0)
        return bool((near(x[0], 0) or near(x[1], 0)) and 
                (not ((near(x[0], 0) and near(x[1], 1)) or 
                        (near(x[0], 1) and near(x[1], 0)))) and on_boundary)

    def map(self, x, y):
        if near(x[0], 1) and near(x[1], 1):
            y[0] = x[0] - 1.
            y[1] = x[1] - 1.
        elif near(x[0], 1):
            y[0] = x[0] - 1.
            y[1] = x[1]
        else:   # near(x[1], 1)
            y[0] = x[0]
            y[1] = x[1] - 1.




if __name__ == "__main__":
    mesh = UnitSquareMesh(32, 32)

    fs = dolfin.FunctionSpace(mesh, "Lagrange", 1, constrained_domain=PeriodicBoundary2D())
    trial_function = dolfin.TrialFunction(fs)
    test_function = dolfin.TestFunction(fs)
    a_K = -1*dolfin.inner(dolfin.nabla_grad(trial_function), dolfin.nabla_grad(test_function)) * dolfin.dx
    C = dolfin.assemble(a_K)
    rows, cols, vals = C.data()
    C = scipy.sparse.csr_matrix((vals, cols, rows))
    C = C.tocsc()

    print C
    plt.spy(C)

