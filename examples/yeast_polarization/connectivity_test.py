#!/usr/bin/env python
from dolfin import *
import dolfin
dolfin.parameters["linear_algebra_backend"] = "uBLAS"
import matplotlib.pyplot as plt
import scipy.sparse



# Sub domain for Periodic boundary condition
class PeriodicBoundary1D(SubDomain):

    # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        # return True if on left or bottom boundary AND NOT on one of the two corners (0, 1) and (1, 0)
        print "PeriodicBoundary1D.inside(x={0}) = {1}".format(x[0], not bool((near(x[0], 0)) and on_boundary))
        return not bool((near(x[0], 0)) and on_boundary)

    def map(self, x, y):
        print "PeriodicBoundary1D.map(x={0}, y={0}) y={0}".format(x, y)
        if near(x[0], 0):
            #y[0] = x[0] + 1.
            y[0] = x[0]+1
            print "\tnew y = {0}".format(y)
        else:
            print "\tnew y = {0}".format(y)



if __name__ == "__main__":
    mesh = UnitIntervalMesh(6)

    fs = dolfin.FunctionSpace(mesh, "Lagrange", 1, constrained_domain=PeriodicBoundary1D())
    #fs = dolfin.FunctionSpace(mesh, "Lagrange", 1)
    trial_function = dolfin.TrialFunction(fs)
    test_function = dolfin.TestFunction(fs)
    a_K = -1*dolfin.inner(dolfin.nabla_grad(trial_function), dolfin.nabla_grad(test_function)) * dolfin.dx
    C = dolfin.assemble(a_K)
    rows, cols, vals = C.data()
    C = scipy.sparse.csr_matrix((vals, cols, rows))
    C = C.tocsc()

    print "C:\n", C
    p = mesh.coordinates()

    # This gives you a mapping from dofs (rows in matrix) to vertices (points in mesh.coordinates)
    dof2vtx =  fs.dofmap().vertex_to_dof_map(mesh)
    # It seems like the number of dofs now is num_coordinates-1, and I guess this makes sense since x[0] = x[1]
    
    
    print "coordinates\n", p

    print "dof to vertex mapping\n", dof2vtx
    print "mapped coordinates\n", p[dof2vtx]

    plt.spy(C)
    plt.title("1D connectivity, with periodic BCs")
    plt.show()

