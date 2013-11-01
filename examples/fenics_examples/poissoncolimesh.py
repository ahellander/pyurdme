from dolfin import *

import scipy.io
parameters["linear_algebra_backend"] = "uBLAS"
# Create mesh and define function space
#mesh = UnitSquare(6, 4)
mesh = Mesh('coli.xml')
V = FunctionSpace(mesh, "Lagrange", 1)

# Define boundary conditions
u0 = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]")
def u0_boundary(x, on_boundary):
	return on_boundary
bc = DirichletBC(V, u0, u0_boundary)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = inner(nabla_grad(u), nabla_grad(v))*dx
L = f*v*dx

b = assemble(L)
A = assemble(a)
AA = A.data()

scipy.io.savemat("Ab.mat", {"A1": AA})

