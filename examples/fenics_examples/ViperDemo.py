#This runs an example problem adapted from the Fenics Book
#on the Coli mesh from Gmsh to look at Viper/Paraview for
#viz needs. It is the Poisson equation with Dirichlet BC.

from dolfin import *
import scipy.io
parameters["linear_algebra_backend"] = "uBLAS"

# Create mesh and define function space
mesh=Mesh('coli.xml')
V = FunctionSpace(mesh, "Lagrange", 1)

# Define boundary conditions
u0 = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]")
def u0_boundary(x, on_boundary):
	return on_boundary
bc = DirichletBC(V, u0, u0_boundary)

#Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = inner(nabla_grad(u), nabla_grad(v))*dx
L = f*v*dx

#Assemble the matrices and output dense/sparse versions
#to MATLAB output
b = assemble(L)
A = assemble(a)

bdense = b.array()
Adense = A.array()

bsparse = b.data()
Asparse = A.data()

scipy.io.savemat("rddense.mat", {"Adense": Adense, "bdense": bdense})
scipy.io.savemat("rdsparse.mat", {"Asparse": Asparse, "bsparse": bsparse})

# Compute solution
u = Function(V)
solve(a == L, u, bc)

# Plot solution and mesh
plot(u, title="Example Problem with E. Coli Mesh")
plot(mesh, title="Mesh for Min Ocsillations in E. Coli", axes=True)

# Dump solution to file in VTK format for ParaView
file = File("poisson.pvd")
file << u

# Hold plot
interactive()

