"""This demo program solves the diffusion equation

     div grad u1 = du1/dt

with only BULK diffusion of one species on the ecoli mesh with homogeneous Neumann boundary conditions and no forcing term (i.e. no reactions).
"""

from dolfin import *
import scipy.io
import numpy
#parameters["linear_algebra_backend"] = "uBLAS"

#Create mesh from Coli problem
mesh = Mesh('coli.xml')
V = FunctionSpace(mesh, "Lagrange", 1)
print mesh

#Diffusion constants for membrane and cytoplasm
Dmem = 1E-14
Dcyt = 2.5E-12

# Define initial conditions
alpha = 2; beta = 1; gamma=8;
u0 = Expression('1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*x[2]*x[2] + gamma*t',
                alpha=alpha, beta=beta, gamma=gamma, t=0)
u0i = interpolate(u0, V)
plot(u0i)

T = 20  # total simulation time
dt = 0.3      # time step

# Define variational problem

# Laplace term
u = TrialFunction(V)
v = TestFunction(V)
a_K = Dcyt*inner(nabla_grad(u), nabla_grad(v))*dx

# "Mass matrix" term
a_M = u*v*dx

M = assemble(a_M)
K = assemble(a_K)


A = M + dt*K


# f term
#f = Expression('0')

# Compute solution
u = Function(V)
t = dt
while t <= T:
    #f_k = interpolate(f, V)
    #F_k = f_k.vector()
    b = M*u0i.vector() #+ dt*M*F_k
    u0.t = t
    solve(A, u.vector(), b)

    t += dt
    u0i.assign(u)

#Try to get solution evaluated at vertex points instead of cell centers
VertexSolution = numpy.zeros(mesh.num_vertices(), dtype='d')
u.compute_vertex_values(VertexSolution, mesh)

Mbulk = M.array()
Kbulk = K.array()

# export to MATLAB as a dense matrix
scipy.io.savemat("rdBULKdense.mat", {"Mbulkdense": Mbulk, "Kbulkdense": Kbulk, "VertexSolution" : VertexSolution})

file = File("BulkDiffusion.pvd")
file << u

plot(u)
interactive()

