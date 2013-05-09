"""This demo program is the first attempt at solving the diffusion part of the MinCDE problem (with a mind to solve the reactions through operator splitting later)

     div grad u1 + f = du1/dt     in Bulk
     div grad u1 + f = du1/dt     in Bulk
     div grad u1 + f = du1/dt     in Bulk
     div grad u1 + f = du1/dt     on Surface
     div grad u1 + f = du1/dt     on Surface

with only BULK diffusion of one species on the ecoli mesh with homogeneous Neumann boundary conditions and an added forcing term (optional). Currently implemented with nonhomogeneous initial conditions so diffusion will indeed occur but can simply be changed to constant when a reaction solver is added.
"""

from dolfin import *
import scipy.io
import numpy
parameters["linear_algebra_backend"] = "uBLAS"

#Create mesh from Coli problem
mesh = Mesh('coli.xml')
v = FunctionSpace(mesh, "Lagrange", 1)
V = MixedFunctionSpace([v,v,v,v,v])
print mesh

#Diffusion coefficients
Dmem = 1
Dcyt = 2.5

# Define initial conditions
alpha = 2; beta = 1; gamma=8;
u01 = Expression('1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*x[2]*x[2] + gamma*t',
                alpha=alpha, beta=beta, gamma=gamma, t=0)
u02 = Expression('1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*x[2]*x[2] + gamma*t',
                alpha=alpha, beta=beta, gamma=gamma, t=0)
u03 = Expression('1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*x[2]*x[2] + gamma*t',
                alpha=alpha, beta=beta, gamma=gamma, t=0)
u04 = Expression('1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*x[2]*x[2] + gamma*t',
                alpha=alpha, beta=beta, gamma=gamma, t=0)
u05 = Expression('1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*x[2]*x[2] + gamma*t',
                alpha=alpha, beta=beta, gamma=gamma, t=0)

u_1a=interpolate(u01,v)
u_1b=interpolate(u02,v)
u_1c=interpolate(u03,v)
u_1d=interpolate(u04,v)
u_1e=interpolate(u05,v)

plot(u_1a)
#plot(u_1b)
#plot(u_1c)
#plot(u_1d)
#plot(u_1e)

#Define Trial and Test functions
u1=TrialFunction(v)
u2=TrialFunction(v)
u3=TrialFunction(v)
u4=TrialFunction(v)
u5=TrialFunction(v)

v1 = TestFunction(v)
v2 = TestFunction(v)
v3 = TestFunction(v)
v4 = TestFunction(v)
v5 = TestFunction(v)

T = 4  # total simulation time
dt = 0.3      # time step

a_K1 = Dcyt*inner(nabla_grad(u1), nabla_grad(v1))*dx
a_K2 = Dcyt*inner(nabla_grad(u2), nabla_grad(v2))*dx
a_K3 = Dcyt*inner(nabla_grad(u3), nabla_grad(v3))*dx
a_K4 = Dcyt*inner(nabla_grad(u4), nabla_grad(v4))*ds
a_K5 = Dcyt*inner(nabla_grad(u5), nabla_grad(v5))*ds

a_M1 = u1*v1*dx
a_M2 = u2*v2*dx
a_M3 = u3*v3*dx
a_M4 = u4*v4*ds
a_M5 = u5*v5*ds

K1=assemble(a_K1)
K2=assemble(a_K2)
K3=assemble(a_K3)
K4=assemble(a_K4)
K5=assemble(a_K5)

M1=assemble(a_M1)
M2=assemble(a_M2)
M3=assemble(a_M3)
M4=assemble(a_M4)
M5=assemble(a_M5)

A1 = M1 + dt*K1
A2 = M2 + dt*K2
A3 = M3 + dt*K3
A4 = M4 + dt*K4
A5 = M5 + dt*K5


#Write a loop to create one large stifness matrix from K1,...,K5 DENSE
K1bulk = K1.array()
K2bulk = K2.array()
K3bulk = K3.array()
K4bulk = K4.array()
K5bulk = K5.array()

size = K1bulk.shape
size = size[0]

K = numpy.zeros(shape=(5*size,size))

j=0;
for i in range(0,size):
	K[j] = K1bulk[i]
	K[j+1] = K2bulk[i]
	K[j+2] = K3bulk[i]
	K[j+3] = K4bulk[i]
	K[j+4] = K5bulk[i]
	j=j+5
	

#Now use row sums to compute the lumped mass matrix and store diagonal elements in a volume vector DENSE
M1bulk = M1.array()
M2bulk = M2.array()
M3bulk = M3.array()
M4bulk = M4.array()
M5bulk = M5.array()


size = M1bulk.shape
size = size[0]

rowsumM1 = [0]*size
rowsumM2 = [0]*size
rowsumM3 = [0]*size
rowsumM4 = [0]*size
rowsumM5 = [0]*size
volume = [0]*size*5

j=0;
for i in range(0,size):
	rowsumM1[i] = M1bulk[i].sum()
	rowsumM2[i] = M2bulk[i].sum()
	rowsumM3[i] = M3bulk[i].sum()
	rowsumM4[i] = M4bulk[i].sum()
	rowsumM5[i] = M5bulk[i].sum()

	volume[j] = rowsumM1[i]
	volume[j+1] = rowsumM2[i]
	volume[j+2] = rowsumM3[i]
	volume[j+3] = rowsumM4[i]
	volume[j+4] = rowsumM5[i]
	j = j+5


#Optional forcing term
#f = Expression(("0","0","0","0","0"))

# Compute solution
u1 = Function(v)
u2 = Function(v)
u3 = Function(v)
u4 = Function(v)
u5 = Function(v)

t = dt
while t <= T:
    #f_k = interpolate(f, V)
    #f_k1 = f_k[0]
    #f_k2 = f_k[1]
    #f_k3 = f_k[2]
    #f_k4 = f_k[3]
    #f_k5 = f_k[4]
	
    #F_k1 = f_k1.vector()
    #F_k2 = f_k2.vector()
    #F_k3 = f_k3.vector()
    #F_k4 = f_k4.vector()
    #F_k5 = f_k5.vector()


    b1 = M1*u_1a.vector() #+ dt*M1*F_k1
    b2 = M2*u_1b.vector()
    b3 = M3*u_1c.vector()
    #b4 = M4*u_1d.vector()
    #b5 = M5*u_1e.vector()

    u01.t = t
    u02.t = t
    u03.t = t
    #u04.t = t
    #u05.t = t

    solve(A1, u1.vector(), b1)
    solve(A2, u2.vector(), b2)
    solve(A3, u3.vector(), b3)
    #solve(A4, u4.vector(), b4)
    #solve(A5, u5.vector(), b5)

    #########
    #########   Here is where the code to solve the reaction step will be then
    #########   the value of u can be updated accordingly
    #########
    t += dt
    u_1a.assign(u1)
    u_1b.assign(u2)
    u_1c.assign(u3)
    #u_1d.assign(u4)
    #u_1e.assign(u5)


scipy.io.savemat("MinCDE_Matrices", {"M1": M1bulk, "M2": M2bulk, "Volume": volume, "Size": size, "rowsumM1": rowsumM1, "M4": M4bulk, "BigK": K, "K1": K1bulk,"K2": K2bulk,"K3": K3bulk})

plot(u1)
#plot(u2)

interactive()
