"""This demo program solves the diffusion equation

     div grad u1 = du1/dt
     div grad u2 = du2/dt
     div grad u3 = du3/dt
     div grad u4 = du4/dt
     div grad u5 = du5/dt

on the ecoli mesh with homogeneous Neumann boundary conditions and no forcing term.
"""

# Copyright (C) 2009-2011 Anders Logg
#
# This file is part of DOLFIN.
#
# DOLFIN is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DOLFIN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
#
# First added:  2009-06-15
# Last changed: 2011-06-28

from dolfin import *
import scipy.io
parameters["linear_algebra_backend"] = "uBLAS"

#constants

Dmem = 1E-14
Dcyt = 2.5E-12

#This is simply a test problem I adapted from the Fenics reference book, it can 
#be changed to check matrices but it must have the same number of components as the
#number of species diffusing (i.e. 5 in this case)
alpha = 3; beta = 1.2
u0 = Expression(("1 + x[0]*x[0] + alpha*x[1]*x[1] + + x[2]*x[2] + beta*t", "1 + x[0]*x[0] + alpha*x[1]*x[1] + + x[2]*x[2] + beta*t", "1 + x[0]*x[0] + alpha*x[1]*x[1] + + x[2]*x[2] + beta*t", "1 + x[0]*x[0] + alpha*x[1]*x[1] + + x[2]*x[2] + beta*t", "1 + x[0]*x[0] + alpha*x[1]*x[1] + + x[2]*x[2] + beta*t"), alpha=alpha, beta=beta, t=0)

# Import mesh and define variational problem
mesh = Mesh('coli.xml')
T1 = FunctionSpace(mesh, "CG", 1)
T2 = FunctionSpace(mesh, "CG", 1)
T3 = FunctionSpace(mesh, "CG", 1)
T4 = FunctionSpace(mesh, "CG", 1)
T5 = FunctionSpace(mesh, "CG", 1)

V = T1 * T2 * T3 * T4 * T5

u_1=project(u0,V)
u_1a=u_1[0]
u_1b=u_1[1]
u_1c=u_1[2]
u_1d=u_1[3]
u_1e=u_1[4]

u = TrialFunction(V)
u1 = u[0]
u2 = u[1]
u3 = u[2]
u4 = u[3]
u5 = u[4]

v = TestFunction(V)
v1 = v[0]
v2 = v[1]
v3 = v[2]
v4 = v[3]
v5 = v[4]

#Define time step
dt = 0.3


a = u1*v1*dx + u2*v2*dx + u3*v3*ds + u4*v4*ds + u5*v5*ds + dt*(Dcyt*dot(grad(u1), grad(v1))*dx + Dcyt*dot(grad(u2), grad(v2))*dx + Dmem*dot(grad(u3), grad(v3))*ds +Dmem*dot(grad(u4), grad(v4))*ds + Dmem*dot(grad(u5), grad(v5))*ds)

L = u_1a*v1*dx + u_1b*v2*dx + u_1c*v3*ds + u_1d*v4*ds + u_1e*v5*ds


# Assemble matrices (loop will yield error as matrices are singular)
A = assemble(a)

u = Function(V)   
T = 2    
t = dt  
while t <= T:
	b = assemble(L)
	u0.t = t
	solve(A, u.vector(), b)

	t += dt
	u_1.assign(u)


bdense = b.array()
Adense = A.array()

bsparse = b.data()
Asparse = A.data()



# export to MATLAB as both a sparse and a dense matrix
scipy.io.savemat("rddense.mat", {"Adense": Adense, "bdense": bdense})
scipy.io.savemat("rdsparse.mat", {"Asparse": Asparse, "bsparse": bsparse})

# Compute solution (will yield an error as matrices are singular)
w = Function(V)
solve(a == L, w)
u1 = w[0]
u2 = w[1]
u3 = w[2]
u4 = w[3]
u5 = w[4]

