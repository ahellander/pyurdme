"""This demo program solves the reaction-diffusion equation

    - div grad u1 + u2 = f
    - div grad u2 + u1 = f

on the unit square with f = sin(x)*sin(y) and homogeneous Neumann
boundary conditions.

This demo is also available in more compact form in short.py,
the world's maybe shortest PDE solver.
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

# Define variational problem
mesh = UnitSquare(8, 8)
T = FunctionSpace(mesh, "CG", 1)
Y = FunctionSpace(mesh, "CG", 1)
V = T * Y
(u1, u2) = TrialFunctions(V)
(v1, v2) = TestFunctions(V)
f = Expression("sin(x[0])*sin(x[1])")
a = dot(grad(u1), grad(v1))*dx + u2*v1*dx + dot(grad(u2), grad(v2))*dx + u1*v2*dx
L = f*(v1+v2)*dx

# Assemble matrices

b = assemble(L)
A = assemble(a)

bb = b.array()
AA = A.array()

# export to MATLAB as dense matrix

scipy.io.savemat("rd.mat", {"A1": AA, "b": bb})

# Compute solution
w = Function(V)
solve(a == L, w)
(u1, u2) = w.split()

