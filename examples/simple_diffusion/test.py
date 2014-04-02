import dolfin
dolfin.parameters["linear_algebra_backend"] = "uBLAS"

import numpy

c = dolfin.Circle(0,0,1)
mesh = dolfin.Mesh(c,20)
function_space = dolfin.FunctionSpace(mesh,"Lagrange",1)
trial_function = dolfin.TrialFunction(function_space)
test_function = dolfin.TestFunction(function_space)

form = dolfin.inner(dolfin.nabla_grad(trial_function), dolfin.nabla_grad(test_function))*dolfin.dx
K = dolfin.assemble(form)
rows, cols, vals = K.data()
print numpy.max(vals)