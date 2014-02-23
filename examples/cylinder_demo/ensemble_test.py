#!/usr/bin/env python
from cylinder_demo3D import cylinderDemo3D
from pyurdme.nsmsolver import NSMSolver
import matplotlib.pyplot as plt
import numpy

# Create the model and solver
model = cylinderDemo3D()
sol = NSMSolver(model, report_level=2)
print "Beginning simulation"
result = sol.run()
results = [result]
#results = sol.run_ensemble(5)

# Plot of the time-average spatial concentration.
x_vals = model.mesh.coordinates()[:, 0]
l = x_vals.shape[0]
for res in results:
    print "result.filename={0} loaded={1}".format(res.filename, res.data_is_loaded)
    plt.clf()
    A_vals = numpy.sum(res['U'], axis=0)[0:2*l-1:2]
    B_vals = numpy.sum(res['U'], axis=0)[1:2*l:2]
    plt.plot(x_vals,A_vals/model.vol,'.r',x_vals,B_vals/model.vol,'.b')
    plt.legend(['A', 'B'])
    plt.show()

