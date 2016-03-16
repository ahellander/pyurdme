import h5py
import numpy

f = h5py.File("result.h5")
T0 = numpy.array(f.get("Trajectories/0/Type_0/positions"))
T1 = numpy.array(f.get("Trajectories/0/Type_1/positions"))

print T0
print numpy.shape(T0)
print numpy.shape(T1)

