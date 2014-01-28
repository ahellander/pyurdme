from pyurdme import pyurdme
from cylinder_demo3D import cylinderDemo3D
import pickle

model = cylinderDemo3D()
model_str = pickle.dumps(model)

model_unpickled = pickle.loads(model_str)
pyurdme.urdme(model_unpickled)
