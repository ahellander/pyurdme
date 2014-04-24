#!/usr/bin/env python
import pyurdme
from examples.cylinder_demo.cylinder_demo3D import cylinderDemo3D

model = cylinderDemo3D()
result1 = model.run(seed=1)
result2 = model.run(seed=1)
if result1 != result2:
    raise Exception("Running the same model with the same seed does not produce the same result.")
