pyurdme
=======

Python API for URDME solvers.

pyurdme is currently under initial and active development, so many things are not fully functional.

Dependencies
=============
You will need to install (at least) all the following external packages

- Fenics (fenicsproject.org/download)
- NumPy (For OSX, http://fonnesbeck.github.io/ScipySuperpack/, provides installers for SciPy,NumPy, Matplotlib)
- Scipy 
- Gmsh (http://geuz.org/gmsh/)
- h5py (pip install h5py)

Quick start
==============

- Install URDME by cloning http://www.github.com/ahellander/urdme, check out the develop branch and follow the install instructions in the README.
- In addition to the install instructions in the README, execute the script urdme/urdme/build.sh, to build the bundeled hdf5 library. 
- Install pyurdme by cloning repository and checking out the dev branch

- Open a FEniCS terminal
- In that terminal, cd into the main pyurdme directory and source the file "pyurdme_init" (simply sets the PYTHONPATH)
- In the same terminal, cd into the examples/simple_diffusion folder and run the file "simple_diffusion.py"
