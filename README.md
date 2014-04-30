pyurdme
=======

pyurdme is a modeling and similation toolkit for spatial stochastic simulations. 


Dependencies
=============
To use pyurdme, you need to install the following dependencies

- Fenics (http://fenicsproject.org/download)
    * Debian: apt-get install fenics
    * Fedora: yum install fenics
    * OSX: We recommend using the binary installer provided at http://fenicsproject.org/download
- numpy, scipy, matplotlib
    Debian: apt-get install python-numpy python-scipy python-matplotlib
    Fedora: yum install python-numpy python-scipy python-matplotlib
    OSX: We recommend using the installer provided by http://fonnesbeck.github.io/ScipySuperpack/


In addition, if you do not use the provided setuptools script, you need to install

- h5py   
- Jinja2 

Installation
=============

pip install https://github.com/ahellander


Quick start
==============

- Open a terminal (On OSX, by opening FeniCS)
- CD into the examples/simple_diffusion folder and run the file "simple_diffusion.py"
