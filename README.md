pyurdme
=======

pyurdme is a modeling and simulation toolkit for spatial stochastic simulations. It makes use of a modified version of the core solver of URDME (www.urdme.org) for mesocopic simulations via the Reaction-Diffusion Master Equation (RDME), and builds on Dolfin/FeniCS (http://fenicsproject.org) for geometric modeling, meshing and Finite Element Assembly.   

Currently supported (tested) platforms are MacOSX >= 10.8 and Ubunutu >= 12.04.   

# Dependencies

To install and use pyurdme, you need to satisfy the following dependencies. Below we point at suggested ways to satisfy them for OSX and Ubuntu. 

- gcc, gfortran (for scipy, matplotlib)

- Fenics
    * Ubuntu:
        + apt-get install fenics
    * OSX: 
        + We recommend using the binary installer provided here: http://fenicsproject.org/download
- numpy, scipy, matplotlib
   * Ubuntu: 
      + apt-get install python-numpy python-scipy python-matplotlib
   * OSX: 
      + We recommend using the installer provided by http://fonnesbeck.github.io/ScipySuperpack/

- h5py:
   * Ubuntu:
     + apt-get intall python-h5py
   * OSX:
      + brew install h5py

In addition, if you do not use the provided setuptools script, you need to install

- Jinja2 

Note: We strongly recommend using the appropriate package manager or binary installer for your platform to satisfy the above dependencies. However, we understand that some users prefer a non-system-wide installation of python packages, such as if using virtualenv. If one of the above listed dependencies is not satiesfied, setuptools will try to install it. For numpy, scipy, matplotlib, h5py this involves building from source. Due to the many non-python dependencies, you will likely need to install development versions of certain libraries (such as freetype and libpng). An easy way to satisfy the dependencies for Ubuntu is

```bash
apt-get build-dep python-numpy python-scipy python-matplotlib
```

If you do not mind system-wide installations, we provide a script to insall all dependecies for Ubunutu, see detailed instrcuctions below. 

# Installation

## Ubuntu
For Ubuntu, we provide a script that will install pyurdme and all dependecies. The install script must be run as root
```bash
git clone https://github.com/ahellander/pyurdme
cd pyurdme
./install_ubuntu.sh
```

If you want to manage the dependencies yourself, after installing them do:

```bash
git clone https://github.com/ahellander/pyurdme.git
cd pyurdme
python setup.py install 

```
or simply 

```bash
pip install https://github.com/ahellander/pyurdme/tarball/master
```

## OSX

```bash 
git clone https://github.com/ahellander/pyurdme.git
cd pyurdme
python setup.py install 
```

or simply 
```bash
pip install https://github.com/ahellander/pyurdme/tarball/master
```

Quick start
==============

