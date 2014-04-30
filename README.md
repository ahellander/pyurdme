pyurdme
=======

pyurdme is a modeling and simulation toolkit for spatial stochastic simulations. 


Dependencies
=============
To install and use pyurdme, you need to install the following dependencies

- gcc, gfortran (for scipy, matplotlib)

- Fenics
    * Debian:
        + apt-get install fenics
    * Fedora: 
        + yum install fenics
    * OSX: 
        + We recommend using the binary installer provided here: http://fenicsproject.org/download
- numpy, scipy, matplotlib
   * Debian: 
      + apt-get install python-numpy python-scipy python-matplotlib
   * Fedora: 
      + yum install python-numpy python-scipy python-matplotlib
   * OSX: 
      + We recommend using the installer provided by http://fonnesbeck.github.io/ScipySuperpack/

- h5py:
   * Debian:
   * Fedora:
   * OSX:
      + brew install h5py

In addition, if you do not use the provided setuptools script, you need to install

- Jinja2 




Installation
=============

After satisfying the above dependencies, you can simply  

```bash
pip install https://github.com/ahellander/pyurdme/tarball/master
```

or if you prefer to download the source

```bash
git clone https://github.com/ahellander/pyurdme.git
cd pyurdme
python setup.py install 

```

Note: We strongly recommend using the appropriate package manager or binary installer for your platform to satisfy the above dependencies. However, we understand that some users prefer a non-system-wide installation of python packages, such as if using virtualenv. If one of the above listed dependencies is not satiesfied, setuptools will try to install it. For numpy, scipy, matplotlib, h5py this involves building from source. Due to the many non-python dependencies, you will likely need to install development versions of certain libraries (such as freetype and libpng). An easy way to satisy the dependencies for Ubuntu is

```bash
apt-get build-dep python-numpy python-scipy python-matplotlib
```

Quick start
==============

- Open a terminal (On OSX, by opening FeniCS)
- CD into the examples/simple_diffusion folder and run the file "simple_diffusion.py"
