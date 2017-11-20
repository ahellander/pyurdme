#!/usr/bin/env bash
apt-get -y update
apt-get -y install git
apt-get -y install build-essential python-dev libhdf5-serial-dev
apt-get -y install python-setuptools
apt-get -y install python-matplotlib python-numpy python-scipy
apt-get -y install make
apt-get -y install python-software-properties
#add-apt-repository ppa:fenics-packages/fenics
apt-get update
#apt-get -y install fenics
apt-get -y install cython python-h5py

make -f pyurdme/urdme/build/Makefile.nsm2 URDME_ROOT=pyurdme/urdme
mkdir -p pyurdme/urdme/bin
cp pyurdme/urdme/{build/nsm2/solver.nsm2,bin}
rm -r pyurdme/urdme/build/nsm2

python setup.py install
