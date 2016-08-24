#!/bin/bash

#Set up PyURDME environment.

apt-get -qq update
apt-get install -y software-properties-common 

add-apt-repository -y ppa:fenics-packages/fenics
apt-get update
apt-get install -y fenics
apt-get dist-upgrade -y

apt-get install -y libgsl0-dev

apt-get install -y python-pip
pip install --upgrade pip

pip install h5py

pip install jupyter
jupyter notebook --generate-config

sed -i "1s/^/c = get_config()\nc.IPKernelApp.pylab = 'inline'\nc.NotebookApp.ip = '*'\nc.NotebookApp.open_browser = False\nc.NotebookApp.port = 8889\n/" /root/.jupyter/jupyter_notebook_config.py

cd /home/pyurdme/pyurdme/mmms/hdf5-1.8.16/
make
make install

cd /home/pyurdme/pyurdme/mmms
make clean
make

source /home/pyurdme/pyurdme_init

cd /home/pyurdme/examples/micro1
jupyter notebook
