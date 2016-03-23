CC=hdf5-1.8.16/c++/src/h5c++
SRC=src/main.cpp mesh.o
EXECUTABLE=bin/mmms
FLAGS=-Wall -pedantic -I$(PWD)/include -g -c -O3 `gsl-config --cflags` --std=c++11 -I/home/fenics/pyurdme/pyurdme/mmms/hdf5-1.8.16/hdf5/include/
CURR_DIR = $(shell pwd)
LFLAGS=`gsl-config --libs`
HDF5FLAGS=-D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -L/home/fenics/pyurdme/pyurdme/mmms/hdf5-1.8.16/hdf5/lib /home/fenics/pyurdme/pyurdme/mmms/hdf5-1.8.16/hdf5/lib/libhdf5_hl_cpp.a /home/fenics/pyurdme/pyurdme/mmms/hdf5-1.8.16/hdf5/lib/libhdf5_cpp.a /home/fenics/pyurdme/pyurdme/mmms/hdf5-1.8.16/hdf5/lib/libhdf5_hl.a /home/fenics/pyurdme/pyurdme/mmms/hdf5-1.8.16/hdf5/lib/libhdf5.a -lz -ldl -lm -Wl,-rpath -Wl,/home/fenics/pyurdme/pyurdme/mmms/hdf5-1.8.16/hdf5/lib
all: mmms

DEPS=main.o mesh.o utils.o urdmemodel.o read_matfile.o coupling.o

mmms: bin/mmms

bin/mmms: $(DEPS)
	$(CC) $(DEPS) $(LFLAGS) $(HDF5FLAGS) -o $(EXECUTABLE)    

main.o: src/main.cpp
	$(CC) src/main.cpp $(FLAGS)


mesh.o:	src/mesh.c
	$(CC) src/mesh.c $(FLAGS)

utils.o: src/utils.c
	$(CC) src/utils.c $(FLAGS)

coupling.o: src/coupling.c
	$(CC) src/coupling.c $(FLAGS)

read_matfile.o: src/read_matfile.c
	$(CC) src/read_matfile.c $(FLAGS)

urdmemodel.o: src/urdmemodel.c
	$(CC) src/urdmemodel.c $(FLAGS)

clean:
	rm -f $(DEPS)
	rm bin/mmms
	rm -rf bin/mmms.dSYM
