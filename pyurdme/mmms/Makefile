CC=clang++
SRC=src/main.cpp mesh.o
EXECUTABLE=bin/mmms
FLAGS=-Wall -pedantic -I$(PWD)/include -g  -c `gsl-config --cflags` --std=c++11 
LFLAGS=`gsl-config --libs` -lhdf5 -lhdf5_cpp -lhdf5_hl
HDF5FLAGS= 
all: mmms

DEPS=main.o mesh.o utils.o urdmemodel.o read_matfile.o coupling.o

mmms: bin/mmms

bin/mmms: $(DEPS)
	$(CC) $(DEPS) $(LFLAGS) -o $(EXECUTABLE)    

main.o: src/main.cpp
	$(CC) src/main.cpp $(FLAGS) $(HDF5FLAGS)


mesh.o:	src/mesh.c
	$(CC) src/mesh.c $(FLAGS) $(HDF5FLAGS)

utils.o: src/utils.c
	$(CC) src/utils.c $(FLAGS) $(HDF5FLAGS)

coupling.o: src/coupling.c
	$(CC) src/coupling.c $(FLAGS) $(HDF5FLAGS)

read_matfile.o: src/read_matfile.c
	$(CC) src/read_matfile.c $(FLAGS) $(HDF5FLAGS)

urdmemodel.o: src/urdmemodel.c
	$(CC) src/urdmemodel.c $(FLAGS) $(HDF5FLAGS)

clean:
	rm -f $(DEPS)
	rm bin/mmms
	rm -rf bin/mmms.dSYM
