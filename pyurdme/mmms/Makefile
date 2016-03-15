CC=g++
SRC=src/main.cpp mesh.o
EXECUTABLE=bin/mmms
FLAGS=-Wall -pedantic -I$(PWD)/include -g -O3 -c
CURR_DIR = $(shell pwd)
LFLAGS=-lhdf5 -lhdf5_hl -lhdf5_cpp `gsl-config --cflags --libs`

all: mmms

DEPS=main.o mesh.o utils.o urdmemodel.o read_matfile.o coupling.o

mmms: bin/mmms

bin/mmms: $(DEPS)
	$(CC) $(DEPS) $(LFLAGS) -o $(EXECUTABLE)    

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