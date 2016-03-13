CC=g++
SRC=src/main.cpp
EXECUTABLE=bin/mmms
FLAGS=-Wall -pedantic `gsl-config --cflags --libs` -g -O3 -lhdf5 -lhdf5_hl -lhdf5_cpp

CURR_DIR = $(shell pwd)

all:
	$(CC) $(SRC) -o $(EXECUTABLE) $(FLAGS)

clean:
	rm bin/mmms
	rm -rf bin/mmms.dSYM