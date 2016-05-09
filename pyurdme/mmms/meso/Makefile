CC=g++
SRC=src/main.cpp
EXECUTABLE=bin/meso
FLAGS=-Wall -Wextra -pedantic-errors -pedantic `gsl-config --cflags --libs` -O3 -fopenmp 

CURR_DIR = $(shell pwd)

main:
	$(CC) $(SRC) -o $(EXECUTABLE) $(FLAGS)

clean:
	rm -rf build/meso2
