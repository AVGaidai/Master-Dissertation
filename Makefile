P=matrix_generation.out
OBJECTS=matrix_generation.c
CC=gcc
CFLAGS=-Wall -g -pedantic -ansi

all:
	$(CC) $(CFLAGS) -o $(P) $(OBJECTS)

clean:
	rm -rf *.out *.txt
