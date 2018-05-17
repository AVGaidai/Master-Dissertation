P=main.out
OBJECTS=main.c
CC=mpicc
CFLAGS=-Wall -g -pedantic 
LIBS=-lm

all:
	$(CC) $(CFLAGS) -o $(P) $(OBJECTS) $(LIBS)

clean:
	rm -rf *.out *.txt
