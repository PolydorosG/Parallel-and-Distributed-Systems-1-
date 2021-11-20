CC=gcc
MPICC=mpicc
CILKCC=/usr/local/OpenCilk-9.0.1-Linux/bin/clang 
CFLAGS=-O3 

default: all

C_cilk:
	$(CILKCC) $(CFLAGS) -o C_cilk C_cilk.c -fcilkplus
	
C_omp:
	$(CC) $(CFLAGS) -o C_omp C_omp.c -fopenmp

C_threads:
	$(CC) $(CFLAGS) -o C_threads C_threads.c -lpthread

.PHONY: clean

all: C_cilk C_omp C_threads

clean:
	rm -f C_cilk C_omp C_threads