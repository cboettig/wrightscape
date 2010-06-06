#-----------------------------------#
#	Makefile for GNU C++ Compiler	#
#-----------------------------------#
CC = gcc
CFLAGS = -Wall -lgsl -lgslcblas -lm -O4 -fpic -funroll-loops --fast-math
EXE = jointprob 
FILES = jointprob.c tree.c multimin.c regimes.c siman.c matrix_method.c linearsoln.c 

build: $(FILES)
	$(CC) $(CFLAGS) -o $(EXE) $(FILES)

profile: $(FILES)
	$(CC) $(CFLAGS) -g -pg -o $(EXE) $(FILES)

debug:  $(FILES)
	$(CC) $(CFLAGS) -g -o $(EXE) $(FILES)
	

