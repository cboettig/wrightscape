#-----------------------------------#
#	Makefile for GNU C++ Compiler	#
#-----------------------------------#
CC = gcc
CFLAGS = -Wall -lgsl -lgslcblas -lm -O4 -fpic -funroll-loops 
EXE = piecewise_regimes 
FILES = piecewise_regimes.c p_multimin.c 

build: $(FILES)
	$(CC) $(CFLAGS) -o $(EXE) $(FILES)

profile: $(FILES)
	$(CC) $(CFLAGS) -g -pg -o $(EXE) $(FILES)

debug:  $(FILES)
	$(CC) $(CFLAGS) -g -o $(EXE) $(FILES)
	

