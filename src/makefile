#-----------------------------------#
#	Makefile for GNU C++ Compiler	#
#-----------------------------------#
CC = gcc
CFLAGS = -Wall -lgsl -lgslcblas -lm -O4 -fpic -funroll-loops 
EXE = wrightscape.exe 
LDLIB = /usr/common/usg/gsl/1.13/lib
INC = /usr/common/usg/gsl/1.13/include
FILES = piecewise_regimes.c likelihood.c multimin.c siman.c mvn.c 

build: $(FILES)
	$(CC) -I$(INC) -L$(LDLIB) $(CFLAGS) -o $(EXE) $(FILES)

profile: $(FILES)
	$(CC) $(CFLAGS) -g -pg -o $(EXE) $(FILES)

debug:  $(FILES)
	$(CC) $(CFLAGS) -g -o $(EXE) $(FILES)
	
clean:
	rm -f *.o *.exe *.bak callgrind.* gmon.out

