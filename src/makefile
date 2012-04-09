#-----------------------------------#
# Makefile for GNU C Compiler	    #
#-----------------------------------#
CC = gcc
CFLAGS = -Wall -O4 -fpic -funroll-loops 
LIB = -lgsl -lgslcblas -lm
EXE = wrightscape.exe 
LDLIB = /usr/common/usg/gsl/1.13/lib
INC = /usr/common/usg/gsl/1.13/include
FILES = piecewise_regimes.c likelihood.c multimin.c siman.c mvn.c 

build: objects shared

objects: $(FILES)
	$(CC) $(CFLAGS) -c piecewise_regimes.c -o piecewise_regimes.o
	$(CC) $(CFLAGS) -c likelihood.c -o likelihood.o
	$(CC) $(CFLAGS) -c siman.c -o siman.o
	$(CC) $(CFLAGS) -c multimin.c -o multimin.o
	$(CC) $(CFLAGS) -c mvn.c -o mvn.o

shared: $(FILES)
	$(CC) $(LIB) -shared -o wrightscape.so *.o  

ibm: $(FILES)
	$(CC) -I$(INC) -L$(LDLIB) $(LIB) $(CFLAGS) -o $(EXE) $(FILES)

profile: $(FILES)
	$(CC)  $(LIB) $(CFLAGS) -g -pg -o $(EXE) $(FILES)

debug:  $(FILES)
	$(CC) $(LIB) $(CFLAGS) -g -o $(EXE) $(FILES)
	
clean:
	rm -f *.o *.exe *.bak callgrind.* gmon.out

