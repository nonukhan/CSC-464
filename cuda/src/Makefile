# Makefile
#
# Generic Linux/Unix Makefile for cuda PHYLIP 3.6.


EXEDIR  = ../exe
TIMEDIR = ../../simple-timing/


ifeq ($(type),cuda)
	CC = nvcc
	CFLAGS = -arch sm_13
else # default to c
	CC = gcc
	CFLAGS = -O3 -fomit-frame-pointer
endif

CC = gcc
#
# An optimized one for gcc
CFLAGS = -O3 -fomit-frame-pointer # -Wall

# math library;
LIBS = -lm


ifeq ($(type),cuda)
all: cuda
else
all: dnaml
endif

clean:
	rm -f *.o dnaml dnaml-cu

phylip.o: phylip.h phylip.c
	$(CC) $(CFLAGS) -c phylip.c

seq.o: seq.h seq.c
	$(CC) $(CFLAGS) -c seq.c

timekeeper.o: $(TIMEDIR)timekeeper.h $(TIMEDIR)timekeeper.c
	$(CC) $(CFLAGS) -c $(TIMEDIR)timekeeper.c

dnaml.o: dnaml.c stdfns.c
	$(CC) $(CFLAGS) -c dnaml.c

dnaml: dnaml.o seq.o phylip.o timekeeper.o
	$(CC) $(CFLAGS) dnaml.o seq.o phylip.o timekeeper.o $(LIBS) -o dnaml

cuda: dnaml.c phylip.c seq.c cudafns.cu $(TIMEDIR)timekeeper.c
	./rename-ext.sh cu
	$(CC) $(CFLAGS) dnaml.cu phylip.cu seq.cu cudafns.cu $(TIMEDIR)timekeeper.cu -o dnaml-cu
	./rename-ext.sh c
