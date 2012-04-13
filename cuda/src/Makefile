# Makefile
#
# Generic Linux/Unix Makefile for cuda PHYLIP 3.6.


EXEDIR  = ../exe
TIMEDIR = ../../simple-timing/


CC = gcc
CFLAGS = -O3 -fomit-frame-pointer
# math library;
LIBS = -lm


all: dnaml

clean:
	rm -f *.o dnaml dnaml-cu $(TIMEDIR)*.o

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

cuda: dnaml.c phylip.c seq.c cudafns.c $(TIMEDIR)timekeeper.c
	./rename-ext.sh cu
	nvcc -arch sm_13 dnaml.cu phylip.cu seq.cu $(TIMEDIR)timekeeper.cu -o dnaml-cu
	./rename-ext.sh c
