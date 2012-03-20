# Makefile
#
# Generic Linux/Unix Makefile for MPI PHYLIP 3.6.
# Derived from the original Makefile for PHYLIP 3.6.
#

PACKAGE=phylip
VERSION=3.67

# You should not need to change anything, though if you want you could
# change the first (noncomment) statements to some directory locations
# that might be more useful.
# The compressed tar archive phylip.tar.gz when uncompressed and extracted
# puts the source code into a directory ./src, and also makes two other
# directories ./exe and ./doc for the final executables and the documentation
# files.  Only change the EXEDIR settings if you want something different
# from that structure.
# If it causes trouble in compiling, the CFLAGS statement below may also need
# to be changed.
#
#  To use the MPI PHYLIP v3.6 Makefile, type
#      make all          to compile the whole package but not install it
#                          or remove the object files
#      make clean        to remove all object files and executables from the
#                          current directory
#      make MPIproml     to compile and link one program, (in this example,
#                          MPIproml) and leave the executable and object files
#                          in the current directory (where the source code is).
# ----------------------------------------------------------------------------
#  (Starting here is the section where you may want to change things)
# ----------------------------------------------------------------------------
#
# In the following statements (the ones that set CFLAGS, DFLAGS, LIBS
# and DLIBS, CC and DC) you should make sure each is set properly.
# Usually this will simply involve making sure that the proper statement
# has no "#" as its first character and that all other possibilities
# have "#" for their first character.
#
# ----------------------------------------------------------------------------
#
# This is the CFLAGS statement:
#
# if these statements say "-g"  and that causes trouble, remove the "-g"
# if the "-lX11" does not work (X compiling not working), you may want to
# remove that switch or fix it.
#
# Here are some possible CFLAGS statements:
#
#
#A minimal one
#CFLAGS =
#
# A basic one for debugging
#CFLAGS  = -g 
#
# An optimized one for gcc
#CFLAGS  = -O3 -fomit-frame-pointer
#
# For some serious debugging using Gnu gcc
#
#CFLAGS=-g -Wall -Wmain -Wmissing-prototypes -Wreturn-type -Wstrict-prototypes  -Wunused -Werror -Wredundant-decls -Waggregate-return -Wcast-align -Wcomment
#
# For Digital Alpha systems with Compaq Tru64 Unix
# (however, be aware that this may cause floating-point problems in programs
#  like Dnaml owing to not using IEEE floating point standards).
#CFLAGS = -fast
#
# ----------------------------------------------------------------------------
#
# These are the libraries for the CC compiles, respectively
#
LIBS    = -lm -lmpi
#
#
# ----------------------------------------------------------------------------
#
# The next two assignments are the invocations of the compiler for the
# ordinary compiles and the tree-drawing programs, CC and DC
#
# This one specifies the "cc" C compiler
CC        = mpicc
#
# This one specifies the "cc" C compiler
#CC        = cc
#
#  To use GCC instead, if it is not the compiler that "cc" invokes
#CC        = gcc
#
# ----------------------------------------------------------------------------
#  (After this point there should not be any reason to change anything)
# ----------------------------------------------------------------------------
#
#
# the list of programs
#
PROGS   =          MPIdnacomp MPIdnadist \
                   MPIdnaml MPIdnamlk MPIdnapars MPIdnapenny \
                   MPIproml MPIpromlk MPIprotdist \
                   MPIprotpars MPIrestdist MPIrestml

#
# general commands
#

#  The first target it executed if you just type "make".  It tells you how to
#  use the Makefile.
#
help:
	@echo ""
	@echo " To use the MPI PHYLIP v3.6 Makefile, type"
	@echo "     make all           to compile the whole package but not install it"
	@echo "                          or remove the object files"
	@echo "     make MPIproml      to compile and link one program, (in this example,"
	@echo "                          MPIproml) and leave the executable and object files"
	@echo "                          in the current directory (where the source code is)."
	@echo " "

introduce:
	@echo "Building MPIPHYLIP version $(VERSION)"

all:        introduce $(PROGS)
	@echo "Finished compiling."
	@echo ""

clean:
	@echo "Removing object files to save space"
	@rm -f *.o
	@echo "Finished removing object files.  Now will remove"
	@echo "executable files from the current directory, but not from the"
	@echo "executables directory.  (If some are not here, the makefile"
	@echo "will terminate with an error message but this is not a problem)"
	@echo ""
	@echo "Removing executables from this directory"
	@rm -f ../bin/$(PROGS)
	@echo "Finished cleanup."
	@echo ""

#
# compile object files shared between programs
# (make's implicit rule for %.o will take care of these)
#

phylip.o:     phylip.h
seq.o:        phylip.h seq.h
mlclock.o:    phylip.h mlclock.h
printree.o:   phylip.h printree.h

#
# compile the individual programs
#

MPIdnacomp.o:      MPIdnacomp.c MPIcopyfile.c seq.h phylip.h
	$(CC) $(CFLAGS) -c MPIdnacomp.c

MPIdnacomp:      MPIdnacomp.o MPIcopyfile.o seq.o phylip.o 
	$(CC) $(CFLAGS) MPIdnacomp.o MPIcopyfile.o seq.o phylip.o $(LIBS) -o ../bin/MPIdnacomp

MPIdnadist.o:      MPIdnadist.c MPIcopyfile.c seq.h phylip.h
	$(CC) $(CFLAGS) -c MPIdnadist.c

MPIdnadist:      MPIdnadist.o MPIcopyfile.c seq.o phylip.o 
	$(CC) $(CFLAGS) MPIdnadist.o MPIcopyfile.o seq.o phylip.o $(LIBS) -o ../bin/MPIdnadist

MPIdnaml.o:      MPIdnaml.c MPIcopyfile.c seq.h phylip.h
	$(CC) $(CFLAGS) -c MPIdnaml.c

MPIdnaml:      MPIdnaml.o MPIcopyfile.o seq.o phylip.o 
	$(CC) $(CFLAGS) MPIdnaml.o MPIcopyfile.o seq.o phylip.o $(LIBS) -o ../bin/MPIdnaml

MPIdnamlk.o:      MPIdnamlk.c MPIcopyfile.c seq.h phylip.h mlclock.h printree.h
	$(CC) $(CFLAGS) -c MPIdnamlk.c

MPIdnamlk:      MPIdnamlk.o MPIcopyfile.o seq.o phylip.o mlclock.o printree.o
	$(CC) $(CFLAGS) $^ $(LIBS) -o ../bin/MPIdnamlk

MPIdnapenny.o:      MPIdnapenny.c MPIcopyfile.c seq.h phylip.h
	$(CC) $(CFLAGS) -c MPIdnapenny.c

MPIdnapenny:      MPIdnapenny.o MPIcopyfile.o seq.o phylip.o 
	$(CC) $(CFLAGS) MPIdnapenny.o MPIcopyfile.o seq.o phylip.o $(LIBS) -o ../bin/MPIdnapenny

MPIdnapars.o:      MPIdnapars.c MPIcopyfile.c seq.h phylip.h
	$(CC) $(CFLAGS) -c MPIdnapars.c 

MPIdnapars:      MPIdnapars.o seq.o phylip.o 
	$(CC) $(CFLAGS) MPIdnapars.o MPIcopyfile.o seq.o phylip.o $(LIBS) -o ../bin/MPIdnapars

MPIproml.o:      MPIproml.c MPIcopyfile.c seq.h phylip.h
	$(CC) $(CFLAGS) -c MPIproml.c

MPIproml:      MPIproml.o MPIcopyfile.o seq.o phylip.o
	$(CC) $(CFLAGS) MPIproml.o MPIcopyfile.c seq.o phylip.o $(LIBS) -o ../bin/MPIproml


MPIpromlk.o:      MPIpromlk.c MPIcopyfile.c seq.h phylip.h mlclock.h printree.h
	$(CC) $(CFLAGS) -c MPIpromlk.c

MPIpromlk:      MPIpromlk.o MPIcopyfile.o seq.o phylip.o mlclock.o printree.o
	$(CC) $(CFLAGS) $^ $(LIBS) -o ../bin/MPIpromlk


MPIprotdist.o:      MPIprotdist.c MPIcopyfile.c seq.h phylip.h
	$(CC) $(CFLAGS) -c MPIprotdist.c MPIcopyfile.c

MPIprotdist:      MPIprotdist.o MPIcopyfile.o seq.o phylip.o 
	$(CC) $(CFLAGS) MPIprotdist.o MPIcopyfile.o seq.o phylip.o $(LIBS) -o ../bin/MPIprotdist

MPIprotpars.o:      MPIprotpars.c MPIcopyfile.c seq.h phylip.h
	$(CC) $(CFLAGS) -c MPIprotpars.c MPIcopyfile.c

MPIprotpars:      MPIprotpars.o MPIcopyfile.o seq.o phylip.o 
	$(CC) $(CFLAGS) MPIprotpars.o MPIcopyfile.o seq.o phylip.o $(LIBS) -o ../bin/MPIprotpars

MPIrestdist.o: MPIrestdist.c MPIcopyfile.c seq.h phylip.h
	 $(CC) $(CFLAGS) -c MPIrestdist.c

MPIrestdist: MPIrestdist.o MPIcopyfile.o seq.o phylip.o 
	$(CC) $(CFLAGS) MPIrestdist.o MPIcopyfile.o seq.o phylip.o $(LIBS) -o ../bin/MPIrestdist

MPIrestml.o: MPIrestml.c MPIcopyfile.c seq.h phylip.h
	 $(CC) $(CFLAGS) -c MPIrestml.c

MPIrestml: MPIrestml.o MPIcopyfile.o seq.o phylip.o 
	$(CC) $(CFLAGS) MPIrestml.o MPIcopyfile.o seq.o phylip.o $(LIBS) -o ../bin/MPIrestml


