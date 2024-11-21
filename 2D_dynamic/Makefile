#
# Makefile for Fd 
#       3d acoustic  calcuration
#

OBJ  =	fault2dPML.o PMLpack.o

FFLAGS = -O8

fault2dPML: arrays.f $(OBJ)
	gfortran $(OBJ) $(FFLAGS) -o fault.x

all: arrays.f  PMLpack.f fault2dPML.f
	make clean; make fault2dPML; make clean

clean:
	rm -f $(OBJ) $(OBJ2)

.SUFFIXES: .f .o

.f.o:
	gfortran -c $(FFLAGS) $<
