#
# Makefile for Fd 
#       3d acoustic  calcuration
#

OBJ  =  triffy.o

FFLAGS = -O8

all: triffy.dec
	make clean_triffy; make triffy; make clean_triffy

triffy:$(OBJ)
	gfortran $(OBJ) $(FFLAGS) -o triffy.x

clean_triffy:
	rm -f $(OBJ) 

.SUFFIXES: .f .o

.f.o:
	gfortran -c $(FFLAGS) $<


