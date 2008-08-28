# Makefile for library and main executable for IBPM code
#
# Clancy Rowley
# Princeton University
#
# $Date$
# $Revision$
# $Author$
# $HeadURL$

all: lib ibpm
	
OBJS = Geometry.o NavierStokesModel.o Grid.o Scalar.o Flux.o \
	VectorOperations.o Regularizer.o ProjectionSolver.o \
	ConjugateGradientSolver.o TimeStepper.o Euler.o RigidBody.o \
	Logger.o OutputRestart.o OutputTecplot.o

CXX = g++
#CXXFLAGS = -Wall -O3 -funroll-loops -fstrict-aliasing -fomit-frame-pointer \
	-ffast-math -DUSE_LIBBLAS -DNDEBUG
# CXXFLAGS = -Wall -g -fast -funroll-loops -DNDEBUG
CXXFLAGS = -Wall -g
LFLAGS = -lfftw3 -lm
AR = ar
ARFLAGS = -r
MAKEDEPEND = gcc -MM

IBPMLIB = ../lib/libibpm.a 
IBPM = ../bin/ibpm

.PHONY: ibpm clean .depend

lib: $(IBPMLIB) 

ibpm: $(IBPM)

$(IBPMLIB): $(OBJS) 
	$(AR) $(ARFLAGS) $(IBPMLIB) $(OBJS)

$(IBPM): ibpm.o $(IBPMLIB)
	$(CXX) -o $@ $^ $(LFLAGS)

.cc.o:
	$(CXX) -c $(CXXFLAGS) $<

clean:
	/bin/rm -rf *.o

.depend:
	$(MAKEDEPEND) $(CXXFLAGS) *.cc > .depend
	
include .depend
