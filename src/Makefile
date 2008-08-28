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
	
include make.inc

OBJS = Geometry.o NavierStokesModel.o Grid.o Scalar.o Flux.o \
	VectorOperations.o Regularizer.o ProjectionSolver.o \
	ConjugateGradientSolver.o TimeStepper.o Euler.o RigidBody.o \
	Logger.o OutputRestart.o OutputTecplot.o

LDFLAGS = -lfftw3 -lm
AR = ar
ARFLAGS = -r
MAKEDEPEND = gcc -MM

LDFLAGS += $(lib_dirs)
CXXFLAGS += $(include_dirs)

IBPMLIB = ../lib/libibpm.a 
IBPM = ../bin/ibpm

.PHONY: ibpm clean depend

lib: $(IBPMLIB) 

ibpm: $(IBPM)

$(IBPMLIB): $(OBJS) 
	$(AR) $(ARFLAGS) $(IBPMLIB) $(OBJS)

$(IBPM): ibpm.o $(IBPMLIB)
	$(CXX) -o $@ $^ $(LDFLAGS)

clean:
	/bin/rm -rf *.o

depend:
	$(MAKEDEPEND) $(CXXFLAGS) *.cc > .depend
	
include .depend
