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

OBJS = \
	CholeskySolver.o \
	ConjugateGradientSolver.o \
	Euler.o \
	Flux.o \
	Geometry.o \
	Grid.o \
	Logger.o \
	NavierStokesModel.o \
	OutputRestart.o \
	OutputTecplot.o \
        OutputForce.o \
	ProjectionSolver.o \
	Regularizer.o \
	RigidBody.o \
	RungeKutta2.o \
	Scalar.o \
	State.o \
	TimeStepper.o \
	VectorOperations.o

LDLIBS = -lfftw3 -lm
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
	$(AR) $(ARFLAGS) $(IBPMLIB) $?

$(IBPM): ibpm.o $(IBPMLIB)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

clean:
	/bin/rm -rf *.o

depend:
	$(MAKEDEPEND) $(CXXFLAGS) *.cc > .depend

-include .depend
