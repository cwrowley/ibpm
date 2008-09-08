# Makefile for library and main executable for IBPM code
#
# Clancy Rowley
# Princeton University
#
# $Date$
# $Revision$
# $Author$
# $HeadURL$

IBPM_VERSION = 0.1

all: lib ibpm
	
include make.inc

OBJS = \
	AdamsBashforth.o \
	BoundaryVector.o \
	CholeskySolver.o \
	ConjugateGradientSolver.o \
	Euler.o \
	Field.o \
	Flux.o \
	Geometry.o \
	Grid.o \
	Logger.o \
	NavierStokesModel.o \
	OutputForce.o \
	OutputRestart.o \
	OutputTecplot.o \
	ParmParser.o \
	ProjectionSolver.o \
	Regularizer.o \
	RigidBody.o \
	RungeKutta2.o \
	RungeKutta3.o \
	Scalar.o \
	State.o \
	TimeStepper.o \
	utils.o \
	VectorOperations.o

LDLIBS = -lfftw3 -lm
AR = ar
ARFLAGS = -r
MAKEDEPEND = gcc -MM

VERSIONDEF = -D'IBPM_VERSION=$(IBPM_VERSION)'
LDFLAGS += $(lib_dirs)
CXXFLAGS += $(VERSIONDEF) $(include_dirs)

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
