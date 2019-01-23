# Makefile for library and main executable for IBPM code
#
# Clancy Rowley
# Princeton University
#
# $Date: 2009-05-12 21:37:53 -0600 (Tue, 12 May 2009) $
# $Revision: 201 $
# $Author: zma $
# $HeadURL: svn+ssh://rainier.princeton.edu/ibpm/trunk/src/Makefile $

EXECUTABLES = ibpm checkgeom

all: libibpm.a $(EXECUTABLES)

include ../config/make.inc

OBJS = \
	BaseFlow.o \
	BC.o \
	BoundaryVector.o \
	CholeskySolver.o \
	ConjugateGradientSolver.o \
	EllipticSolver.o \
	EllipticSolver2d.o \
	Field.o \
	Flux.o \
	Geometry.o \
	Grid.o \
	IBSolver.o \
	Logger.o \
	NavierStokesModel.o \
	OutputEnergy.o \
	OutputForce.o \
	OutputRestart.o \
	OutputTecplot.o \
	OutputProbes.o\
	ParmParser.o \
	ProjectionSolver.o \
	Regularizer.o \
	RigidBody.o \
	Scalar.o \
	ScalarToTecplot.o \
	State.o \
    StateVector.o \
	utils.o \
	VectorOperations.o

VPATH = ../src

AR = ar
ARFLAGS = -r
MAKEDEPEND = gcc -MM

LDLIBS = -lfftw3 -lm
LDFLAGS += $(lib_dirs)
CXXFLAGS += $(include_dirs)

.PHONY: clean distclean depend

libibpm.a: $(OBJS)
	$(AR) $(ARFLAGS) $@ $?

$(EXECUTABLES) : % : %.o libibpm.a
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

clean:
	-/bin/rm -rf *.o .depend

distclean: clean
	-/bin/rm -rf libibpm.a $(EXECUTABLES)

depend:
	$(MAKEDEPEND) $(CXXFLAGS) ../src/*.cc > .depend

-include .depend
