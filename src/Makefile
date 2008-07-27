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
	
OBJS = Geometry.o NavierStokesModel.o Grid.o Scalar.o Flux.o VectorOperations.o

CXX = g++
CXXFLAGS = -g
AR = ar
ARFLAGS = -r
MAKEDEPEND = gcc -MM

IBPMLIB = ../lib/libibpm.a
IBPM = ../bin/ibpm

.PHONY: ibpm clean

lib: $(IBPMLIB)

ibpm: $(IBPM)

$(IBPMLIB): $(OBJS)
	$(AR) $(ARFLAGS) $(IBPMLIB) $(OBJS)

$(IBPM): ibpm.o $(IBPMLIB)
	$(CXX) -o $@ $^

.cc.o:
	$(CXX) -c $(CXXFLAGS) $<

clean:
	/bin/rm -rf *.o

.depend:
	$(MAKEDEPEND) $(CXXFLAGS) *.cc > .depend
	
include .depend