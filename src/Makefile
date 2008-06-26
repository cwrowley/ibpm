# Makefile for library and main executable for IBPM code
#
# Clancy Rowley
# Princeton University
#
# $Date: $
# $Revision: $
# $Author: $
# $HeadURL: $

all: lib ibpm
	
OBJS = Geometry.o NavierStokes.o
CXX = g++
CXXFLAGS = -g
AR = ar
ARFLAGS = -r

IBPMLIB = ../lib/libibpm.a
IBPM = ../bin/ibpm

OBJS = Geometry.o NavierStokes.o

.PHONY: ibpm clean

lib: $(IBPMLIB)

ibpm: $(IBPM)

$(IBPMLIB): $(OBJS)
	$(AR) $(ARFLAGS) $(IBPMLIB) $(OBJS)

$(IBPM): ibpm.o $(IBPMLIB)
	$(CXX) -o $@ $^

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $<

clean:
	/bin/rm -rf *.o
