# Makefile for examples using IBPM library
#
# Clancy Rowley
# Princeton University
#
# $Date$
# $Revision$
# $Author$
# $HeadURL$

include ../src/make.inc

LDFLAGS = -lfftw3 -lm
MAKEDEPEND = gcc -MM

LDFLAGS += $(lib_dirs)
CXXFLAGS += -I../src $(include_dirs)

IBPMLIB = ../lib/libibpm.a 

.PHONY:  clean depend

plunging: plunging.o $(IBPMLIB)
	$(CXX) -o $@ $^ $(LDFLAGS)

pitching: pitching.o $(IBPMLIB)
	$(CXX) -o $@ $^ $(LDFLAGS)

clean:
	/bin/rm -rf *.o

depend:
	$(MAKEDEPEND) $(CXXFLAGS) *.cc > .depend

-include .depend
