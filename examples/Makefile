# Makefile for examples using IBPM library
#
# Clancy Rowley
# Princeton University
#
# $Date: 2008-08-29 06:01:14 -0700 (Fri, 29 Aug 2008) $
# $Revision: 103 $
# $Author: sbrunton $
# $HeadURL: svn+ssh://zion.princeton.edu/ibpm/trunk/src/Makefile $

include ../src/make.inc

LDFLAGS = -lfftw3 -lm
MAKEDEPEND = gcc -MM

LDFLAGS += $(lib_dirs)
CXXFLAGS += -I../src $(include_dirs)

IBPMLIB = ../lib/libibpm.a 

.PHONY:  clean depend

plunging: plunging.o $(IBPMLIB)
	$(CXX) -o $@ $^ $(LDFLAGS)

clean:
	/bin/rm -rf *.o

depend:
	$(MAKEDEPEND) $(CXXFLAGS) *.cc > .depend

-include .depend
