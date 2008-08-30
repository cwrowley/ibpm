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

TARGETS = pitching plunging Oseen

LDLIBS = -lfftw3 -lm
MAKEDEPEND = gcc -MM

LDFLAGS += $(lib_dirs)
CXXFLAGS += -I../src $(include_dirs)

IBPMLIB = ../lib/libibpm.a 

.PHONY:  clean depend

all: $(TARGETS)

% : %.cc $(IBPMLIB)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)
	
clean:
	-$(RM) -r *.o $(TARGETS)
	-$(RM) tecplot/*.plt

depend:
	$(MAKEDEPEND) $(CXXFLAGS) *.cc > .depend

-include .depend
