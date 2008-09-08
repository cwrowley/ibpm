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

TARGETS = pitching plunging Oseen RigidBodyLoad bin2plt bininfo

LDLIBS = -lfftw3 -lm
MAKEDEPEND = gcc -MM

LDFLAGS += $(lib_dirs)
CXXFLAGS += -I../src $(include_dirs)

IBPMLIB = ../lib/libibpm.a 

.PHONY:  clean depend

all: $(TARGETS)

% : %.o $(IBPMLIB)
	cd ../src && $(MAKE) lib
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

%.o : %.cc
	$(CXX) $(CXXFLAGS) -c $^

% : %.cc

clean:
	-$(RM) -r *.o $(TARGETS)
	-$(RM) tecplot/*.plt

depend:
	$(MAKEDEPEND) $(CXXFLAGS) *.cc > .depend

-include .depend
