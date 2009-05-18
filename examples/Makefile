# Makefile for examples using IBPM library
#
# Clancy Rowley
# Princeton University
#
# $Date$
# $Revision$
# $Author$
# $HeadURL$

include ../config/make.inc

TARGETS = pitching plunging Oseen RigidBodyLoad bin2plt bininfo

LDLIBS = -lfftw3 -lm
MAKEDEPEND = gcc -MM

LDFLAGS += $(lib_dirs)
CXXFLAGS += -I../src $(include_dirs)

IBPMLIB = libibpm.a 

.PHONY: lib clean depend

all: $(TARGETS)

% : %.o lib
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $< ../build/$(IBPMLIB) $(LDLIBS)

lib:
	cd ../build && $(MAKE)

%.o : %.cc
	$(CXX) $(CXXFLAGS) -c $^

% : %.cc

clean:
	-$(RM) -r *.o $(TARGETS)
	-$(RM) tecplot/*.plt

depend:
	$(MAKEDEPEND) $(CXXFLAGS) *.cc > .depend

-include .depend
