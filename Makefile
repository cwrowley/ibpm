# Main makefile for IBPM code
#
# Clancy Rowley
# Princeton University
#
# $Date$
# $Revision$
# $Author$
# $HeadURL$

all: lib test ibpm

.PHONY: lib test ibpm doc clean
DIRS = src doc test

lib:
	(cd src; make lib)

test:
	(cd test; make run_tests)

ibpm:
	(cd src; make ibpm)

doc:
	(cd doc; make doc)

clean:
	for dir in $(DIRS); do (cd $$dir; make clean); done

distclean: clean
	for dir in $(DIRS); do (cd $$dir; make distclean); done
	/bin/rm -f bin/ibpm
	/bin/rm -f lib/libibpm.a
