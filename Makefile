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
	cd src && $(MAKE) lib

test:
	cd test && $(MAKE) run_tests

ibpm:
	cd src && $(MAKE) ibpm

doc:
	cd doc && $(MAKE) doc

clean:
	for dir in $(DIRS); do ( cd $$dir && $(MAKE) clean; ) done

distclean: clean
	for dir in $(DIRS); do \
	  ( cd $$dir && $(MAKE) distclean; )\
	done
	/bin/rm -f bin/ibpm
	/bin/rm -f lib/libibpm.a
