# Main makefile for IBPM code
#
# Clancy Rowley
# Princeton University
#
# $Date$
# $Revision$
# $Author$
# $HeadURL$

.PHONY: ibpm test doc clean distclean
DIRS = build test doc

ibpm:
	cd build && $(MAKE)

test:
	cd test && $(MAKE)

doc:
	cd doc && $(MAKE)

all: ibpm test doc

clean:
	for dir in $(DIRS); do ( cd $$dir && $(MAKE) clean; ) done

distclean: clean
	for dir in $(DIRS); do \
	  ( cd $$dir && $(MAKE) distclean; )\
	done
