# Makefile for documentation for IBPM code
#
# Clancy Rowley
# Princeton University
#
# $Date$
# $Revision$
# $Author$
# $HeadURL$

all: doc ibpm_design.pdf

.PHONY: doc clean distclean
.SUFFIXES: .tex .pdf

doc:
	doxygen

.tex.pdf:
	(pdflatex $<; pdflatex $<)

clean:
	/bin/rm -f *.aux *.fdb_latexmk *.log *.out *.pdfsync

distclean: clean
	/bin/rm -rf html latex ibpm_design.pdf
