# Makefile for documentation for IBPM code
#
# Clancy Rowley
# Princeton University
#
# $Date$
# $Revision$
# $Author$
# $HeadURL$

all: doc ibpm_manual.pdf

SUBDIRS = examples snippets

.PHONY: doc clean distclean $(SUBDIRS)
.SUFFIXES: .tex .pdf

include ../config/make.inc

doc: $(SUBDIRS)
	doxygen

$(SUBDIRS):
	$(MAKE) -C $@

.tex.pdf:
	(pdflatex $<; pdflatex $<; pdflatex $<)

clean:
	-$(RM) -f *.aux *.toc *.blg *.fdb_latexmk *.log *.out *.pdfsync
	for dir in $(SUBDIRS); do \
	  $(MAKE) -C $$dir clean; \
	done

distclean: clean
	-$(RM) -rf html latex ibpm_design.pdf
