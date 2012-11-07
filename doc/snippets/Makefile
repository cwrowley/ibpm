# Makefile for compiling and running snippets
# CWR, 23 Sep 2009

.PHONY: all clean lib

output_files := $(patsubst %.cc, %.out, $(wildcard *.cc))
IBPM_DIR = ../..
SRC_DIR = $(IBPM_DIR)/src
IBPM_LIB = $(IBPM_DIR)/build/libibpm.a

include $(IBPM_DIR)/config/make.inc

all: $(output_files)

%.out: %.x
	./$< > $@

%.o: %.cc
	$(CXX) -c -o $@ -I$(SRC_DIR) -x c++ -DSNIPPET_FNAME="\"$<\"" template.cc.tpl

%.x: %.o
	make lib
	$(CXX) -o $@ $< $(IBPM_LIB)

lib:
	$(MAKE) -C $(IBPM_DIR)/build libibpm.a

clean:
	-$(RM) *.out *.x
