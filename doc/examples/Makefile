# Makefile for compiling and running examples
# CWR, 23 Sep 2009

.PHONY: all clean lib

IBPM_DIR = ../..
SRC_DIR = $(IBPM_DIR)/src
IBPM_LIB = $(IBPM_DIR)/build/libibpm.a

include $(IBPM_DIR)/config/make.inc

output_files := $(patsubst %.cc, %.out, $(wildcard *.cc))

all: $(output_files)

%.out: %.x
	./$< > $@

%.x: %.cc
	make lib
	$(CXX) -o $@ -I$(SRC_DIR) $< $(IBPM_LIB) -lfftw3 -lm

cylinder.png: cylinder.out cyl010.plt cylinder.lay plot_cylinder.mcr
	tec360 -b -p plot_cylinder.mcr cylinder.lay

lib:
	$(MAKE) -C $(IBPM_DIR)/build libibpm.a

clean:
	-$(RM) *.out *.x *.plt *.log *.png

