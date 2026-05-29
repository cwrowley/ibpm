SHELL = bash

mode = default # Choose "default" or "debug".
compiler = g++ # Choose "g++", "clang++", or "icpc".

extension := .cc
depExtension := .d

sources := $(wildcard *$(extension))
objects := $(sources:$(extension)=.o)
depend := $(sources:$(extension)=$(depExtension))
executables := $(sources:$(extension)=)

# Default options.

CXX := $(compiler) -std=c++11
OPT := -O3

ifeq ($(strip $(compiler)), g++)
WARN := -Wall -Wextra -pedantic -Wno-comment -Wno-unused-local-typedefs \
	-Wno-deprecated-declarations
else ifeq ($(strip $(compiler)), clang++)
CXX += -stdlib=libc++
WARN := -Wall -Wextra -pedantic
else ifeq ($(strip $(compiler)), icpc)
WARN := -Wall -diag-disable 279,2196
endif

CPPFLAGS := -DNDEBUG
CXXFLAGS = $(DEBUG) $(INCLUDE) $(OPT) $(WARN)
DEBUG :=

include include.mk
INCLUDE += -I.. # For Newton_GMRES headers.

# Other options.

ifeq ($(strip $(mode)), debug)
DEBUG := -g
OPT := -Og
endif

.PHONY : all clean

all : $(executables)
# Create a .d dependency Makefile from each source file.
%$(depExtension) : %$(extension)
	@set -e; \
	$(RM) $@; \
	$(CXX) -MM $(CPPFLAGS) $(INCLUDE) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	$(RM) $@.$$$$
% : %.o
	$(CXX) -o $@ $<
clean :
	$(RM) $(executables) $(objects) $(depend)

include $(depend) # Include .d dependency Makefiles.
