all: run_tests

TEST_FILES= \
	BCTest.o \
	BoundaryVectorTest.o \
	EllipticSolver2dTest.o \
	EllipticSolverTest.o \
	FluxTest.o \
	GeometryTest.o \
	GridTest.o \
	MotionTest.o \
	NavierStokesModelTest.o \
	OutputProbesTest.o\
	ParmParserTest.o \
	ProjectionSolverTest.o \
	RegularizerTest.o \
	RigidBodyTest.o \
	ScalarTest.o \
	StateTest.o \
	TangentSE2Test.o \
	VectorOperationsTest.o \

include ../config/make.inc

GTEST_DIR = gtest-1.6.0
CPPFLAGS += -I$(GTEST_DIR)/include
CXXFLAGS += $(include_dirs) -I../src
LDFLAGS += $(lib_dirs)
BUILDDIR = ../build
IBPMLIB = libibpm.a
LIBS = $(BUILDDIR)/$(IBPMLIB) -lfftw3 -lm

# All Google Test headers.  Usually you shouldn't change this
# definition.
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h

MAKEDEPEND = gcc -MM

.PHONY: lib clean distclean depend

run_tests: runner
	./runner 2> runner.err

runner: lib $(TEST_FILES) gtest_main.a
	$(CXX) $(LDFLAGS) -o $@ $(TEST_FILES) gtest_main.a $(LIBS)

lib:
	cd $(BUILDDIR) && make $(IBPMLIB)

% : %.o lib
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)

%.o : %.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $(INCPATH) $<

# Usually you shouldn't tweak such internal variables, indicated by a
# trailing _.
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

# For simplicity and to avoid depending on Google Test's
# implementation details, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Test
# compiles fast and for ordinary users its source rarely changes.
gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest_main.cc

gtest.a : gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

gtest_main.a : gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

clean:
	-$(RM) -r *.o

distclean: clean
	-$(RM) runner

depend:
	$(MAKEDEPEND) $(CXXFLAGS) $(INCPATH) *.cc > .depend

-include .depend
