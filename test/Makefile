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

CXXFLAGS += $(include_dirs) -I../src
LDFLAGS += $(lib_dirs)
BUILDDIR = ../build
IBPMLIB = libibpm.a
LIBS = $(BUILDDIR)/$(IBPMLIB) -lgtest -lgtest_main -lfftw3 -lm

MAKEDEPEND = gcc -MM

.PHONY: lib clean distclean depend

run_tests: runner
	./runner 2> runner.err

runner: lib $(TEST_FILES)
	$(CXX) $(LDFLAGS) -o $@ $(TEST_FILES) $(LIBS)

lib:
	cd $(BUILDDIR) && make $(IBPMLIB)

% : %.o lib
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)

%.o : %.cc
	$(CXX) -c $(CXXFLAGS) $(INCPATH) $<

clean:
	-$(RM) -r *.o

distclean: clean
	-$(RM) runner

depend:
	$(MAKEDEPEND) $(CXXFLAGS) $(INCPATH) *.cc > .depend

-include .depend
