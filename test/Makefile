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
IBPMLIB = ../build/libibpm.a
MAKEDEPEND = gcc -MM

.PHONY: $(IBPMLIB) clean distclean depend

run_tests: runner
	./runner 2> runner.err

% : %.o $(IBPMLIB)
	$(CXX) $(LDFLAGS) -o $@ $^ -lgtest -lgtest_main -lfftw3 -lm

runner: $(TEST_FILES) $(IBPMLIB)
	$(CXX) $(LDFLAGS) -o $@ $^ -lgtest -lgtest_main -lfftw3 -lm

$(IBPMLIB):
	cd ../build && $(MAKE)

%.o : %.cc
	$(CXX) -c $(CXXFLAGS) $(INCPATH) $<

clean:
	-$(RM) -r *.o

distclean: clean
	-$(RM) runner

depend:
	$(MAKEDEPEND) $(CXXFLAGS) $(INCPATH) *.cc > .depend

-include .depend
