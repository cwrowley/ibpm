all: run_tests

TEST_FILES=GridTest.o GeometryTest.o BoundaryVectorTest.o FluxTest.o \
	ScalarTest.o RegularizerTest.o VectorOperationsTest.o \
	ProjectionSolverTest.o RigidBodyTest.o TangentSE2Test.o \
	MotionTest.o NavierStokesModelTest.o StateTest.o

include ../src/make.inc

CXXFLAGS += $(include_dirs) -I../src
LDFLAGS += $(lib_dirs)
IBPMLIB = ../lib/libibpm.a
MAKEDEPEND = gcc -MM

.PHONY: $(IBPMLIB) clean distclean depend

run_tests: runner
	./runner

$(IBPMLIB):
	cd ../src && $(MAKE) lib

% : %.o $(IBPMLIB)
	$(CXX) $(LDFLAGS) -o $@ $^ -lgtest -lgtest_main -lfftw3 -lm

runner: $(TEST_FILES) $(IBPMLIB)
	$(CXX) $(LDFLAGS) -o $@ $^ -lgtest -lgtest_main -lfftw3 -lm

%.o : %.cc
	$(CXX) -c $(CXXFLAGS) $(INCPATH) $<

clean:
	-$(RM) -r *.o

distclean: clean
	-$(RM) runner

depend:
	$(MAKEDEPEND) $(CXXFLAGS) $(INCPATH) *.cc > .depend

-include .depend
