all: run_tests

TEST_FILES=GridTest.o GeometryTest.o BoundaryVectorTest.o FluxTest.o \
	ScalarTest.o RegularizerTest.o VectorOperationsTest.o \
	ProjectionSolverTest.o RigidBodyTest.o TangentSE2Test.o \
	MotionTest.o 

CXX = g++
CXXFLAGS = -g
MAKEDEPEND = gcc -MM
INCPATH = -I ../src
IBPMLIB = ../lib/libibpm.a

.PHONY: $(IBPMLIB) clean distclean .depend

run_tests: runner
	./runner

$(IBPMLIB):
	(cd ../src; make lib)

runner: $(TEST_FILES) $(IBPMLIB)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lgtest -lgtest_main -lfftw3 -lm

oseen: Oseen.o $(IBPMLIB)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lfftw3 -lm	

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) $<

clean:
	/bin/rm -rf *.o

distclean: clean
	/bin/rm -f runner

.depend:
	$(MAKEDEPEND) $(CXXFLAGS) $(INCPATH) *.cc > .depend

include .depend
