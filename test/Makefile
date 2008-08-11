all: run_tests

TEST_FILES=GridTest.cc GeometryTest.cc BoundaryVectorTest.cc FluxTest.cc \
	ScalarTest.cc RegularizerTest.cc VectorOperationsTest.cc \
	ProjectionSolverTest.cc

CXX = g++
CXXFLAGS = -g
INCPATH = -I ../src
IBPMLIB = ../lib/libibpm.a

.PHONY: $(IBPMLIB) clean distclean

run_tests: runner
	./runner

$(IBPMLIB):
	(cd ../src; make lib)

runner: $(TEST_FILES) $(IBPMLIB)
	$(CXX) $(CXXFLAGS) -o $@ $(INCPATH) $^ -lgtest -lgtest_main -lfftw3 -lm

clean:
	/bin/rm -rf *.o

distclean: clean
	/bin/rm -f runner
