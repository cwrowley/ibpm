all: run_tests

#TEST_FILES=SimpleTest.h ComplicatedTest.h
TEST_FILES=TestGrid.h TestGeometry.h
CXX = g++
INCPATH = -I ../src
IBPMLIB = ../lib/libibpm.a

.PHONY: $(IBPMLIB) clean distclean

run_tests: runner
	./runner

$(IBPMLIB):
	(cd ../src; make lib)

runner: runner.cc $(IBPMLIB)
	$(CXX) -o $@ $(INCPATH) $^

runner.cc: $(TEST_FILES) $(IBPMLIB)
	cxxtestgen.pl -o $@ --error-printer $^

clean:
	/bin/rm -rf *.o runner.cc

distclean: clean
	/bin/rm -f runner