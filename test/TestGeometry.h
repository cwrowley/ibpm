#include <cxxtest/TestSuite.h>
#include "Geometry.h"

class TestGeometry : public CxxTest::TestSuite
{
public:
	void testPlusOne(void)
	{
		TS_ASSERT_EQUALS(plus_one(5),6);
	}

	void testPlusOneAgain(void) {
		TS_ASSERT_EQUALS(plus_one(-1),0);
	}
};