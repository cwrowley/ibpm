#include "Geometry.h"
#include <gtest/gtest.h>

namespace {

TEST(GeometryTest, PlusOne) {
	EXPECT_EQ(plus_one(5),6);
}

TEST(GeometryTest, PlusOneAgain) {
	EXPECT_EQ(plus_one(-1),0);
}

TEST(GeometryTest, PlusOneYetAgain) {
    EXPECT_EQ(plus_one(0),1);
}

}  // namespace

