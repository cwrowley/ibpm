#include "Grid.h"
#include <gtest/gtest.h>

using namespace ibpm;

namespace {

double tol = 1e-14;
    
class GridTest : public testing::Test {
protected:
    GridTest() : 
      _nx(10),
      _ny(20),
      _length(2),
      _xOffset(-1),
      _yOffset(-3),
      _grid(_nx, _ny, _length, _xOffset, _yOffset) {}

	virtual ~GridTest() {}

    int _nx;
    int _ny;
    double _length;
    double _xOffset;
    double _yOffset;
    Grid _grid;
};

TEST_F(GridTest, TestDx) {
    EXPECT_DOUBLE_EQ(_grid.getDx(), 0.2);
}

TEST_F(GridTest, TestXCenter) {
     EXPECT_DOUBLE_EQ(_grid.getXCenter(0),-0.9);
     EXPECT_DOUBLE_EQ(_grid.getXCenter(9), 0.9); 
     EXPECT_NEAR(_grid.getXCenter(5), 0.1, tol);
}
// 
TEST_F(GridTest, testYCenter) {
     EXPECT_DOUBLE_EQ(_grid.getYCenter(0), -2.9);
     EXPECT_DOUBLE_EQ(_grid.getYCenter(10),-0.9);    
     EXPECT_DOUBLE_EQ(_grid.getYCenter(19), 0.9);
}

TEST_F(GridTest, testXEdge) {
    EXPECT_DOUBLE_EQ(_grid.getXEdge(0),-1);
    EXPECT_DOUBLE_EQ(_grid.getXEdge(5), 0);
    EXPECT_DOUBLE_EQ(_grid.getXEdge(10),1);     
}

TEST_F(GridTest, testYEdge) {
    EXPECT_DOUBLE_EQ(_grid.getYEdge(0), -3);
    EXPECT_DOUBLE_EQ(_grid.getYEdge(10),-1);
    EXPECT_DOUBLE_EQ(_grid.getYEdge(20), 1);        
}

// NOTE: For Google Test, Death tests available only under Linux.
//
// TEST(GridTest, testBoundaries) {
//     EXPECT_DEATH(_grid.getXCenter(-1),"");
// }

//  void testBoundaries() {
//      TS_ASSERT_THROWS_ANYTHING(_grid->getXCenter(-1));
//      TS_ASSERT_THROWS_ANYTHING(_grid->getYCenter(-1));
//      TS_ASSERT_THROWS_ANYTHING(_grid->getXCenter(10));
//      TS_ASSERT_THROWS_ANYTHING(_grid->getYCenter(20));
//      TS_ASSERT_THROWS_ANYTHING(_grid->getXEdge(-1));
//      TS_ASSERT_THROWS_ANYTHING(_grid->getYEdge(-1));
//      TS_ASSERT_THROWS_ANYTHING(_grid->getXEdge(11));
//      TS_ASSERT_THROWS_ANYTHING(_grid->getYEdge(21));
//  }
// 

}  // namespace

