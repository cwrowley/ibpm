#include "ProjectionSolver.h"
#include <gtest/gtest.h>

namespace {

class ProjectionSolverTest : public testing::Test {
protected:
    ProjectionSolverTest() {}

    // data
    int _nx;
    int _ny;
};

TEST_F(ProjectionSolverTest, True) {
    EXPECT_DOUBLE_EQ(1., 1.);
}

} // namespace
