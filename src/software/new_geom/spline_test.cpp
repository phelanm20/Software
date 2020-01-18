#include "software/new_geom/spline.h"

#include <gtest/gtest.h>

TEST(TestSpline, test_eigen)
{
    Spline s = Spline();
    EXPECT_EQ(s.mariaTest(6), 3);

}
