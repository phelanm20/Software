#include "software/new_geom/spline.h"

#include <gtest/gtest.h>

TEST(mariaTest, maria_test)
{
    Spline splin = Spline();
    EXPECT_EQ(splin.mariaTest(3), 3);
}

TEST(mariaTest, maria_testTWO)
{
    Spline splin2 = Spline();
    EXPECT_EQ(splin2.mariaTest(9), 9);
}

/*
TEST(mariaTest, maria_testTHREE)
{
    Spline splin2 = Spline();
    EXPECT_EQ(splin2.mariaTest(9), (1.0/3.0));
}
*/
