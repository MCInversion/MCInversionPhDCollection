
// Copyright 2017-2021 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#include "gtest/gtest.h"

#include <pmp/Types.h>

using namespace pmp;

TEST(MatVecTest, allFinite)
{
    vec3 v(0, 0, 0);
    EXPECT_TRUE(v.allFinite());
}

TEST(MatVecTest, nonFiniteInf)
{
    vec3 v(std::numeric_limits<Scalar>::infinity());
    EXPECT_FALSE(v.allFinite());
}

TEST(MatVecTest, nonFiniteNaN)
{
    vec3 v(std::numeric_limits<Scalar>::quiet_NaN());
    EXPECT_FALSE(v.allFinite());
}
