#include "Common.h"

#include <gtest/gtest.h>

std::vector<Point> segmentation(const std::vector<Point>& points);

void check(const std::vector<Point>& input, const std::vector<Point>& answer)
{
  ASSERT_TRUE(segmentation(input) == answer);
}

TEST(Segmentation, RightLeft)
{
  EXPECT_NO_FATAL_FAILURE(check(
  { { 1, 1 }, { 1, 2 }, { -1, 1 } },
  { { 1, 1 }, { 1, 2 } }));
}

TEST(Segmentation, LeftRight)
{  
  EXPECT_NO_FATAL_FAILURE(check(
  { { -1, 1 }, { 1, 1 }, { 1, 2 } },
  { { 1, 1 }, { 1, 2 } }));
}

TEST(Segmentation, LeftRightLeft)
{
  EXPECT_NO_FATAL_FAILURE(check(
  { { -1, 1 }, { 1, 2 }, { 1, 3 }, { -1, 4 } },
  { { 1, 2 }, { 1, 3 } }));
}

TEST(Segmentation, RightLeftRight)
{
  EXPECT_NO_FATAL_FAILURE(check(
  { { 1, 1 }, { -1, 2 }, { 1, 3 } },
  { { 1, 3 }, { 1, 1 } }));
}
