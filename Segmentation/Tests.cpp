#include "Common.h"

#include <gtest/gtest.h>

std::vector<Point> segmentationNaive(const std::vector<Point>& points);

void checkAnswer(const std::vector<Point>& input, const std::vector<Point>& answer)
{
  ASSERT_TRUE(segmentationNaive(input) == answer);
}

void checkFailure(const std::vector<Point>& input)
{
  auto answer = segmentationNaive(input);
  ASSERT_TRUE(answer.size() == 1);
  ASSERT_TRUE(std::isnan(answer.front().x) && std::isnan(answer.front().y));
}

TEST(Segmentation, RightLeft)
{
  EXPECT_NO_FATAL_FAILURE(checkAnswer(
  { { 1, 1 }, { 1, 2 }, { -1, 1 } },
  { { 1, 1 }, { 1, 2 } }));
}

TEST(Segmentation, LeftRight)
{
  EXPECT_NO_FATAL_FAILURE(checkAnswer(
  { { -1, 1 }, { 1, 1 }, { 1, 2 } },
  { { 1, 1 }, { 1, 2 } }));
}

TEST(Segmentation, LeftRightLeft)
{
  EXPECT_NO_FATAL_FAILURE(checkAnswer(
  { { -1, 1 }, { 1, 2 }, { 1, 3 }, { -1, 4 } },
  { { 1, 2 }, { 1, 3 } }));
}

TEST(Segmentation, RightLeftRight)
{
  EXPECT_NO_FATAL_FAILURE(checkAnswer(
  { { 1, 1 }, { -1, 2 }, { 1, 3 } },
  { { 1, 3 }, { 1, 1 } }));
}

TEST(Segmentation, OnlyLeft)
{
  EXPECT_NO_FATAL_FAILURE(checkAnswer(
  { { -1, 1 }, { -1, 2 }, { -1, 3 } },
  {}));
}

TEST(Segmentation, OnlyRight)
{
  EXPECT_NO_FATAL_FAILURE(checkAnswer(
  { { 1, 1 }, { 1, 2 }, { 1, 3 } },
  { { 1, 1 }, { 1, 2 }, { 1, 3 } }));
}

TEST(Segmentation, Empty)
{
  EXPECT_NO_FATAL_FAILURE(checkAnswer(
  {},
  {}));
}

TEST(Segmentation, Incorrect1)
{
  EXPECT_NO_FATAL_FAILURE(checkFailure(
  { { -1, 1 }, { 1, 1 }, { -1, 1 }, { 1, 2 } }));
}

TEST(Segmentation, Incorrect2)
{
  EXPECT_NO_FATAL_FAILURE(checkFailure(
  { { 1, 2 }, { -1, 1 }, { 1, 1 }, { -1, 1 } }));
}
