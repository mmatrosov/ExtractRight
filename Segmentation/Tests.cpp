#include "Common.h"

#include <gtest/gtest.h>

std::vector<Point> segmentationNaive(const std::vector<Point>& points);
std::vector<Point> segmentationNaiveRefactored(const std::vector<Point>& points);
std::vector<Point> segmentation(std::vector<Point> points);

void checkAnswer(const std::vector<Point>& input, const std::vector<Point>& answer)
{
  EXPECT_TRUE(segmentationNaive(input) == answer);
  EXPECT_TRUE(segmentationNaiveRefactored(input) == answer);
  EXPECT_TRUE(segmentation(input) == answer);
}

void checkFailure(const std::vector<Point>& input)
{
  auto answer = segmentationNaive(input);
  EXPECT_TRUE(answer.size() == 1);
  EXPECT_TRUE(std::isnan(answer.front().x) && std::isnan(answer.front().y));

  EXPECT_THROW(segmentationNaiveRefactored(input), std::runtime_error);
  EXPECT_THROW(segmentation(input), std::runtime_error);
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
