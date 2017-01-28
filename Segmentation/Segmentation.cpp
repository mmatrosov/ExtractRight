#include <boost/range/algorithm.hpp>
#include <boost/algorithm/cxx11/is_partitioned.hpp>
#include <boost/algorithm/cxx11/partition_point.hpp>
#include <boost/iterator/iterator_adaptor.hpp>

#include <gtest/gtest.h>

#include <vector>

struct Point
{
  double x;
  double y;
};

inline bool operator==(const Point& a, const Point& b)
{
  return a.x == b.x && a.y == b.y;
}

std::vector<Point> segmentationNaive(const std::vector<Point>& points)
{
  std::vector<Point> result;
  result.clear();

  if (points.size() == 0)
    return result;

  int p = 0;
  for (int i = 1; i < points.size(); ++i)
    if (points[i - 1].x < 0 && points[i].x >= 0)
    {
      p = i;
      break;
    }

  int q = 0;
  for (int i = 1; i < points.size(); ++i)
    if (points[i - 1].x >= 0 && points[i].x < 0)
    {
      q = i;
      break;
    }

  if (p == q)
  {
    if (points[0].x >= 0)
      return points;
    else
      return result;
  }

  int i = p;
  while (i != q)
  {
    if (points[i].x < 0)
    {
      result.clear();
      Point nan;
      nan.x = sqrt(-1);
      nan.y = sqrt(-1);
      result.push_back(nan);
      return result;
    }
    result.push_back(points[i]);
    if (++i >= points.size())
      i = 0;
  }

  i = q;
  while (i != p)
  {
    if (points[i].x >= 0)
    {
      result.clear();
      Point nan;
      nan.x = sqrt(-1);
      nan.y = sqrt(-1);
      result.push_back(nan);
      return result;
    }
    if (++i >= points.size())
      i = 0;
  }

  return std::move(result);
}

std::vector<Point> segmentationNaiveRefactored(const std::vector<Point>& points)
{
  auto isRight = [](const Point& pt) { return pt.x >= 0; };

  int p = 0;
  int q = 0;
  for (int i = 1; i < points.size(); ++i)
  {
    if (!isRight(points[i - 1]) && isRight(points[i]))
      p = i;
    if (isRight(points[i - 1]) && !isRight(points[i]))
      q = i;
  }

  std::vector<Point> result;

  if (p == q)
  {
    if (!points.empty() && isRight(points[0]))
      result = points;
    return result;
  }

  auto iterate = [&](int from, int to, bool shouldBeRight)
  {
    int i = from;
    while (i != to)
    {
      if (isRight(points[i]) != shouldBeRight)
        throw std::runtime_error("Unexpected order");
      if (shouldBeRight)
        result.push_back(points[i]);
      if (++i >= points.size())
        i = 0;
    }
  };

  iterate(p, q, true);
  iterate(q, p, false);

  return result;
}

std::vector<Point> segmentationMirrada(const std::vector<Point>& points)
{
  using PointsIterator = std::vector<Point>::const_iterator;

  auto getLastNegative = [](PointsIterator i, PointsIterator end)
  {
    while (i != end && i->x < 0) { ++i; }
    return i;
  };

  auto getLastPositive = [](PointsIterator i, PointsIterator end)
  {
    while (i != end && i->x >= 0) { ++i; }
    return i;
  };

  std::vector<Point> segment;
  PointsIterator secondStart = points.begin();
  PointsIterator secondEnd = getLastPositive(secondStart, points.end());
  PointsIterator firstStart = getLastNegative(secondEnd, points.end());
  PointsIterator firstEnd = getLastPositive(firstStart, points.end());
  PointsIterator lastIt = (secondStart == secondEnd) ? getLastNegative(firstEnd, points.end()) : firstEnd;

  if (lastIt != points.end())
  {
    throw std::runtime_error("Unexpected order");
  }

  segment.insert(segment.end(), firstStart, firstEnd);
  segment.insert(segment.end(), secondStart, secondEnd);

  return segment;
}

std::vector<Point> segmentation(std::vector<Point> points)
{
  using namespace boost::range;
  using namespace boost::algorithm;

  auto isRight = [](const Point& pt) { return pt.x >= 0; };

  auto middle = adjacent_find(points,
    [&](auto&& pt1, auto&& pt2) { return !isRight(pt1) && isRight(pt2); });

  if (middle != points.end())
    rotate(points, std::next(middle));

  if (!is_partitioned(points, isRight))
    throw std::runtime_error("Unexpected order");

  points.erase(partition_point(points, isRight), points.end());

  return points;
}

template<class It>
class WrappingIterator : public boost::iterator_adaptor<WrappingIterator<It>, It>
{
public:
  WrappingIterator() = default;
  WrappingIterator(It it, It begin, It end) : iterator_adaptor_(it), m_begin(begin), m_end(end) {}

private:
  friend class boost::iterator_core_access;

  auto dereference() const
  {
    return *(m_begin + (base_reference() - m_begin) % (m_end - m_begin));
  }

  It m_begin;
  It m_end;
};

template<class It>
auto makeWrappingIterator(It it, It begin, It end)
{
  return WrappingIterator<It>(it, begin, end);
}

auto segmentationIter(const std::vector<Point>& points)
{
  using namespace boost::range;
  using namespace boost::algorithm;

  auto isRight = [](const Point& pt) { return pt.x >= 0; };

  auto middle = adjacent_find(points,
    [&](auto&& pt1, auto&& pt2) { return !isRight(pt1) && isRight(pt2); });

  middle = middle != points.end() ? std::next(middle) : points.begin();

  auto begin = makeWrappingIterator(middle, points.begin(), points.end());
  auto end = begin + points.size();

  if (!std::is_partitioned(begin, end, isRight))
    throw std::runtime_error("Unexpected order");

  end = std::partition_point(begin, end, isRight);

  return boost::make_iterator_range(begin, end);
}

void checkAnswer(const std::vector<Point>& input, const std::vector<Point>& answer)
{
  EXPECT_TRUE(segmentationNaive(input) == answer);
  EXPECT_TRUE(segmentationNaiveRefactored(input) == answer);
  EXPECT_TRUE(segmentation(input) == answer);
  EXPECT_TRUE(segmentationIter(input) == answer);
}

void checkFailure(const std::vector<Point>& input)
{
  auto answer = segmentationNaive(input);
  EXPECT_TRUE(answer.size() == 1);
  EXPECT_TRUE(std::isnan(answer.front().x) && std::isnan(answer.front().y));

  EXPECT_THROW(segmentationNaiveRefactored(input), std::runtime_error);
  EXPECT_THROW(segmentation(input), std::runtime_error);
  EXPECT_THROW(segmentationIter(input), std::runtime_error);
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

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}