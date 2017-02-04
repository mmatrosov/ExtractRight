#include <boost/range/algorithm.hpp>
#include <boost/range/join.hpp>
#include <boost/algorithm/cxx11/is_partitioned.hpp>
#include <boost/algorithm/cxx11/partition_point.hpp>
#include <boost/iterator/iterator_adaptor.hpp>

#include <iostream>
#include <vector>
#include <list>

struct Point
{
  double x;
  double y;

  Point() = default;
  Point(double x, double y) : x(x), y(y) {}
};

inline bool operator==(const Point& a, const Point& b)
{
  return a.x == b.x && a.y == b.y;
}
inline bool operator!=(const Point& a, const Point& b)
{
  return !(a == b);
}

const std::vector<Point> segmentationNaive(const std::vector<Point>& points)
{
  std::vector<Point> result;
  result.clear();

  if (points.size() == 0)
    return result;

  int p = 0;
  bool found = false;
  for (int i = 1; i < points.size() && ~found; ++i)
    if (points[i - 1].x < 0 && points[i].x >= 0)
    {
      p = i;
      found = true;
    }

  int q = 0;
  found = false;
  for (int i = 1; i < points.size() && ~found; ++i)
    if (points[i - 1].x >= 0 && points[i].x < 0)
    {
      q = i;
      found = true;
    }

  if (p == q)
  {
    if ((*points.begin()).x >= 0)
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
  std::vector<Point> result;

  if (points.empty())
    return result;

  auto isRight = [](const Point& pt) { return pt.x >= 0; };

  auto findBoundary = [&](bool rightToLeft)
  {
    for (int i = 1; i < points.size(); ++i)
      if (isRight(points[i - 1]) == rightToLeft && 
          isRight(points[i]) != rightToLeft)
        return i;
    return 0;
  };

  int p = findBoundary(false);
  int q = findBoundary(true);

  if (p == q)
    return isRight(points[0]) ? points : result;

  auto appendResult = [&](int from, int to, bool shouldBeRight)
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

  appendResult(p, q, true);
  appendResult(q, p, false);

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

// Provide pointers instead of iterators to avoid debug range checks in VS
template<class It>
class WrappingIterator : public boost::iterator_adaptor<WrappingIterator<It>, It>
{
  using Base = boost::iterator_adaptor<WrappingIterator<It>, It>;

public:
  WrappingIterator() = default;
  WrappingIterator(It it, It begin, It end) : 
    Base(it), m_begin(begin), m_size(end - begin) {}

private:
  friend class boost::iterator_core_access;

  typename Base::reference dereference() const  // Can use decltype(auto) instead, but no just auto!
  {
    return *(m_begin + (this->base_reference() - m_begin) % m_size);
  }

  It m_begin;
  size_t m_size;
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

  if (!std::is_partitioned(begin, end, isRight))  // Note how std:: is necessary
    throw std::runtime_error("Unexpected order");

  end = std::partition_point(begin, end, isRight);

  return boost::make_iterator_range(begin, end);
}

// gsl::span?

template<class It, class Predicate>
auto segmentationRange(It first, It last, Predicate p)
{
  using namespace boost::range;
  using namespace boost::algorithm;

  auto middle = adjacent_find(first, last,
    [&](auto&& a, auto&& b) { return !p(a) && p(b); });

  middle = middle != last ? std::next(middle) : first;

  auto rotated = boost::join(
    boost::make_iterator_range(middle, last), 
    boost::make_iterator_range(first, middle));

  if (!is_partitioned(rotated, p))
    throw std::runtime_error("Unexpected order");

  auto end = partition_point(rotated, p);

  return boost::make_iterator_range(rotated.begin(), end);
}

#define EXPECT_TRUE(x) if(!(x)) throw std::runtime_error("test failed")
#define EXPECT_THROW(x, E) do { try { x; } catch (const E&) { break; } throw std::runtime_error("test failed"); } while (false)

bool isRight(const Point& pt)
{
  return pt.x >= 0;
};


void checkAnswer(const std::vector<Point>& input, const std::vector<Point>& answer)
{
  EXPECT_TRUE(segmentationNaive(input) == answer);
  EXPECT_TRUE(segmentationNaiveRefactored(input) == answer);
  EXPECT_TRUE(segmentation(input) == answer);
  EXPECT_TRUE(segmentationIter(input) == answer);

  EXPECT_TRUE(segmentationRange(input.begin(), input.end(), isRight) == answer);

  auto inputList = std::list<Point>(input.begin(), input.end());
  auto outputRange = segmentationRange(inputList.begin(), inputList.end(), isRight);
  EXPECT_TRUE(outputRange == answer);

  Point p{};
  p = *outputRange.begin();
  *outputRange.begin() = p;
}

void checkFailure(const std::vector<Point>& input)
{
  auto answer = segmentationNaive(input);
  EXPECT_TRUE(answer.size() == 1);
  EXPECT_TRUE(std::isnan(answer.front().x) && std::isnan(answer.front().y));

  EXPECT_THROW(segmentationNaiveRefactored(input), std::runtime_error);
  EXPECT_THROW(segmentation(input), std::runtime_error);
  EXPECT_THROW(segmentationIter(input), std::runtime_error);
  EXPECT_THROW(segmentationRange(input.begin(), input.end(), isRight), std::runtime_error);
}

void testRightLeft()
{
  checkAnswer(
  { { 1, 1 }, { 1, 2 }, { -1, 3 } },
  { { 1, 1 }, { 1, 2 } });
}

void testLeftRight()
{
  checkAnswer(
  { { -1, 1 }, { 1, 2 }, { 1, 3 } },
  { { 1, 2 }, { 1, 3 } });
}

void testLeftRightLeft()
{
  checkAnswer(
  { { -1, 1 }, { 1, 2 }, { 1, 3 }, { -1, 4 } },
  { { 1, 2 }, { 1, 3 } });
}

void testRightLeftRight()
{
  checkAnswer(
  { { 1, 1 }, { -1, 2 }, { 1, 3 } },
  { { 1, 3 }, { 1, 1 } });
}

void testOnlyLeft()
{
  checkAnswer(
  { { -1, 1 }, { -1, 2 }, { -1, 3 } },
  {});
}

void testOnlyRight()
{
  checkAnswer(
  { { 1, 1 }, { 1, 2 }, { 1, 3 } },
  { { 1, 1 }, { 1, 2 }, { 1, 3 } });
}

void testEmpty()
{
  checkAnswer(
  {},
  {});
}

void testIncorrect1()
{
  checkFailure(
  { { -1, 1 }, { 1, 2 }, { -1, 3 }, { 1, 4 } });
}

void testIncorrect2()
{
  checkFailure(
  { { 1, 1 }, { -1, 2 }, { 1, 3 }, { -1, 4 } });
}

int main(int argc, char* argv[])
{
  try
  {
    testRightLeft();
    testLeftRight();
    testLeftRightLeft();
    testRightLeftRight();
    testOnlyLeft();
    testOnlyRight();
    testEmpty();
    testIncorrect1();
    testIncorrect2();
  }
  catch (const std::runtime_error&)
  {
    std::cout << "Failure!" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Success." << std::endl;

  return EXIT_SUCCESS;
}