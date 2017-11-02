#include <boost/range/algorithm.hpp>
#include <boost/range/join.hpp>
#include <boost/algorithm/cxx11/is_partitioned.hpp>
#include <boost/algorithm/cxx11/partition_point.hpp>

#pragma warning (push)
#pragma warning (disable : 4141)
#include <benchmark/benchmark.h>
#pragma warning (pop)

#include <iostream>
#include <vector>
#include <list>

struct Point
{
  using T = double;

  T x;
  T y;

  Point() = default;
  Point(T x, T y) : x(x), y(y) {}
};

inline bool operator==(const Point& a, const Point& b)
{
  return a.x == b.x && a.y == b.y;
}
inline bool operator!=(const Point& a, const Point& b)
{
  return !(a == b);
}

#pragma warning (push)
#pragma warning (disable : 4804)
const std::vector<Point> extractNaive(const std::vector<Point>& points)
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
#pragma warning (pop)

std::vector<Point> extractRefactored(const std::vector<Point>& points)
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

std::vector<Point> extractManual(const std::vector<Point>& points)
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

std::vector<Point> extractAlgoBasic(std::vector<Point> points)
{
  using namespace boost::range;
  using namespace boost::algorithm;

  auto isRight = [](const Point& pt) { return pt.x >= 0; };

  auto middle = adjacent_find(points,
    [&](auto&& pt1, auto&& pt2) { return !isRight(pt1) && isRight(pt2); });
  middle = middle != points.end() ? std::next(middle) : points.begin();

  rotate(points, middle);

  if (!is_partitioned(points, isRight))
    throw std::runtime_error("Unexpected order");

  points.erase(partition_point(points, isRight), points.end());

  return points;
}

template<class It>
class WrappingIterator : public boost::iterator_facade<WrappingIterator<It>, 
                                                       typename It::value_type, 
                                                       boost::random_access_traversal_tag, 
                                                       typename It::reference>
{
public:
  WrappingIterator() = default;
  WrappingIterator(It it, It begin, It end) : 
    m_begin(begin), m_size(end - begin), m_offset(it - begin) {}

  template <class OtherIt>
  WrappingIterator(const WrappingIterator<OtherIt>& other) :
    m_begin(other.m_begin), m_size(other.m_size), m_offset(other.m_offset) {}

private:
  friend class boost::iterator_core_access;

  template<class> friend class WrappingIterator;

  using Base = boost::iterator_facade<WrappingIterator<It>, 
                                      typename It::value_type, 
                                      boost::random_access_traversal_tag, 
                                      typename It::reference>;

  typename Base::reference dereference() const
  {
    return *(m_begin + (m_offset < m_size ? m_offset : m_offset - m_size));
  }

  template <class OtherIt>
  bool equal(const WrappingIterator<OtherIt>& other) const
  {
    assert(other.m_begin == m_begin && other.m_size == m_size);
    return other.m_offset == m_offset;
  }

  void advance(typename Base::difference_type n)
  {
    m_offset += n;
  }

  void increment()
  {
    ++m_offset;
  }

  void decrement()
  {
    --m_offset;
  }

  template <class OtherIt>
  typename Base::difference_type distance_to(const WrappingIterator<OtherIt>& other) const
  {
    assert(other.m_begin == m_begin && other.m_size == m_size);
    return other.m_offset - m_offset;
  }

  It m_begin;
  size_t m_size;
  size_t m_offset;
};

template<class It>
auto makeWrappingIterator(It it, It begin, It end)
{
  return WrappingIterator<It>(it, begin, end);
}

auto extractAlgoWrappingIterator(const std::vector<Point>& points)
{
  using namespace boost::range;
  using namespace boost::algorithm;

  auto isRight = [](const Point& pt) { return pt.x >= 0; };

  auto middle = adjacent_find(points,
    [&](auto&& pt1, auto&& pt2) { return !isRight(pt1) && isRight(pt2); });

  middle = middle != points.end() ? std::next(middle) : points.begin();

  auto begin = makeWrappingIterator(middle, points.begin(), points.end());
  auto rotated = boost::make_iterator_range(begin, begin + points.size());

  if (!is_partitioned(rotated, isRight))
    throw std::runtime_error("Unexpected order");

  auto end = partition_point(rotated, isRight);

  return boost::make_iterator_range(begin, end);
}

// gsl::span?

template<class It, class Predicate>
auto extractAlgoGeneric(It first, It last, Predicate p)
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

void used()
{
#if 0
  using namespace boost::range;
  using namespace boost::algorithm;

  // lambda functions (including generic)
  // ternary operator
  // exceptions
  // transient parameters
  std::vector<int>::empty();
  adjacent_find();
  rotate();
  is_partitioned();
  partition_point();
  std::next();

  // custom make-function
  // template parameters for iterators
  // template parameters for predicates
  // function return type deduction 
  boost::range;
  boost::algorithm;
  boost::iterator_adaptor<int> a;
  boost::make_iterator_range();
  boost::join();
#endif
}

void checkAnswer(const std::vector<Point>& input, const std::vector<Point>& answer)
{
  EXPECT_TRUE(extractNaive(input) == answer);
  EXPECT_TRUE(extractRefactored(input) == answer);
  EXPECT_TRUE(extractManual(input) == answer);
  EXPECT_TRUE(extractAlgoBasic(input) == answer);
  EXPECT_TRUE(extractAlgoWrappingIterator(input) == answer);

  EXPECT_TRUE(extractAlgoGeneric(input.begin(), input.end(), isRight) == answer);

  auto copy = input;
  auto cit = makeWrappingIterator(input.begin(), input.begin(), input.end());
  auto it = makeWrappingIterator(copy.begin(), copy.begin(), copy.end());
  cit = it;
  auto i2(it);

  auto inputList = std::list<Point>(input.begin(), input.end());
  auto outputRange = extractAlgoGeneric(inputList.begin(), inputList.end(), isRight);
  EXPECT_TRUE(outputRange == answer);

  Point p{};
  p = *outputRange.begin();
  *outputRange.begin() = p;
}

void checkFailure(const std::vector<Point>& input)
{
  auto answer = extractNaive(input);
  EXPECT_TRUE(answer.size() == 1);
  EXPECT_TRUE(std::isnan(answer.front().x) && std::isnan(answer.front().y));

  EXPECT_THROW(extractRefactored(input), std::runtime_error);
  EXPECT_THROW(extractManual(input), std::runtime_error);
  EXPECT_THROW(extractAlgoBasic(input), std::runtime_error);
  EXPECT_THROW(extractAlgoWrappingIterator(input), std::runtime_error);
  EXPECT_THROW(extractAlgoGeneric(input.begin(), input.end(), isRight), std::runtime_error);
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

std::vector<Point> getTestArray()
{
  static const int count = 50'000'000;
  std::vector<Point> points(count, { -1, 1 });
  std::fill_n(points.begin(), count / 4, Point{ 1, 1 });
  std::fill_n(points.rbegin(), count / 4, Point{ 1, 1 });
  return points;
}

constexpr int numIterations = 3;

void BM_testNaive(benchmark::State& state)
{
  const auto points = getTestArray();

  for (auto _ : state)
  {
    benchmark::DoNotOptimize(extractNaive(points));
  }
}
BENCHMARK(BM_testNaive)->Unit(benchmark::kMillisecond)->Iterations(numIterations);

void BM_testManual(benchmark::State& state)
{
  const auto points = getTestArray();

  for (auto _ : state)
  {
    benchmark::DoNotOptimize(extractManual(points));
  }
}
BENCHMARK(BM_testManual)->Unit(benchmark::kMillisecond)->Iterations(numIterations);

void BM_testAlgoBasic(benchmark::State& state)
{
  const auto points = getTestArray();

  for (auto _ : state)
  {
    benchmark::DoNotOptimize(extractAlgoBasic(points));
  }
}
BENCHMARK(BM_testAlgoBasic)->Unit(benchmark::kMillisecond)->Iterations(numIterations);

void BM_testAlgoWrappingIterator(benchmark::State& state)
{
  const auto points = getTestArray();

  for (auto _ : state)
  {
    benchmark::DoNotOptimize(boost::copy_range<std::vector<Point>>(extractAlgoWrappingIterator(points)));
  }
}
BENCHMARK(BM_testAlgoWrappingIterator)->Unit(benchmark::kMillisecond)->Iterations(numIterations);

void BM_testAlgoGeneric(benchmark::State& state)
{
  const auto points = getTestArray();

  for (auto _ : state)
  {
    benchmark::DoNotOptimize(boost::copy_range<std::vector<Point>>(extractAlgoGeneric(points.begin(), points.end(), isRight)));
  }
}
BENCHMARK(BM_testAlgoGeneric)->Unit(benchmark::kMillisecond)->Iterations(numIterations);

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

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();

  return EXIT_SUCCESS;
}