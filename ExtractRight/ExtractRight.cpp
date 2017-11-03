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

using namespace boost::range;
using namespace boost::algorithm;

struct Point
{
  using T = int;

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

bool isRight(const Point& pt)
{
  return pt.x >= 0;
};

class ExtractNaive
{
public: 
  #pragma warning (push)
  #pragma warning (disable : 4804)
  const std::vector<Point> operator()(const std::vector<Point>& points) const
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
        Point left;
        left.x = -1;
        left.y = -1;
        result.push_back(left);
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
        Point left;
        left.x = -1;
        left.y = -1;
        result.push_back(left);
        return result;
      }
      if (++i >= points.size())
        i = 0;
    }

    return std::move(result);
  }
#pragma warning (pop)
};

class ExtractRefactored
{
public:
  std::vector<Point> operator()(const std::vector<Point>& points) const
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
};

class ExtractCopyReference
{
public:
  std::vector<Point> operator()(const std::vector<Point>& points) const
  {
    auto begin1 = points.begin();
    auto end1 = find_if_not(begin1, points.end(), isRight);
    auto begin2 = find_if(end1, points.end(), isRight);
    auto end2 = find_if_not(begin2, points.end(), isRight);
    auto endPts = begin1 == end1 ? find_if(end2, points.end(), isRight) : end2;

    if (endPts != points.end())
      throw std::runtime_error("Unexpected order");

    std::vector<Point> result;
    result.insert(result.end(), begin2, end2);
    result.insert(result.end(), begin1, end1);

    return result;
  }
};

class ExtractCopyRotate
{
public:
  std::vector<Point> operator()(const std::vector<Point>& points) const
  {
    auto isRight = [](const Point& pt) { return pt.x >= 0; };

    auto middle = adjacent_find(points,
      [&](auto&& pt1, auto&& pt2) { return !isRight(pt1) && isRight(pt2); });
    middle = middle != points.end() ? std::next(middle) : points.begin();

    std::vector<Point> result(points.size());
    rotate_copy(points, middle, result.begin());

    if (!is_partitioned(result, isRight))
      throw std::runtime_error("Unexpected order");

    result.erase(partition_point(result, isRight), result.end());

    return result;
  }
};

class ExtractInplace
{
public:
  void operator()(std::vector<Point>& points) const
  {
    auto isRight = [](const Point& pt) { return pt.x >= 0; };

    auto middle = adjacent_find(points,
      [&](auto&& pt1, auto&& pt2) { return !isRight(pt1) && isRight(pt2); });
    middle = middle != points.end() ? std::next(middle) : points.begin();

    rotate(points, middle);

    if (!is_partitioned(points, isRight))
      throw std::runtime_error("Unexpected order");

    points.erase(partition_point(points, isRight), points.end());
  }
};

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

class ExtractViewWrappingIterator
{
public:
  auto operator()(const std::vector<Point>& points) const
  {
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
};

// gsl::span?

class ExtractViewGeneric
{
public:
  auto operator()(const std::vector<Point>& points) const
  {
    return operator()(points.begin(), points.end(), isRight);
  }

  template<class It, class Predicate>
  auto operator()(It first, It last, Predicate p) const
  {
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
};

#define EXPECT_TRUE(x) if(!(x)) throw std::runtime_error("test failed")
#define EXPECT_THROW(x, E) do { try { x; } catch (const E&) { break; } throw std::runtime_error("test failed"); } while (false)

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
  EXPECT_TRUE(ExtractNaive()(input) == answer);
  EXPECT_TRUE(ExtractRefactored()(input) == answer);
  EXPECT_TRUE(ExtractCopyReference()(input) == answer);
  EXPECT_TRUE(ExtractCopyRotate()(input) == answer);
  EXPECT_TRUE(ExtractViewWrappingIterator()(input) == answer);
  EXPECT_TRUE(ExtractViewGeneric()(input) == answer);

  auto temp = input;
  ExtractInplace()(temp);
  EXPECT_TRUE(temp == answer);


  auto copy = input;
  auto cit = makeWrappingIterator(input.begin(), input.begin(), input.end());
  auto it = makeWrappingIterator(copy.begin(), copy.begin(), copy.end());
  cit = it;
  auto i2(it);

  auto inputList = std::list<Point>(input.begin(), input.end());
  auto outputRange = ExtractViewGeneric()(inputList.begin(), inputList.end(), isRight);
  EXPECT_TRUE(outputRange == answer);

  Point p{};
  p = *outputRange.begin();
  *outputRange.begin() = p;
}

void checkFailure(const std::vector<Point>& input)
{
  auto answer = ExtractNaive()(input);
  EXPECT_TRUE(answer.size() == 1);
  EXPECT_TRUE(answer.front().x < 0 && answer.front().y < 0);

  EXPECT_THROW(ExtractRefactored()(input), std::runtime_error);
  EXPECT_THROW(ExtractCopyReference()(input), std::runtime_error);
  EXPECT_THROW(ExtractCopyRotate()(input), std::runtime_error);
  EXPECT_THROW(ExtractViewWrappingIterator()(input), std::runtime_error);
  EXPECT_THROW(ExtractViewGeneric()(input), std::runtime_error);

  auto temp = input;
  EXPECT_THROW(ExtractInplace()(temp), std::runtime_error);
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
  static const int count = 1'000'000;
  std::vector<Point> points(count, { -1, 1 });
  std::fill_n(points.begin(), count / 4, Point{ 1, 1 });
  std::fill_n(points.rbegin(), count / 4, Point{ 1, 1 });
  return points;
}

void setupExtractBenchmark(benchmark::internal::Benchmark* benchmark)
{
  benchmark->Unit(benchmark::kMicrosecond);
}

template<class T>
void extractCopy(benchmark::State& state)
{
  const auto points = getTestArray();

  for (auto _ : state)
  {
    benchmark::DoNotOptimize(T()(points));
  }
}
BENCHMARK_TEMPLATE(extractCopy, ExtractNaive)->Apply(setupExtractBenchmark);
BENCHMARK_TEMPLATE(extractCopy, ExtractCopyReference)->Apply(setupExtractBenchmark);
BENCHMARK_TEMPLATE(extractCopy, ExtractCopyRotate)->Apply(setupExtractBenchmark);

template<class T>
void extractInplace(benchmark::State& state)
{
  const auto points = getTestArray();

  for (auto _ : state)
  {
    state.PauseTiming();
    auto temp = points;
    state.ResumeTiming();
    T()(temp);
    benchmark::DoNotOptimize(temp);
  }
}
BENCHMARK_TEMPLATE(extractInplace, ExtractInplace)->Apply(setupExtractBenchmark);

template<class T>
void extractView(benchmark::State& state)
{
  const auto points = getTestArray();

  for (auto _ : state)
  {
    benchmark::DoNotOptimize(T()(points));
  }
}
BENCHMARK_TEMPLATE(extractView, ExtractViewGeneric)->Apply(setupExtractBenchmark);
BENCHMARK_TEMPLATE(extractView, ExtractViewWrappingIterator)->Apply(setupExtractBenchmark);

void setupTraverseBenchmark(benchmark::internal::Benchmark* benchmark)
{
  benchmark->Unit(benchmark::kMicrosecond);
}

template<class T>
void traverse(benchmark::State& state)
{
  const auto points = getTestArray();
  auto range = T()(points);
  std::vector<Point> result(range.size());

  for (auto _ : state)
  {
    copy(range, result.begin());
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK_TEMPLATE(traverse, ExtractCopyRotate)->Apply(setupTraverseBenchmark);
BENCHMARK_TEMPLATE(traverse, ExtractViewWrappingIterator)->Apply(setupTraverseBenchmark);
BENCHMARK_TEMPLATE(traverse, ExtractViewGeneric)->Apply(setupTraverseBenchmark);

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
