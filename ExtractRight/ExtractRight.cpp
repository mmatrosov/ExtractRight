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

class ExtractCopy
{
public:
  std::vector<Point> operator()(const std::vector<Point>& points) const
  {
    auto begin1 = std::find_if    (points.begin(), points.end(), isRight);
    auto end1   = std::find_if_not(begin1,         points.end(), isRight);
    auto begin2 = std::find_if    (end1,           points.end(), isRight);
    auto end2   = std::find_if_not(begin2,         points.end(), isRight);

    if (!(begin2 == end2 || begin1 == points.begin() && end2 == points.end()))
      throw std::runtime_error("Unexpected order");

    std::vector<Point> result;
    result.reserve(end1 - begin1 + end2 - begin2);
    result.insert(result.end(), begin2, end2);
    result.insert(result.end(), begin1, end1);

    return result;
  }
};

class ExtractView
{
public:
  auto operator()(const std::vector<Point>& points) const
  {
    auto begin1 = std::find_if    (points.begin(), points.end(), isRight);
    auto end1   = std::find_if_not(begin1,         points.end(), isRight);
    auto begin2 = std::find_if    (end1,           points.end(), isRight);
    auto end2   = std::find_if_not(begin2,         points.end(), isRight);

    if (!(begin2 == end2 || begin1 == points.begin() && end2 == points.end()))
      throw std::runtime_error("Unexpected order");

    return boost::join(boost::make_iterator_range(begin2, end2),
                       boost::make_iterator_range(begin1, end1));
  }
};

class ExtractViewGeneric
{
public:
  template<class It, class Predicate>
  auto operator()(It first, It last, Predicate p) const
  {
    It begin1 = std::find_if    (first,  last, p);
    It end1   = std::find_if_not(begin1, last, p);
    It begin2 = std::find_if    (end1,   last, p);
    It end2   = std::find_if_not(begin2, last, p);

    if (!(begin2 == end2 || begin1 == first && end2 == last))
      throw std::runtime_error("Unexpected order");

    return boost::join(boost::make_iterator_range(begin2, end2),
                       boost::make_iterator_range(begin1, end1));
  }

  auto operator()(const std::vector<Point>& points) const
  {
    return operator()(points.begin(), points.end(), isRight);
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

  template<class OtherIt>
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
    auto begin1 = std::find_if    (points.begin(), points.end(), isRight);
    auto end1   = std::find_if_not(begin1,         points.end(), isRight);
    auto begin2 = std::find_if    (end1,           points.end(), isRight);
    auto end2   = std::find_if_not(begin2,         points.end(), isRight);

    if (!(begin2 == end2 || begin1 == points.begin() && end2 == points.end()))
      throw std::runtime_error("Unexpected order");

    auto middle = begin2 == end2 ? begin1 : begin2;
    auto beginRes = makeWrappingIterator(middle, points.begin(), points.end());

    size_t count = end1 - begin1 + end2 - begin2;
    return boost::make_iterator_range(beginRes, beginRes + count);
  }
};

class GatherNaive
{
public:
  template<class It>
  void operator()(It first, It last, It begin1, It end1, It begin2, It end2)
  {
    auto middle = begin2 == end2 ? begin1 : begin2;
    std::rotate(first, middle, last);
  }
};

class GatherSmart
{
public:
  template<class It>
  void operator()(It first, It last, It begin1, It end1, It begin2, It end2)
  {
    assert(begin2 == end2 || begin1 == first && end2 == last);
    if (begin2 == end2)
    {
      if (begin1 != first)
        std::move(begin1, end1, first);
    }
    else
    {
      auto len2 = std::distance(begin2, end2);
      auto lenFree = std::distance(end1, begin2);
      if (len2 <= lenFree)
      {
        auto len1 = std::distance(begin1, end1);
        std::move_backward(begin1, end1, first + len1 + len2);
        std::move(begin2, end2, first);
      }
      else
        std::rotate(first, begin2, last);
    }
  }
};

template<class Gather>
class ExtractInplace
{
public:
  void operator()(std::vector<Point>& points) const
  {
    auto begin1 = std::find_if    (points.begin(), points.end(), isRight);
    auto end1   = std::find_if_not(begin1,         points.end(), isRight);
    auto begin2 = std::find_if    (end1,           points.end(), isRight);
    auto end2   = std::find_if_not(begin2,         points.end(), isRight);

    if (!(begin2 == end2 || begin1 == points.begin() && end2 == points.end()))
      throw std::runtime_error("Unexpected order");

    Gather()(points.begin(), points.end(), begin1, end1, begin2, end2);

    size_t count = end1 - begin1 + end2 - begin2;
    points.erase(points.begin() + count, points.end());
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
  EXPECT_TRUE(ExtractCopy()(input) == answer);
  EXPECT_TRUE(ExtractViewWrappingIterator()(input) == answer);
  EXPECT_TRUE(ExtractView()(input) == answer);
  EXPECT_TRUE(ExtractViewGeneric()(input) == answer);

  auto temp = input;
  ExtractInplace<GatherSmart>()(temp);
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
  EXPECT_THROW(ExtractCopy()(input), std::runtime_error);
  EXPECT_THROW(ExtractViewWrappingIterator()(input), std::runtime_error);
  EXPECT_THROW(ExtractView()(input), std::runtime_error);
  EXPECT_THROW(ExtractViewGeneric()(input), std::runtime_error);

  auto temp = input;
  EXPECT_THROW(ExtractInplace<GatherSmart>()(temp), std::runtime_error);
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
void testCopy(benchmark::State& state)
{
  const auto points = getTestArray();

  for (auto _ : state)
  {
    benchmark::DoNotOptimize(T()(points));
  }
}
BENCHMARK_TEMPLATE(testCopy, ExtractNaive)->Apply(setupExtractBenchmark);
BENCHMARK_TEMPLATE(testCopy, ExtractCopy)->Apply(setupExtractBenchmark);

template<class T>
void testView(benchmark::State& state)
{
  const auto points = getTestArray();

  for (auto _ : state)
  {
    benchmark::DoNotOptimize(T()(points));
  }
}
BENCHMARK_TEMPLATE(testView, ExtractView)->Apply(setupExtractBenchmark);
BENCHMARK_TEMPLATE(testView, ExtractViewGeneric)->Apply(setupExtractBenchmark);
BENCHMARK_TEMPLATE(testView, ExtractViewWrappingIterator)->Apply(setupExtractBenchmark);

template<template<class> class T, class F>
void testInplace(benchmark::State& state)
{
  const auto points = getTestArray();

  for (auto _ : state)
  {
    state.PauseTiming();
    auto temp = points;
    state.ResumeTiming();
    T<F>()(temp);
    benchmark::DoNotOptimize(temp);
  }
}
BENCHMARK_TEMPLATE(testInplace, ExtractInplace, GatherNaive)->Apply(setupExtractBenchmark);
BENCHMARK_TEMPLATE(testInplace, ExtractInplace, GatherSmart)->Apply(setupExtractBenchmark);

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
    boost::range::copy(range, result.begin());
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK_TEMPLATE(traverse, ExtractCopy)->Apply(setupTraverseBenchmark);
BENCHMARK_TEMPLATE(traverse, ExtractView)->Apply(setupTraverseBenchmark);
BENCHMARK_TEMPLATE(traverse, ExtractViewGeneric)->Apply(setupTraverseBenchmark);
BENCHMARK_TEMPLATE(traverse, ExtractViewWrappingIterator)->Apply(setupTraverseBenchmark);

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
