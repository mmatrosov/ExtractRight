#define _SILENCE_CXX17_RESULT_OF_DEPRECATION_WARNING

#include <boost/range/algorithm.hpp>
#include <boost/range/join.hpp>
#include <boost/algorithm/cxx11/is_partitioned.hpp>
#include <boost/algorithm/cxx11/partition_point.hpp>

#pragma warning (push)
#pragma warning (disable : 4141)
#include <benchmark/benchmark.h>
#pragma warning (pop)

#include <gtest/gtest.h>

#include <range/v3/core.hpp>
#include <range/v3/view/concat.hpp>

#include <experimental/generator>

#include <iostream>
#include <vector>
#include <list>

namespace std
{
  using std::experimental::generator;
}

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
    result.reserve((end1 - begin1) + (end2 - begin2));
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

template<class It>
struct Bounds
{
  It begin1, end1, begin2, end2;
};

template<class It, class Predicate>
Bounds<It> findBounds(It first, It last, Predicate p)
{
  auto begin1 = std::find_if    (first,  last, p);
  auto end1   = std::find_if_not(begin1, last, p);
  auto begin2 = std::find_if    (end1,   last, p);
  auto end2   = std::find_if_not(begin2, last, p);
  return { begin1, end1, begin2, end2 };
}

class ExtractNoCheck
{
public:
  template<class It, class Predicate>
  auto operator()(It first, It last, Predicate p) const
  {
    auto bounds = findBounds(first, last, p);

    return boost::join(boost::make_iterator_range(bounds.begin2, bounds.end2),
                       boost::make_iterator_range(bounds.begin1, bounds.end1));
  }

  auto operator()(const std::vector<Point>& points) const
  {
    return operator()(points.begin(), points.end(), isRight);
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

class ExtractViewGenericRanges
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

    return ranges::view::concat(ranges::range<It>(begin2, end2),
                                ranges::range<It>(begin1, end1));
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

    size_t count = (end1 - begin1) + (end2 - begin2);
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
    if (begin2 == end2) {
      if (begin1 != first)
        std::move(begin1, end1, first);
      return;
    }
    auto len2 = std::distance(begin2, end2);
    auto lenFree = std::distance(end1, begin2);
    if (len2 <= lenFree) {
      auto len1 = std::distance(begin1, end1);
      std::move_backward(begin1, end1, first + len1 + len2);
      std::move(begin2, end2, first);
      return;
    }
    std::rotate(first, begin2, last);
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

    size_t count = (end1 - begin1) + (end2 - begin2);
    points.erase(points.begin() + count, points.end());
  }
};

class ExtractViewCoroutine
{
public:
  // Generic version is slow for some reason, due to templated predicate

  std::generator<Point> operator()(const std::vector<Point>& points) const
  {
    auto begin1 = std::find_if    (points.begin(), points.end(), isRight);
    auto end1   = std::find_if_not(begin1,         points.end(), isRight);
    auto begin2 = std::find_if    (end1,           points.end(), isRight);
    auto end2   = std::find_if_not(begin2,         points.end(), isRight);

    if (!(begin2 == end2 || begin1 == points.begin() && end2 == points.end()))
      throw std::runtime_error("Unexpected order");

    for (auto it = begin2; it != end2; ++it) co_yield *it;
    for (auto it = begin1; it != end1; ++it) co_yield *it;  
  }
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

std::vector<Point> maskToVector(const std::string& mask)
{
  auto nx = boost::range::count(mask, 'X');
  auto ns = boost::range::count(mask, '*');
  if (!(nx == 0 && ns == 0 || ns == 1))
    throw std::logic_error("Incorrect mask!");

  auto result = std::vector<Point>(mask.size(), { -1, 0 });

  int index = -1;  // Not known by default

  for (int pass : { 0, 1 })
  {
    for (size_t i = 0; i < mask.size(); ++i)
    {
      const auto c = mask[i];
      if (c == '*')
      {
        if (index < 0)
          index = 0;  // First pass
        else
          break;  // Second pass
      }
      if (c == 'X' || c == '*')
      {
        if (index >= 0)
          result[i] = Point(1, index++);
      }
    }
    
  }

  return result;
}

void checkAnswer(const std::string& inputMask, const std::string& answerMask)
{
  const auto input = maskToVector(inputMask);
  const auto answer = maskToVector(answerMask);

  EXPECT_EQ(answer, ExtractCopy()(input));
  EXPECT_EQ(answer, ExtractViewWrappingIterator()(input));
  EXPECT_EQ(answer, ExtractView()(input));
  EXPECT_EQ(answer, ExtractViewGeneric()(input));
  EXPECT_EQ(answer, ExtractViewGenericRanges()(input) | ranges::to_vector);
  EXPECT_EQ(answer, ExtractNoCheck()(input));

  auto temp = input;
  ExtractInplace<GatherSmart>()(temp);
  EXPECT_EQ(answer, temp);

  auto copy = input;
  auto cit = makeWrappingIterator(input.begin(), input.begin(), input.end());
  auto it = makeWrappingIterator(copy.begin(), copy.begin(), copy.end());
  cit = it;
  auto i2(it);

  auto inputList = std::list<Point>(input.begin(), input.end());
  auto outputRange = ExtractViewGeneric()(inputList.begin(), inputList.end(), isRight);
  EXPECT_EQ(answer, outputRange);

  if (!outputRange.empty())
  {
    Point p{};
    p = *outputRange.begin();
    *outputRange.begin() = p;
  }

  auto&& gen = ExtractViewCoroutine()(input);
  auto yielded = std::vector<Point>(gen.begin(), gen.end());
  EXPECT_EQ(answer, yielded);
}

void checkFailure(const std::string& inputMask)
{
  const auto input = maskToVector(inputMask);

  EXPECT_THROW(ExtractCopy()(input), std::runtime_error);
  EXPECT_THROW(ExtractViewWrappingIterator()(input), std::runtime_error);
  EXPECT_THROW(ExtractView()(input), std::runtime_error);
  EXPECT_THROW(ExtractViewGeneric()(input), std::runtime_error);
  EXPECT_THROW(ExtractViewGenericRanges()(input), std::runtime_error);

  auto temp = input;
  EXPECT_THROW(ExtractInplace<GatherSmart>()(temp), std::runtime_error);
}

TEST(ExtractTest, RightLeft)
{
  checkAnswer("*X.", "*X");
}

TEST(ExtractTest, LeftRight)
{
  checkAnswer(".*X", "*X");
}

TEST(ExtractTest, LeftRightLeft)
{
  checkAnswer(".*X.", "*X");
}

TEST(ExtractTest, RightLeftRight)
{
  checkAnswer("X.*", "*X");
}

TEST(ExtractTest, OnlyLeft)
{
  checkAnswer("...", "");
}

TEST(ExtractTest, OnlyRight)
{
  checkAnswer("*XX", "*XX");
}

TEST(ExtractTest, Empty)
{
  checkAnswer("", "");
}

TEST(ExtractTest, Lonely)
{
  checkAnswer(".............*.....", "*");
}

TEST(ExtractTest, MostlyAtTheEnd)
{
  checkAnswer("........*XXXXXXXX.", "*XXXXXXXX");
}

TEST(ExtractTest, MostlyAtTheBeginning)
{
  checkAnswer(".*XXXXXXXX........", "*XXXXXXXX");
}

TEST(ExtractTest, OneMissingNearEnd)
{
  checkAnswer("XXXXXXXXXXXXXX.*", "*XXXXXXXXXXXXXX");
}

TEST(ExtractTest, OneMissingNearBeginning)
{
  checkAnswer("X.*XXXXXXXXXXXXX", "*XXXXXXXXXXXXXX");
}

TEST(ExtractTest, Incorrect1)
{
  checkFailure(".*.X");
}

TEST(ExtractTest, Incorrect2)
{
  checkFailure("*.X.");
}

TEST(ExtractTest, Incorrect3)
{
  checkFailure("*XXXX.XXXXXXXX.XXXX");
}

std::vector<Point> getBenchmarkArray()
{
  static const int count = 1'000'000;
  std::vector<Point> points(count, { -1, 1 });
  std::fill_n(points.begin(), count / 4, Point{ 1, 1 });
  std::fill_n(points.rbegin(), count / 4, Point{ 1, 1 });
  return points;
}

void setupBenchmark(benchmark::internal::Benchmark* benchmark)
{
  benchmark->Unit(benchmark::kMicrosecond);
}

template<template<class> class T, class F>
void runInplace(benchmark::State& state)
{
  const auto points = getBenchmarkArray();

  for (auto _ : state)
  {
    state.PauseTiming();
    auto temp = points;
    state.ResumeTiming();
    T<F>()(temp);
    benchmark::DoNotOptimize(temp);
  }
}
BENCHMARK_TEMPLATE(runInplace, ExtractInplace, GatherNaive)->Apply(setupBenchmark);
BENCHMARK_TEMPLATE(runInplace, ExtractInplace, GatherSmart)->Apply(setupBenchmark);

template<class T>
void run(benchmark::State& state)
{
  const auto points = getBenchmarkArray();
  std::vector<Point> result(points.size());

  for (auto _ : state)
  {
    state.PauseTiming();
    auto range = T()(points);
    state.ResumeTiming();
    std::copy(range.begin(), range.end(), result.begin());
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK_TEMPLATE(run, ExtractCopy)->Apply(setupBenchmark);
BENCHMARK_TEMPLATE(run, ExtractView)->Apply(setupBenchmark);
BENCHMARK_TEMPLATE(run, ExtractViewGeneric)->Apply(setupBenchmark);
BENCHMARK_TEMPLATE(run, ExtractViewGenericRanges)->Apply(setupBenchmark);
BENCHMARK_TEMPLATE(run, ExtractViewWrappingIterator)->Apply(setupBenchmark);
BENCHMARK_TEMPLATE(run, ExtractViewCoroutine)->Apply(setupBenchmark);

int main(int argc, char* argv[])
{
  testing::InitGoogleTest(&argc, argv);
  if (auto status = RUN_ALL_TESTS(); status != 0)
    return status;

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();

  return EXIT_SUCCESS;
}
