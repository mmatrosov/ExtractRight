#include "ExtractNoCheck.hpp"

#pragma warning (push)
#pragma warning (disable : 4141)
#include <benchmark/benchmark.h>
#pragma warning (pop)

#include <boost/range/algorithm/count.hpp>
#include <boost/range/algorithm/count_if.hpp>

#include <gtest/gtest.h>

#ifdef _MSC_VER
#include <range/v3/core.hpp>
#define iterator_range range<It>
#else
#include <range/v3/iterator_range.hpp>
#endif
#include <range/v3/view/concat.hpp>

#include <cctype>
#include <list>

class ExtractCopy
{
public:
  std::vector<Point> extract(const std::vector<Point>& points)
  {
    auto begin1 = std::find_if    (points.begin(), points.end(), isPositive);
    auto end1   = std::find_if_not(begin1,         points.end(), isPositive);
    auto begin2 = std::find_if    (end1,           points.end(), isPositive);
    auto end2   = std::find_if_not(begin2,         points.end(), isPositive);

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
  auto extract(const std::vector<Point>& points)
  {
    auto begin1 = std::find_if    (points.begin(), points.end(), isPositive);
    auto end1   = std::find_if_not(begin1,         points.end(), isPositive);
    auto begin2 = std::find_if    (end1,           points.end(), isPositive);
    auto end2   = std::find_if_not(begin2,         points.end(), isPositive);

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
  auto extract(It first, It last, Predicate p)
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

  auto extract(const std::vector<Point>& points)
  {
    return extract(points.begin(), points.end(), isPositive);
  }
};

class ExtractViewRanges
{
public:
  template<class It, class Predicate>
  auto extract(It first, It last, Predicate p)
  {
    It begin1 = std::find_if    (first,  last, p);
    It end1   = std::find_if_not(begin1, last, p);
    It begin2 = std::find_if    (end1,   last, p);
    It end2   = std::find_if_not(begin2, last, p);

    if (!(begin2 == end2 || begin1 == first && end2 == last))
      throw std::runtime_error("Unexpected order");

    return ranges::view::concat(ranges::iterator_range(begin2, end2),
                                ranges::iterator_range(begin1, end1));
  }

  auto extract(const std::vector<Point>& points)
  {
    return extract(points.begin(), points.end(), isPositive);
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

#ifdef _MSC_VER
template<class It>
auto makeWrappingIterator(It it, It begin, It end)
{
  return WrappingIterator<It>(it, begin, end);
}
#define WrappingIterator makeWrappingIterator
#endif

class ExtractViewWrappingIterator
{
public:
  auto extract(const std::vector<Point>& points)
  {
    auto begin1 = std::find_if    (points.begin(), points.end(), isPositive);
    auto end1   = std::find_if_not(begin1,         points.end(), isPositive);
    auto begin2 = std::find_if    (end1,           points.end(), isPositive);
    auto end2   = std::find_if_not(begin2,         points.end(), isPositive);

    if (!(begin2 == end2 || begin1 == points.begin() && end2 == points.end()))
      throw std::runtime_error("Unexpected order");

    auto begin = WrappingIterator(begin2 == end2 ? begin1 : begin2,
                                  points.begin(), points.end());

    size_t count = (end1 - begin1) + (end2 - begin2);
    return boost::make_iterator_range(begin, begin + count);
  }
};

class GatherNaive
{
public:
  template<class It>
  void gather(It first, It last, It begin1, It end1, It begin2, It end2)
  {
    assert(begin2 == end2 || begin1 == first && end2 == last);
    auto middle = begin2 == end2 ? begin1 : begin2;
    std::rotate(first, middle, last);
  }
};

class GatherSmart
{
public:
  template<class It>
  void gather(It first, It last, It begin1, It end1, It begin2, It end2)
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
class ExtractMove
{
public:
  std::vector<Point> extract(std::vector<Point>&& points)
  {
    auto begin1 = std::find_if    (points.begin(), points.end(), isPositive);
    auto end1   = std::find_if_not(begin1,         points.end(), isPositive);
    auto begin2 = std::find_if    (end1,           points.end(), isPositive);
    auto end2   = std::find_if_not(begin2,         points.end(), isPositive);

    if (!(begin2 == end2 || begin1 == points.begin() && end2 == points.end()))
      throw std::runtime_error("Unexpected order");

    Gather().gather(points.begin(), points.end(), begin1, end1, begin2, end2);

    points.resize((end1 - begin1) + (end2 - begin2));

    return std::move(points);  // Note move! https://stackoverflow.com/a/29128776/261217
  }
};

class ExtractList
{
public:
  std::list<Point> extract(std::list<Point>&& points)
  {
    auto begin1 = std::find_if    (points.begin(), points.end(), isPositive);
    auto end1   = std::find_if_not(begin1,         points.end(), isPositive);
    auto begin2 = std::find_if    (end1,           points.end(), isPositive);
    auto end2   = std::find_if_not(begin2,         points.end(), isPositive);

    if (!(begin2 == end2 || begin1 == points.begin() && end2 == points.end()))
      throw std::runtime_error("Unexpected order");

    points.erase(points.begin(), begin1);
    points.splice(points.begin(), points, begin2, end2);
    points.erase(end1, points.end());

    return std::move(points);
  }
};

#ifdef __clang__
#include <cppcoro/generator.hpp>
namespace std
{
  using cppcoro::generator;
}
#endif

#ifdef _MSC_VER
#include <experimental/generator>
namespace std
{
  using std::experimental::generator;
}
#endif

#if defined(__clang__) || defined(_MSC_VER)
class ExtractViewCoroutine
{
public:
  // Generic version is slow for some reason, due to templated predicate

  std::generator<const Point> extract(const std::vector<Point>& points)
  {
    auto begin1 = std::find_if    (points.begin(), points.end(), isPositive);
    auto end1   = std::find_if_not(begin1,         points.end(), isPositive);
    auto begin2 = std::find_if    (end1,           points.end(), isPositive);
    auto end2   = std::find_if_not(begin2,         points.end(), isPositive);

    if (!(begin2 == end2 || begin1 == points.begin() && end2 == points.end()))
      throw std::runtime_error("Unexpected order");

    for (auto it = begin2; it != end2; ++it) co_yield *it;
    for (auto it = begin1; it != end1; ++it) co_yield *it;  
  }
};
#else
#define ExtractViewCoroutine ExtractViewGeneric
#endif

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

  EXPECT_EQ(answer, ExtractCopy().extract(input));
  EXPECT_EQ(answer, ExtractViewWrappingIterator().extract(input));
  EXPECT_EQ(answer, ExtractView().extract(input));
  EXPECT_EQ(answer, ExtractViewGeneric().extract(input));
  EXPECT_EQ(answer, ExtractViewRanges().extract(input) | ranges::to_vector);
  EXPECT_EQ(answer, ExtractNoCheck().extract(input));
  EXPECT_EQ(answer, ExtractNoCheckSimple().extract(input));
  EXPECT_EQ(answer, ExtractMove<GatherSmart>().extract(std::vector<Point>(input)));

  auto copy = input;
  auto cit = WrappingIterator(input.begin(), input.begin(), input.end());
  auto it = WrappingIterator(copy.begin(), copy.begin(), copy.end());
  cit = it;
  auto i2(it);

  auto inputList = std::list<Point>(input.begin(), input.end());
  auto outputRange = ExtractViewGeneric().extract(inputList.begin(), inputList.end(), isPositive);
  EXPECT_EQ(answer, outputRange);

  if (!outputRange.empty())
  {
    Point p{};
    p = *outputRange.begin();
    *outputRange.begin() = p;
  }

  auto outputList = ExtractList().extract(std::move(inputList));
  EXPECT_EQ(answer, std::vector<Point>(outputList.begin(), outputList.end()));

  auto&& gen = ExtractViewCoroutine().extract(input);
  auto yielded = std::vector<Point>(gen.begin(), gen.end());
  EXPECT_EQ(answer, yielded);
}

void checkFailure(const std::string& inputMask)
{
  const auto input = maskToVector(inputMask);

  EXPECT_THROW(ExtractCopy().extract(input), std::runtime_error);
  EXPECT_THROW(ExtractViewWrappingIterator().extract(input), std::runtime_error);
  EXPECT_THROW(ExtractView().extract(input), std::runtime_error);
  EXPECT_THROW(ExtractViewGeneric().extract(input), std::runtime_error);
  EXPECT_THROW(ExtractViewRanges().extract(input), std::runtime_error);
  EXPECT_THROW(ExtractMove<GatherSmart>().extract(std::vector<Point>(input)), std::runtime_error);
  EXPECT_THROW(ExtractList().extract(std::list<Point>(input.begin(), input.end())), std::runtime_error);
}

TEST(ExtractTest, PosNeg)
{
  checkAnswer("*X.", "*X");
}

TEST(ExtractTest, NegPos)
{
  checkAnswer(".*X", "*X");
}

TEST(ExtractTest, NegPosNeg)
{
  checkAnswer(".*X.", "*X");
}

TEST(ExtractTest, PosNegPos)
{
  checkAnswer("X.*", "*X");
}

TEST(ExtractTest, OnlyNeg)
{
  checkAnswer("...", "");
}

TEST(ExtractTest, OnlyPos)
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

void checkFindAny(const std::string& str)
{
  auto isalpha = [](unsigned char c) { return std::isalpha(c); };
  auto it = findAny(str.begin(), str.end(), isalpha);

  if (it != str.end())
    ASSERT_TRUE(std::isalpha(*it));
  else
    ASSERT_EQ(0, boost::range::count_if(str, isalpha));
}

TEST(FindAnyTest, Empty)
{
  checkFindAny("");
}

TEST(FindAnyTest, Absent)
{
  checkFindAny(".");
  checkFindAny("..");
  checkFindAny("...");
  checkFindAny("...............");
  checkFindAny("................");  // 16 chars
  checkFindAny(".................");
}

TEST(FindAnyTest, Existing)
{
  checkFindAny("X");
  checkFindAny("X.");
  checkFindAny(".X");
  checkFindAny("X..");
  checkFindAny(".X.");
  checkFindAny("..X");
  checkFindAny("X..............");
  checkFindAny("X...............");  // 16 chars
  checkFindAny("X................");
  checkFindAny("..............X");
  checkFindAny("...............X");  // 16 chars
  checkFindAny("................X");
  checkFindAny("......XXX......");
  checkFindAny("......XXX.......");  // 16 chars
  checkFindAny("......XXX........");
}

TEST(FindAnyTest, Batch)
{
  for (int len = 0; len < 128; ++len)
  {
    const auto baseStr = std::string(len, '.');
    ASSERT_NO_FATAL_FAILURE(checkFindAny(baseStr));

    for (int i = 0; i < len; ++i)
    {
      auto str = baseStr;
      str[i] = 'X';
      ASSERT_NO_FATAL_FAILURE(checkFindAny(str));
    }
  }
}

void checkFindAnyOptimalSampling(const int count)
{
  int numSamples = 0;
  auto p = [&](char c) { ++numSamples; return std::isalpha(c); };

  auto str = std::string(count, '.');
  auto it = findAny(str.begin(), str.end(), p);

  ASSERT_EQ(str.end(), it);
  ASSERT_EQ(count, numSamples);
}

TEST(FindAnyTest, OptimalSampling)
{
  for (int count = 0; count < 128; ++count)
    ASSERT_NO_FATAL_FAILURE(checkFindAnyOptimalSampling(count));
}

static const int BenchDataSize = 1 << 20;

enum Distribution
{
  BothEnds, Beginning
};

template<class Container>
Container getBenchmarkArray(Distribution dist, int fraction)
{
  static const auto sample = Point{ 1, 1 };
  static const auto hole = Point{ -1, 1 };

  Container points(BenchDataSize, hole);

  const int len = BenchDataSize / fraction;

  if (dist == BothEnds)
  {
    std::fill_n(points.begin(), len / 2, sample);
    std::fill_n(points.rbegin(), len / 2, sample);
  }
  else
  {
    std::fill_n(points.begin(), len, sample);   
  }

  return points;
}

void setupBenchmark(benchmark::internal::Benchmark* benchmark)
{
  benchmark->Unit(benchmark::kMicrosecond)->ArgName("Fraction");
}

enum Mode { Full, NoTraverse, OnlyInit };

template<class T>
struct ImplTraits
{
  using ContainerT = std::vector<Point>;
  using BufferT = const ContainerT&;
};
template<class Gather>
struct ImplTraits<ExtractMove<Gather>>
{
  using ContainerT = std::vector<Point>;
  using BufferT = ContainerT;
};
template<>
struct ImplTraits<ExtractList>
{
  using ContainerT = std::list<Point>;
  using BufferT = ContainerT;
};

template<class T, Mode mode = Full, Distribution dist = BothEnds>
void run(benchmark::State& state)
{
  const auto points = getBenchmarkArray<typename ImplTraits<T>::ContainerT>(dist, state.range(0));

  for (auto _ : state)
  {
    typename ImplTraits<T>::BufferT buffer = points;

    if constexpr (mode == OnlyInit)
    {
      benchmark::DoNotOptimize(buffer);
      continue;
    }

    auto range = T().extract(std::move(buffer));  // Most of implementations will not actually move from buffer

    if constexpr (mode == NoTraverse)
    {
      benchmark::DoNotOptimize(range.begin());
      continue;
    }

    volatile int sink;
    for (const auto& p : range) 
      sink = p.x;
    benchmark::DoNotOptimize(sink);
  }
}
BENCHMARK_TEMPLATE(run, ExtractMove<GatherNaive>)->Apply(setupBenchmark)->Arg(2);
BENCHMARK_TEMPLATE(run, ExtractMove<GatherSmart>)->Apply(setupBenchmark)->Arg(2);
BENCHMARK_TEMPLATE(run, ExtractList)->Apply(setupBenchmark)->Arg(2);
BENCHMARK_TEMPLATE(run, ExtractMove<GatherNaive>, OnlyInit)->Apply(setupBenchmark)->Arg(2);
BENCHMARK_TEMPLATE(run, ExtractMove<GatherSmart>, OnlyInit)->Apply(setupBenchmark)->Arg(2);
BENCHMARK_TEMPLATE(run, ExtractList, OnlyInit)->Apply(setupBenchmark)->Arg(2);
BENCHMARK_TEMPLATE(run, ExtractMove<GatherNaive>, NoTraverse)->Apply(setupBenchmark)->Arg(2);
BENCHMARK_TEMPLATE(run, ExtractMove<GatherSmart>, NoTraverse)->Apply(setupBenchmark)->Arg(2);
BENCHMARK_TEMPLATE(run, ExtractList, NoTraverse)->Apply(setupBenchmark)->Arg(2);

BENCHMARK_TEMPLATE(run, ExtractCopy)->Apply(setupBenchmark)->Arg(2);
BENCHMARK_TEMPLATE(run, ExtractView)->Apply(setupBenchmark)->Arg(2);
BENCHMARK_TEMPLATE(run, ExtractViewWrappingIterator)->Apply(setupBenchmark)->Arg(2);
BENCHMARK_TEMPLATE(run, ExtractViewGeneric)->Apply(setupBenchmark)->Arg(2);
BENCHMARK_TEMPLATE(run, ExtractViewRanges)->Apply(setupBenchmark)->Arg(2);
BENCHMARK_TEMPLATE(run, ExtractViewCoroutine)->Apply(setupBenchmark)->Arg(2);
BENCHMARK_TEMPLATE(run, ExtractNoCheck)->Apply(setupBenchmark)->Arg(2);
BENCHMARK_TEMPLATE(run, ExtractNoCheckSimple)->Apply(setupBenchmark)->Arg(2);
BENCHMARK_TEMPLATE(run, ExtractCopy, NoTraverse)->Apply(setupBenchmark)->Arg(2);
BENCHMARK_TEMPLATE(run, ExtractView, NoTraverse)->Apply(setupBenchmark)->Arg(2);
BENCHMARK_TEMPLATE(run, ExtractViewWrappingIterator, NoTraverse)->Apply(setupBenchmark)->Arg(2);
BENCHMARK_TEMPLATE(run, ExtractViewGeneric, NoTraverse)->Apply(setupBenchmark)->Arg(2);
BENCHMARK_TEMPLATE(run, ExtractViewRanges, NoTraverse)->Apply(setupBenchmark)->Arg(2);
BENCHMARK_TEMPLATE(run, ExtractViewCoroutine, NoTraverse)->Apply(setupBenchmark)->Arg(2);
BENCHMARK_TEMPLATE(run, ExtractNoCheck, NoTraverse)->Apply(setupBenchmark)->Arg(2);
BENCHMARK_TEMPLATE(run, ExtractNoCheckSimple, NoTraverse)->Apply(setupBenchmark)->Arg(2);

BENCHMARK_TEMPLATE(run, ExtractViewGeneric, Full, Beginning)->Apply(setupBenchmark)->RangeMultiplier(2)->Range(1, BenchDataSize * 2);
BENCHMARK_TEMPLATE(run, ExtractNoCheck, Full, Beginning)->Apply(setupBenchmark)->RangeMultiplier(2)->Range(1, BenchDataSize * 2);
BENCHMARK_TEMPLATE(run, ExtractViewGeneric, NoTraverse, Beginning)->Apply(setupBenchmark)->RangeMultiplier(2)->Range(1, BenchDataSize * 2);
BENCHMARK_TEMPLATE(run, ExtractNoCheck, NoTraverse, Beginning)->Apply(setupBenchmark)->RangeMultiplier(2)->Range(1, BenchDataSize * 2);

int main(int argc, char* argv[])
{
  testing::InitGoogleTest(&argc, argv);
  if (auto status = RUN_ALL_TESTS(); status != 0)
    return status;

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();

  return EXIT_SUCCESS;
}
