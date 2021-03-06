#pragma once

#include "Common.hpp"

#include <boost/range/join.hpp>

#include <algorithm>
#include <functional>
#include <vector>

template<class It>
struct Bounds
{
  It begin1, end1, begin2, end2;
};

template<class It, class Predicate>
Bounds<It> findBounds(It first, It last, Predicate p, std::forward_iterator_tag)
{
  It begin1 = std::find_if    (first,  last, p);
  It end1   = std::find_if_not(begin1, last, p);
  It begin2 = std::find_if    (end1,   last, p);
  return { begin1, end1, begin2, last };
}

template<class It, class Predicate>  // It is RandomAccess
It findAny(It first, It last, Predicate p)
{
  using diff_t = typename std::iterator_traits<It>::difference_type;

  diff_t n = last - first;

  diff_t step = 1;
  while (step <= n)
    step *= 2;

  while (step > 1)
  {
    for (diff_t i = step / 2 - 1; i < n; i += step)
      if (p(first[i]))
        return first + i;
    step /= 2;
  }

  return last;
}

template<class It, class Predicate>
Bounds<It> findBounds(It first, It last, Predicate p, std::random_access_iterator_tag)
{
  Bounds<It> bounds = { first, last, last, last };  // Whole sequence by default

  if (first == last)
    return bounds;

  if (!p(*first) || !p(last[-1]))
  {
    // One segment, or empty
    It sample = findAny(first, last, p);
    bounds.begin1 = std::partition_point(first, sample, std::not_fn(p));
    bounds.end1 = std::partition_point(sample, last, p);
  }
  else if (It hole = findAny(first, last, std::not_fn(p)); hole != last)
  {
    // Two segments
    bounds.end1 = std::partition_point(first, hole, p);
    bounds.begin2 = std::partition_point(hole, last, std::not_fn(p));
  }

  return bounds;
}

template<class It, class Predicate>
Bounds<It> findBounds(It first, It last, Predicate p)
{
  return findBounds(first, last, p, typename std::iterator_traits<It>::iterator_category{});
}

template<class It, class Predicate>
Bounds<It> findBoundsForSlides(It first, It last, Predicate p)
{
  if constexpr (std::is_same_v<typename std::iterator_traits<It>::iterator_category, 
                               std::random_access_iterator_tag>)
    return findBoundsRandomAccess(first, last, p);
  else
    return findBoundsDefault(first, last, p);
}

class ExtractNoCheck
{
public:
  template<class It, class Predicate>
  auto extract(It first, It last, Predicate p)
  {
    auto bounds = findBounds(first, last, p);

    return boost::join(boost::make_iterator_range(bounds.begin2, bounds.end2),
                       boost::make_iterator_range(bounds.begin1, bounds.end1));
  }

  auto extract(const std::vector<Point>& points)
  {
    return extract(points.begin(), points.end(), isPositive);
  }
};

class ExtractNoCheckSimple
{
public:
  template<class It, class Predicate>
  auto extract(It first, It last, Predicate p)
  {
    It begin1 = std::find_if    (first,  last, p);
    It end1   = std::find_if_not(begin1, last, p);
    It begin2 = std::find_if    (end1,   last, p);
    It end2   = last;

    assert(end2 == std::find_if_not(begin2, last, p) &&
      (begin2 == end2 || begin1 == first && end2 == last));

    return boost::join(boost::make_iterator_range(begin2, end2),
                       boost::make_iterator_range(begin1, end1));
  }

  auto extract(const std::vector<Point>& points)
  {
    return extract(points.begin(), points.end(), isPositive);
  }
};
