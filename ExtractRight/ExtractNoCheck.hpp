#pragma once

#include "Common.hpp"

template<class It>
struct Bounds
{
  It begin1, end1, begin2, end2;
};

template<class It, class Predicate>
Bounds<It> findBounds(It first, It last, Predicate p, std::forward_iterator_tag)
{
  auto begin1 = std::find_if    (first,  last, p);
  auto end1   = std::find_if_not(begin1, last, p);
  auto begin2 = std::find_if    (end1,   last, p);
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
    {
      if (p(first[i]))
        return first + i;
    }
    step /= 2;
  }

  return last;
}

template<class It, class Predicate>
Bounds<It> findBounds(It first, It last, Predicate p, std::random_access_iterator_tag)
{
  return findBounds(first, last, p, std::forward_iterator_tag{});
}

template<class It, class Predicate>
Bounds<It> findBounds(It first, It last, Predicate p)
{
  return findBounds(first, last, p, std::iterator_traits<It>::iterator_category{});
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