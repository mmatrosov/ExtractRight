#pragma once

#include <boost/range/algorithm.hpp>
#include <boost/range/join.hpp>
#include <boost/algorithm/cxx11/is_partitioned.hpp>
#include <boost/algorithm/cxx11/partition_point.hpp>

#include <range/v3/core.hpp>
#include <range/v3/view/concat.hpp>

#include <experimental/generator>

#include <vector>
#include <cctype>

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

inline bool isRight(const Point& pt)
{
  return pt.x >= 0;
};
