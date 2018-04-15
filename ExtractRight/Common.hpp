#pragma once

#include <iostream>

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

inline bool isPositive(const Point& pt)
{
  return pt.x >= 0;
};

// Convenient for unit tests
inline std::ostream& operator<<(std::ostream& os, const Point& p)
{
  return os << '(' << p.x << ", " << p.y << ')';
}