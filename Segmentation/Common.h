#pragma once

struct Point
{
  double x;
  double y;
};

inline bool operator==(const Point& a, const Point& b)
{
  return a.x == b.x && a.y == b.y;
}