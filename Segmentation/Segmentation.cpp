#include "Common.h"

#include <boost/range/algorithm.hpp>
#include <boost/algorithm/cxx11/is_partitioned.hpp>
#include <boost/algorithm/cxx11/find_if_not.hpp>

#include <vector>

std::vector<Point> segmentationNaive(const std::vector<Point>& points)
{
  std::vector<Point> result;
  result.clear();

  if (points.size() == 0)
    return result;

  int p = 0;
  for (int i = 1; i < points.size(); ++i)
    if (points[i - 1].x < 0 && points[i].x >= 0)
    {
      p = i;
      break;
    }

  int q = 0;
  for (int i = 1; i < points.size(); ++i)
    if (points[i - 1].x >= 0 && points[i].x < 0)
    {
      q = i;
      break;
    }

  if (p == q)
  {
    if (points[0].x >= 0)
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

std::vector<Point> segmentationNaiveRefactored(const std::vector<Point>& points)
{
  auto isRight = [](const Point& pt) { return pt.x >= 0; };

  int p = 0;
  int q = 0;
  for (int i = 1; i < points.size(); ++i)
  {
    if (!isRight(points[i - 1]) && isRight(points[i]))
      p = i;
    if (isRight(points[i - 1]) && !isRight(points[i]))
      q = i;
  }

  std::vector<Point> result;

  if (p == q)
  {
    if (!points.empty() && isRight(points[0]))
      result = points;
    return result;
  }

  auto iterate = [&](int from, int to, bool shouldBeRight)
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

  iterate(p, q, true);
  iterate(q, p, false);

  return result;
}


std::vector<Point> segmentation(std::vector<Point> points)
{
  using namespace boost::range;
  using namespace boost::algorithm;

  auto isRight = [](const Point& pt) { return pt.x >= 0; };

  auto middle = adjacent_find(points,
    [&](auto&& pt1, auto&& pt2) { return !isRight(pt1) && isRight(pt2); });

  if (middle != points.end())
    rotate(points, std::next(middle));

  if (!is_partitioned(points, isRight))
    throw std::runtime_error("Unexpected order");

  points.erase(find_if_not(points, isRight), points.end());

  return points;
}