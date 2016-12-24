#include <vector>
#include "Common.h"

std::vector<Point> segmentation(const std::vector<Point>& points)
{
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

  std::vector<Point> result;

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

  return result;
}
