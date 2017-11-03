#include <iostream>
#include <vector>

struct Point
{
  int x;
  int y;
};

std::istream& operator>>(std::istream& stream, Point& p)
{
  return stream >> p.x >> p.y;
}
std::ostream& operator<<(std::ostream& stream, const Point& p)
{
  return stream << p.x << ' ' << p.y;
}

void extractRight(std::vector<Point>& points)
{
  // Write everything you need inside this function, or in global namespace
}

// Do not modify main function
int main()
{
  std::vector<Point> points;

  Point p;
  while (std::cin >> p)
    points.push_back(p);

  try
  {
    extractRight(points);
  }
  catch (const std::runtime_error& e)
  {
    std::cout << e.what();
    return 0;
  }

  for (const auto& p : points)
    std::cout << p << '\n';
}
