// Disc

#include "libs/circle.h"

/*!
\class Disc2 circle.h
\brief A disc in two dimensions.

Discs in the plane derive from Circle2, and only reimplement a few functions, such as Disc2::R() and Disc2::Normal().

\ingroup PlanarGroup
*/

/*!
\brief Creates a disc.
\param r Radius.
*/
Disc2::Disc2(const double& r) :Circle2(r)
{
}

/*!
\brief Creates a disc.
\param c Center.
\param r Radius.
*/
Disc2::Disc2(const Vector2& c, const double& r) :Circle2(c, r)
{
}

/*!
\brief Creates a disc.
\param circle %Circle.
*/
Disc2::Disc2(const Circle2& circle) :Circle2(circle)
{
}

/*!
\brief Creates a disc given three vertices.

Basically set the center to the center of the triangle,
and compute the radius.

\param x, y, z Argument vertices
*/
Disc2::Disc2(const Vector2& x, const Vector2& y, const Vector2& z) :Circle2(x, y, z)
{
}

/*!
\brief Computes the distance between a point and the disc.

\param p Point.
*/
double Disc2::R(const Vector2& p) const
{
  Vector2 n = p - c;

  // Distance to center
  double y = n * n;

  // Inside
  if (y < r * r)
  {
    return 0.0;
  }

  // Distance to circle
  y = sqrt(y) - r;
  return y;
}

/*!
\brief Computes the signed distance between a point and the disc.

\param p Point.
*/
double Disc2::Signed(const Vector2& p) const
{
  return Norm(p - c) - r;
}

/*!
\brief Computes the distance vector between a disc
and a point.

\param p Point.
*/
Vector2 Disc2::Normal(const Vector2& p) const
{
  Vector2 n = p - c;

  double y = n * n;

  // Inside
  if (y < r * r)
  {
    return Vector2::Null;
  }
  // Outside
  else
  {
    return (1.0 - r / sqrt(y)) * n;
  }
}

/*!
\brief Overloaded.
\param s Stream.
\param disc The disc.
*/
std::ostream& operator<<(std::ostream& s, const Disc2& disc)
{
  s << "Disc2(" << disc.c << ',' << disc.r << ')';
  return s;
}