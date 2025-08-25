// Lines

#include "libs/segment.h"
#include "libs/circle.h"

/*!
\class Line2 segment.h
\brief Lines in the plane.

\sa Segment2
\ingroup PlanarGroup
*/
const Line2 X(Vector2::Null, Vector2::X);
const Line2 Y(Vector2::Null, Vector2::Y);


/*!
\brief Creates a line in the plane.
\param a,b Two points.
*/
Line2::Line2(const Vector2& a, const Vector2& b) :Segment2(a, b)
{
}

/*!
\brief Creates a line from a segment.
\param s %Segment.
*/
Line2::Line2(const Segment2& s) :Segment2(s.Vertex(0), s.Vertex(1))
{
}


/*!
\brief Compute the intersection between two lines.

Note that the algorithm is the same as for segments, except that the
range of the intersection parameters is not checked.

\param l Other line.
\param p Intersection point.
*/
bool Line2::Intersection(const Line2& l, Vector2& p) const
{
  Vector2 u = b - a;
  Vector2 v = l.b - l.a;
  double d = u / v;

  // Segments are parallel
  if (fabs(d) < Segment2::epsilon) return false;

  Vector2 w = a - l.a;
  double s = v / w;
  //  double t = u/w;	

  p = a + u * (s / d);
  return true;
}

/*!
\brief Compute the intersection between a line and a segment.

Note that the algorithm is the same as for segments,
except that the range of the intersection parameters
is checked only for the segment.

\param e %Segment.
\param p Intersection point.
*/
bool Line2::Intersection(const Segment2& e, Vector2& p) const
{
  Vector2 u = b - a;
  Vector2 v = e.Vertex(1) - e.Vertex(0);
  double d = u / v;

  // Segments are parallel
  if (fabs(d) < Segment2::epsilon) return false;

  Vector2 w = a - e.Vertex(0);
  double s = (v / w) / d;

  // Test if intersection point lies outside range
  //if (s < 0.0 || s>1.0) return false;	
  double t = (u / w) / d;
  if (t < 0.0 || t>1.0) return false;

  p = a + u * s;
  return true;
}

/*!
\brief Compute the position of a point with respect to a line.
\param p Point
\param epsilon Precision.
*/
bool Line2::IsLeftOrOn(const Vector2& p, const double& epsilon) const
{
  return (b - a) / (p - a) >= 0.0 - epsilon;
}

/*!
\brief Compute the position of a point with respect to a line.
\param p Point
\param epsilon Precision.
*/
bool Line2::IsRightOrOn(const Vector2& p, const double& epsilon) const
{
  return (b - a) / (p - a) <= 0.0 + epsilon;
}

/*!
\brief Compute the squared distance to the line.
\param p Point.
*/
double Line2::R(const Vector2& p) const
{
  // Orthogonal
  Vector2 abo = (b - a).Orthogonal();

  Vector2 ap = p - a;

  return Math::Sqr(abo * ap) / (abo * abo);
}

/*!
\brief Compute the point symmetric to the line.
\param p Point.
*/
Vector2 Line2::Symmetry(const Vector2& p) const
{
  Vector2 u = Normalized(b - a);
  Vector2 v = u.Orthogonal();
  return p - 2.0 * ((p - a) - ((p - a) * v) * v);
}

/*!
\brief Compute the box bounding the box symmetric to the line.
\param box The box.
*/
Box2 Line2::Symmetry(const Box2& box) const
{
  Box2 s(Symmetry(box.Vertex(0)));
  for (int i = 1; i < 4; i++)
  {
    s.Extend(Symmetry(box.Vertex(i)));
  }
  return s;
}

/*!
\brief Compute the circle symmetric to the line.
\param c %Circle.
*/
Circle2 Line2::Symmetry(const Circle2& c) const
{
  return Circle2(Symmetry(c.Center()), c.Radius());
}
