// Segments

#include <QtWidgets/QGraphicsScene>

#include "libs/segment.h"
#include "libs/triangle.h"

/*!
\class Segment2 segment.h
\brief Segments in the plane.

A segment stores its end vertices.

\ingroup PlanarGroup
*/

/*!
\brief Overloaded.
\param s Stream.
\param segment The segment.
*/
std::ostream& operator<<(std::ostream& s, const Segment2& segment)
{
  s << "Segment2(" << segment.a << ',' << segment.b << ")";
  return s;
}

/*!
\brief Draw a segment.
\param scene Graphics scene.
\param pen The pen.
*/
void Segment2::Draw(QGraphicsScene& scene, const QPen& pen) const
{
  scene.addLine(a[0], a[1], b[0], b[1], pen);
}

/*!
\brief Draw a segment as an arrow.
\param scene Graphics scene.
\param pen The pen.
\param brush The brush.
\param s Size of the arrow head.
*/
void Segment2::DrawArrow(QGraphicsScene& scene, const double& s, const QPen& pen, const QBrush& brush) const
{
  Vector2 u = Normalized(b - a);
  Vector2 v = u.Orthogonal();

  Vector2 c = b - u * s;

  // Arrow segment
  scene.addLine(a[0], a[1], c[0], c[1], pen);

  // Arrow head
  Triangle2(c - 0.5 * v * s, b, c + 0.5 * v * s).Draw(scene, pen, brush);
}

/*!
\brief Compute the squared distance to the segment.
\param p Point.
*/
double Segment2::R(const Vector2& p) const
{
  Vector2 n = p - b;
  // Do not normalize there
  Vector2 axis = b - a;
  double s = n * axis;

  // Squared distance 
  double r = n * n;

  // Not vertex b
  if (s < 0.0)
  {
    n = p - a;
    s = n * axis;
    r = n * n;
    // Not vertex a
    if (s > 0.0)
    {
      s = (s * s) / (axis * axis);
      r -= s;
      // Added an absolute value as imprecisions may cause rounding errors
      r = fabs(r);
    }
  }
  return r;
}

/*!
\brief Compute the normal vector between a point and
its projection onto the segment.
\param p Point.
*/
Vector2 Segment2::Normal(const Vector2& p) const
{
  Vector2 n = p - b;

  // Do not normalize there
  Vector2 axis = b - a;

  // Not that vertex
  if (axis * n <= 0.0)
  {
    n = p - a;
    double s = axis * n;
    // Not that one either
    if (s >= 0.0)
    {
      // Normalize only there : compute length and scale scalar product
      double length = axis * axis;

      // Dividing s by length is more efficient than normalizing the axis
      s /= length;

      n -= s * axis;
    }
  }
  return n;
}

/*!
\brief Compute the squared distance to the axis edge characterized by its end points.
\param p Point.
\param s Linear coefficient corresponding to the projection of argument point p onto the edge.
*/
double Segment2::R(const Vector2& p, double& s) const
{
  Vector2 n = p - b;
  // Do not normalize there
  Vector2 axis = b - a;
  s = n * axis;

  // Squared distance 
  double r = n * n;

  // Not vertex b
  if (s < 0.0)
  {
    n = p - a;
    s = n * axis;
    r = n * n;
    // Not vertex a
    if (s > 0.0)
    {
      double length = axis * axis;
      s = s * s / length;
      r -= s; // Pythagore
      s = sqrt(s / length);
      // Added an absolute value as imprecisions may cause rounding errors
      r = fabs(r);
    }
    else
    {
      s = 0.0;
    }
  }
  else
  {
    s = 1.0;
  }
  return r;
}

/*!
\brief Compute the intersection between two segments.
\param e Other segment.
\param p Intersection point.
*/
bool Segment2::Intersection(const Segment2& e, Vector2& p) const
{
  Vector2 u = b - a;
  Vector2 v = e.b - e.a;
  double d = u / v;

  // Segments are parallel
  if (fabs(d) < Segment2::epsilon) return false;

  Vector2 w = a - e.a;
  double s = (v / w) / d;

  // Test if intersection point lies outside range
  if (s < 0.0 || s>1.0) return false;
  double t = (u / w) / d;
  if (t < 0.0 || t>1.0) return false;

  p = a + u * s;
  return true;
}

/*!
\brief Test if two segments intersect.

This function computes the intersection and test the intersection parameters.
Two segments that share an end point are detected as intersecting.

\param segment Other segment.
*/
bool Segment2::Intersect(const Segment2& segment) const
{
  Vector2 u = b - a;
  Vector2 v = segment.b - segment.a;
  double d = u / v;

  // Segments are almost parallel
  if (fabs(d) < Segment2::epsilon) return false;

  Vector2 w = a - segment.a;
  double s = (v / w) / d;

  // Test if intersection point lies outside range
  if ((s < 0.0) || (s > 1.0)) return false;
  double t = (u / w) / d;
  if ((t < 0.0) || (t > 1.0)) return false;

  return true;
}

/*!
\brief This functions tests if two segments intersect.

The algorithm avoids the computation of the intersection.
Instead, it checks whether the two points of one segment
lie on the same half space defined by the line of the other segment.

Contrary to bool Intersect(const Segment2&),
two segments that share one point are not detected as intersecting.

\param segment Other segment.
*/
bool Segment2::IntersectOpen(const Segment2& segment) const
{
  Vector2 x = b - a;
  Vector2 y = segment.b - segment.a;

  double sbax = x / (segment.a - a);
  double sbbx = x / (segment.b - b);
  if ((sbax >= 0.0) && (sbbx >= 0.0)) return false;
  if ((sbax <= 0.0) && (sbbx <= 0.0)) return false;

  double saay = y / (a - segment.a);
  double saby = y / (b - segment.a);
  if ((saay >= 0.0) && (saby >= 0.0)) return false;
  if ((saay <= 0.0) && (saby <= 0.0)) return false;

  return true;
}

/*!
\brief Compute the bounding box of the segment.
*/
Box2 Segment2::GetBox() const
{
  return Box2(Vector2::Min(a, b), Vector2::Max(a, b));
}

/*!
\brief Compute an orthogonal vector to the segment.

The returned vector is not normalized.
*/
Vector2 Segment2::Orthogonal() const
{
  return (b - a).Orthogonal();
}

/*!
\copydoc Segment::Intersect(const Vector& , const Vector& , const double , const double , const double )
*/
Vector2 Segment2::Intersect(const Vector2& a, const Vector2& b, const double va, const double vb, const double y)
{
  double t = (y - va) / (vb - va);
  return a + t * (b - a);
}


/*!
\brief Translates a segment.

\param t Translation vector.
*/
void Segment2::Translate(const Vector2& t)
{
  a += t;
  b += t;
}

/*!
\brief Translates a segment.

\param t Translation vector.
*/
Segment2 Segment2::Translated(const Vector2& t) const
{
  return Segment2(a + t, b + t);
}

/*!
\brief Translates a segment.

\param r Rotation matrix.
*/
void Segment2::Rotate(const Matrix2& r)
{
  a = r * a;
  b = r * b;
}

/*!
\brief Uniformly scales a segment.

\param s Scaling factor.
*/
void Segment2::Scale(const double& s)
{
  a *= s;
  b *= s;
}

/*!
\brief Extend the length of a segment by a given distance.

A negative value shrinks the segment.

\sa Box2::Extended
\param r Distance.
*/
Segment2 Segment2::Extended(const double& r) const
{
  Vector2 ab = b - a;
  double length = Norm(ab);
  if (length == 0.0)
  {
    return Segment2(a, b);
  }
  Vector2 u = (r / length) * ab;
  return Segment2(a - u, b + u);
}

