// Curves 

#include "libs/curvepoint.h"
#include "libs/segment.h"
#include "libs/curveset.h"

/*!
\class PointCurve2 curvepoint.h
\brief Piecewise point curve in the plane.

\ingroup KernelGroup
*/

/*!
\brief Creates a curve with points.

\param p Set of points.
\param closed Flag defining whether the curve is closed (true) or open (false).
*/
PointCurve2::PointCurve2(const VectorSet2& p, bool closed) :VectorSet2(p), closed(closed)
{
}

/*!
\brief Creates a curve with points.

\param p Set of points.
\param closed Flag defining whether the curve is closed (true) or open (false).
*/
PointCurve2::PointCurve2(const QVector<Vector2>& p, bool closed) :VectorSet2(p), closed(closed)
{
}

/*!
\brief Create a curve from a set of vertexes and a subset of indexes.
\param p Array of points.
\param indexes Set of indexes.
\param closed Flag defining whether the curve is closed (true) or open (false).
*/
PointCurve2::PointCurve2(const QVector<Vector2>& p, const QVector<int>& indexes, bool closed) :closed(closed)
{
  int n = indexes.size();
  v.resize(n);

  for (int i = 0; i < n; i++)
  {
    v[i] = p.at(indexes.at(i));
  }
}

/*!
\brief Creates a point curve from a cubic curve.

\param c %Cubic curve.
\param n Number of points.
\param a,b Internval.
*/
PointCurve2::PointCurve2(const CubicCurve2& c, int n, const double& a, const double& b) : closed(false)
{
  v.reserve(n);
  for (int i = 0; i < n; i++)
  {
    v.append(c.Eval(Math::Lerp(a, b, Math::Unit(i, n))));
  }
}

/*!
\brief Creates a point curve from a quadric curve.

\param c %Quadric curve.
\param n Number of points.
\param a,b Internval.
*/
PointCurve2::PointCurve2(const QuadricCurve2& c, int n, const double& a, const double& b) : closed(false)
{
  v.reserve(n);
  for (int i = 0; i < n; i++)
  {
    v.append(c.Eval(Math::Lerp(a, b, Math::Unit(i, n))));
  }
}

/*!
\brief Compute the length of the curve.
*/
double PointCurve2::Length() const
{
  double length = 0.0;
  for (int i = 0; i < v.size() - 1; i++)
  {
    length += Norm(v.at(i + 1) - v.at(i));
  }

  // Add last segment if closed
  if (closed)
  {
    length += Norm(v.at(v.size() - 1) - v.at(0));
  }
  return length;
}

/*!
\brief Compute the length of a discrete curve.
\param p %Array of points.
\param closed Boolean.
*/
double PointCurve2::Length(const QVector<Vector2>& p, bool closed)
{
  double length = 0.0;
  for (int k = 0; k < p.size() - 1; k++)
  {
    length += Norm(p.at(k) - p.at(k + 1));
  }
  // Add last segment if closed
  if (closed)
  {
    length += Norm(p.at(p.size() - 1) - p.at(0));
  }
  return length;
}

/*!
\brief Compute the tangent at a given point.

\image html tangent.png

\param i Index.
 */
Vector2 PointCurve2::Tangent(int i) const
{
  const int k = v.size() - 1;

  if (!closed)
  {
    if (i == 0)
    {
      return v.at(1) - v.at(0);
    }
    else if (i == k)
    {
      return v.at(k) - v.at(k - 1);
    }
    else
    {
      return (v.at(i + 1) - v.at(i - 1)) / 2.0;
    }
  }
  else
  {
    if (i == 0)
    {
      return v.at(1) - v.at(k);
    }
    else if (i == k)
    {
      return v.at(0) - v.at(k - 1);
    }
    else
    {
      return (v.at(i + 1) - v.at(i - 1)) / 2.0;
    }
  }
}

/*!
\brief Transform the planar point curve into a three dimensional point curve.

Transformation if performed as follows: (x,y) are mapped to (0,x,y).
\param f %Frame.
*/
PointCurve PointCurve2::Transform(const Frame& f) const
{
  QVector<Vector> t(v.size());
  for (int i = 0; i < v.size(); i++)
  {

    t[i] = f.Transform(Vector(0.0, v[i][0], v[i][1]));
  }

  return PointCurve(t, closed);
}

/*!
\brief Computes the projection of p onto the piecewise linear curve and returns the squared distance.
\param p Point.
\param u Parameter defining the coordinate of the projection of the argument vertex onto the curve.
\param k Index of the curve for which the minimum distance was found.
*/
double PointCurve2::R(const Vector2& p, double& u, int& k) const
{
  if (v.size() == 0)
    return 0.0;

  k = 0;
  Vector2 a = v.at(0);
  double d = SquaredNorm(a - p);
  for (int i = 1; i < v.size(); i++)
  {
    const Vector2 b = v.at(i);
    const Segment2 ab(a, b);
    double s;
    double di = ab.R(p, s);
    if (di < d)
    {
      d = di;
      u = s;
      k = i - 1;
    }
    a = b;
  }
  return d;
}

/*!
\brief Computes the projection of p onto the piecewise linear curve, returns the squared distance, and the side with respect to the curve
\param p Point.
\param u Parameter defining the coordinate of the projection of the argument vertex onto the curve.
\param k Index of the curve for which the minimum distance was found.
\param side Side of the point with respect to the closest segment (1 or -1)
*/
double PointCurve2::R(const Vector2& p, double& u, int& k, int& side) const
{
  if (v.size() == 0)
    return 0.0;

  k = 0;
  Vector2 a = v.at(0);
  double d = SquaredNorm(a - p);
  for (int i = 1; i < v.size(); i++)
  {
    const Vector2 b = v.at(i);
    const Segment2 ab(a, b);
    double s;
    double di = ab.R(p, s);
    if (di < d)
    {
      d = di;
      u = s;
      k = i - 1;
      if ((p - a) / (b - a) > 0.0)
        side = 1;
      else
        side = -1;
    }
    a = b;
  }
  return d;
}

/*!
\brief Compute the sinuosity.

Sinusoity is defined as the ratio between the length of the curve and the distance between end points.
*/
double PointCurve2::Sinuosity() const
{
  return Length() / Norm(v.first() - v.last());
}

/*!
\brief Compute the closed bounding polygon.
*/
Polygon2 PointCurve2::GetPolygon() const
{
  return Polygon2(v);
}

/*!
\brief Compute the cubic curve passing through all the points of the polyline.
*/
CubicCurve2Set PointCurve2::ToCubicCurve() const
{
  return CubicCurve2Set(v);
}

/*!
\brief Compute the quadric curve passing through all the points of the polyline.
*/
QuadricCurve2Set PointCurve2::ToQuadricCurve() const
{
  return QuadricCurve2Set(v);
}

/*!
\brief Compute the curvature at a given point.
\param i Index of the point on the curve.
*/
double PointCurve2::Curvature(int i) const
{
  if (i == 0 || i == v.size() - 1)
    return 0.0;

  // First derivative
  Vector2 dxy = Tangent(i);

  // Second derivative
  Vector2 dxy0 = Tangent(i - 1);
  Vector2 dxy1 = Tangent(i + 1);
  Vector2 ddxy = (dxy1 - dxy0) / 2.0;

  // Compute curvature
  double dx = dxy[0];
  double dy = dxy[1];
  double ddx = ddxy[0], ddy = ddxy[1];
  return (dx * ddy - dy * ddx) / Math::Pow(Math::Sqr(dx) + Math::Sqr(dy), 1.5);
}

/*!
\brief Draw the point curve.
\param scene Graphics scene.
\param pen The pen.
*/
void PointCurve2::Draw(QGraphicsScene& scene, const QPen& pen) const
{
  for (int i = 0; i < v.size() - 1; i++)
  {
    Segment2(v.at(i), v.at(i + 1)).Draw(scene, pen);
  }
}

/*!
\brief Return a translated point curve.
\param t Translation vector.
*/
PointCurve2 PointCurve2::Translated(const Vector2& t) const
{
  return PointCurve2(VectorSet2::Translated(t), closed);
}

/*!
\brief Return a transformed point curve.
\param frame The frame.
*/
PointCurve2 PointCurve2::Transformed(const Frame2& frame) const
{
  return PointCurve2(VectorSet2::Transformed(frame), closed);
}

/*!
\brief Return a rotated point curve.
\param a Rotation angle.
*/
PointCurve2 PointCurve2::Rotated(const double& a) const
{
  return PointCurve2(VectorSet2::Rotated(a), closed);

}

/*!
\brief Return a scaled point curve.
\param s Scaling vector.
*/
PointCurve2 PointCurve2::Scaled(const Vector2& s) const
{
  return PointCurve2(VectorSet2::Scaled(s), closed);
}