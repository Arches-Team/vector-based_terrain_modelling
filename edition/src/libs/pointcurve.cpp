// Curves 

#include "libs/curvepoint.h"
#include "libs/segment.h"

/*!
\class PointCurve curvepoint.h
\brief Piecewise point curve.

\ingroup KernelGroup
*/

/*!
\brief Creates a curve with points.

\param p Set of points.
\param closed Flag defining whether the curve is closed (true) or open (false).
*/
PointCurve::PointCurve(const QVector<Vector>& p, bool closed) :VectorSet(p), closed(closed)
{
}

/*!
\brief Creates a point curve from a cubic curve.

Samples are obtained by a simple discretization of the interval

\param c %Cubic curve.
\param n Number of points.
\param a,b Interval.
*/
PointCurve::PointCurve(const CubicCurve& c, int n, const double& a, const double& b) : closed(false)
{
  v.reserve(n);
  for (int i = 0; i < n; i++)
  {
    v.append(c.Eval(Math::Lerp(a, b, Math::Unit(i, n))));
  }
}

/*!
\brief Creates a point curve from a quadric curve.

Samples are obtained by a simple discretization of the interval

\param c %Quadric curve.
\param n Number of points.
\param a,b Internval.
*/
PointCurve::PointCurve(const QuadricCurve& c, int n, const double& a, const double& b) : closed(false)
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
double PointCurve::Length() const
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
double PointCurve::Length(const QVector<Vector>& p, bool closed)
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
\param i Index.
 */
Vector PointCurve::Tangent(int i) const
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
      return v.at(i + 1) - v.at(i - 1);
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
      return v.at(i + 1) - v.at(i - 1);
    }
  }
}


/*!
\brief Computes the projection of p onto the piecewise linear curve and returns the squared distance.
\param p Point.
\param u Parameter defining the coordinate of the projection of the argument vertex onto the curve.
\param k Index of the segment for which the minimum distance was found.
*/
double PointCurve::R(const Vector& p, double& u, int& k) const
{
  if (v.size() == 0)
    return 0.0;

  k = 0;
  Vector a = v.at(0);
  double d = SquaredNorm(a - p);
  for (int i = 1; i < v.size(); i++)
  {
    const Vector b = v.at(i);
    const Segment ab(a, b);
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
\brief Rotate the curve.
\param r Rotation matrix.
*/
void PointCurve::Rotate(const Matrix& r)
{
  for (Vector& p : v)
  {
    p = r * p;
  }
}

/*!
\brief Rotate the curve.
\param t Translation vector.
*/
void PointCurve::Translate(const Vector& t)
{
  for (Vector& p : v)
  {
    p += t;
  }
}

/*!
\brief Scale the curve.
\param s Scaling.
*/
void PointCurve::Scale(const double& s)
{
  for (Vector& p : v)
  {
    p *= s;
  }
}

/*!
\brief Tranform the curve.
\param frame The frame.
*/
void PointCurve::Transform(const Frame& frame)
{
  for (Vector& p : v)
  {
    p = frame.Transform(p);
  }
}

/*!
\brief Tranform the curve.
\param frame The frame.
*/
PointCurve PointCurve::Transformed(const Frame& frame) const
{
  QVector<Vector> q = v;

  for (Vector& p : q)
  {
    p = frame.Transform(p);
  }
  return PointCurve(q, closed);
}