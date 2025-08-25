// Cubic vurves 

#include "libs/curve.h"
#include "libs/quartic.h"
#include "libs/polynomial.h"

/*!
\class CubicCurve curve.h
\brief %Cubic curves.

%Cubic curves are implemented with three cubic polynomials defining the
coordinates of the points on the curve. It is a minimum storage class in
the sense that it does not keep extra constant data which may be useful
for computing the bounding box of the quadric curve parameterized over
the unit interval [0,1] or the distance between a point an the curve.

\image html cubic.png

\ingroup KernelGroup

*/

/*!
\brief Creates a cubic curve.

\param x, y, z Cubic equations of the coordinates.
*/
CubicCurve::CubicCurve(const Cubic& x, const Cubic& y, const Cubic& z) :x(x), y(y), z(z)
{
}

/*!
\brief Translate a curve.

\param t %Translation vector.
*/
CubicCurve CubicCurve::Translated(const Vector& t) const
{
  Cubic tx = x;
  tx[0] += t[0];
  Cubic ty = y;
  ty[0] += t[1];
  Cubic tz = z;
  tz[0] += t[2];
  return CubicCurve(tx, ty, tz);
}

/*!
\brief Computes the quadric curve whose components are defined as the derivative of the cubic coordinates of the curve.

This function is useful for preprocessing tangent computations.
*/
QuadricCurve CubicCurve::Tangent() const
{
  return QuadricCurve(x.Prime(), y.Prime(), z.Prime());
}

/*!
\brief Computes the tangent to the cubic curve at a given position on the curve.

The tangent is the obtained by evaluating the first derivative.
\sa CubicCurve::Tangent()
\param t Position.
*/
Vector CubicCurve::Tangent(const double& t) const
{
  return Vector(x.Prime()(t), y.Prime()(t), z.Prime()(t));
}

/*!
\brief Computes the normal vector to the cubic curve at a given position on the curve.

The normal is the obtained by evaluating the second derivative.
\param t Position.
*/
Vector CubicCurve::Normal(const double& t) const
{
  return Vector(x.Second()(t), y.Second()(t), z.Second()(t));
}

/*!
\brief Computes the length of the cubic curve over interval [a,b].

Since we have the parametric equations of the coordinates of the points on the curve,
then the length is the integral of the square root of the sum of the squared derivatives
where the limits of integration are [a,b].

s = Integral Sqrt((dx/dt)<SUP>2</SUP> + (dy/dt)<SUP>2</SUP>) dt,

In the general case, the function we end up with is not integrable in closed
form, leading to an elliptic integral or some such. Thus we evaluate the integral
numerically.

\param a,b Interval.
\param n Discretization of the integration interval.
*/
double CubicCurve::Length(const double& a, const double& b, int n) const
{
  double s = 0.0;

  // Compute tangent quadric curve
  QuadricCurve tangent = Tangent();

  // Evaluate quartic defined as the squared norm of the tangent coordinates
  Quartic quartic = (tangent[0]) * (tangent[0]) + (tangent[1]) * (tangent[1]) + (tangent[2]) * (tangent[2]);

  double epsilon = (b - a) / n;

  double u = sqrt(fabs(quartic(a)));
  for (double t = a + epsilon; t <= b + epsilon / 2.0; t += epsilon)
  {
    double v = sqrt(fabs(quartic(t)));
    s += (u + v) / 2.0;
    u = v;
  }
  s *= epsilon;

  return s;
}

/*!
\brief Compute the parameter t corresponding to a given input arc length.
\param s Arc length, which should be positive.
\param n Discretization of the curve.
*/
double CubicCurve::InverseArcLength(const double& s, int n) const
{
  if (s <= 0.0)
    return 0.0;

  double length = 0.0;

  // Compute tangent quadric curve
  QuadricCurve tangent = Tangent();

  // Evaluate quartic defined as the squared norm of the tangent coordinates
  Quartic quartic = tangent[0] * tangent[0] + tangent[1] * tangent[1] + tangent[2] * tangent[2];

  double epsilon = 1.0 / n;

  double u = sqrt(quartic(0.0));
  double t = epsilon;
  while (length < s)
  {
    double v = sqrt(quartic(t));
    length += epsilon * (u + v) / 2.0;
    u = v;
    t += epsilon;
  }
  return t - epsilon / 2.0;
}

/*!
\brief Computes the length of the cubic curve. Integration is performed over interval [0,1].

This function simply calls CubicCurve::Length(0.0,1.0,n) with n as the passed argument parameter.

\param n Discretization of the curve.
*/
double CubicCurve::Length(int n) const
{
  return Length(0.0, 1.0, n);
}

/*!
\brief Negative curve unary operator.
*/
CubicCurve CubicCurve::operator- () const
{
  return CubicCurve(-x, -y, -z);
}

/*!
\brief Evaluates curve point at a given location.
\param t Parameter.
*/
Vector CubicCurve::operator() (const double& t) const
{
  return Vector(x(t), y(t), z(t));
}

/*!
\brief Evaluates curve point at a given location.
\param t Parameter.
*/
Vector CubicCurve::Eval(const double& t) const
{
  return Vector(x(t), y(t), z(t));
}
/*!
\brief Creates an Hermite cubic curve on interval [0,1] given vertex locations and tangent vectors (in that order).
\image html hermite.png

\param a,b End vertices of the Hermite cubic curve.
\param ta,tb Tangent vectors at the end vertices.
*/
CubicCurve CubicCurve::Hermite(const Vector& a, const Vector& b, const Vector& ta, const Vector& tb)
{
  return CubicCurve(Cubic::Hermite(a[0], b[0], ta[0], tb[0]), Cubic::Hermite(a[1], b[1], ta[1], tb[1]), Cubic::Hermite(a[2], b[2], ta[2], tb[2]));
}

/*!
\brief Creates a Bezier cubic curve on interval [0,1] given four vertex locations.
\image html bezier.png

\param a,b,c,d Vertexes of the Bezier curve.
*/
CubicCurve CubicCurve::Bezier(const Vector& a, const Vector& b, const Vector& c, const Vector& d)
{
  return CubicCurve(Cubic::Bezier(a[0], b[0], c[0], d[0]), Cubic::Bezier(a[1], b[1], c[1], d[1]), Cubic::Bezier(a[2], b[2], c[2], d[2]));
}

/*!
\brief Overloaded.
\param s Stream.
\param curve The cubic curve.
*/
std::ostream& operator<<(std::ostream& s, const CubicCurve& curve)
{
  s << "CubicCurve(" << curve.x << ',' << curve.y << ',' << curve.z << ')';
  return s;
}

/*!
\brief Compute the bounding box of the cubic curve with parameter interval [0,1].

This function is computationnally intensive as the bounds are found by
evaluating the range of every coordinate within [0,1].
*/
Box CubicCurve::GetBox() const
{
  Box box;

  // Define bounding box
  x.Range(box[0][0], box[1][0]);
  y.Range(box[0][1], box[1][1]);
  z.Range(box[0][2], box[1][2]);

  return box;
}

/*!
\brief Returns the squared distance between
a point and the cubic curve.

This function computes the projection of a point p onto
the curve by minimizing the distance function along the curve,
expressed as a degree 5 polynomial.

Most of the coefficients of this polynomial are constant,
and could be pre-processed in the constructor.

Note that this function is computationnally intensive, and
~25 times slower than QuadricCurve::R(). Therefore, cubic
curves should be used with care as distance queries are
called many times.

\sa QuadricCurve::R()

\param p Point.
\param u Parameter defining the coordinate of the projection of the argument vertex onto the curve.
\return The distance.
*/
double CubicCurve::R(const Vector& p, double& u) const
{
  double xyz[5];

  // Quartic constant coefficients
  xyz[4] = 3.0 * (x[3] * x[3] + y[3] * y[3] + z[3] * z[3]);
  xyz[3] = 5.0 * (x[3] * x[2] + y[3] * y[2] + z[3] * z[2]);
  xyz[2] = 4.0 * (x[3] * x[1] + y[3] * y[1] + z[3] * z[1]) + 2.0 * (x[2] * x[2] + y[2] * y[2] + z[2] * z[2]);
  xyz[1] = 3.0 * (x[2] * x[1] + y[2] * y[1] + z[2] * z[1]);
  xyz[0] = x[1] * x[1] + y[1] * y[1] + z[1] * z[1];

  double t[5];

  Polynomial polynomial(
    xyz[4], xyz[3],
    xyz[2],
    xyz[1] + 3.0 * (x[3] * (x[0] - p[0]) + y[3] * (y[0] - p[1]) + z[3] * (z[0] - p[2])),
    xyz[0] + 2.0 * (x[2] * (x[0] - p[0]) + y[2] * (y[0] - p[1]) + z[2] * (z[0] - p[2])),
    x[1] * (x[0] - p[0]) + y[1] * (y[0] - p[1]) + z[1] * (z[0] - p[2]));

  int k = polynomial.Solve(t);

  // Search at caps
  Vector e = p - (*this)(0.0);
  double r = e * e;

  Vector f = p - (*this)(1.0);
  double s = f * f;

  u = 0.0;
  if (s < r)
  {
    r = s;
    u = 1.0;
  }

  for (int i = 0; i < k; i++)
  {
    if ((t[i] > 0.0) && (t[i] < 1.0))
    {
      f = p - Vector(x(t[i]), y(t[i]), z(t[i]));
      s = f * f;
      if (s < r)
      {
        u = t[i];
        r = s;
      }
    }
  }
  return r;
}

/*!
\brief Compute the parameter corresponding to the argument curvilign absisca of a point on the curve.

\param s Argument curvilign absisca (should be between 0.0 and the length of the curve).
\param n Discretization of the integration interval.
*/
double CubicCurve::U(const double& s, int n) const
{
  if (s == 0.0) return 0.0;

  const double epsilon = 1.0 / n;
  double ss = 0.0;
  double old_ss = 0.0;

  Vector a = Vector(x(0.0), y(0.0), z(0.0));

  for (double u = epsilon; u < 1.0; u += epsilon)
  {
    Vector b = Vector(x(u), y(u), z(u));
    ss += Norm(b - a);
    if (ss > s)
    {
      // Local linear interpolation, should not degrade the performances too much, increases precision
      return u - epsilon + epsilon * ((s - old_ss) / (ss - old_ss));
    }
    a = b;
    old_ss = ss;
  }
  return 1.0;
}

/*!
\brief Compute the curvilign absisca of a point on the curve.

\param t Parameter that should be within [0,1] interval.
\param n Discretization of the integration interval.
*/
double CubicCurve::S(const double& t, int n) const
{
  const double epsilon = t / n;
  double s = 0.0;
  Vector a = Vector(x(0.0), y(0.0), z(0.0));

  double u = epsilon;

  for (int i = 1; i <= n; i++)
  {
    Vector b = Vector(x(u), y(u), z(u));
    s += Norm(b - a);
    a = b;
    u += epsilon;
  }
  return s;
}

/*!
\brief Compute the Frenet vectors at a given point on the curve.

The columns of the returned matrix contain the tangent, normal and binormal vectors.

\param u Parameter.
\param vertical Boolean, set to true to compute frame along the horizontal curve first, false (default) otherwise.
*/
Matrix CubicCurve::GetMatrix(const double& u, bool vertical) const
{
  if (vertical == false)
  {
    Vector t = Tangent(u);
    Vector n = Normal(u);
    Vector b = t / n;
    n = b / t;
    Normalize(t);
    Normalize(b);
    Normalize(n);
    return Matrix(t, n, b);
  }
  else
  {
    Vector t = Tangent(u);
    Vector n = Vector::Z / t;
    Vector b = t / n;
    Normalize(t);
    Normalize(b);
    Normalize(n);
    return Matrix(t, n, b);
  }
}

/*!
\brief Compute the Frenet frame at a given point on the curve.

\param u Parameter.
\param vertical Boolean, set to true to compute frame along the horizontal curve first, false (default) otherwise.
*/
Frame CubicCurve::GetFrame(const double& u, bool vertical) const
{
  return Frame(GetMatrix(u, vertical), Eval(u));
}

/*!
\brief Computes the curvature of the curve at a given positon.
\param t Parameter.
*/
double CubicCurve::Curvature(const double& t) const
{
  // First and second derivatives
  double xp = x.Prime()(t);
  double yp = y.Prime()(t);
  double zp = z.Prime()(t);
  double xpp = x.Second()(t);
  double ypp = y.Second()(t);
  double zpp = z.Second()(t);

  return sqrt((zpp * yp - ypp * zp) * (zpp * yp - ypp * zp) + (xpp * zp - xp * zpp) * (xpp * zp - xp * zpp) + (ypp * xp - yp * xpp) * (ypp * xp - yp * xpp)) / Math::Sqrt32(xp * xp + yp * yp + zp * zp);
}

/*!
\brief Compute the i-th BÃ©zier control point of the curve.

Computations have been optimized.
\param i Index.
*/
inline Vector CubicCurve::BezierControl(int i) const
{
  if (i == 0)
  {
    // Equivalent to c(0.0)
    return Vector(x[0], y[0], z[0]);
  }
  else if (i == 1)
  {
    return Vector(x[0] + x[1] / 3.0, y[0] + y[1] / 3.0, z[0] + z[1] / 3.0);
  }
  else if (i == 2)
  {
    return Vector(x[0] + (2.0 * x[1] + x[2]) / 3.0, y[0] + (2.0 * y[1] + y[2]) / 3.0, z[0] + (2.0 * z[1] + z[2]) / 3.0);
  }
  else
  {
    // Equivalent to c(1.0)
    return Vector(x[0] + x[1] + x[2] + x[3], y[0] + y[1] + y[2] + y[3], z[0] + z[1] + z[2] + z[3]);
  }
}

/*!
\brief Approximate a cubic curve by two quadric curves.

N. Truong, C. Yuksel and L. Seiler. Quadratic Approximation of Cubic Curves. <I> Proceedings of the ACM on Computer Graphics and Interactive Techniques.</I> 2020.
\param gamma Subdision parameter, in [0,1].
\param a,b Returned quadric curves.
*/
void CubicCurve::Approximate(QuadricCurve& a, QuadricCurve& b, const double& gamma) const
{
  Vector c[4];
  for (int i = 0; i < 4; i++)
  {
    c[i] = BezierControl(i);
  }

  // Control points
  Vector q[5];

  q[0] = c[0];
  q[4] = c[3];
  q[1] = c[0] + 1.5 * gamma * (c[1] - c[0]);
  q[3] = c[3] - 1.5 * (1.0 - gamma) * (c[3] - c[2]);
  q[2] = (1.0 - gamma) * q[1] + gamma * q[3];

  a = QuadricCurve::Bezier(q[0], q[1], q[2]);
  b = QuadricCurve::Bezier(q[2], q[3], q[4]);
}

/*!
\brief Compute the set of frames along the curve.
\param step Distance beween samples.
\param vertical Boolean, set to true to compute frame along the horizontal curve first, false (default) otherwise.
*/
QVector<Frame> CubicCurve::GetFrames(const double& step, bool vertical) const
{
  QVector<Frame> frames;

  double length = Length(512);

  int n = int(length / step) + 1;

  frames.append(GetFrame(0.0, vertical));
  for (int i = 1; i < n - 1; i++)
  {
    double s = length * Math::Unit(i, n);
    double u = U(s);

    frames.append(GetFrame(u, vertical));

  }
  frames.append(GetFrame(1.0, vertical));

  return frames;
}