// Quadric curves 

#include "libs/curve.h"
#include "libs/quartic.h"

/*!
\class QuadricCurve curve.h
\brief %Quadric curves.

Quadric curves are implemented with three quadric polynomials defining the
coordinates of the points on the curve. It is a minimum storage class in
the sense that it does not keep extra constant data
for accelerating the computation of the bounding box of the quadric curve parameterized over
the unit interval [0,1], or the distance between a point an the curve.

\ingroup KernelGroup

*/

/*!
\brief Creates a quadric curve.
\param x, y, z Quadric parametric equations of the corresponding x, y, z coordinates.
*/
QuadricCurve::QuadricCurve(const Quadric& x, const Quadric& y, const Quadric& z) :x(x), y(y), z(z)
{
}

/*!
\brief Negates the coefficients of a quadric curve.
*/
QuadricCurve QuadricCurve::operator- () const
{
  return QuadricCurve(-x, -y, -z);
}

/*!
\brief Computes the point on the curve.
\param t Parameter.
*/
Vector QuadricCurve::operator() (const double& t) const
{
  return Vector(x(t), y(t), z(t));
}

/*!
\brief Computes the point on the curve.

This function is redundant with the operator().
\param t Parameter.
*/
Vector QuadricCurve::Eval(const double& t) const
{
  return Vector(x(t), y(t), z(t));
}

/*!
\brief Computes the curvature of the curve at a given positon.
\param t Parameter.
*/
double QuadricCurve::Curvature(const double& t) const
{
  // First and second derivatives, could have used Prime() and Second() instead of inlining computations
  double xp = 2.0 * x[2] * t + x[1];
  double yp = 2.0 * y[2] * t + y[1];
  double zp = 2.0 * z[2] * t + z[1];
  double xpp = 2.0 * x[2];
  double ypp = 2.0 * y[2];
  double zpp = 2.0 * z[2];

  return sqrt((zpp * yp - ypp * zp) * (zpp * yp - ypp * zp) + (xpp * zp - xp * zpp) * (xpp * zp - xp * zpp) + (ypp * xp - yp * xpp) * (ypp * xp - yp * xpp)) / Math::Sqrt32(xp * xp + yp * yp + zp * zp);
}

/*!
\brief Computes the tangent to the curve.
\param t Parameter.
*/
Vector QuadricCurve::Tangent(const double& t) const
{
  return Vector(2.0 * x[2] * t + x[1], 2.0 * y[2] * t + y[1], 2.0 * z[2] * t + z[1]);
}

/*!
\brief Computes the secondary derivative to the curve.

Note that although the parameter has no influence over the secondary derivative for quadric curves, it
is still there out of closed form consistency.
*/
Vector QuadricCurve::Normal(const double&) const
{
  return Vector(2.0 * x[2], 2.0 * y[2], 2.0 * z[2]);
}

/*!
\brief Create a quadratic Bezier curve on interval [0,1].
\param a,b,c The control points of the Bezier curve.
*/
QuadricCurve QuadricCurve::Bezier(const Vector& a, const Vector& b, const Vector& c)
{
  return QuadricCurve(Quadric::Bezier(a[0], b[0], c[0]), Quadric::Bezier(a[1], b[1], c[1]), Quadric::Bezier(a[2], b[2], c[2]));
}

/*!
\brief Create a set of quadratic Bezier curves.
\param c Control points of the Bezier curve.
*/
QVector<QuadricCurve> QuadricCurve::Bezier(const QVector<Vector>& c)
{
  const int n = c.size();

  QVector<Vector> knot;
  knot.clear();

  // Starting vertices
  knot.append(c[0]);
  knot.append(c[1]);

  // Mid vertices
  for (int i = 1; i < n - 2; i++)
  {
    knot.append(0.5 * (c[i] + c[i + 1]));
    knot.append(c[i + 1]);
  }

  // End vertex
  knot.append(c[n - 1]);

  QVector<QuadricCurve> array(n-2);

  for (int i = 0; i < n - 2; i++)
  {
    array[i] = QuadricCurve::Bezier(knot[2 * i], knot[2 * i + 1], knot[2 * i + 2]);
  }
  return array;
}

/*!
\brief Overloaded.
\param s Stream.
\param curve The quadric curve.
*/
std::ostream& operator<<(std::ostream& s, const QuadricCurve& curve)
{
  s << "QuadricCurve(" << curve.x << ',' << curve.y << ',' << curve.z << ')';
  return s;
}

/*!
\brief Compute the bounding box of the quadric curve with parameter interval [0,1].
*/
Box QuadricCurve::GetBox() const
{
  Box box;

  // Define bounding box
  x.Range(box[0][0], box[1][0]);
  y.Range(box[0][1], box[1][1]);
  z.Range(box[0][2], box[1][2]);

  return box;
}

/*!
\brief Compute signed squared distance between a point and the quadric curve.

Simply compute the square root of the squared distance.

\sa double QuadricCurve::R(const Vector& , double& ) const

\param p Point.
*/
double QuadricCurve::Signed(const Vector& p) const
{
  double u;
  return sqrt(R(p, u));
}

/*!
\brief Compute the squared distance between a point and the quadric curve.

This function computes the projection of a point p onto
the curve by minimizing the distance function along the curve,
expressed as a cubic polynomial.

Most of the coefficients of this cubic polynomial are constant,
and could be pre-processed in the constructor.

\warn This function does not work if the quadric degenerates to a linear curve.

\param p Point.
\param u Parameter defining the coordinate of the projection of the argument vertex onto the curve.
*/
double QuadricCurve::R(const Vector& p, double& u) const
{
  double xyz[3];
  double aa[3];

  // Quadric constant coefficients
  xyz[2] = 2.0 * (x[2] * x[2] + y[2] * y[2] + z[2] * z[2]);
  xyz[1] = 3.0 * (x[2] * x[1] + y[2] * y[1] + z[2] * z[1]);
  xyz[0] = x[1] * x[1] + y[1] * y[1] + z[1] * z[1];

  // Normalize by highest coefficient
  xyz[1] /= xyz[2];
  xyz[0] /= xyz[2];

  aa[0] = xyz[1] * xyz[1];
  aa[1] = 2.0 / xyz[2];
  aa[2] = 1.0 / xyz[2];

  double t[3];

  // Normalized form
  Cubic c(
    1.0,
    xyz[1],
    xyz[0] + aa[1] * (x[2] * (x[0] - p[0]) + y[2] * (y[0] - p[1]) + z[2] * (z[0] - p[2])),
    aa[2] * ((x[0] - p[0]) * x[1] + (y[0] - p[1]) * y[1] + (z[0] - p[2]) * z[1]));

  int k = c.SolveNormalized(t);

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
\brief Compute the curvilign absisca of a point on the curve.

\param t Parameter.
\param n Discretization of the integration interval.
*/
double QuadricCurve::S(const double& t, int n) const
{
  const double epsilon = 1.0 / n;
  double s = 0.0;
  Vector a = Vector(x(0.0), y(0.0), z(0.0));
  for (double u = epsilon; u < t; u += epsilon)
  {
    Vector b = Vector(x(u), y(u), z(u));
    s += Norm(b - a);
    a = b;
  }
  return s;
}

/*!
\brief Compute the parameter corresponding to the argument curvilign absisca of a point on the curve.

\param s Argument curvilign absisca (should be between 0.0 and the length of the curve).
\param n Discretization of the integration interval.
*/
double QuadricCurve::U(const double& s, int n) const
{
  if (s == 0.0) return 0.0;
  const double epsilon = 1.0 / n;
  double ss = 0.0;
  Vector a = Vector(x(0.0), y(0.0), z(0.0));
  for (double u = epsilon; u < 1.0; u += epsilon)
  {
    Vector b = Vector(x(u), y(u), z(u));
    ss += Norm(b - a);
    if (ss > s)
    {
      return u - epsilon / 2.0;
    }
    a = b;
  }
  return 1.0;
}

/*!
\brief Computes the length of the quadric curve over interval [a,b].

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
double QuadricCurve::Length(const double& a, const double& b, int n) const
{
  double s = 0.0;

  // Evaluate quadric defined as the squared norm of the tangent coordinates
  Quadric quadric = Quadric(4.0 * (x[2] * x[2] + y[2] * y[2] + z[2] * z[2]), 4.0 * (x[2] * x[1] + y[2] * y[1] + z[2] * z[1]), x[1] * x[1] + y[1] * y[1] + z[1] * z[1]);

  double epsilon = (b - a) / n;

  double u = sqrt(quadric(a));
  for (double t = a + epsilon; t <= b + epsilon / 2.0; t += epsilon)
  {
    double v = sqrt(quadric(t));
    s += (u + v) / 2.0;
    u = v;
  }
  s *= epsilon;

  return s;
}

/*!
\brief Computes the length of the cubic curve.

Integration is performed over interval [0,1] using a uniform discretization.
\param n Discretization.
*/
double QuadricCurve::Length(int n) const
{
  return Length(0.0, 1.0, n);
}

/*!
\brief Computes the length of the quadric curve over interval [a,b] using a closed form integral.

\param a,b Interval.
*/
double QuadricCurve::LengthIntegral(const double& a, const double& b) const
{
  // Integral
  auto Mu = [](double u, double k)-> double {
    double v = sqrt(u * u + k);
    return u * v + k * log(u + v);
    };

  const Vector a2(x[2], y[2], z[2]);
  const Vector a1(x[1], y[1], z[1]);

  double ab = a2 * a1;
  double aa = a2 * a2;

  double d = ab / (2.0 * aa);
  double k = (aa * (a1 * a1) - Math::Sqr(ab)) / (4.0 * aa * aa);
  double l = sqrt(aa) * (Mu(b + d, k) - Mu(a + d, k));

  return l;
}

/*!
\brief Compute the Frenet vectors at a given point on the curve.

The columns of the returned matrix contain the tangent, normal and binormal vectors.

\param u Parameter.
\param vertical Boolean, set to true to compute frame along the horizontal curve first, false (default) otherwise.
*/
Matrix QuadricCurve::GetMatrix(const double& u, bool vertical) const
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
    Vector n = Vector::Z /t;
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
Frame QuadricCurve::GetFrame(const double& u, bool vertical) const
{
  return Frame(GetMatrix(u,vertical), Eval(u));
}

/*!
\brief Compute the set of frames along the curve.
\param step Distance beween samples.
\param vertical Boolean, set to true to compute frame along the horizontal curve first, false (default) otherwise.
*/
QVector<Frame> QuadricCurve::GetFrames(const double& step, bool vertical) const
{
  QVector<Frame> frames;

  double length = Length(512);

  int n = int(length / step) + 1;

  frames.append(GetFrame(0.0,vertical));
  for (int i = 0; i < n - 1; i++)
  {
    double s = length * Math::Unit(i, n);
    double u = U(s);

    frames.append(GetFrame(u,vertical));

  }
  frames.append(GetFrame(1.0,vertical));

  return frames;
}