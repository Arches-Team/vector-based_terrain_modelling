// Curves 

#include "libs/curve.h"
#include "libs/quartic.h"

#include <QtGui/QPainter>
#include <QtGui/QPainterPath>
#include <QtWidgets/QGraphicsScene>

/*!
\class QuadricCurve2 curve.h
\brief %Quadric curves in the plane.

Quadric curves are implemented with two quadric polynomials defining the
coordinates of the points on the curve. It is a minimum storage class in
the sense that it does not keep extra constant data which may be useful
for computing the bounding box of the quadric curve parameterized over
the unit interval [0,1] or the distance between a point an the curve.

\ingroup PlanarGroup

*/

/*!
\brief Creates a quadric curve.
\param x, y Quadric parametric equations of the corresponding x, y coordinates.
*/
QuadricCurve2::QuadricCurve2(const Quadric& x, const Quadric& y) :x(x), y(y)
{
}

/*!
\brief Creates a quadric curve in the plane.

\param c %Quadric curve.
*/
QuadricCurve2::QuadricCurve2(const QuadricCurve& c) :x(c[0]), y(c[1])
{
}

/*!
\brief Negates the coefficients of a quadric curve.
*/
QuadricCurve2 QuadricCurve2::operator- () const
{
  return QuadricCurve2(-x, -y);
}

/*!
\brief Computes the point on the curve.
\param t Parameter.
*/
Vector2 QuadricCurve2::operator() (const double& t) const
{
  return Vector2(x(t), y(t));
}

/*!
\brief Computes the point on the curve.
\param t Parameter.
*/
Vector2 QuadricCurve2::Eval(const double& t) const
{
  return Vector2(x(t), y(t));
}

/*!
\brief Computes the signed curvature of the curve at a given point <B>p</B>(t).
\param t Parameter.
*/
double QuadricCurve2::Curvature(const double& t) const
{
  double xp = 2.0 * x[2] * t + x[1];
  double yp = 2.0 * y[2] * t + y[1];
  double xpp = 2.0 * x[2];
  double ypp = 2.0 * y[2];

  return (ypp * xp - yp * xpp) / Math::Sqrt32(xp * xp + yp * yp);
}

/*!
\brief Computes the tangent to the curve.
\param t Parameter.
*/
Vector2 QuadricCurve2::Tangent(const double& t) const
{
  return Vector2(2.0 * x[2] * t + x[1], 2.0 * y[2] * t + y[1]);
}

/*!
\brief Computes the normal vector to the quadric curve at a given position on the curve.

The normal is the obtained by evaluating the second derivative.
\param t Position (not used).
*/
#pragma warning(push)
#pragma warning(disable: 4100)  
Vector2 QuadricCurve2::Normal(const double& t) const
{
  return Vector2(2.0 * x[2], 2.0 * y[2]);
}
#pragma warning(pop)

/*!
\brief Translate the curve.
\param t Translation vector.
*/
void QuadricCurve2::Translate(const Vector2& t)
{
  x[0] += t[0];
  y[0] += t[1];
}

/*!
\brief Creates a quadratic Bezier curve on interval [0,1].
\param a,b,c The control points of the Bezier curve.
*/
QuadricCurve2 QuadricCurve2::Bezier(const Vector2& a, const Vector2& b, const Vector2& c)
{
  return QuadricCurve2(Quadric::Bezier(a[0], b[0], c[0]), Quadric::Bezier(a[1], b[1], c[1]));
}

/*!
\brief Creates a quadratic Bezier curve on interval [0,1] that passes through argument points a,b and c for parameters 0, t, and 1.
\param a,b,c Points of the quadric curve.
\param t Parameter for second point.

\author Lea Godot, Eric Galin.
*/
QuadricCurve2 QuadricCurve2::ThroughPoints(const Vector2& a, const Vector2& b, const Vector2& c, const double& t)
{
  // Matrix with parameters 0.0, t and 1.0
  Matrix A(0.0 * 0.0, t * t, 1.0 * 1.0, 0.0, t, 1.0, 1.0, 1.0, 1.0);

  // Solve
  Matrix IA = Inverse(A);

  Vector cx = IA * Vector(a[0], b[0], c[0]);
  Vector cy = IA * Vector(a[1], b[1], c[1]);

  return QuadricCurve2(Quadric(cx[0], cx[1], cx[2]), Quadric(cy[0], cy[1], cy[2]));
}


/*!
\brief Overloaded.
\param s Stream.
\param curve The quadric curve.
*/
std::ostream& operator<<(std::ostream& s, const QuadricCurve2& curve)
{
  s << "QuadricCurve2(" << curve.x << ',' << curve.y << ')';
  return s;
}

/*!
\brief Compute the bounding box of the quadric curve with parameter interval [0,1].
*/
Box2 QuadricCurve2::GetBox() const
{
  Box2 box;

  // Define bounding box
  x.Range(box[0][0], box[1][0]);
  y.Range(box[0][1], box[1][1]);

  return box;
}

/*!
\brief Compute the squared distance between a point and the quadric curve.

This function computes the projection of the point p onto
the curve by minimizing the distance function along the curve,
expressed as a cubic polynomial.

Most of the coefficients of this cubic polynomial are constant,
and could be pre-processed in the constructor.

\warn This function does not work if the quadric degenerates to a linear curve.

\param p Vertex.
\param u Parameter defining the coordinate of the projection of the argument vertex onto the curve.
*/
double QuadricCurve2::R(const Vector2& p, double& u) const
{
  double xy[3];
  double aa[3];

  // Quadric constant coefficients
  xy[2] = 2.0 * (x[2] * x[2] + y[2] * y[2]);
  xy[1] = 3.0 * (x[2] * x[1] + y[2] * y[1]);
  xy[0] = x[1] * x[1] + y[1] * y[1];

  // Normalize by highest coefficient
  xy[1] /= xy[2];
  xy[0] /= xy[2];

  aa[0] = xy[1] * xy[1];
  aa[1] = 2.0 / xy[2];
  aa[2] = 1.0 / xy[2];

  double t[3];

  // Normalized form
  Cubic c(
    1.0,
    xy[1],
    xy[0] + aa[1] * (x[2] * (x[0] - p[0]) + y[2] * (y[0] - p[1])),
    aa[2] * ((x[0] - p[0]) * x[1] + (y[0] - p[1]) * y[1]));

  int k = c.SolveNormalized(t);

  // Search at caps
  Vector2 e = p - Eval(0.0);
  double r = e * e;

  Vector2 f = p - Eval(1.0);
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
      f = p - Vector2(x(t[i]), y(t[i]));
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
double QuadricCurve2::S(const double& t, int n) const
{
  const double epsilon = 1.0 / n;
  double s = 0.0;
  Vector2 a = Vector2(x(0.0), y(0.0));
  for (double u = epsilon; u < t; u += epsilon)
  {
    Vector2 b = Vector2(x(u), y(u));
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
double QuadricCurve2::U(const double& s, int n) const
{
  if (s == 0.0) return 0.0;
  const double epsilon = 1.0 / n;
  double ss = 0.0;
  Vector2 a = Vector2(x(0.0), y(0.0));
  for (double u = epsilon; u < 1.0; u += epsilon)
  {
    Vector2 b = Vector2(x(u), y(u));
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
double QuadricCurve2::Length(const double& a, const double& b, int n) const
{
  double s = 0.0;

  // Evaluate quadric defined as the squared norm of the tangent coordinates
  Quadric quadric = Quadric(4.0 * (x[2] * x[2] + y[2] * y[2]), 4.0 * (x[2] * x[1] + y[2] * y[1]), x[1] * x[1] + y[1] * y[1]);

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
\brief Computes the length of the quadric curve over interval [a,b] using a closed form integral.

\param a,b Interval.
*/
double QuadricCurve2::LengthIntegral(const double& a, const double& b) const
{
  // Integral
  auto Mu = [](double u, double k)-> double {
    double v = sqrt(u * u + k);
    return u * v + k * log(u + v);
    };

  const Vector2 a2(x[2], y[2]);
  const Vector2 a1(x[1], y[1]);

  double ab = a2 * a1;
  double aa = a2 * a2;

  double d = ab / (2.0 * aa);
  double k = (aa * (a1 * a1) - Math::Sqr(ab)) / (4.0 * aa * aa);
  double l = sqrt(aa) * (Mu(b + d, k) - Mu(a + d, k));

  return l;
}

/*!
\brief Computes the length of the cubic curve. Integration is performed over interval [0,1].
*/
double QuadricCurve2::Length(int n) const
{
  return Length(0.0, 1.0, n);
}

/*!
\brief Compute the inverse mapping for a given input point.
\param p Point.
\param u, v Inverse Coordinates, u will denote the position along the curve and v the signed distance to the curve.
*/
double QuadricCurve2::UV(const Vector2& p, double& u, double& v) const
{
  Vector2 n = Normal(p, u);

  if (u == 0.0)
  {
    Vector2 na = Normalized(Vector2(-y.Prime()(u), x.Prime()(u)));
    v = n * na;
  }
  else if (u == 1.0)
  {
    Vector2 nb = Normalized(Vector2(-y.Prime()(u), x.Prime()(u)));
    v = n * nb;
  }
  else
  {
    v = Norm(n);
    Vector2 nc = Vector2(-y.Prime()(u), x.Prime()(u));
    if (n * nc < 0.0)
    {
      v = -v;
    }
  }
  return n * n;
}

/*!
\brief Compute the projection of a point on the curve.
\param p Point.
\param u Parameter defining the projected point.
*/
Vector2 QuadricCurve2::Normal(const Vector2& p, double& u) const
{
  double xy[3];
  double aa[3];

  // Quadric constant coefficients
  xy[2] = 2.0 * (x[2] * x[2] + y[2] * y[2]);
  xy[1] = 3.0 * (x[2] * x[1] + y[2] * y[1]);
  xy[0] = x[1] * x[1] + y[1] * y[1];

  // Normalize by highest coefficient
  xy[1] /= xy[2];
  xy[0] /= xy[2];

  aa[0] = xy[1] * xy[1];
  aa[1] = 2.0 / xy[2];
  aa[2] = 1.0 / xy[2];

  double t[3];

  // Normalized form
  Cubic c(
    1.0,
    xy[1],
    xy[0] + aa[1] * (x[2] * (x[0] - p[0]) + y[2] * (y[0] - p[1])),
    aa[2] * ((x[0] - p[0]) * x[1] + (y[0] - p[1]) * y[1]));

  int k = c.SolveNormalized(t);

  // Search at caps
  Vector2 e = p - (*this)(0.0);
  double r = e * e;

  Vector2 f = p - (*this)(1.0);
  double s = f * f;

  u = 0.0;
  if (s < r)
  {
    e = f;
    r = s;
    u = 1.0;
  }

  for (int i = 0; i < k; i++)
  {
    if ((t[i] > 0.0) && (t[i] < 1.0))
    {
      f = p - Vector2(x(t[i]), y(t[i]));
      s = f * f;
      if (s < r)
      {
        e = f;
        u = t[i];
        r = s;
      }
    }
  }
  return e;
}

/*!
\brief Draw a quadric curve.
\param scene Graphics scene.
\param pen The pen.
\param brush The brush, should the triangle be filled.
*/
#pragma warning(push)
#pragma warning(disable: 4100)  
void QuadricCurve2::Draw(QGraphicsScene& scene, const QPen& pen, const QBrush& brush) const
{
  Vector2 h0 = (*this)(0.0);
  Vector2 h1 = (*this)(1.0);
  Vector2 h2 = Tangent(0.0);
  Vector2 h3 = Tangent(1.0);

  Vector2 b0 = h0;
  Vector2 b1 = (3.0 * h0 + h2) / 3.0;
  Vector2 b2 = (3.0 * h1 - h3) / 3.0;
  Vector2 b3 = h1;

  QPointF qb0 = QPointF(b0[0], b0[1]);
  QPointF qb1 = QPointF(b1[0], b1[1]);
  QPointF qb2 = QPointF(b2[0], b2[1]);
  QPointF qb3 = QPointF(b3[0], b3[1]);

  QPainterPath path(qb0);

  path.cubicTo(qb1, qb2, qb3);

  scene.addPath(path, pen);
}
#pragma warning(pop)

/*!
\brief Create a set of quadratic Bezier curves.
\param c Control points of the Bezier curve.
*/
QVector<QuadricCurve2> QuadricCurve2::Bezier(const QVector<Vector2>& c)
{
  const int n = c.size();

  QVector<Vector2> knot;
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

  QVector<QuadricCurve2> array(n - 2);

  for (int i = 0; i < n - 2; i++)
  {
    array[i] = QuadricCurve2::Bezier(knot[2 * i], knot[2 * i + 1], knot[2 * i + 2]);
  }

  return array;
}

/*!
\brief Compute the intersection between the curve and a ray.

\param ray The ray.
\param t Intersection depths.
*/
int QuadricCurve2::Intersect(const Ray2& ray, double t[2]) const
{
  Quadric q = ray.Direction()[0] * y - ray.Direction()[1] * x;
  q[0] += ray.Direction()[1] * ray.Origin()[0] - ray.Direction()[0] * ray.Origin()[1];

  double u[2];
  // Solve quadric equation over unit interval
  int n = q.Solve(u, 0.0, 1.0);

  if (n == 0)
  {
    return 0;
  }

  // Compute intersection depths
  if (fabs(ray.Direction()[0]) > fabs(ray.Direction()[1]))
  {
    for (int i = 0; i < n; i++)
    {
      t[i] = (x(u[i]) - ray.Origin()[0]) / ray.Direction()[0];
    }
  }
  else
  {
    for (int i = 0; i < n; i++)
    {
      t[i] = (y(u[i]) - ray.Origin()[1]) / ray.Direction()[1];
    }
  }

  if (n == 1)
  {
    if (t[0] < 0.0)
    {
      return 0;
    }
    else
    {
      return 1;
    }
  }

  // Intersection along ray only yields positive solutions
  if (t[0] < 0.0)
  {
    if (t[1] < 0.0)
    {
      return 0;
    }
    else
    {
      t[0] = t[1];
      return 1;
    }
  }
  else
  {
    if (t[1] < 0.0)
    {
      return 1;
    }
    else
    {
      return 2;
    }
  }
}


/*!
\brief Compute the Frenet vectors (tangent and normal) at a given point on the curve.

\param u Parameter.
*/
Matrix2 QuadricCurve2::GetMatrix(const double& u) const
{
  Vector2 t = Tangent(u);
  Normalize(t);
  return Matrix2(t, t.Orthogonal());
}

/*!
\brief Compute the frame at a given point on the curve.

\param u Parameter.
*/
Frame2 QuadricCurve2::GetFrame(const double& u) const
{
  return Frame2(GetMatrix(u), Eval(u));
}


/*!
\brief Compute the set of frames along the curve.
\param step Distance beween samples.
*/
QVector<Frame2> QuadricCurve2::GetFrames(const double& step) const
{
  QVector<Frame2> frames;

  double length = Length(512);

  int n = int(length / step) + 1;

  frames.append(GetFrame(0.0));
  for (int i = 0; i < n - 1; i++)
  {
    double s = length * Math::Unit(i, n);
    double u = U(s);

    frames.append(GetFrame(u));

  }
  frames.append(GetFrame(1.0));

  return frames;
}

/*!
\brief Compute the closest point on the curve.

Given a set of curvilear abscica u, possible outside [0,1], compute which point on the curve with abscica taken amont u and {0,1} is closest to the point p.

\param p Point.
\param t Set of curvilinear abscica.
\param n Size.
*/
double QuadricCurve2::R(const Vector2& p, double* t, int n) const
{
  // Search at caps
  double u = SquaredNorm(p - Eval(0.0));
  double t0 = 0.0;

  double v = SquaredNorm(p - Eval(1.0));
  if (v < u)
  {
    t0 = 1.0;
    u = v;
  }

  for (int i = 0; i < n; i++)
  {
    if ((t[i] > 0.0) && (t[i] < 1.0))
    {
      v = SquaredNorm(p - Vector2(x(t[i]), y(t[i])));

      if (v < u)
      {
        t0 = t[i];
        u = v;
      }
    }
  }

  return u;
}
