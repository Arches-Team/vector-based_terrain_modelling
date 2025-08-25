// Curves 

#include "libs/curve.h"
#include "libs/quartic.h"
#include "libs/polynomial.h"

#include <QtGui/QPainter>
#include <QtGui/QPainterPath>
#include <QtWidgets/QGraphicsScene>

/*!
\class CubicCurve2 curve.h
\brief %Cubic curves in the plane.

\ingroup PlanarGroup
*/

/*!
\brief Creates a cubic curve.

\param x, y %Cubic parametric equations of the corresponding x, y coordinates.
*/
CubicCurve2::CubicCurve2(const Cubic& x, const Cubic& y) :x(x), y(y)
{
}

/*!
\brief Creates a cubic curve.

\param c %Cubic curve.
*/
CubicCurve2::CubicCurve2(const CubicCurve& c)
{
  CubicCurve2::x = c[0];
  CubicCurve2::y = c[1];
}

/*!
\brief Computes the quadric curve whose components are the derivatives of the curve.

This function is useful for preprocessing tangent computations.
*/
QuadricCurve2 CubicCurve2::Tangent() const
{
  return QuadricCurve2(x.Prime(), y.Prime());
}

/*!
\brief Computes the tangent to the cubic curve at a given position on the curve.

The tangent is the obtained by evaluating the first derivative.
\param t Position.
*/
Vector2 CubicCurve2::Tangent(const double& t) const
{
  return Vector2(x.Prime()(t), y.Prime()(t));
}

/*!
\brief Computes the normal vector to the cubic curve at a given position on the curve.

The normal is the obtained by evaluating the second derivative.
\param t Position.
*/
Vector2 CubicCurve2::Normal(const double& t) const
{
  return Vector2(x.Second()(t), y.Second()(t));
}

/*!
\brief Computes the length of the cubic curve over interval [a,b].

Since we have the parametric equations of the coordinates of the points on the curve,
then the length is the integral of the square root of the sum of the squared derivatives
where the limits of integration are [a,b].

In the general case, the function we end up with is not integrable in closed
form, leading to an elliptic integral or some such. Thus we evaluate the integral
numerically.

\param a,b Interval.
\param n Discretization of the integration interval.
*/
double CubicCurve2::Length(const double& a, const double& b, int n) const
{
  double s = 0.0;

  // Compute tangent quadric curve
  QuadricCurve2 tangent = Tangent();

  // Evaluate quartic defined as the squared norm of the tangent coordinates
  Quartic quartic = (tangent[0]) * (tangent[0]) + (tangent[1]) * (tangent[1]);

  double epsilon = (b - a) / n;

  double u = sqrt(quartic(a));
  for (double t = a + epsilon; t <= b + epsilon / 2.0; t += epsilon)
  {
    double v = sqrt(quartic(t));
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
double CubicCurve2::InverseArcLength(const double& s, int n) const
{
  double length = 0.0;

  // Compute tangent quadric curve
  QuadricCurve2 tangent = Tangent();

  // Evaluate quartic defined as the squared norm of the tangent coordinates
  Quartic quartic = tangent[0] * tangent[0] + tangent[1] * tangent[1];

  double epsilon = 1.0 / n;

  double u = sqrt(quartic(0.0));
  double t = epsilon;
  while (length < s)
  {
    double v = sqrt(quartic(t));
    length += epsilon * (u + v) / 2.0;
    u = v;
  }
  return t - epsilon / 2.0;
}

/*!
\brief Computes the length of the cubic curve. Integration is performed over interval [0,1].

Simply call CubicCurve::Length(0.0,1.0,n) with n as the passed argument parameter.
\param n Discretization of the curve.
*/
double CubicCurve2::Length(int n) const
{
  return Length(0.0, 1.0, n);
}

/*!
\brief Negative curve unary operator.
*/
CubicCurve2 CubicCurve2::operator- () const
{
  return CubicCurve2(-x, -y);
}

/*!
\brief Evaluates curve point at a given location.
\param t Parameter.
*/
Vector2 CubicCurve2::operator() (const double& t) const
{
  return Vector2(x(t), y(t));
}

/*!
\brief Evaluates curve point at a given location.
\param t Parameter.
*/
Vector2 CubicCurve2::Eval(const double& t) const
{
  return Vector2(x(t), y(t));
}

/*!
\brief Creates an Hermite cubic curve on interval [0,1] given
vertex locations and tangent vectors (in that order).
\param a,b End vertices of the Hermite cubic curve.
\param ta,tb Tangent vectors at the end vertices.
*/
CubicCurve2 CubicCurve2::Hermite(const Vector2& a, const Vector2& b, const Vector2& ta, const Vector2& tb)
{
  return CubicCurve2(Cubic::Hermite(a[0], b[0], ta[0], tb[0]), Cubic::Hermite(a[1], b[1], ta[1], tb[1]));
}

/*!
\brief Creates a Bezier cubic curve on interval [0,1] given
four vertex locations.
\param a,b,c,d Vertexes of the Bezier curve.
*/
CubicCurve2 CubicCurve2::Bezier(const Vector2& a, const Vector2& b, const Vector2& c, const Vector2& d)
{
  return CubicCurve2(Cubic::Bezier(a[0], b[0], c[0], d[0]), Cubic::Bezier(a[1], b[1], c[1], d[1]));
}

/*!
\brief Overloaded.
\param s Stream.
\param curve The cubic curve.
*/
std::ostream& operator<<(std::ostream& s, const CubicCurve2& curve)
{
  s << "CubicCurve(" << curve.x << ',' << curve.y << ')';
  return s;
}

/*!
\brief Compute the bounding box of the cubic curve with parameter interval [0,1].
*/
Box2 CubicCurve2::GetBox() const
{
  Box2 box;

  // Define bounding box
  x.Range(box[0][0], box[1][0]);
  y.Range(box[0][1], box[1][1]);

  return box;
}

/*!
\brief Computes the signed curvature of the curve at a given point <B>p</B>(t).
\param t Parameter.
*/
double CubicCurve2::Curvature(const double& t) const
{
  double xp = x.Prime()(t);
  double yp = y.Prime()(t);
  double xpp = x.Second()(t);
  double ypp = y.Second()(t);

  return (ypp * xp - yp * xpp) / Math::Sqrt32(xp * xp + yp * yp);
}

/*!
\brief Compute the squared distance between
a point in the plane and the cubic curve.

This function computes the projection of a point p onto
the curve by minimizing the distance function along the curve,
expressed as a degree 5 polynomial.

Most of the coefficients of this polynomial are constant,
and could be pre-processed in the constructor.
\param p Point.
\param u Parameter defining the coordinate of the projection of the argument vertex onto the curve.
*/
double CubicCurve2::R(const Vector2& p, double& u) const
{
  Vector2 e = Normal(p, u);
  return e * e;
}

/*!
\brief Compute the parameter corresponding to the argument curvilign absisca of a point on the curve.

\param s Argument curvilign absisca (should be between 0.0 and the length of the curve).
\param n Discretization of the integration interval.
*/
double CubicCurve2::U(const double& s, int n) const
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
\brief Compute the curvilign absisca of a point on the curve.

\param t Parameter that should be within [0,1] interval.
\param n Discretization of the integration interval.
*/
double CubicCurve2::S(const double& t, int n) const
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
\brief Compute the projection of a point on the curve.
\param p Point.
\param u Parameter defining the projected point.
*/
Vector2 CubicCurve2::Normal(const Vector2& p, double& u) const
{
  Quartic xy;

  // Quartic constant coefficients
  xy[4] = 3.0 * (x[3] * x[3] + y[3] * y[3]);
  xy[3] = 5.0 * (x[3] * x[2] + y[3] * y[2]);
  xy[2] = 4.0 * (x[3] * x[1] + y[3] * y[1]) + 2.0 * (x[2] * x[2] + y[2] * y[2]);
  xy[1] = 3.0 * (x[2] * x[1] + y[2] * y[1]);
  xy[0] = x[1] * x[1] + y[1] * y[1];

  double t[5];

  Polynomial polynomial(
    xy[4], xy[3],
    xy[2],
    xy[1] + 3.0 * (x[3] * (x[0] - p[0]) + y[3] * (y[0] - p[1])),
    xy[0] + 2.0 * (x[2] * (x[0] - p[0]) + y[2] * (y[0] - p[1])),
    x[1] * (x[0] - p[0]) + y[1] * (y[0] - p[1]));

  int k = polynomial.Solve(t);

  // Search at caps
  Vector2 e = p - Vector2(x(0.0), y(0.0));
  double r = e * e;

  Vector2 f = p - Vector2(x(1.0), y(1.0));
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
\brief Compute the inverse mapping for a given input point.
\param p Point.
\param u, v Inverse Coordinates, u will denote the position along the curve and v the signed distance to the curve.
\return Squared distance to the curve.
*/
double CubicCurve2::UV(const Vector2& p, double& u, double& v) const
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
\brief Compute the inverse mapping for a given input point.

\sa CubicCurve2::UV(const Vector2&, double&, double&)

\param p Point.
\return Inverse uv coordinates compacted into a Vector2.
*/
Vector2 CubicCurve2::UV(const Vector2& p) const
{
  double u, v;
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
  return Vector2(u, v);
}

/*!
\brief Compute the Frenet vectors (tangent and normal) at a given point on the curve.

\param u Parameter.
*/
Matrix2 CubicCurve2::GetMatrix(const double& u) const
{
  Vector2 t = Tangent(u);
  Normalize(t);
  return Matrix2(t, t.Orthogonal());
}

/*!
\brief Compute the frame at a given point on the curve.

\param u Parameter.
*/
Frame2 CubicCurve2::GetFrame(const double& u) const
{
  return Frame2(GetMatrix(u), Eval(u));
}


/*!
\brief Compute the set of frames along the curve.
\param step Distance beween samples.
*/
QVector<Frame2> CubicCurve2::GetFrames(const double& step) const
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
\brief Draw a cubic curve.

The brush is not used.
\param scene Graphics scene.
\param pen The pen.
*/
void CubicCurve2::Draw(QGraphicsScene& scene, const QPen& pen, const QBrush&) const
{
  Vector2 pa = BezierControl(0);
  Vector2 pb = BezierControl(1);
  Vector2 pc = BezierControl(2);
  Vector2 pd = BezierControl(3);

  // Path
  QPainterPath path(QPointF(pa[0], pa[1]));

  QPen joining = pen;
  joining.setCapStyle(Qt::RoundCap); // Could have been FlatCap, however this yields discontinuities
  joining.setJoinStyle(Qt::RoundJoin);

  path.cubicTo(
    QPointF(pb[0], pb[1]),
    QPointF(pc[0], pc[1]),
    QPointF(pd[0], pd[1]));

  scene.addPath(path, joining);
}

/*!
\brief Compute the maximum curvature over a given interval.

Compute the analytic equation of the numerator of the curvature expression, search the roots, find extrema and returns  the maximum.
\param a,b Interval
*/
double CubicCurve2::Curvature(const double& a, const double& b) const
{
  double c = fabs(Curvature(a));
  c = Math::Max(c, fabs(Curvature(b)));

  Quadric xp = x.Prime();
  Quadric yp = y.Prime();
  Linear xpp = x.Second();
  Linear ypp = y.Second();

  Cubic n = xp * ypp - yp * xpp;
  Quartic d = xp * xp + yp * yp;

  // Derivative of numerator of curvature
  Polynomial np((-8.0) * n[2] * d[4], (-5.0) * n[2] * d[3] - 10.0 * n[1] * d[4], (-2.0) * n[2] * d[2] - 7.0 * n[1] * d[3] - 12.0 * n[0] * d[4], n[2] * d[1] - 4.0 * n[1] * d[2] - 9.0 * n[0] * d[3],
    4.0 * n[2] * d[0] - n[1] * d[1] - 6.0 * n[0] * d[2],
    2.0 * n[1] * d[0] - 3.0 * n[0] * d[1]
  );

  double t[6];

  // Check roots and compute extrema of curvature according to the previously computed polynomials
  int s = np.Solve(t, a, b);
  for (int i = 0; i < s; i++)
  {
    c = Math::Max(c, n(t[i]) / Math::Sqrt32(d(t[i])));

  }
  return c;
}

/*!
\brief Compute the intersection between the curve and a ray.

\param ray The ray.
\param t Intersection depths.
*/
int CubicCurve2::Intersect(const Ray2& ray, double t[3]) const
{
  Cubic q = ray.Direction()[0] * y - ray.Direction()[1] * x;
  q[0] += ray.Direction()[1] * ray.Origin()[0] - ray.Direction()[0] * ray.Origin()[1];

  double u[3];
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

  if (n == 2) {
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

  if (t[0] < 0.0)
  {
    if (t[1] < 0.0)
    {
      if (t[2] < 0.0)
      {
        return 0;
      }
      else
      {
        t[0] = t[2];
        return 1;
      }
    }
    else
    {
      if (t[2] < 0.0)
      {
        t[0] = t[1];
        return 1;
      }
      else
      {
        t[0] = t[1];
        t[1] = t[2];
        return 2;
      }
    }
  }
  else
  {
    if (t[1] < 0.0)
    {
      if (t[2] < 0.0)
      {
        return 1;
      }
      else
      {
        t[1] = t[2];
        return 2;
      }
    }
    else
    {
      if (t[2] < 0.0)
      {
        return 2;
      }
      else
      {
        return 3;
      }
    }
  }

  /*
  // Only one root
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
  else if (n == 2)
  {
    Math::Sort(t[0], t[1]);
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
  }
  else
  {
    Math::Sort(t[0], t[1], t[2]);
    if (t[0] < 0.0)
    {
      if (t[1] < 0.0)
      {
        if (t[2] < 0.0)
        {
          return 0;
        }
        else
        {
          t[0] = t[2];
          return 1;
        }
      }
      else
      {
        t[0] = t[1];
        t[1] = t[2];
        return 2;
      }
    }
    else
    {
      return 3;
    }
  }
  */
}
