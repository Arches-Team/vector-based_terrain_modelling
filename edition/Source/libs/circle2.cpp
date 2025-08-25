// Circle

#include <QtWidgets/QGraphicsScene>

#include "libs/circle.h"
#include "libs/polygon.h"
#include "libs/curve.h"

/*!
\class Circle2 circle.h
\brief Circles in the plane.

\sa Circle
\ingroup PlanarGroup
*/

const Circle2 Circle2::Unit(1.0);
const Circle2 Circle2::Infinite(Math::Infinity);

/*!
\brief Create a circle.
\param c Center.
\param r Radius.
*/
Circle2::Circle2(const Vector2& c, const double& r) :c(c), r(r)
{
}

/*!
\brief Create a circle centered at origin.
\param r Radius.
*/
Circle2::Circle2(const double& r) :Circle2(Vector2::Null, r)
{
}

/*!
\brief Create the circle passing through 3 points.

This is the circumscribed circle of the triangle.

\param a, b, c Points.
*/
Circle2::Circle2(const Vector2& a, const Vector2& b, const Vector2& c)
{
  Vector2 ab = b - a;
  Vector2 ac = c - a;
  Vector2 bc = c - b;

  double E = ab[0] * (a[0] + b[0]) + ab[1] * (a[1] + b[1]);
  double F = ac[0] * (a[0] + c[0]) + ac[1] * (a[1] + c[1]);

  double G = 2.0 * (ab[0] * bc[1] - ab[1] * bc[0]);

  double x = (ac[1] * E - ab[1] * F) / G;
  double y = (ab[0] * F - ac[0] * E) / G;

  Circle2::c = Vector2(x, y);
  Circle2::r = Norm(a - Circle2::c);
}

/*!
\brief Vertex on the circle.
\param a Angle.
*/
Vector2 Circle2::Vertex(const double& a) const
{
  return c + r * Vector2::Polar(a);
}

/*!
\brief Check if a point is inside the circle.
\param p Point.
*/
bool Circle2::Inside(const Vector2& p) const
{
  if (SquaredNorm(p - c) <= r * r)
  {
    return true;
  }
  else
  {
    return false;
  }
}

/*!
\brief Test if a circle is entirely embedded in another one.

\param circle The other circle.
*/
bool Circle2::Inside(const Circle2& circle) const
{
  // Check if circles do not intersect
  double cc = SquaredNorm(c - circle.c);
  double rr = Math::Sqr(r + circle.r);
  if (cc >= rr)
  {
    return false;
  }

  // Sort radii so that ra is the large radius and rb is the small one
  double ra = r;
  double rb = circle.r;
  Math::Sort(rb, ra);

  // Distance between centers
  double scc = sqrt(cc);

  // Small circle is fully inside the large one
  if (scc + rb <= ra)
  {
    return true;
  }

  return false;
}

/*!
\brief Test if a box is inside the circle.

\param box The box.
*/
bool Circle2::Inside(const Box2& box) const
{
  for (int i = 0; i < 4; i++)
  {
    if (!Inside(box.Vertex(i)))
    {
      return false;
    }
  }

  return true;
}

/*!
\brief Check if a point is within a given range of the circle.
\param p Point.
\param e Range.
*/
bool Circle2::InsideRange(const Vector2& p, const double& e) const
{
  if (SquaredNorm(p - c) <= (r + e) * (r + e))
  {
    return true;
  }
  else
  {
    return false;
  }
}

/*!
\brief Check intersection with a box.

Simply compute the distance between the center of the
circle and the box and compare to radius of the circle.

\sa Box::R(const Vector&) const
\param box The box.
*/
bool Circle2::Intersect(const Box2& box) const
{
  return (box.R(c) < r * r);
}

/*!
\brief Draw a circle.
\param scene Graphics scene.
\param pen The pen.
\param brush The brush, should the circle be filled.
*/
void Circle2::Draw(QGraphicsScene& scene, const QPen& pen, const QBrush& brush) const
{
  scene.addEllipse(c[0] - r, c[1] - r, r + r, r + r, pen, brush);
}

/*!
\brief Generate a random vector inside the circle.
\param random %Random number generator.
*/
Vector2 Circle2::RandomInside(Random& random) const
{
  // Unit radius
  double u;

  // Random samples
  double rx, ry;
  do
  {
    rx = 2.0 * random.Uniform() - 1.0;
    ry = 2.0 * random.Uniform() - 1.0;
    u = rx * rx + ry * ry;
  } while (u >= 1.0);

  return c + r * Vector2(rx, ry);
}

/*!
\brief Compute a random vertex inside a unit disc.
\param random %Random number generator.
*/
Vector2 Circle2::RandomUnit(Random& random)
{
  // Unit radius
  double u;

  // Random samples
  double rx, ry;
  do
  {
    rx = 2.0 * random.Uniform() - 1.0;
    ry = 2.0 * random.Uniform() - 1.0;
    u = rx * rx + ry * ry;
  } while (u >= 1.0);

  return Vector2(rx, ry);
}

/*!
\brief Generate a random vector on the circle.
\param random %Random number generator.
*/
Vector2 Circle2::RandomOn(Random& random) const
{
  Vector2 uv = RandomUnit(random);

  double s = SquaredNorm(uv);

  double rx = (uv[0] * uv[0] - uv[1] * uv[1]) / s;
  double ry = 2.0 * uv[0] * uv[1] / s;

  return c + r * Vector2(rx, ry);
}

/*!
\brief Generate a spiraling vector in the disc.
\param i Integer, should be <n.
\param n Number of points.
*/
Vector2 Circle2::Vogel(int i, int n) const
{
  const double goldenAngle = Math::Pi * (3.0 - sqrt(5.0));

  double r = sqrt((i + 0.5) / n);
  double a = i * goldenAngle; // Could add a small offset

  return r * Vector2(cos(a), sin(a));
}

/*!
\brief Creates a circle enclosing a set of points.

The algorithm computes a fast approximation of the
bounding ball for a point set based on the algorithm
given by [Jack Ritter, 1990].

\param p Array of vertices.
\param n Number of vertices.
*/
Circle2::Circle2(Vector2* p, int n)
{
  Box2 box = Box2(p[0], p[0]);

  int Pxmin = 0;
  int Pxmax = 0;
  int Pymin = 0;
  int Pymax = 0;

  for (int i = 1; i < n; i++)
  {
    if (p[i][0] < box[0][0])
    {
      box[0][0] = p[i][1];
      Pxmin = i;
    }
    else if (p[i][0] > box[1][0])
    {
      box[1][0] = p[i][1];
      Pxmax = i;
    }
    if (p[i][1] < box[0][1])
    {
      box[0][1] = p[i][2];
      Pymin = i;
    }
    else if (p[i][1] > box[1][1])
    {
      box[1][1] = p[i][2];
      Pymax = i;
    }
  }

  // Select the largest extent as an initial diameter for the ball
  Vector2 x = p[Pxmax] - p[Pxmin];
  Vector2 y = p[Pymax] - p[Pymin];

  // Direction is largest extent
  if (x * x >= y * y)
  {
    c = (p[Pxmin] + p[Pxmax]) * 0.5;
    r = x * x * 0.25;
  }
  else
  {
    c = (p[Pymin] + p[Pymax]) * 0.5;
    r = y * y * 0.25;
  }

  // Check that all points p[i] are in the circle
  // If not, expand the circle just enough to include them
  for (int j = 0; j < n; j++)
  {
    Vector2 cp = p[j] - c;
    double s = cp * cp;
    // p[i] is outside the circle, so expand to include it
    if (s > r)
    {
      // Compute radii
      s = sqrt(s);
      r = sqrt(r);

      // Shift center towards p[i]
      c += (0.5 * (s - r) / s) * cp;

      // Enlarge radius just enough
      r = 0.5 * (r + s);

      // Restore squared radius
      r *= r;
    }
  }

  // Radius
  r = sqrt(r);
}

/*!
\brief Extend the radius of the circle.

\param re Radius extension.
*/
Circle2 Circle2::Extended(const double& re) const
{
  return Circle2(c, Circle2::r + re);
}

/*!
\brief Extend the circle so that the argument point should be embedded in the new circle.

This is the same as computing the smallest enclosing circle of a circle and a point.

\param p Point.
*/
void Circle2::Extend(const Vector2& p)
{
  Vector2 pc = p - c;

  // Squared distance
  double re = pc * pc;

  // Escape if point lies inside the circle
  if (re < r * r)
    return;

  // Distance
  re = sqrt(re);

  // Normalized direction
  pc /= re;

  // Compute hal difference
  re = 0.5 * (re - r);

  // Move center by half difference
  c += pc * re;

  r += re;
}

/*!
\brief Compute the area of the surface between two circles.

If the circles do not intersect, the function returns 0; if the smallest
is included in the largest, the function returns the area of the smallest,
otherwise complex computations are performed.

\param circle The other circle.
*/
double Circle2::Area(const Circle2& circle) const
{
  // Check if circles do not intersect
  double cc = SquaredNorm(c - circle.c);
  double rr = Math::Sqr(r + circle.r);
  if (cc >= rr)
  {
    return 0.0;
  }

  // Sort radii so that ra is the large radius and rb is the small one
  double ra = r;
  double rb = circle.r;
  Math::Sort(rb, ra);

  // Distance between centers
  double scc = sqrt(cc);

  // Small circle is fully inside the large one
  if (scc + rb <= ra)
  {
    return Math::Pi * rb * rb;
  }
  double rra = ra * ra;
  double rrb = rb * rb;
  double s = rrb * acos((cc + rrb - rra) / (2.0 * cc * rb)) + rra * acos((cc + rra - rrb) / (2.0 * cc * ra)) - 0.5 * sqrt((-cc + rb + ra) * (cc + rb - ra) * (cc - rb + ra) * (cc + rb + ra));
  return s;
}

/*!
\brief Check if two circles intersect each other.

\param circle The other circle.
*/
bool Circle2::Intersect(const Circle2& circle) const
{
  double cc = SquaredNorm(c - circle.c);
  double rr = Math::Sqr(r + circle.r);
  if (cc > rr)
  {
    return false;
  }
  else
  {
    return true;
  }
}

/*!
\brief Check the intersection between a circle and a ray.

Note that intersections are sorted.

This function assumes that the ray is normalized, i.e. has a unit direction vector.

\param ray The (normalized) ray.
\param ta, tb Intersection depths.
*/
bool Circle2::Intersect(const Ray2& ray, double& ta, double& tb) const
{
  // Center
  Vector2 n = c - ray.Origin();

  // Distance to center
  double k = n * n;
  double t = n * ray.Direction();
  if (!(k < r * r) && (t < Epsilon))
    return false;

  double h = r * r - k + t * t;

  if (h < Epsilon)
    return false;

  // Intersection depths
  h = sqrt(h);
  ta = t - h;
  tb = t + h;

  return true;
}

/*!
\brief Check the intersection between a circle and a ray.

\param ray The (normalized) ray.
*/
bool Circle2::Intersect(const Ray2& ray) const
{
  // Center
  Vector2 n = c - ray.Origin();

  // Distance to center
  double k = n * n;
  double t = n * ray.Direction();

  if (!(k < r * r) && (t < Epsilon))
    return false;

  double h = r * r - k + t * t;

  if (h < Epsilon)
    return false;

  return true;
}

/*!
\brief Check the intersection between a circle and a segment.

\param segment The segment.
*/
bool Circle2::Intersect(const Segment2& segment) const
{
  Ray2 ray(segment.Vertex(0), segment.GetAxis());

  double ta, tb;
  if (!Intersect(ray, ta, tb))
  {
    return false;
  }

  if ((tb < 0.0) || (ta > segment.Length()))
  {
    return false;
  }

  return true;
}


/*!
\brief Compute the distance vector between a ciecle.
\param p Point.
*/
Vector2 Circle2::Normal(const Vector2& p) const
{
  Vector2 n = p - c;

  // Distance to center
  double y = sqrt(n * n);

  // Distance to circle
  y = sqrt(y) - r;

  n *= (y - r) / y;
  return n;
}

/*!
\brief Compute the squared distance between a point and the circle.

\param p Point.
*/
double Circle2::R(const Vector2& p) const
{
  Vector2 n = p - c;

  // Distance to center
  double y = n * n;

  // Distance to circle
  y = sqrt(y) - r;

  return y * y;
}

/*
\brief Creates a circle enclosing a set of points.

The algorithm computes a fast approximation of the
bounding ball for a point set based on the algorithm
given by [Jack Ritter, 1990].

\param p Array of vertices.
\param n Number of vertices.
*/
Circle2::Circle2(const QVector<Vector2>& s)
{
  Vector2 a = s.at(0);
  Vector2 b = a;

  int Pxmin = 0;
  int Pxmax = 0;
  int Pymin = 0;
  int Pymax = 0;

  // Compute largest extent along major axes
  for (int i = 1; i < s.size(); i++)
  {
    const Vector2& q = s.at(i);
    if (q[0] < a[0])
    {
      a[0] = q[0];
      Pxmin = i;
    }
    else if (q[0] > b[0])
    {
      b[0] = q[0];
      Pxmax = i;
    }
    if (q[1] < a[1])
    {
      a[1] = q[1];
      Pymin = i;
    }
    else if (q[1] > b[1])
    {
      b[1] = q[1];
      Pymax = i;
    }
  }

  // Select the largest extent as an initial diameter for the ball
  Vector2 x = s.at(Pxmax) - s.at(Pxmin);
  Vector2 y = s.at(Pymax) - s.at(Pymin);

  // Compute center according to the largest extent
  if (x * x >= y * y)
  {
    c = (s.at(Pxmin) + s.at(Pxmax)) * 0.5;
    r = x * x * 0.25;
  }
  else
  {
    c = (s.at(Pymin) + s.at(Pymax)) * 0.5;
    r = y * y * 0.25;
  }

  // Check that all points p[i] are in the circle
  // If not, expand the circle just enough to include them
  for (int j = 0; j < s.size(); j++)
  {
    Vector2 cp = s.at(j) - c;
    double ncp = cp * cp;
    // p[i] is outside the circle, so expand to include it
    if (ncp > r)
    {
      // Compute radii
      ncp = sqrt(ncp);
      r = sqrt(r);

      // Shift center towards p[i]
      c += (0.5 * (ncp - r) / ncp) * cp;

      // Enlarge radius just enough
      r = 0.5 * (r + ncp);

      // Restore squared radius
      r *= r;
    }
  }

  // Radius
  r = sqrt(r);
}

/*!
\brief Overloaded.
\param s Stream.
\param circle The circle.
*/
std::ostream& operator<<(std::ostream& s, const Circle2& circle)
{
  s << "Circle2(" << circle.c << ',' << circle.r << ')';
  return s;
}

/*!
\brief Compute a Poisson sphere distribution inside a circle.

This function uses a simple O(n<SUP>3</SUP>) dart throwing algorithm.

\sa Box2::Poisson

\param ra Radius of the sphere.
\param n Number of candidate points.
\param border True to sample the border
\param random %Random number generator.
*/
QVector<Vector2> Circle2::Poisson(const double& ra, int n, bool border, Random& random) const
{
  QVector<Vector2> p;

  if (border)
  {
    // Put points on border
    int b = int(Math::TwoPi * r / (2.0 * ra));
    for (int i = 0; i < b; i++)
    {
      p.append(Vertex(Math::Angle(i, b)));
    }
  }

  // Collision radius
  double r24 = 4.0 * ra * ra;

  // Create instances
  for (int i = 0; i < n; i++)
  {
    Vector2 t = RandomInside(random);
    bool hit = false;
    for (int j = 0; j < p.size(); j++)
    {
      if (SquaredNorm(t - p.at(j)) < r24)
      {
        hit = true;
        break;
      }
    }
    if (hit == false)
      p.append(t);
  }
  return p;
}


/*!
\brief Compute the quadratic BÃ©zier approximation of a circle arc.

Approximation yields best results for 0 &le; t &le; \pi;/4.
\param t Angle.
*/
QuadricCurve2 Circle2::QuadricBezierArc(const double& t) const
{
  Vector2 a = c + Vector2(r, 0.0);
  Vector2 b = c + r * Vector2(cos(t), sin(t));
  Vector2 m = c + Vector2(r, r * (1.0 - cos(t)) / sin(t));
  return QuadricCurve2::Bezier(a, m, b);
}

/*!
\brief Compute the cubic BÃ©zier approximation of a circle arc.

\param t Angle.
*/
CubicCurve2 Circle2::CubicBezierArc(const double& t) const
{
  const double k = (4.0 / 3.0) * tan(t / 4.0);
  Vector2 b[4];
  b[0] = c + r * Vector2(1.0, 0.0);
  b[1] = c + r * Vector2(1.0, k);

  b[3] = c + r * Vector2(cos(t), sin(t));

  b[2] = b[3] + r * k * Vector2(sin(t), -cos(t));
  return CubicCurve2::Bezier(b[0], b[1], b[2], b[3]);
}
