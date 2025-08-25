// Circle

#include "libs/circle.h"
#include "libs/polygon.h"

/*!
\class Circle circle.h
\brief Circles.

Circles are characterized by their center, axis and radius.
They are sometimes interpreted as discs, as for ray intersection.

\sa Disc

\ingroup KernelGroup
*/

/*!
\brief Creates a circle.
\param r Radius.
*/
Circle::Circle(const double& r) :Circle(Vector::Null, Vector::Z, r)
{
}

/*!
\brief Creates a circle.
\param c Center.
\param axis %Axis, which should be normalized.
\param r Radius.
*/
Circle::Circle(const Vector& c, const Vector& axis, const double& r) :c(c), axis(axis), r(r)
{
}

/*!
\brief Compute the axis aligned bounding box of a circle.

Some math will demonstrate that the box can be computed easily
directly from the axis of the circle. Computations involve three
square roots. For a larger but even faster bounding box, use
the bounding box of the corresponding sphere:
\code
Circle c;
Box b(c.Center(),c.Radius());
\endcode
*/
Box Circle::GetBox() const
{
  Vector s = Axis::BoxVector(axis);

  return Box(c - r * s, c + r * s);
}

/*!
\brief Compute the ray-disc intersection.
\param ray The ray.
\param t Intersection depth.
*/
bool Circle::Intersect(const Ray& ray, double& t) const
{
  // Normal ray-plane intersection
  double e = axis * ray.Direction();
  if (fabs(e) < epsilon)
    return false;

  t = ((c - ray.Origin()) * axis) / e;
  if (t < 0.0)
    return false;

  Vector p = ray(t) - c;
  if (p * p > r * r)
    return false;

  return true;
}

/*!
\brief Rotates a circle.

\param r Rotation matrix.
*/
void Circle::Rotate(const Matrix& r)
{
  Vector d = c + axis;
  c = r * c;
  d = r * d;
  axis = d - c;
}

/*!
\brief Translates a circle.

\param t Translation vector.
*/
void Circle::Translate(const Vector& t)
{
  c += t;
}

/*!
\brief Uniformly scales a circle.

\param s Scaling factor.
*/
void Circle::Scale(const double& s)
{
  c *= s;

  r *= fabs(s);
}

/*!
\brief Compute the distance between a point and a circle.

Cost: 10+, 8*, 1sqrt

\param p Point.
*/
double Circle::R(const Vector& p) const
{
  Vector n = p - c;
  double z = n * axis;
  z *= z;

  // Distance to center
  double d = n * n;

  // Squared radial distance
  double y = d - z;

  // Distance to circle
  y = r - sqrt(y);

  return y * y + z;
}

/*!
\brief Computes the distance vector between a circle and a point.

This function projects the point
in the plane of the circle and computes the distance to the circle
in its plane, before adding the extra term along the normal.
\param p Point.
*/
Vector Circle::Normal(const Vector& p) const
{
  Vector n = p - c;
  Vector z = (n * axis) * axis;
  Vector x = n - z;

  double y = x * x;

  // Circle
  return z + (1.0 - r / sqrt(y)) * x;
}

/*!
\brief Creates a circle given three vertices.

Basically set the center to the center of the triangle,
and compute the radius.

\param x, y, z Argument vertices
*/
Circle::Circle(const Vector& x, const Vector& y, const Vector& z)
{
  Vector a = x - z;
  Vector b = y - z;

  double na = a * a;
  double nb = b * b;
  double ab = a * b;

  double d = na * nb - ab * ab;

  if (fabs(d) > 1e-06)
  {
    double e = 0.5 / d;
    double u = e * nb * (na - ab);
    double v = e * na * (nb - ab);

    Circle::c = u * x + v * y + (1.0 - u - v) * z;
    Circle::r = Norm(c - x);
    Circle::axis = Normalized(a / b);
  }
  else
  {
    Circle::c = Vector::Null;
    Circle::r = 0.0;
    Circle::axis = Vector::Null;
  }
}

/*!
\brief Overloaded.
\param s Stream.
\param circle The circle.
*/
std::ostream& operator<<(std::ostream& s, const Circle& circle)
{
  s << "Circle(" << circle.c << ',' << circle.axis << ',' << circle.r << ')';
  return s;
}

/*!
\brief Generate a random vector inside the circle.
\param random %Random number generator.
*/
Vector Circle::RandomInside(Random& random) const
{
  // Random point inside unit circle
  Vector2 uv = Circle2::RandomUnit(random);

  Vector x, y;
  axis.Orthonormal(x, y);

  return c + r * (uv[0] * x + uv[1] * y);
}

/*!
\brief Generate a spiraling vector in the disc.
\param i Integer, should be <n.
\param n Number of points.
*/
Vector Circle::Vogel(int i, int n) const
{
  // Random point inside unit circle
  Vector2 uv = Circle2::Unit.Vogel(i, n);

  Vector x, y;
  axis.Orthonormal(x, y);

  return c + r * (uv[0] * x + uv[1] * y);
}

/*!
\brief Generate a random vector on the circle.
\param random %Random number generator.
*/
Vector Circle::RandomOn(Random& random) const
{
  Vector2 uv = Circle2::RandomUnit(random);

  // Random point inside unit circle
  double s = SquaredNorm(uv);

  // Random point on the circle
  double rx = (uv[0] * uv[0] - uv[1] * uv[1]) / s;
  double ry = 2.0 * uv[0] * uv[1] / s;

  Vector x, y;
  axis.Orthonormal(x, y);

  return c + r * (rx * x + ry * y);
}
