// Axis

#include "libs/axis.h"
#include "libs/ray.h"
#include "libs/sphere.h"

/*!
\class Axis axis.h
\brief A core axis class.

The class includes the two end vertices Axis::a and Axis::b, the normalized direction Axis::axis
and the length Axis::length. It is used in several geometric classes such as Cylinder, Capsule,
Cone and deformations such as Taper or Twist.

\ingroup KernelGroup
*/

const Axis Axis::X(Vector::Null, Vector::X);
const Axis Axis::Y(Vector::Null, Vector::Y);
const Axis Axis::Z(Vector::Null, Vector::Z);

/*!
\brief Creates an axis given end vertices.
\param a, b End vertices of the axis.
*/
Axis::Axis(const Vector& a, const Vector& b) :a(a), b(b)
{
  axis = b - a;
  length = Norm(axis);
  axis /= length;
}

/*!
\brief Compute the polynomial equation of the Euclidean distance between a ray and the axis.
\param ray The ray.
*/
Quadric Axis::Equation(const Ray& ray) const
{
  Vector pa = a - ray.Origin();

  double dx = ray.Direction() * axis;
  double pax = pa * axis;
  double dpa = ray.Direction() * pa;

  return Quadric(ray.Direction() * ray.Direction() - dx * dx, 2.0 * (dx * pax - dpa), pa * pa - pax * pax);
}

/*!
\brief Rotates an axis.

Only the vertices and the axis vector are modified, whereas the length is perserved.

\param r Rotation matrix.
*/
void Axis::Rotate(const Matrix& r)
{
  a = r * a;
  b = r * b;
  axis = r * axis;
}

/*!
\brief Translates an axis.

\param t Translation vector.
*/
void Axis::Translate(const Vector& t)
{
  a += t;
  b += t;
}

/*!
\brief Uniformly scales an axis.

\param s Scaling factor.
*/
void Axis::Scale(const double& s)
{
  a *= s;
  b *= s;
  length *= s;
}

/*!
\brief Scales an axis.

\param s Scaling vector.
*/
void Axis::Scale(const Vector& s)
{
  a *= s;
  b *= s;
  axis = b - a;
  length = Norm(axis);
  axis /= length;
}

/*!
\brief Overloaded.
\param s Stream.
\param axis The axis.
*/
std::ostream& operator<<(std::ostream& s, const Axis& axis)
{
  s << "Axis (" << axis.a << ',' << axis.b << ")";
  return s;
}

/*!
\brief Compute the normal vector between a point and
its projection onto the edge.
\param p Point.
*/
Vector Axis::Normal(const Vector& p) const
{
  Vector n = p - b;
  // Not that vertex
  if (axis * n <= 0.0)
  {
    n = p - a;
    double s = axis * n;
    // Not that one either
    if (s >= 0.0)
    {
      n -= s * axis;
    }
  }
  return n;
}

/*!
\brief Compute the squared distance to the segment.

This function is more efficient than Segment::R(const Vector&) const
as Axis stores the unit vector.
\param p Point.
*/
double Axis::R(const Vector& p) const
{
  Vector n = p - a;
  double s = n * axis;

  // Squared distance 
  double r = n * n;

  // Not vertex a
  if (s > 0.0)
  {
    r -= s * s;

    // Vertex b
    if (s > length)
    {
      s -= length;
      r += s * s;
    }
  }
  return r;
}

/*!
\brief Compute the squared distance to the axis.
\param p Point.
\param s Signed distance of the projection of argument point p onto the edge.
*/
double Axis::R(const Vector& p, double& s) const
{
  Vector n = p - a;
  s = n * axis;

  // Squared distance 
  double r = n * n;

  // Not vertex a
  if (s > 0.0)
  {
    r -= s * s;

    // Vertex b
    if (s > length)
    {
      s -= length;
      r += s * s;
      s = 1.0;
    }
    else
    {
      s /= length;
    }
  }
  // Vertex a
  else
  {
    s = 0.0;
  }
  return r;
}

/*!
\brief Compute an orthonormal frame attached to the axis.

\sa Axis::GetFrame(const Vector&)
*/
Matrix Axis::GetFrame() const
{
  Vector x = Normalized(axis.Orthogonal());
  Vector y = axis / x;

  return Matrix(x, y, axis);
}

/*!
\brief Compute an orthonormal frame attached to the axis.

The axes of the frame are defined as the column vectors of the returned matrix.

\code
Matrix m = Axis::GetFrame(Vector(1.0,2.0,-1.0));
for (int i=0;i<16;i++)
{
Vector p=m*Vector(cos(i/16.0),sin(i/16.0),0.0);
}
\endcode

\param z %Vector, which should be unit.

\sa Axis::GetFrame()
*/
Matrix Axis::GetFrame(const Vector& z)
{
  // Orthonormal vector
  Vector x, y;
  z.Orthonormal(x, y);

  return Matrix(x, y, z);
}

/*!
\brief Compute the squared distance between two axes.

\param axis %The other axis.
*/
double Axis::R(const Axis& axis) const
{
  Vector u = b - a;
  Vector v = axis.b - axis.a;
  Vector w = a - axis.a;
  double uu = Math::Sqr(length);
  double uv = u * v;
  double vv = Math::Sqr(axis.length);
  double uw = u * w;
  double vw = v * w;
  double e = uu * vv - uv * uv;

  double sN, sD = e;       // sc = sN / sD, default sD = D >= 0
  double tN, tD = e;       // tc = tN / tD, default tD = D >= 0

  // Compute the line parameters of the two closest points
  // Lines are almost parallel
  if (e < 1e-6)
  {
    // force using point P0 on segment S1 to prevent possible division by 0.0 later
    sN = 0.0;
    sD = 1.0;
    tN = vw;
    tD = vv;
  }
  // Get the closest points on the infinite lines
  else
  {
    sN = uv * vw - vv * uw;
    tN = uu * vw - uv * uw;
    if (sN < 0.0)
    {        // sc < 0 => the s=0 edge is visible
      sN = 0.0;
      tN = vw;
      tD = vv;
    }
    else if (sN > sD)
    {  // sc > 1  => the s=1 edge is visible
      sN = sD;
      tN = vw + uv;
      tD = vv;
    }
  }

  if (tN < 0.0)
  {            // tc < 0 => the t=0 edge is visible
    tN = 0.0;
    // recompute sc for this edge
    if (-uw < 0.0)
    {
      sN = 0.0;
    }
    else if (-uw > uu)
    {
      sN = sD;
    }
    else
    {
      sN = -uw;
      sD = uu;
    }
  }
  else if (tN > tD)
  {      // tc > 1  => the t=1 edge is visible
    tN = tD;
    // recompute sc for this edge
    if ((-uw + uv) < 0.0)
    {
      sN = 0;
    }
    else if ((-uw + uv) > uu)
    {
      sN = sD;
    }
    else
    {
      sN = (-uw + uv);
      sD = uu;
    }
  }

  double s = (fabs(sN) < 1e-6 ? 0.0 : sN / sD);
  double t = (fabs(tN) < 1e-6 ? 0.0 : tN / tD);

  // Difference of the two closest points
  Vector st = w + (s * u) - (t * v);

  return st * st;
}

/*!
\brief Symmetric point.
\param p The point.
*/
Vector Axis::Symmetric(const Vector& p) const
{
  Vector ap = p - a;
  double s = ap * axis;
  Vector n = ap - s * axis;
  return p - 2 * n;
}

/*!
\brief Symmetric sphere.
\param sphere The sphere.
*/
Sphere Axis::Symmetric(const Sphere& sphere) const
{
  return Sphere(Symmetric(sphere.Center()), sphere.Radius());
}
