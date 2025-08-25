// Hemisphere

// Self include
#include "libs/hemisphere.h"

#include "libs/sphere.h"

/*!
\class HemiSphere hemisphere.h
\brief %Hemispheres.

\ingroup ExtendedKernelGroup

*/

/*!
\brief Create a hemi-sphere.
\param c Center.
\param r Radius.
\param a %Axis, should be normalized.
*/
HemiSphere::HemiSphere(const Vector& c, const Vector& a, const double& r) :c(c), r(r), axis(a)
{
}

/*!
\brief Create a hemi-sphere centered at the origin and with vertical direction.
\param r Radius.
*/
HemiSphere::HemiSphere(const double& r) :HemiSphere(Vector::Null, Vector::Z, r)
{
}

/*!
\brief Compute the distance vector between the hemisphere and a point.
\param p Point.
*/
Vector HemiSphere::Normal(const Vector& p) const
{
  Vector n = p - c;

  // Side
  double s = n * axis;

  // Sphere side
  if (s > 0.0)
  {
    double t = n * n;

    if (t < r * r)
      return Vector::Null;

    t = sqrt(t);
    n *= (1.0 - r / t);

    return n;
  }

  // Otherwise, use disc side
  Vector z = s * axis;
  Vector x = n - z;

  double y = x * x;

  // Inside
  if (y < r * r)
  {
    return z;
  }
  // Outside
  else
  {
    return z + (1.0 - r / sqrt(y)) * x;
  }
}

/*!
\brief Compute the distance between the hemisphere
and a point.
\param p Point.
*/
double HemiSphere::R(const Vector& p) const
{
  Vector n = p - c;

  // Side
  double s = n * axis;

  // Sphere side
  if (s > 0.0)
  {
    double t = n * n;

    if (t < r * r)
      return 0.0;

    t = sqrt(t) - r;

    return t * t;
  }

  // Otherwise, use disc side
  Vector z = s * axis;
  Vector x = n - z;

  double y = x * x;

  // Inside
  if (y < r * r)
  {
    return z * z;
  }
  // Outside
  else
  {
    y = sqrt(y) - r;
    return z * z + y * y;
  }
}

/*!
\brief Compute the Euclidean signed distance between the hemisphere and a point.
\param p Point.
*/
double HemiSphere::Signed(const Vector& p) const
{
  Vector n = p - c;

  // Side, which is also the distance to the plane
  double s = n * axis;

  // Sphere side
  if (s > 0.0)
  {
    double t = n * n;
    t = sqrt(t) - r;

    // Inside the hemisphere
    if (t < 0.0)
    {
      // Distance to sphere
      return Math::Max(-s, t);
    }

    return t;
  }

  // Otherwise, use disc side
  Vector z = s * axis;
  Vector x = n - z;

  double y = x * x;

  // Inside
  if (y < r * r)
  {
    return -s;
  }
  // Outside
  else
  {
    y = sqrt(y) - r;
    return sqrt(z * z + y * y);
  }
}

/*!
\brief Compute the bouding box.
*/
Box HemiSphere::GetBox() const
{
  return Box(c, r);
}

/*!
\brief Generate a random direction.
\param axis %Axis.
\param random %Random number generator.
*/
Vector HemiSphere::RandomDirection(const Vector& axis, Random& random)
{
  Vector d = Sphere::RandomNormal(random);

  if (axis * d < 0.0)
  {
    d = -d;
  }
  return d;
}

/*!
\brief Compute the i-th Fibonacci direction on the vertical unit hemisphere.
\param i Point index.
\param n Number of samples on the hemisphere.
*/
Vector HemiSphere::FibonacciUnit(int i, int n)
{
  double phiI = Math::TwoPi * ((i / Math::Golden) - floor(i / Math::Golden));
  double cosThetaI = 1.0 - ((2.0 * i + 1.0) / (2.0 * n));
  double sinThetaI = sqrt(1.0 - (cosThetaI * cosThetaI));
  return Vector(cos(phiI) * sinThetaI, sin(phiI) * sinThetaI, cosThetaI);
}

/*!
\brief Compute the i-th Fibonacci direction on the hemisphere.
\param i Point index.
\param n Number of samples on the hemisphere.
*/
Vector HemiSphere::Fibonacci(int i, int n)
{
  // Get Fibonacci direction on unit hemisphere
  Vector d = FibonacciUnit(i, n);

  // Transform into world space direction using an orthonormal basis
  Vector b1, b2;
  axis.Orthonormal(b1, b2);
  Vector dFib = d[0] * b1 + d[1] * b2 + d[2] * axis;

  // Project on surface
  return c + r * dFib;
}

/*!
\brief Overloaded output-stream operator.
\param s Stream.
\param hemisphere The hemisphere.
*/
std::ostream& operator<<(std::ostream& s, const HemiSphere& hemisphere)
{
  s << "HemiSphere(" << hemisphere.c << "," << hemisphere.axis << "," << hemisphere.r << ")";

  return s;
}