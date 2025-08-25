// Torus

#include "libs/torus.h"
#include "libs/quartic.h"
#include "libs/cylinder.h"

/*!
\class Torus torus.h
\brief A simple torus.

The structure does not store the squared radius of the torus so as to save memory.

Rotation and translation are inherited from Circle.

\ingroup KernelGroup
*/

const Torus Torus::Unit(Vector::Null, Vector::Z, 1.0, 1.0);

/*!
\brief Creates a torus.
\param c Center.
\param axis %Axis, which should be normalized.
\param r Major radius.
\param s Minor (small) radius.
*/
Torus::Torus(const Vector& c, const Vector& axis, const double& r, const double& s) :Circle(c, axis, r), s(s)
{
}

/*!
\brief Create a torus.
\param c Circle.
\param s Minor (small) radius.
*/
Torus::Torus(const Circle& c, const double& s) :Circle(c), s(s)
{
}

/*!
\brief Create a horizontal torus centered at origin.
\param r,s Major and minor radii.
*/
Torus::Torus(const double& r, const double& s) :Torus(Vector::Null, Vector::Z, r, s)
{
}

/*!
\brief Compute the squared Euclidean distance between a point and the torus.

\param p Point.
*/
double Torus::R(const Vector& p) const
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

  y = sqrt(y * y + z);

  // Inside torus
  if (y < s)
    return 0.0;

  y -= s;

  return y * y;
}

/*!
\brief Compute the signed Euclidean distance between a point and the torus.

\param p Point.
*/
double Torus::Signed(const Vector& p) const
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

  y = sqrt(y * y + z);

  y -= s;

  return y;
}

/*!
\brief Computes the distance vector between a torus and a given point.

This function projects the point
in the plane of the circle and computes the distance to the circle
in its plane, before adding the extra term along the normal.

\param p Point.
*/
Vector Torus::Normal(const Vector& p) const
{
  Vector n = p - c;

  double y = n * axis;
  Vector z = y * axis;
  Vector x = n - z;

  double u = x * x;

  // Normal to circle
  n = z + (1.0 - r / sqrt(u)) * x;

  u = Norm(n);

  // Inside torus
  if (u < s)
    return Vector::Null;

  n *= (1.0 - s / u);
  return n;
}

//! May be used to sort roots
inline void Sort(double* a, int n)
{
  for (int i = 0; i < n; i++)
  {
    for (int j = i + 1; i < n; i++)
    {
      if (a[j] < a[i])
      {
        double t = a[i];
        a[i] = a[j];
        a[j] = t;
      }
    }
  }
}

/*!
\brief Computes the intersection depths and the normals
between a ray and a torus.

The quartic polynomial in t is derived from the
canonical algebraic expression. Roots are found using the
fast analytical Quartic::Solve() function of the Quartic class.
This function returns an integer which is the number of
intersections found. Note that intersections are sorted.

\param ray The ray.
\param depth Array of intersection depths.
\param normal Array of surface normals at intersection points.
*/
int Torus::Intersect(const Ray& ray, double* depth, Vector* normal) const
{
  Vector pa = c - ray.Origin();

  double dx = ray.Direction() * axis;
  double pax = pa * axis;
  double dpa = ray.Direction() * pa;

  // Squared distance to axis
  Quadric er(ray.Direction() * ray.Direction() - dx * dx, 2.0 * (dx * pax - dpa), pa * pa - pax * pax);

  // Squared distance along axis
  Quadric ea(dx * dx, -2.0 * dx * pax, pax * pax);

  // Create quartic equation
  er += ea;

  er[0] -= s * s + r * r;

  Quartic e(er[2] * er[2], 2.0 * er[2] * er[1], er[1] * er[1] + 2.0 * er[0] * er[2], 2.0 * er[1] * er[0], er[0] * er[0]);

  ea[0] -= s * s;

  ea *= r * r * 4.0;

  e[2] += ea[2];
  e[1] += ea[1];
  e[0] += ea[0];

  // Roots
  double t[4];

  int k = e.Solve(t);
  Sort(t, k);
  int i = 0;
  while (k--)
  {
    // Compute the intersection point out of the torus space
    Vector p = ray(t[k]);

    Vector n = p - c;
    Vector z = (n * axis) * axis;
    Vector x = n - z;

    double y = x * x;

    // Vector normal to circle
    n = z + (1.0 - r / sqrt(y)) * x;

    // Depth
    depth[i] = t[k];

    // Get normal from derivatives   
    normal[i] = Normalized(n);
    i++;
  }
  return i;
}

/*!
\brief Check if a point is inside or outside the torus.
\param p Point.
*/
bool Torus::Inside(const Vector& p) const
{
  Vector n = p - c;

  // Axial distance
  double axial = n * axis;
  axial *= axial;

  // Radial distance
  double radial = sqrt(n * n - axial);

  radial -= r;
  radial *= radial;

  if (axial + radial <= s * s)
    return true;

  return false;
}

/*!
\brief Computes the bounding box of a torus.
*/
Box Torus::GetBox() const
{
  return Circle::GetBox().Extended(s);
}

/*!
\brief Uniformly scales a torus.

\param s Scaling factor.
*/
void Torus::Scale(const double& s)
{
  Circle::Scale(s);

  Torus::s *= Math::Abs(s);
}

/*!
\brief Compute the volume of a torus.
*/
double Torus::Volume() const
{
  return 2.0 * Math::Pi * Math::Pi * r * s * s;
}

/*!
\brief Compute the surface area of a torus.
*/
double Torus::Area() const
{
  return 4.0 * Math::Pi * Math::Pi * r * s;
}

/*!
\brief Overloaded.

\param s Stream.
\param torus The torus.
*/
std::ostream& operator<<(std::ostream& s, const Torus& torus)
{
  s << "Torus(" << torus.c << ',' << torus.axis << ',' << torus.r << ',' << torus.s << ')';
  return s;
}

/*!
\brief Generate a random vector inside the dodecahedron.

Uses a simple rejection method.

\param random %Random number generator.
*/
Vector Torus::RandomInside(Random& random) const
{
  Box box = GetBox();

  while (true)
  {
    Vector p = box.RandomInside(random);
    if (Inside(p))
      return p;
  }

  // Should never happen
  return Vector::Null;
}
