// Dodecahedron 

// Self include
#include "libs/dodecahedron.h"
#include "libs/plane.h"

/*!
\class Dodecahedron dodecahedron.h
\brief A dodecahedron.

The following Cartesian coordinates define the vertices
of a dodecahedron with edge-length 2/\phi, centered at the origin.

(\sign 1, \sign 1, \sign 1)

(0, \sign 1/\phi, \sign\phi)

(\sign 1/\phi, \sign\phi, 0)

(\sign\phi, 0, \sign 1/\phi)

where \phi=(1+sqrt(5))/2 is the golden ratio.

\ingroup ExtendedKernelGroup
*/
const double Dodecahedron::Phi = (1 + sqrt(5.0)) / 2.0;

const Vector Dodecahedron::vertex[20] = {
  Vector(0, 1.0 / Phi, Phi), Vector(0, -1 / Phi, Phi), Vector(0, -1 / Phi, -Phi), Vector(0, 1 / Phi, -Phi),
  Vector(Phi, 0, 1 / Phi), Vector(-Phi, 0, 1 / Phi), Vector(-Phi, 0, -1 / Phi), Vector(Phi, 0, -1 / Phi),
  Vector(1.0 / Phi, Phi, 0), Vector(-1.0 / Phi, Phi, 0), Vector(-1.0 / Phi, -Phi, 0), Vector(1.0 / Phi, -Phi, 0),
  Vector(1, 1, 1), Vector(-1, 1, 1), Vector(-1, -1, 1), Vector(1, -1, 1),
  Vector(1, -1, -1), Vector(1, 1, -1), Vector(-1, 1, -1), Vector(-1, -1, -1)
};

const Vector Dodecahedron::normal[12] = {
  Vector(1.0 / Phi, 0.0, 1.0) / sqrt(1.0 + 1.0 / (Phi * Phi)),  // Normal n0
  Vector(0, 1.0, 1.0 / Phi) / sqrt(1.0 + 1.0 / (Phi * Phi)),  // n1
  Vector(-1.0 / Phi, 0, 1.0) / sqrt(1.0 + 1.0 / (Phi * Phi)), // -n3
  Vector(0.0, -1.0, 1.0 / Phi) / sqrt(1.0 + 1.0 / (Phi * Phi)), // n4
  Vector(1.0 / Phi, 0.0, -1.0) / sqrt(1.0 + 1.0 / (Phi * Phi)), // n3
  Vector(0.0, -1.0, -1.0 / Phi) / sqrt(1.0 + 1.0 / (Phi * Phi)), // -n1
  Vector(-1.0 / Phi, 0.0, -1.0) / sqrt(1.0 + 1.0 / (Phi * Phi)), // -n0
  Vector(0.0, 1.0, -1.0 / Phi) / sqrt(1.0 + 1.0 / (Phi * Phi)),  // -n4
  Vector(1.0, -1.0 / Phi, 0.0) / sqrt(1.0 + 1.0 / (Phi * Phi)),  // n5
  Vector(1.0, 1.0 / Phi, 0.0) / sqrt(1.0 + 1.0 / (Phi * Phi)),   // n2
  Vector(-1.0, 1.0 / Phi, 0.0) / sqrt(1.0 + 1.0 / (Phi * Phi)),  // -n5
  Vector(-1.0, -1.0 / Phi, 0.0) / sqrt(1.0 + 1.0 / (Phi * Phi))  // -n2
};

const int Dodecahedron::face[12][5] = {
  { 0, 1, 15, 4, 12 }, { 0, 12, 8, 9, 13 }, { 0, 13, 5, 14, 1 }, { 1, 14, 10, 11, 15 }, { 2, 3, 17, 7, 16 },
  { 2, 16, 11, 10, 19 }, { 2, 19, 6, 18, 3 }, { 18, 9, 8, 17, 3 }, { 15, 11, 16, 7, 4 }, { 4, 7, 17, 8, 12 },
  { 13, 9, 18, 6, 5 }, { 5, 6, 19, 10, 14 }
};

const double Dodecahedron::R2 = sqrt(3.0);

/*!
\brief Create a dodecahedron.
\param r Radius.
*/
Dodecahedron::Dodecahedron(const double& r) :Dodecahedron(Vector::Null, r)
{
}

/*!
\brief Create a dodecahedron.
\param c Center.
\param r Radius.
*/
Dodecahedron::Dodecahedron(const Vector& c, const double& r)
{
  Dodecahedron::c = c;
  Dodecahedron::r = r;
}

/*!
\brief Compute the volume of the dodecahedron.
*/
double Dodecahedron::Volume() const
{
  double u = r / R2;
  return (1.0 / 4.0) * (15.0 + 7.0 * sqrt(5.0)) * u * u * u;
}

/*!
\brief Compute the area of the dodecahedron.
*/
double Dodecahedron::Area() const
{
  double u = r / R2;
  return 3.0 * sqrt(25.0 + 10.0 * sqrt(5.0)) * u * u;
}


/*!
\brief Compute the intersection between a dodecahedron and a ray.

Intersect the ray with the planes corresponding to the faces of the convex dodecahedron.

\param ray The ray.
\param ta, tb Intersection depths.
\param fa, fb Index of the intersected faces.
*/
bool Dodecahedron::Intersect(const Ray& ray, double& ta, double& tb, int& fa, int& fb) const
{
  // Intersection depths with planes
  ta = -Math::Infinity;
  tb = Math::Infinity;

  // Indexes
  fa = -1;
  fb = -1;

  for (int i = 0; i < 12; i++)
  {
    double t;
    if (Plane(Normal(i), Vertex(i, 0)).Intersect(ray, t))
    {
      if ((Normal(i) * ray.Direction()) < 0.0)
      {
        if (ta > t)
        {
          t = ta;
          fa = i;
        }
      }
      else
      {
        if (tb < t)
        {
          t = tb;
          fb = i;
        }
      }
    }
    if (ta > tb) return false;
  }
  return true;
}

/*!
\brief Compute the approximate signed distance to the dodecahedron.
\param p Point.
*/
double Dodecahedron::Signed(const Vector& p) const
{
  double s = GetPlane(0).Signed(p);
  for (int i = 1; i < 12; i++)
  {
    s = Math::Max(s, GetPlane(i).Signed(p));
  }
  return s;
}

/*!
\brief Check if a point is inside.

\param p The point.
*/
bool Dodecahedron::Inside(const Vector& p) const
{
  for (int i = 0; i < 12; i++)
  {
    if (!Plane(Normal(i), Vertex(i, 0)).Inside(p))
    {
      return false;
    }
  }
  return true;
}

/*!
\brief Return the plane of the k-th face.
\param k Index.
*/
Plane Dodecahedron::GetPlane(int k) const
{
  return Plane(normal[k], c + r * vertex[face[k][0]]);
}

/*!
\brief Generate a random vector inside the dodecahedron.

Uses a simple rejection method.

\param random %Random number generator.
*/
Vector Dodecahedron::RandomInside(Random& random) const
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
