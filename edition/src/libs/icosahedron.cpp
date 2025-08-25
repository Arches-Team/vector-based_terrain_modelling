// Icosahedron 

#include "libs/icosahedron.h"
#include "libs/plane.h"

/*!
\defgroup ExtendedKernelGroup Extended geometric classes

\brief Extended classes group several more specific classes such as Tetrahedra, Octahedron, Dodecahedron, Icosahedron.
*/

/*!
\class Icosahedron icosahedron.h

\brief An icosahedron.

The following Cartesian coordinates define the vertices
of an icosahedron with edge-length 2, centered at the origin.

(0,\sign 1,\sign\phi)

(\sign 1,\sign\phi,0)

(\sign\phi,0,\sign 1)

where \phi=(1+sqrt(5))/2 is the golden ratio.

\ingroup ExtendedKernelGroup
*/
const double Icosahedron::Phi = (1 + sqrt(5.0)) / 2.0;

//! Radius of the circumscribed sphere of an icosahedron with edge length 2.
const double Icosahedron::R2 = (1.0 / 4.0) * sqrt(10 + 2.0 * sqrt(5.0)) * 2.0;

// Array of vertices
const Vector Icosahedron::vertex[12] = {
  Vector(0.0, Phi, 1), Vector(0.0, Phi, -1), Vector(0.0, -Phi, 1), Vector(0.0, -Phi, -1),
  Vector(Phi, 1, 0.0), Vector(-Phi, 1, 0.0), Vector(Phi, -1, 0.0), Vector(-Phi, -1, 0.0),
  Vector(-1, 0.0, Phi), Vector(1, 0.0, Phi), Vector(-1, 0.0, -Phi), Vector(1, 0.0, -Phi)
};

// Array of normals
const Vector Icosahedron::normal[20] = {
  Vector(1.0 / Phi, Phi, 0) / sqrt(3.0), Vector(-1.0 / Phi, Phi, 0) / sqrt(3.0), Vector(1, 1, 1) / sqrt(3.0), Vector(-1, 1, 1) / sqrt(3.0),
  Vector(0, 1.0 / Phi, Phi) / sqrt(3.0), Vector(1, 1,-1) / sqrt(3.0), Vector(-1, 1,-1) / sqrt(3.0), Vector(0, 1.0 / Phi,-Phi) / sqrt(3.0),
  Vector(1.0 / Phi,-Phi, 0) / sqrt(3.0), Vector(-1.0 / Phi,-Phi, 0) / sqrt(3.0), Vector(1,-1, 1) / sqrt(3.0), Vector(-1,-1,1) / sqrt(3.0),
  Vector(0,-1.0 / Phi,Phi) / sqrt(3.0), Vector(1,-1,-1) / sqrt(3.0), Vector(-1,-1,-1) / sqrt(3.0), Vector(0,-1.0 / Phi,-Phi) / sqrt(3.0),
  Vector(Phi, 0, 1.0 / Phi) / sqrt(3.0), Vector(Phi, 0,-1.0 / Phi) / sqrt(3.0), Vector(-Phi, 0, 1.0 / Phi) / sqrt(3.0), Vector(-Phi,0,-1.0 / Phi) / sqrt(3.0)
};

// Array of vertex indexes
const int Icosahedron::face[20][3] = {
  { 1,0,4 },{ 0,1,5 },{ 4,0,9 },{ 0,5,8 },{ 0,8,9 },
  { 1,4,11 },{ 5,1,10 },{ 10,1,11 },{ 2,3,6 },{ 3,2,7 },
  { 2,6,9 },{ 7,2,8 },{ 8,2,9 },{ 6,3,11 },{ 3,7,10 },
  { 3,10,11 },{ 6,4,9 },{ 4,6,11 },{ 5,7,8 },{ 7,5,10 }
};

// Array of edges indexes
const int Icosahedron::edge[30][2] = {
  { 1, 0}, { 0, 4}, { 4, 1},
  { 1, 5}, { 5, 0},
  { 0, 9}, { 9, 4},
  { 5, 8}, { 8, 0},
  { 8, 9},
  { 4,11}, {11, 1},
  { 1,10}, {10, 5},
  {11,10},
  { 2, 3}, { 3, 6}, { 6, 2},
  { 2, 7}, { 7, 3},
  { 6, 9}, { 9, 2},
  { 2, 8}, { 8, 7},
  { 3,11}, {11, 6},
  { 7,10}, {10, 3},
  { 4, 6},
  { 5, 7}
};

/*!
\brief Creates an icosahedron.
\param r Radius.
*/
Icosahedron::Icosahedron(const double& r) :Icosahedron(Vector::Null, r)
{
}

/*!
\brief Creates an icosahedron given its center and radius.
\param c Center.
\param r Radius.
*/
Icosahedron::Icosahedron(const Vector& c, const double& r)
{
  Icosahedron::c = c;
  Icosahedron::r = r;
}

/*!
\brief Compute the volume of the icosahedron.
*/
double Icosahedron::Volume() const
{
  double u = r / R2;
  return (5.0 / 12.0) * (3.0 + sqrt(5.0)) * u * u * u;
}

/*!
\brief Compute the area of the icosahedron.
*/
double Icosahedron::Area() const
{
  double u = r / R2;
  return 5.0 * sqrt(3.0) * u * u;
}


/*!
\brief Compute the intersection between an icosahedron and a ray.

Intersect the ray with the planes corresponding to the faces of the convex icosahedron.

\param ray The ray.
\param ta, tb Intersection depths.
\param fa, fb Index of the intersected faces.
*/
bool Icosahedron::Intersect(const Ray& ray, double& ta, double& tb, int& fa, int& fb) const
{
  // Intersection depths with planes
  ta = -Math::Infinity;
  tb = Math::Infinity;

  // Indexes
  fa = -1;
  fb = -1;

  for (int i = 0; i < 20; i++)
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
\brief Compute the approximate signed distance to the icosahedron.
\param p Point.
*/
double Icosahedron::Signed(const Vector& p) const
{
  double s = GetPlane(0).Signed(p);
  for (int i = 1; i < 20; i++)
  {
    s = Math::Max(s, GetPlane(i).Signed(p));
  }
  return s;
}

/*!
\brief Check if a point is inside.

\param p The point.
*/
bool Icosahedron::Inside(const Vector& p) const
{
  for (int i = 0; i < 20; i++)
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
*/
Plane Icosahedron::GetPlane(int k) const
{
  return Plane(normal[k], c + r * vertex[face[k][0]]);
}

/*!
\brief Generate a random vector inside the icosahedron.

Uses a simple rejection method.

\param random %Random number generator.
*/
Vector Icosahedron::RandomInside(Random& random) const
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

/*!
\brief Overloaded output-stream operator.
\param s Stream.
\param icosahedron The icosahedron.
*/
std::ostream& operator<<(std::ostream& s, const Icosahedron& icosahedron)
{
  s << "Icosahedron(" << icosahedron.c << "," << icosahedron.r << ")";

  return s;
}