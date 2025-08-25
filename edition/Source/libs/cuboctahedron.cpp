// Cuboctahedron 

#include "libs/cuboctahedron.h"

/*!
\class Cuboctahedron cuboctahedron.h

\brief A cuboctahedron is a polyhedron with 8 triangular faces and 6 square faces.

A cuboctahedron has 12 identical vertices, with 2 triangles and 2 squares meeting at each,
and 24 identical edges, each separating a triangle from a square.

\image html cuboctahedron.png

\ingroup ExtendedKernelGroup
*/

// Array of vertices
const Vector Cuboctahedron::vertex[12] = {
  Vector(-1.0,-1.0, 0.0), Vector(1.0,-1.0, 0.0), Vector(-1.0, 1.0, 0.0), Vector(1.0, 1.0, 0.0),
  Vector(-1.0, 0.0,-1.0), Vector(1.0, 0.0,-1.0), Vector(-1.0, 0.0, 1.0), Vector(1.0, 0.0, 1.0),
  Vector(0.0,-1.0,-1.0), Vector(0.0, 1.0,-1.0), Vector(0.0,-1.0, 1.0), Vector(0.0, 1.0, 1.0)
};

const Vector Cuboctahedron::normal[14] = {
  Vector(-1.0, -1.0, -1.0) / Math::Sqrt3,
  Vector(1.0,  1.0, -1.0) / Math::Sqrt3,
  Vector(-1.0, -1.0, -1.0) / Math::Sqrt3,
  Vector(1.0,  1.0, -1.0) / Math::Sqrt3,
  Vector(-1.0, -1.0, 1.0) / Math::Sqrt3,
  Vector(1.0,  1.0, 1.0) / Math::Sqrt3,
  Vector(-1.0, -1.0, 1.0) / Math::Sqrt3,
  Vector(1.0,  1.0, 1.0) / Math::Sqrt3,
  Vector(-1.0, 0.0, 0.0),
  Vector(1.0, 0.0, 0.0),
  Vector(0.0, -1.0, 0.0),
  Vector(0.0, 1.0, 0.0),
  Vector(0.0, 0.0,-1.0),
  Vector(0.0, 0.0, 1.0)
};

/*!
\brief Creates a cuboctahedron.
\param c Center.
\param r Radius.
*/
Cuboctahedron::Cuboctahedron(const Vector& c, const double& r) :c(c), r(r)
{
}

/*!
\brief Creates a cuboctahedron.
\param r Radius.
*/
Cuboctahedron::Cuboctahedron(const double& r) :Cuboctahedron(Vector::Null, r)
{
}

/*!
\brief Compute the volume of the cuboctahedron.
*/
double Cuboctahedron::Volume() const
{
  return (5.0 / 3.0) * sqrt(2.0) * Math::Cube(EdgeLength());
}

/*!
\brief Compute the area of the cuboctahedron.
*/
double Cuboctahedron::Area() const
{
  return (6.0 + 2.0 * sqrt(3.0)) * Math::Sqr(EdgeLength());
}

/*!
\brief Check if a point is inside a cuboctahedron.
\param p Point.
*/
bool Cuboctahedron::Inside(const Vector& p) const
{
  Vector q = p - c;

  // Triple planar symetry
  q = Abs(q);

  // Cube
  if ((q[0] > r) || (q[1] > r) || (q[2] > r))
  {
    return false;
  }
  // Octahedron
  if (q[0] + q[1] + q[2] > 2.0 * r)
  {
    return false;
  }

  return true;
}

/*!
\brief Compute the signed distance approximation.
\param p Point.
*/
double Cuboctahedron::Signed(const Vector& p) const
{
  Vector q = p - c;

  // Triple planar symetry
  q = Abs(q);

  // Cube
  double s = Math::Max(q[0] - r, q[1] - r, q[2] - r);

  // Intersection with octahedron
  s = Math::Max(s, (q[0] + q[1] + q[2] - 2.0 * r) / sqrt(3.0));

  return s;
}

/*!
\brief Compute the squared distance approximation.
\param p Point.
*/
double Cuboctahedron::R(const Vector& p) const
{
  double s = Signed(p);
  return (s < 0.0) ? 0.0 : s * s;
}