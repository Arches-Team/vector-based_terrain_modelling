// Octahedron 

#include "libs/octahedron.h"
#include "libs/plane.h"

/*!
\class Octahedron octahedron.h

\brief An axis aligned octahedron.

\ingroup ExtendedKernelGroup
*/

// Array of vertices
const Vector Octahedron::vertex[6] = {
  Vector(0.0, 0.0, -1.0), Vector(1.0, 0.0, 0.0), Vector(0.0, 1.0, 0.0), Vector(-1.0, 0.0, 0.0), Vector(0.0, -1.0, 0.0), Vector(0.0, 0.0, 1.0)
};

// Array of vertex indexes
const int Octahedron::face[8][3] = {
  { 0, 2, 1 }, { 0, 3, 2 }, { 0, 4, 3 }, { 0, 1, 4 },
  { 5, 1, 2 }, { 5, 2, 3 }, { 5, 3, 4 }, { 5, 4, 1 }
};

/*!
\brief Creates an octahedron.
\param r Radius.
*/
Octahedron::Octahedron(const double& r) :Octahedron(Vector::Null, r)
{
}

/*!
\brief Creates an octahedron given its center and radius.
\param c Center.
\param r Radius.
*/
Octahedron::Octahedron(const Vector& c, const double& r) :center(c), half(r)
{
  planenormal = Normalized(Vector(half[1] * half[2], half[0] * half[2], half[0] * half[1]));

  Vector exy(-half[0], half[1], 0.0);
  Vector eyz(0.0, -half[1], half[2]);
  Vector ezx(half[0], 0.0, -half[2]);

  exyn = exy / planenormal;
  eyzn = eyz / planenormal;
  ezxn = ezx / planenormal;
}

/*!
\brief Creates an octahedron given the vertices of its bounding box.
\param a,b Vectors.
*/
Octahedron::Octahedron(const Vector& a, const Vector& b) :center(0.5 * (a + b)), half(0.5 * (b - a))
{
  planenormal = Normalized(Vector(half[1] * half[2], half[0] * half[2], half[0] * half[1]));

  Vector exy(-half[0], half[1], 0.0);
  Vector eyz(0.0, -half[1], half[2]);
  Vector ezx(half[0], 0.0, -half[2]);

  exyn = exy / planenormal;
  eyzn = eyz / planenormal;
  ezxn = ezx / planenormal;
}

/*!
\brief Check if a point is inside the octahedron.
\param p Point.
*/
bool Octahedron::Inside(const Vector& p) const
{
  Vector q = p - center;

  // Take advantage of symmetry
  q = Abs(q);

  Vector qa = q - Vector(half[0], 0.0, 0.0);

  if (qa * planenormal > 0.0)
  {
    return false;
  }

  return true;
}

/*!
\brief Return the plane of the k-th face.
\param k Index.
*/
Plane Octahedron::GetPlane(int k) const
{
  return Plane(Normal(k), center + half.Scaled(vertex[face[k][0]]));
}

/*!
\brief Compute the normal vector between a point and the octahedron.
\param p Point.
*/
Vector Octahedron::Normal(const Vector& p) const
{
  Vector q = p - center;

  // Take advantage of symmetry
  Vector symmetry = Vector(q[0] < 0.0 ? -1.0 : 1.0, q[1] < 0.0 ? -1.0 : 1.0, q[2] < 0.0 ? -1.0 : 1.0);
  q = Abs(q);

  // Position against one point
  Vector qa = q - Vector(half[0], 0.0, 0.0);

  // Inside
  if (qa * planenormal < 0.0)
  {
    return Vector::Null;
  }

  // Other points
  Vector qb = q - Vector(0.0, half[1], 0.0);
  Vector qc = q - Vector(0.0, 0.0, half[2]);

  Vector exy(-half[0], half[1], 0.0);
  Vector eyz(0.0, -half[1], half[2]);
  Vector ezx(half[0], 0.0, -half[2]);

  double qaezx = (qa * ezx);
  double qaexy = (qa * exy);
  double qbexy = (qb * exy);
  double qbeyz = (qb * eyz);
  double qceyz = (qc * eyz);
  double qcezx = (qc * ezx);

  // Booleans for vertexes
  bool taezx = qaezx >= 0.0;
  bool taexy = qaexy <= 0.0;
  bool tbexy = qbexy >= 0.0;
  bool tbeyz = qbeyz <= 0.0;
  bool tceyz = qceyz >= 0.0;
  bool tcezx = qcezx <= 0.0;

  // Vertexes
  if (taezx && taexy)
  {
    return qa.Scaled(symmetry);
  }
  if (tbexy && tbeyz)
  {
    return qb.Scaled(symmetry);
  }
  if (tceyz && tcezx)
  {
    return qc.Scaled(symmetry);
  }

  // Booleans for edges
  bool texyn = (qa * exyn) < 0.0;
  bool teyzn = (qb * eyzn) < 0.0;
  bool tezxn = (qc * ezxn) < 0.0;

  // Plane
  if (texyn && teyzn && tezxn)
  {
    return ((qa * planenormal) * planenormal).Scaled(symmetry);
  }
  // Edges
  if (!texyn && !taexy && !tbexy)
  {
    double u = qaexy / (exy * exy);
    return (qa - u * exy).Scaled(symmetry);
  }
  if (!teyzn && !tbeyz && !tceyz)
  {
    double u = qbeyz / (eyz * eyz);
    return (qb - u * eyz).Scaled(symmetry);
  }
  if (!tezxn && !tcezx && !taezx)
  {
    double u = qcezx / (ezx * ezx);
    return (qc - u * ezx).Scaled(symmetry);
  }
  // Should never happen
  return Vector::Null;
}

/*!
\brief Compute the squared Euclidean distance between a point and the octahedron.
\param p The point.
*/
double Octahedron::R(const Vector& p) const
{
  Vector q = p - center;

  // Take advantage of symmetry
  q = Abs(q);

  // Position against one point
  Vector qa = q - Vector(half[0], 0.0, 0.0);

  // Inside
  if (qa * planenormal < 0.0)
  {
    return 0.0;
  }

  // Other points
  Vector qb = q - Vector(0.0, half[1], 0.0);
  Vector qc = q - Vector(0.0, 0.0, half[2]);

  Vector exy(-half[0], half[1], 0.0);
  Vector eyz(0.0, -half[1], half[2]);
  Vector ezx(half[0], 0.0, -half[2]);

  double qaezx = (qa * ezx);
  double qaexy = (qa * exy);
  double qbexy = (qb * exy);
  double qbeyz = (qb * eyz);
  double qceyz = (qc * eyz);
  double qcezx = (qc * ezx);

  // Booleans for vertexes
  bool taezx = qaezx >= 0.0;
  bool taexy = qaexy <= 0.0;
  bool tbexy = qbexy >= 0.0;
  bool tbeyz = qbeyz <= 0.0;
  bool tceyz = qceyz >= 0.0;
  bool tcezx = qcezx <= 0.0;

  // Vertexes
  if (taezx && taexy)
  {
    return qa * qa;
  }
  if (tbexy && tbeyz)
  {
    return qb * qb;
  }
  if (tceyz && tcezx)
  {
    return qc * qc;
  }

  // Booleans for edges
  bool texyn = (qa * exyn) < 0.0;
  bool teyzn = (qb * eyzn) < 0.0;
  bool tezxn = (qc * ezxn) < 0.0;

  // Plane
  if (texyn && teyzn && tezxn)
  {
    return Math::Sqr(qa * planenormal);
  }
  // Edges
  if (!texyn && !taexy && !tbexy)
  {
    double u = qaexy / (exy * exy);
    return SquaredNorm(qa - u * exy);
  }
  if (!teyzn && !tbeyz && !tceyz)
  {
    double u = qbeyz / (eyz * eyz);
    return SquaredNorm(qb - u * eyz);
  }
  if (!tezxn && !tcezx && !taezx)
  {
    double u = qcezx / (ezx * ezx);
    return SquaredNorm(qc - u * ezx);
  }
  // Should never happen
  return 0.0;
}

/*!
\brief Compute the signed Euclidean distance between a point and the octahedron.
\param p Point.
*/
double Octahedron::Signed(const Vector& p) const
{
  Vector q = p - center;

  // Take advantage of symmetry
  q = Abs(q);

  // Position against one point
  Vector qa = q - Vector(half[0], 0.0, 0.0);

  // Inside
  double d = qa * planenormal;
  if (d < 0.0)
  {
    return d;
  }

  // Other points
  Vector qb = q - Vector(0.0, half[1], 0.0);
  Vector qc = q - Vector(0.0, 0.0, half[2]);

  Vector exy(-half[0], half[1], 0.0);
  Vector eyz(0.0, -half[1], half[2]);
  Vector ezx(half[0], 0.0, -half[2]);

  double qaezx = (qa * ezx);
  double qaexy = (qa * exy);
  double qbexy = (qb * exy);
  double qbeyz = (qb * eyz);
  double qceyz = (qc * eyz);
  double qcezx = (qc * ezx);

  // Booleans for vertexes
  bool taezx = qaezx >= 0.0;
  bool taexy = qaexy <= 0.0;
  bool tbexy = qbexy >= 0.0;
  bool tbeyz = qbeyz <= 0.0;
  bool tceyz = qceyz >= 0.0;
  bool tcezx = qcezx <= 0.0;

  // Vertexes
  if (taezx && taexy)
  {
    return Norm(qa);
  }
  if (tbexy && tbeyz)
  {
    return Norm(qb);
  }
  if (tceyz && tcezx)
  {
    return Norm(qc);
  }

  // Booleans for edges
  bool texyn = (qa * exyn) < 0.0;
  bool teyzn = (qb * eyzn) < 0.0;
  bool tezxn = (qc * ezxn) < 0.0;

  // Plane
  if (texyn && teyzn && tezxn)
  {
    return (qa * planenormal);
  }

  // Edges
  if (!texyn && !taexy && !tbexy)
  {
    double u = qaexy / (exy * exy);
    return Norm(qa - u * exy);
  }
  if (!teyzn && !tbeyz && !tceyz)
  {
    double u = qbeyz / (eyz * eyz);
    return Norm(qb - u * eyz);
  }
  if (!tezxn && !tcezx && !taezx)
  {
    double u = qcezx / (ezx * ezx);
    return Norm(qc - u * ezx);
  }

  // Should never happen
  return 0.0;
}
