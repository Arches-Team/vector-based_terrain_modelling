// Octahedron 

#pragma once

#include "libs/box.h"
#include "libs/triangle.h"

class Plane;

class Octahedron {
protected:
  Vector center = Vector::Null; //!< Center.
  Vector half; //!< Half diagonal of the box.
  Vector planenormal; //!< Normal to the plane in the first octant.
  Vector exyn, eyzn, ezxn; //!< Edge normals.
public:
  explicit Octahedron(const Vector&, const Vector&);
  explicit Octahedron(const Vector&, const double&);
  explicit Octahedron(const double&);

  double Volume() const;
  double Area() const;

  Vector Vertex(int) const;
  Vector Vertex(int, int) const;
  Vector Normal(int) const;

  Vector Normal(const Vector&) const;
  double R(const Vector&) const;
  double Signed(const Vector&) const;

  Plane GetPlane(int) const;

  Box GetBox() const;
  bool Inside(const Vector&) const;
public:
  static int VertexIndex(int, int);
protected:
  static const int face[8][3]; //!< Array of triangles.
  static const Vector vertex[6]; //!< Array of internal vertexes for computing vertices.
};

/*!
\brief Compute the volume of the octahedron.

It is 8 times the volume of a tetrahedron.
Since three vectors are orthogonal, the volume
is 1/6 the volume of a half box (thus 1/48 the volume of the box).
*/
inline double Octahedron::Volume() const
{
  return half[0] * half[1] * half[2] / 6.0;
}

/*!
\brief Compute the area of the octahedron.

It is 8 times the area of one of its triangular faces.

The code has been optimized, instead of computing the
normal of the triangle using Triangle::Normal(), the
cross product has been inlined.
*/
inline double Octahedron::Area() const
{
  return 8.0 * Norm(Vector(half[1] * half[2], half[0] * half[2], half[0] * half[1]));
}

//! Return the i-th vertex as a vector. 
inline Vector Octahedron::Vertex(int i) const
{
  return center + half.Scaled(vertex[i]);
}

/*!
\brief Compute the normal to the i-th face as a normalized vector.
*/
inline Vector Octahedron::Normal(int i) const
{
  return Triangle(Vertex(i, 0), Vertex(i, 1), Vertex(i, 2)).Normal();
}

//! Return the index of the j-th vertex of the i-th face. 
inline int Octahedron::VertexIndex(int i, int j)
{
  return face[i][j];
}

//! Return the j-th vertex of the i-th face. 
inline Vector Octahedron::Vertex(int i, int j) const
{
  return Vertex(face[i][j]);
}

//! Return the bounding box of the octahedron.
inline Box Octahedron::GetBox() const
{
  return Box(center - half, center + half);
}
