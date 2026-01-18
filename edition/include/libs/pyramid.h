// Pyramid

#pragma once

#include "libs/box.h"

class Pyramid {
protected:
  Vector center = Vector::Null; //!< Center.
  double height = 1.0; //!< Height.
  double a = 1.0; //!< Half base-side size.

  Vector nt; //!< Normal to triangular face.
  Vector neaa; //!< Unit diagonal vector of the triangular faces.
public:
  explicit Pyramid(const Vector&, const double&);
  explicit Pyramid(const Vector&, const double&, const double&);
  explicit Pyramid(const double&, const double&);

  double Volume() const;
  double Area() const;

  Vector Vertex(int) const;
  Vector Vertex(int, int) const;
  Vector Normal(int) const;

  Vector Normal(const Vector&) const;
  double R(const Vector&) const;
  double Signed(const Vector&) const;

  Box GetBox() const;
  bool Inside(const Vector&) const;

  friend std::ostream& operator<<(std::ostream&, const Pyramid&);

public:
  static int VertexIndex(int, int); // removed constexpr for gcc compatibility
protected:
  static const int face[6][3]; //!< Array of triangles.
  static const Vector vertex[5]; //!< Array of vertices.
};

//! Return the i-th vertex as a vector. 
inline Vector Pyramid::Vertex(int i) const
{
  return center + vertex[i].Scaled(Vector(a, a, height));
}

//! Return the index of the j-th vertex of the i-th face. 
inline int Pyramid::VertexIndex(int i, int j)
{
  return face[i][j];
}

//! Return the j-th vertex of the i-th face. 
inline Vector Pyramid::Vertex(int i, int j) const
{
  return Vertex(face[i][j]);
}

//! Return the bounding box of the pyramid.
inline Box Pyramid::GetBox() const
{
  return Box(center - Vector(a, a, 0.0), center + Vector(a, a, height));
}