// Dodecahedron 

#pragma once

#include "libs/box.h"

class Plane;

class Dodecahedron {
protected:
  Vector c = Vector::Null; //!< Center of dodecahedron.
  double r = 1.0; //!< Radius.
public:
  explicit Dodecahedron(const double& = 1.0);
  explicit Dodecahedron(const Vector&, const double&);
  double Volume() const;
  double Area() const;
  double Radius() const;

  Vector Vertex(int) const;
  Vector Vertex(int, int) const;
  Vector Normal(int) const;

  Box GetBox() const;

  Plane GetPlane(int) const;

  bool Inside(const Vector&) const;
  Vector RandomInside(Random & = Random::R239) const;

  double Signed(const Vector&) const;

  bool Intersect(const Ray&, double&, double&, int&, int&) const;
public:
  static int VertexIndex(int, int); // removed constexpr for gcc compatibility
protected:
  static const double R2; //!< Radius of the circumscribed sphere of an dodecahedron with edge length 2.0/Phi.
  static const double Phi; //!< Golden ratio.
  static const Vector vertex[20]; //!< Array of vertices.
  static const Vector normal[12]; //!< Array of normals.
  static const int face[12][5]; //!< Array of pentagons.
};

//! Compute the radius of the circmscribed sphere of the dodecahedron.
inline double Dodecahedron::Radius() const
{
  return r;
}

/*!
\brief Return the i-th vertex.
\param i Index.
*/
inline Vector Dodecahedron::Vertex(int i) const
{
  return c + r * vertex[i];
}

/*!
\brief Return the index of the j-th vertex of the i-th face.
\param i Face (pentagon).
\param j Index representing the j-th vertex of the face.
*/
inline int Dodecahedron::VertexIndex(int i, int j) // for gcc compatibility constexpr has been removed
{
  return face[i][j];
}

/*!
\brief Return the normal to the i-th face as a normalized vector.
\param i Index.
*/
inline Vector Dodecahedron::Normal(int i) const
{
  return normal[i];
}

/*!
\brief Return the j-th vertex of the i-th face.
\param i Face (pentagon).
\param j Index representing the j-th vertex of the face.
*/
inline Vector Dodecahedron::Vertex(int i, int j) const
{
  return Vertex(face[i][j]);
}

//! Return the bounding box of the shape.
inline Box Dodecahedron::GetBox() const
{
  return Box(c, r * Phi);
}

