// Icosahedron 
#pragma once


#include "libs/box.h"

class Icosahedron {
protected:
  Vector c = Vector::Null; //!< Center.
  double r = 1.0; //!< Radius.
public:
  explicit Icosahedron(const Vector&, const double&);
  explicit Icosahedron(const double& = 1.0);

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

  friend std::ostream& operator<<(std::ostream&, const Icosahedron&);

public:
  static int VertexIndex(int, int); // removed the constexpr for gcc compatibility
  static int VertexEdgeIndex(int, int);// removed the constexpr for gcc compatibility

protected:
  static const double R2; //!< Radius of the circumscribed sphere of an icosahedron with edge length 2.
  static const double Phi; //!< Phi.
  static const Vector vertex[12]; //!< Array of vertices.
  static const Vector normal[20]; //!< Array of normals.
  static const int edge[30][2]; //!< Array of edges.
  static const int face[20][3]; //!< Array of pentagons.
};

/*!
\brief Return the radius of the circmscribed sphere of the icosi-dodecahedron.
*/
inline double Icosahedron::Radius() const
{
  return r;
}

//! Return the i-th vertex as a vector. 
inline Vector Icosahedron::Vertex(int i) const
{
  return c + r * vertex[i];
}

//! Return the normal to the i-th face as a normalized vector. 
inline Vector Icosahedron::Normal(int i) const
{
  return normal[i];
}

//! Return the index of the j-th vertex of the i-th face. 
inline int Icosahedron::VertexIndex(int i, int j)
{
  return face[i][j];
}

//! Return the j-th vertex of the i-th face. 
inline Vector Icosahedron::Vertex(int i, int j) const
{
  return Vertex(face[i][j]);
}

/*!
\brief Return the index of the j-th vertex of the i-th edge.
\param i Edge.
\param j Vertex (should be O or 1).
*/
inline int Icosahedron::VertexEdgeIndex(int i, int j)
{
  return edge[i][j];
}

//! Return the bounding box.
inline Box Icosahedron::GetBox() const
{
  return Box(c, r * Phi);
}