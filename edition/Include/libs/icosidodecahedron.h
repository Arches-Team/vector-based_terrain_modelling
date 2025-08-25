// Icosidodecahedron 

#pragma once

#include "libs/box.h"

class Plane;

class Icosidodecahedron {
protected:
  Vector c = Vector::Null; //!< Center of icosi-dodecahedron.
  double r = 1.0; //!< Radius.
public:
  explicit Icosidodecahedron(const Vector&, const double&);
  explicit Icosidodecahedron(const double&);
  double Volume() const;
  double Area() const;
  double Radius() const;

  Vector Vertex(int) const;
  Vector Vertex(int, int) const;
  Vector Normal(int) const;

  QVector<Vector> Vertexes() const;

  Box GetBox() const;
  Plane GetPlane(int) const;

  bool Intersect(const Ray&, double&, double&, int&, int&) const;
  bool Inside(const Vector&) const;
  double Signed(const Vector&) const;

  friend std::ostream& operator<<(std::ostream&, const Icosidodecahedron&);

public:
  static int VertexTriangleIndex(int, int);
  static int VertexPentagonIndex(int, int);
protected:
  static const double R2; //!< Radius of the circumscribed sphere of an icosidodecahedron with edge length 2.0/Phi.
  static const double Phi; //!< Golden ratio.
  static const Vector vertex[30]; //!< Array of vertices.
  static const Vector normal[32]; //!< Array of normals.
  static const int pentagon[12][5]; //!< Array of pentagons.
  static const int triangle[20][3]; //!< Array of triangles.
};

/*!
\brief Return the radius of the circmscribed sphere of the icosi-dodecahedron.
*/
inline double Icosidodecahedron::Radius() const
{
  return r;
}

/*!
\brief Return the i-th vertex.
\param i Index.
*/
inline Vector Icosidodecahedron::Vertex(int i) const
{
  return c + r * vertex[i];
}

/*!
\brief Return the index of the j-th vertex of the i-th triangle face.
\param i Face (pentagon).
\param j Index representing the j-th vertex of the face.
*/
inline int Icosidodecahedron::VertexTriangleIndex(int i, int j)
{
  return triangle[i][j];
}

/*!
\brief Return the index of the j-th vertex of the i-th pentagon face.
\param i Face (pentagon).
\param j Index representing the j-th vertex of the face.
*/
inline int Icosidodecahedron::VertexPentagonIndex(int i, int j)
{
  return pentagon[i][j];
}

/*!
\brief Return the normal to the i-th face as a normalized vector.
\param i Index.
*/
inline Vector Icosidodecahedron::Normal(int i) const
{
  return normal[i];
}

/*!
\brief Return the j-th vertex of the i-th face.
\param i Face number, starting with the 12 pentagons, and following with the 20 triangles.
\param j Index representing the j-th vertex of the face.
*/
inline Vector Icosidodecahedron::Vertex(int i, int j) const
{
  if (i < 12)
  {
    return Vertex(pentagon[i][j]);
  }
  else
  {
    return Vertex(triangle[i][j]);
  }
}

//! Return the bounding box of the shape.
inline Box Icosidodecahedron::GetBox() const
{
  return Box(c, r * Phi);
}
