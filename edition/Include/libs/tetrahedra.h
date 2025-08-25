// Tetrahedras

#pragma once

#include "libs/triangle.h"
#include "libs/box.h"

class Tetrahedra
{
protected:
  Vector p[4]; //!< Array of vertices.
  Vector n[4]; //!< Normals to the faces.
public:
  //! Empty
  Tetrahedra() {}
  explicit Tetrahedra(const Vector&, const Vector&, const Vector&, const Vector&);
  explicit Tetrahedra(const Box&, int);

  //! Empty.
  ~Tetrahedra() {}

  // Access to vertices, normals, triangles
  Vector Vertex(int) const;
  Triangle GetTriangle(int) const;
  Plane GetPlane(int) const;
  Vector Normal(int) const;

  constexpr int VertexIndex(int, int) const;

  Vector Normal(const Vector&) const;
  Vector Normal2(const Vector&) const;
  double R(const Vector&) const;
  double Signed(const Vector&) const;

  double Volume() const;
  double Area() const;

  Box GetBox() const;
  Sphere GetSphere() const;

  Tetrahedra Scaled(const Vector&) const;
  Tetrahedra Scaled(const double&) const;
  Tetrahedra Translated(const Vector&) const;
  Tetrahedra Rotated(const Matrix&) const;

  bool Intersect(const Ray&, double&, double&) const;
  bool Intersect(const Ray&, double&, double&, Vector&, Vector&) const;
  bool Inside(const Vector&) const;

  bool Intersect(const Tetrahedra&) const;
  bool Inside(const Plane&) const;

  friend std::ostream& operator<<(std::ostream&, const Tetrahedra&);
  Vector RandomInside(Random & = Random::R239) const;
  QVector<Vector> Poisson(const double&, int, Random & = Random::R239) const;
public:
  static constexpr double epsilon = 1e-7; //!< Internal epsilon constant.
protected:
  static constexpr int vertex[4][3] = {
  { 0, 1, 2 },
  { 0, 3, 1 },
  { 0, 2, 3 },
  { 1, 3, 2 }
  }; //!< Array of indexes defining the faces of the tetrahedron.
  static constexpr int split[6][4] = {
  { 0, 1, 5, 7 },
  { 0, 1, 3, 7 },
  { 0, 3, 7, 2 },
  { 1, 6, 7, 2 },
  { 1, 6, 4, 7 },
  { 1, 4, 7, 5 }
  }; //!< Indexes of the vertexes of the sub-tetrahedron of a cube.
};


/*!
\brief Returns the i-thvertex of the tetrahedron.
\param i Index.
*/
inline Vector Tetrahedra::Vertex(int i) const
{
  return p[i];
}

/*!
\brief Returns the i-th triangle of the tetrahedron.
\param i Index.
*/
inline Triangle Tetrahedra::GetTriangle(int i) const
{
  return Triangle(Vertex(vertex[i][0]), Vertex(vertex[i][1]), Vertex(vertex[i][2]));
}

/*!
\brief Returns the normal of the i-th triangle of the tetrahedron.
\param i Index.
*/
inline Vector Tetrahedra::Normal(int i) const
{
  return n[i];
}

/*!
\brief Return the index of the j-th vertex of the i-th triangle.
\param i %Triangle index.
\param j Vertex index.
*/
inline constexpr int Tetrahedra::VertexIndex(int i, int j) const
{
  return vertex[i][j];
}
