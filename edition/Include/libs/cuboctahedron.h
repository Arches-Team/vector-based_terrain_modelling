// Cuboctahedron 

#pragma once

#include "libs/box.h"

class Cuboctahedron
{
protected:
  Vector c = Vector::Null; //!< Center.
  double r = 1.0; //!< Radius.
public:
  explicit Cuboctahedron(const Vector&, const double&);
  explicit Cuboctahedron(const double&);
  double Volume() const;
  double Area() const;

  Vector Vertex(int) const;
  Vector Vertex(int, int) const;

  Vector Normal(int) const;

  bool Inside(const Vector&) const;
  
  double Signed(const Vector&) const;
  double R(const Vector&) const;

  Box GetBox() const;
public:
  static constexpr int VertexIndex(int, int);
protected:
  static const Vector vertex[12]; //!< Array of vertices.
  static const Vector normal[14]; //!< Array of normals.

  // Array of vertex indexes
  static constexpr int face[20][3] = {
      // 8 Triangles
      { 0, 4, 8 }, { 5, 1, 8 }, { 1, 7, 10 }, { 0, 10, 6 },
      { 7, 11, 3 }, { 6, 11, 2}, { 4, 2, 9 }, { 9, 3, 5 },
      // 6 Squares
      { 0, 8, 1 },{ 0, 1, 10 },{ 1, 5, 3 },{ 1, 3, 7 },
      { 6, 10, 7 },{ 6, 7, 11 },{ 0, 2, 4 },{ 0, 6, 2 },
      { 4, 9, 5 },{ 4, 5, 8 },{ 2, 3, 9 },{ 2, 11, 3 }
  };


protected:
  double EdgeLength() const;
};

/*!
\brief Compute the length of an edge.
*/
inline double Cuboctahedron::EdgeLength() const
{
  return r * Math::Sqrt2;
}

//! Return the i-th vertex as a vector. 
inline Vector Cuboctahedron::Vertex(int i) const
{
  return c + r * vertex[i];
}

//! Return the index of the j-th vertex of the i-th face. 
inline constexpr int Cuboctahedron::VertexIndex(int i, int j)
{
  return face[i][j];
}

//! Return the j-th vertex of the i-th face. 
inline Vector Cuboctahedron::Vertex(int i, int j) const
{
  return Vertex(face[i][j]);
}

/*!
\brief Return the normal of the i-th face.

Note that the first 8 faces are  triangles, and the 6 last are squares.
\param i Index, should be 0 &#8804; i &#8804; 14.
*/

inline Vector Cuboctahedron::Normal(int i) const
{
  return normal[i];
}

/*!
\brief Return the bounding box of the cuboctahedron.
*/
inline Box Cuboctahedron::GetBox() const
{
  return Box(c, r);
}
