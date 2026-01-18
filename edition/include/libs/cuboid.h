// Cuboid
#pragma once


#include "libs/box.h"
#include "libs/quadrangle.h"

class Cuboid
{
protected:
  Vector a[8]; //!< Vertexes.
public:
  explicit Cuboid(const Box&);
  explicit Cuboid(const Vector&, const Vector&, const Vector&, const Vector&, const Vector&, const Vector&, const Vector&, const Vector&);

  //! Empty.
  ~Cuboid() {}

  // Access vertexes
  Vector& operator[] (int);
  Vector operator[] (int) const;

  Vector Vertex(int) const;
  Vector Vertex(int, int) const;
  Vector Normal(int) const;
  
  Quadrangle GetQuadrangle(int) const;

  void Translate(const Vector&);
  void Rotate(const Matrix&);
  void Scale(const Vector&);

  double Area() const;
public:
  static const int edge[24]; //!< Edge vertices.
  static const int quadrangle[24]; //!< Quadrangle faces.
};

//! Returns either end vertex of the box.
inline Vector& Cuboid::operator[] (int i)
{
  return a[i];
}

//! Overloaded.
inline Vector Cuboid::operator[] (int i) const
{
  return a[i];
}

/*!
\brief Returns the i-th vertex of the cuboid.

Same as operator[].

\param i Vertex.
*/
inline Vector Cuboid::Vertex(int i) const
{
  return a[i];
}

/*!
\brief Returns the j-th vertex of i-th quadrangle face of the cuboid.
\param i %Quadrangle.
\param j Vertex.
*/
inline Vector Cuboid::Vertex(int i, int j) const
{
  return a[quadrangle[i * 4 + j]];
}

/*!
\brief Returns the normal to the i-th quadrangle face of the cuboid.
\param i %Quadrangle index.
*/
inline Vector Cuboid::Normal(int i) const
{
  return Normalized((Vertex(i, 1) - Vertex(i, 0)) / (Vertex(i, 2) - Vertex(i, 0)));
}

/*!
\brief Returns tthe i-th quadrangle.
\param i %Quadrangle index.
*/
inline Quadrangle Cuboid::GetQuadrangle(int i) const
{
  return Quadrangle(a[quadrangle[i * 4]], a[quadrangle[i * 4 + 1]], a[quadrangle[i * 4 + 2]], a[quadrangle[i * 4 + 3]]);
}

