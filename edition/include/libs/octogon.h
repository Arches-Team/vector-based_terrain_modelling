// Octogon

#pragma once

#include "libs/box.h"

class Octogon2 {
protected:
  Vector2 c = Vector::Null; //!< Center.
  double r = 1.0; //!< Radius.
public:
  //! \brief Empty.
  Octogon2() {}
  explicit Octogon2(const Vector2&, const double&);
  explicit Octogon2(const double&);

  void Translate(const Vector2&);
  void Scale(const double&);

  bool Inside(const Vector2&) const;
  double R(const Vector2&) const;
  double Signed(const Vector2&) const;
  Vector2 Normal(const Vector2&) const;

  Vector2 Center() const;
  double Radius() const;
  double InscribedRadius() const;

  Vector2 Radial(int) const;
  Vector2 Edge(int) const;
  Vector2 Vertex(int) const;
  Box2 GetBox() const;

  double Area() const;
  double Perimeter() const;
  double Edge() const;

  Vector2 RandomInside(Random&) const;
  QVector<Vector2> Poisson(const double&, int, const QVector<Vector2>&, bool = true, Random & = Random::R239) const;
  QVector<Vector2> Poisson(const double&, int, Random & = Random::R239) const;

  friend std::ostream& operator<<(std::ostream&, const Octogon2&);
  void Draw(QGraphicsScene&, const QPen & = QPen(), const QBrush & = QBrush()) const;

protected:
  static const Vector2 vertex[8]; //!< Array of vertexes.
  static const Vector2 normal[8]; //!< Array of normal vectors to the edges.
  static const Vector2 edge[8]; //!< Unit edge vectors.
};

/*!
\brief Return the coordinates of the k-th vertex.
\param k Index, should be in [0,7].
*/
inline Vector2 Octogon2::Vertex(int k) const
{
  return c + vertex[k] * r;
}

/*!
\brief Return the radial vector of the k-th vertex.
\param k Index, should be in [0,7].
*/
inline Vector2 Octogon2::Radial(int k) const
{
  return vertex[k] * r;
}

/*!
\brief Return the edge vector connecting vertexes k and k+1.
\param k Index, should be in [0,7].
*/
inline Vector2 Octogon2::Edge(int k) const
{
  return edge[k] * r;
}

/*!
\brief Return the center of the octogon.
*/
inline Vector2 Octogon2::Center() const
{
  return c;
}

/*!
\brief Return the radius of the octogon.
*/
inline double Octogon2::Radius() const
{
  return r;
}

/*!
\brief Edge length.
*/
inline double Octogon2::Edge() const
{
  return sqrt(2.0 - sqrt(2.0)) * r;
}

/*!
\brief Inscribed radius.
*/
inline double Octogon2::InscribedRadius() const
{
  return ((1 + sqrt(2.0)) / 2.0) * Edge();
}

class IrregularOctogon2
{
protected:
  Vector2 c = Vector::Null; //!< Center.
  double s = 1.0; //!< Side length.
  double t = 1.0; //!< Side length of diagonal edges.
public:
  //! \brief Empty.
  IrregularOctogon2() {}
  explicit IrregularOctogon2(const double&, const double&);
  explicit IrregularOctogon2(const Vector2&, const double&, const double&);

  void Translate(const Vector2&);
  void Scale(const double&);

  bool Inside(const Vector2&) const;
  double Signed(const Vector2&) const;
  double R(const Vector2&) const;

  Vector2 Center() const;
  double Radius() const;

  Vector2 Vertex(int) const;
  Box2 GetBox() const;

  double Area() const;
  double Perimeter() const;

  friend std::ostream& operator<<(std::ostream&, const IrregularOctogon2&);
  void Draw(QGraphicsScene&, const QPen & = QPen(), const QBrush & = QBrush()) const;

  Vector2 RandomInside(Random&) const;
  QVector<Vector2> Poisson(const double&, int, const QVector<Vector2>&, bool = true, Random & = Random::R239) const;
  QVector<Vector2> Poisson(const double&, int, Random & = Random::R239) const;

protected:
};

/*!
\brief Return the coordinates of the k-th vertex.
\param k Index, should be in [0,7].
*/
inline Vector2 IrregularOctogon2::Vertex(int k) const
{
  Vector2 v(t + s, s);

  // Other direction
  if (((k + 1) % 8) & 2)
  {
    v = Vector2(s, t + s); 
  }

  // Symmetry 
  if (((k + 2) % 8) >> 2)
  { 
    v[0] = -v[0];
  }

  // Symmetry
  if (k >> 2) { v[1] = -v[1]; }

  return c + v;
}


/*!
\brief Return the center of the octogon.
*/
inline Vector2 IrregularOctogon2::Center() const
{
  return c;
}

/*!
\brief Return the radius of the octogon.
*/
inline double IrregularOctogon2::Radius() const
{
  return (s + t) * Math::Sqrt2;
}