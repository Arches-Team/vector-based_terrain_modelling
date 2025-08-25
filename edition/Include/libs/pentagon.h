// Pentagon

#pragma once

#include "libs/box.h"

class Pentagon2 {
protected:
  Vector2 c = Vector::Null; //!< Center.
  double r = 1.0; //!< Radius.
public:
  //! \brief Empty.
  Pentagon2() {}
  explicit Pentagon2(const Vector2&, const double&);
  explicit Pentagon2(const double&);

  void Translate(const Vector2&);
  void Scale(const double&);

  bool Inside(const Vector2&) const;
  double R(const Vector2&) const;
  double Signed(const Vector2&) const;
  Vector2 Normal(const Vector2&) const;

  bool Intersect(const Ray2&, double&, double&) const;

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

  Vector2 RandomInside(Random & = Random::R239) const;
  QVector<Vector2> Poisson(const double&, int, Random & = Random::R239) const;

  friend std::ostream& operator<<(std::ostream&, const Pentagon2&);
  void Draw(QGraphicsScene&, const QPen & = QPen(), const QBrush & = QBrush()) const;

protected:
  static const Vector2 vertex[5]; //!< Array of vertexes.
  static const Vector2 normal[5]; //!< Array of normal vectors to the edges.
  static const Vector2 edge[5]; //!< Unit edge vectors.
  static const double Alpha; //!< Constant.
protected:
  static int Sector(const Vector2&);
private:
  static const double epsilon; //!< Internal epsilon constant.
};

/*!
\brief Return the coordinates of the k-th vertex.
\param k Index, should be in [0,4].
*/
inline Vector2 Pentagon2::Vertex(int k) const
{
  return c + vertex[k] * r;
}

/*!
\brief Return the radial vector of the k-th vertex.
\param k Index, should be in [0,4].
*/
inline Vector2 Pentagon2::Radial(int k) const
{
  return vertex[k] * r;
}

/*!
\brief Return the edge vector connecting vertexes k and k+1.
\param k Index, should be in [0,4].
*/
inline Vector2 Pentagon2::Edge(int k) const
{
  return edge[k] * r;
}

/*!
\brief Return the center of the pentagon.
*/
inline Vector2 Pentagon2::Center() const
{
  return c;
}

/*!
\brief Return the radius of the pentagon.
*/
inline double Pentagon2::Radius() const
{
  return r;
}

/*!
\brief Edge length.
*/
inline double Pentagon2::Edge() const
{
  return sqrt((5.0 - sqrt(5.0)) / 2.0) * r;
}

/*!
\brief Inscribed radius.
*/
inline double Pentagon2::InscribedRadius() const
{
  return Alpha * Edge();
}
