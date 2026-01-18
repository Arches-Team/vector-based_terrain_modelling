// Hexagon

#pragma once

#include "libs/box.h"

class Circle2;

class Hexagon2 {
protected:
  Vector2 c = Vector2::Null; //!< Center.
  double r = 1.0; //!< Radius.
public:
  //! Empty
  Hexagon2() {}
  explicit Hexagon2(const Vector2&, const double&);
  explicit Hexagon2(const double&);

  void Translate(const Vector2&);
  void Scale(const double&);

  double Width() const;
  double Height() const;

  bool Inside(const Vector2&) const;
  bool Intersect(const Hexagon2&) const;

  double R(const Vector2&) const;
  double Signed(const Vector2&) const;
  Vector2 Normal(const Vector2&) const;

  Vector2 Center() const;
  double Radius() const;

  Vector2 Radial(int) const;
  Vector2 Edge(int) const;
  Vector2 Vertex(int) const;
  Box2 GetBox() const;

  double Area() const;
  double Perimeter() const;

  Vector2 RandomInside(Random & = Random::R239) const;
  QVector<Vector2> Poisson(const double&, int, Random & = Random::R239) const;

  bool Intersect(const Circle2&) const;

  friend std::ostream& operator<<(std::ostream&, const Hexagon2&);
  void Draw(QGraphicsScene&, const QPen & = QPen(), const QBrush & = QBrush()) const;

protected:
  static const Vector2 vertex[6]; //!< Array of vertices.
  static const Vector2 normal[6]; //!< Array of normal vectors to the edges.
  static const Vector2 edge[6]; //!< Unit edge vectors.
  static const double Alpha; //!< Constant sin(&pi;/6)=&radic;3/2.
protected:
  static int Sector(const Vector2&);
};

/*!
\brief Return the coordinates of the k-th vertex.
\param k Index, should be in [0,5].
*/
inline Vector2 Hexagon2::Vertex(int k) const
{
  return c + vertex[k] * r;
}

/*!
\brief Return the radial vector of the k-th vertex.
\param k Index, should be in [0,5].
*/
inline Vector2 Hexagon2::Radial(int k) const
{
  return vertex[k] * r;
}

/*!
\brief Return the edge vector connecting vertexes k and k+1.
\param k Index, should be in [0,5].
*/
inline Vector2 Hexagon2::Edge(int k) const
{
  return edge[k] * r;
}

/*!
\brief Return the center of the hexagonal cell.
*/
inline Vector2 Hexagon2::Center() const
{
  return c;
}

/*!
\brief Return the radius of the hexagonal cell.
*/
inline double Hexagon2::Radius() const
{
  return r;
}

/*!
\brief Width of the hexagon.
*/
inline double Hexagon2::Width() const
{
  return r * 2.0;
}

/*!
\brief Height of the hexagon.
*/
inline double Hexagon2::Height() const
{
  return r * Math::Sqrt3;
}

class HexagonAlpha2 :public Hexagon2
{
protected:
  double a; //!< Angle.
  Vector2 x; //!< Base vector.
public:
  //! Empty
  HexagonAlpha2();
  explicit HexagonAlpha2(const Vector2&, const double&, const double& = 0.0);
  explicit HexagonAlpha2(const Hexagon2&, const double& = 0.0);

  void Rotate(const double&);

  Vector2 Vertex(int) const;

  Box2 GetBox() const;
  Box2 GetTightBox() const;

  // Drawing
  void Draw(QGraphicsScene&, const QPen & = QPen(), const QBrush & = QBrush()) const;
  friend std::ostream& operator<<(std::ostream&, const HexagonAlpha2&);
};

/*!
\brief Constructor.
*/
inline HexagonAlpha2::HexagonAlpha2() : a(0.0), x(Vector2(1.0, 0.0))
{
}

class Hexagonal :protected Hexagon2
{
protected:
  double a = 0.0, b = 0.0; //! Heights.
public:
  //! Empty
  Hexagonal() {}
  explicit  Hexagonal(const Vector2&, const double&, const double&, const double&);
  explicit Hexagonal(const Hexagon2&, const double&, const double&);
  //! Empty
  ~Hexagonal() {}

  Vector Center(int) const;

  double Radius() const;
  Vector Vertex(int) const;
  Vector Radial(int) const;
  Vector Edge(int) const;
  Vector Normal(int) const;

  double Area() const;
  double Volume() const;

  Box GetBox() const;

  double Signed(const Vector&) const;
  double R(const Vector&) const;

  friend std::ostream& operator<<(std::ostream&, const Hexagonal&);
public:
  static const Hexagonal Unit; //!< Unit hexagonal prism.
};

/*!
\brief Return the base or apex vertex of the hexagonal prism.
\param k Index.
*/
inline Vector Hexagonal::Center(int k) const
{
  return (k == 0) ? c.ToVector(a) : c.ToVector(b);
}

/*!
\brief Return the radius of the hexagonal cell.
*/
inline double Hexagonal::Radius() const
{
  return r;
}

/*!
\brief Return the coordinates of the k-th vertex.
\param k Index, should be in [0,5] for lower cap, in [6,11] for upper cap.
*/
inline Vector Hexagonal::Vertex(int k) const
{
  return Hexagon2::Vertex(k % 6).ToVector(k < 6 ? a : b);
}

/*!
\brief Return the radial vector of the k-th vertex.
\param k Index, should be in [0,5].
*/
inline Vector Hexagonal::Radial(int k) const
{
  return (vertex[k] * r).ToVector();
}

/*!
\brief Return the edge vector connecting vertexes k and k+1.
\param k Index, should be in [0,5].
*/
inline Vector Hexagonal::Edge(int k) const
{
  return (edge[k] * r).ToVector();
}

/*!
\brief Return the normal vector to the k-th edge connecting vertexes k and k+1.
\param k Index, should be in [0,5].
*/
inline Vector Hexagonal::Normal(int k) const
{
  return normal[k].ToVector();
}

/*!
\brief Area of the hexagonal prism.
*/
inline double Hexagonal::Area() const
{
  return Hexagon2::Area() * 2.0 + Hexagon2::Perimeter() * (b - a);
}

/*!
\brief Volume of the hexagonal prism.
*/
inline double Hexagonal::Volume() const
{
  return Hexagon2::Area() * (b - a);
}
