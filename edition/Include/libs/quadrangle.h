// Quadrangle

#pragma once

#include "libs/box.h"
#include "libs/rectangle.h"

// Quadrangle class
class Quadrangle
{
protected:
  Vector p[4] = { Vector::Null,Vector::X,Vector::X + Vector::Y,Vector::Y }; //!< Vertices, stored in counter-clockwize order.
public:
  //! Empty.
  Quadrangle() {}
  explicit Quadrangle(const Vector&, const Vector&, const Vector&, const Vector&);
  explicit Quadrangle(double, double);
  explicit Quadrangle(double);

  //! Empty.
  ~Quadrangle() {}

  Vector Vertex(const double&, const double&) const;
  Vector Normal(const double&, const double&) const;
  Vector Vertex(int) const;

  Vector Normal() const;

  Box GetBox() const;

  void Rotate(const Matrix&);
  void Translate(const Vector&);
  void Scale(const double&);
  void Transform(const FrameScaled&);

  Quadrangle Transformed(const FrameScaled&) const;
  Quadrangle Translated(const Vector&) const;

  double Area() const;

  // Stream
  friend std::ostream& operator<<(std::ostream&, const Quadrangle&);
protected:
  static double epsilon; //!< Internal \htmlonly\epsilon;\endhtmlonly constant.
};

/*
\brief Return the i-th vertex.
\param i Index.
*/
inline Vector Quadrangle::Vertex(int i) const
{
  return p[i];
}

// Quadrangle class
class Quadrangle2
{
protected:
  Vector2 p[4] = { Vector2::Null,Vector2::X,Vector2::X + Vector2::Y,Vector2::Y }; //!< Vertices, stored in counter-clockwize order.
public:
  //! Empty.
  Quadrangle2() {}
  explicit Quadrangle2(const Vector2&, const Vector2&, const Vector2&, const Vector2&);
  explicit Quadrangle2(const Quadrangle&);
  explicit Quadrangle2(const double&);

  //! Empty.
  ~Quadrangle2() {}

  Vector2 Vertex(const double&, const double&) const;
  Vector2 Vertex(int) const;

  bool InverseBilinear(const Vector2&, double&, double&) const;

  Box2 GetBox() const;
  bool Inside(const Vector2&) const;
  double R(const Vector2&) const;
  Vector2 Normal(const Vector2&) const;

  void Translate(const Vector2&);
  void Scale(const double&);
  void Rotate(const Matrix2&);
  void Transform(const Frame2&);

  // Mathematics
  double Perimeter() const;
  double Area() const;

  // Stream
  friend std::ostream& operator<<(std::ostream&, const Quadrangle2&);

  void Draw(QGraphicsScene&, const QPen & = QPen(), const QBrush & = QBrush()) const;
protected:
  static double epsilon; //!< Internal \htmlonly\epsilon;\endhtmlonly constant.
};

/*
\brief Return the i-th vertex.
\param i Index.
*/
inline Vector2 Quadrangle2::Vertex(int i) const
{
  return p[i];
}