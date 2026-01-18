// Rectangle

#pragma once

#include "libs/box.h"

class Rectangles2;

// Rectangle class
class Rectangles
{
protected:
  Vector x = Vector::X, y = Vector::Y, z = Vector::Z;    //!< Axes.
  Vector c = Vector::Null;    //!< Center.
  double a = 1.0, b = 1.0;    //!< Side lengths.
public:
  //! Empty.
  Rectangles() {}
  explicit Rectangles(const double&, const double&);
  explicit Rectangles(const Vector&, const Vector&, const Vector&);
  explicit Rectangles(const Vector&, const Vector&, const Vector&, const double&, const double&);

  //! Empty.
  ~Rectangles() {}

  Vector Center() const;

  Vector Vertex(const double&, const double&) const;
  Vector Normal() const;

  // Distance
  double R(const Vector&) const;
  Vector Normal(const Vector&) const;

  void Rotate(const Matrix&);
  void Translate(const Vector&);
  void Scale(const double&);

  Box GetBox() const;

  double Area() const;
  double Width() const;
  double Height() const;

  // Stream
  friend std::ostream& operator<<(std::ostream&, const Rectangles&);
protected:
  static double epsilon;//!< Internal epsilon constant.
};

//! Gets the center of the rectangle.
inline Vector Rectangles::Center() const
{
  return c;
}

//! Compute the area of the rectangle.
inline double Rectangles::Area() const
{
  return a * b;
}

//! Return the half width.
inline double Rectangles::Width() const
{
  return a;
}

//! Return the half height.
inline double Rectangles::Height() const
{
  return b;
}

/*!
\brief Return the normal.
*/
inline Vector Rectangles::Normal() const
{
  return z;
}

class Circle2;

class Rectangles2 {
protected:
  Vector2 x = Vector2::X; //!< Base horizontal vector.
  Vector2 c = Vector2::Null; //!< Center.
  double a = 1.0, b = 1.0; //!< Width and height.
public:
  //! Empty.
  Rectangles2() {}
  explicit Rectangles2(const Vector2&, const Vector2&, const double&, const double&);
  explicit Rectangles2(const Box2&);

  //! Empty.
  ~Rectangles2() {}

  Vector2 Center() const;
  Vector2 Vertex(int) const;
  double Width() const;
  double Height() const;

  // Distance
  double R(const Vector2&) const;
  Vector2 Normal(const Vector2&) const;

  Box2 GetBox() const;

  Rectangles ToRectangles() const;

  bool Intersect(const Circle2&) const;
  bool Intersect(const Rectangles2&) const;

  void Rotate(const Matrix2&);
  void Translate(const Vector2&);
  void Scale(const double&);
  Rectangles2 Transformed(const Frame2&) const;

  QPolygonF GetQt() const;
  void Draw(QGraphicsScene&, const QPen & = QPen(), const QBrush & = QBrush()) const;
protected:
  static double epsilon;//!< Internal \htmlonly\epsilon;\endhtmlonly constant.
protected:
  bool Overlap(const Rectangles2&) const;
};

/*
\brief Return the center of the rectangle.
*/
inline Vector2 Rectangles2::Center() const
{
  return c;
}

/*!
\brief Returns the k-th vertex of the rectangle.
\param k Index.
\sa Box2::Vertex(int) const
*/
inline Vector2 Rectangles2::Vertex(int k) const
{
  return c + x * ((k & 1) ? a : -a) + x.Orthogonal() * ((k & 2) ? b : -b);
}

/*!
\brief Return the width.
*/
inline double Rectangles2::Width() const
{
  return a;
}

/*!
\brief Return the height.
*/
inline double Rectangles2::Height() const
{
  return b;
}
