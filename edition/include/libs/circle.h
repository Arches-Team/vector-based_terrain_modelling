// Circle

#pragma once

#include <QtCore/QVector> 

#include "libs/box.h"

// Circle class
class Circle
{
protected:
  Vector c = Vector::Null;    //!< Center of the circle.
  Vector axis = Vector::Z; //!< %Axis.
  double r = 1.0;    //!< Radius.
public:
  //! Empty.
  Circle() {}
  explicit Circle(const double&);
  explicit Circle(const Vector&, const Vector&, const double&);
  explicit Circle(const Vector&, const Vector&, const Vector&);

  //! Empty.
  ~Circle() {}

  bool Intersect(const Ray&, double&) const;

  Vector Center() const;
  Vector Axis() const;
  double Radius() const;

  // Distance
  double R(const Vector&) const;
  Vector Normal(const Vector&) const;

  void Rotate(const Matrix&);
  void Translate(const Vector&);
  void Scale(const double&);

  Box GetBox() const;

  double Area() const;
  static double Area(const double&);
  // Stream
  friend std::ostream& operator<<(std::ostream&, const Circle&);

  // Random
  Vector RandomInside(Random & = Random::R239) const;
  Vector RandomOn(Random & = Random::R239) const;
  Vector Vogel(int, int) const;
protected:
  static constexpr const double epsilon = 1.0e-10; //!< Internal epsilon constant.
};

//! Gets the center of a circle.
inline Vector Circle::Center() const
{
  return c;
}

//! Radius of the circle.
inline double Circle::Radius() const
{
  return r;
}

//! %Axis of the circle.
inline Vector Circle::Axis() const
{
  return axis;
}

/*!
\brief Area of a circle.

\sa Circle::Area(const double&)
*/
inline double Circle::Area() const
{
  return Math::Pi * r * r;
}

/*!
\brief Area of a circle given its radius.

This static function avoids the creation of an instance of a circle.
\sa Circle::Area()
\param r Radius.
*/
inline double Circle::Area(const double& r)
{
  return Math::Pi * r * r;
}

class Disc :public Circle
{
protected:
public:
  //! Empty.
  Disc() {}
  explicit Disc(const double&);
  explicit Disc(const Vector&, const Vector&, const double&);
  explicit Disc(const Vector&, const Vector&, const Vector&);
  //! Empty.
  ~Disc() {}

  double R(const Vector&) const;
  Vector Normal(const Vector&) const;
};

class QuadricCurve2;
class CubicCurve2;

class Circle2
{
protected:
  Vector2 c = Vector2::Null; //!< Center.
  double r = 1.0; //!< Radius.
public:
  //! Empty
  Circle2() {}
  explicit Circle2(const double&);
  explicit Circle2(const Vector2&, const double&);
  explicit Circle2(const Vector2&, const Vector2&, const Vector2&);
  explicit Circle2(Vector2*, int);
  explicit Circle2(const QVector<Vector2>&);

  bool Inside(const Vector2&) const;
  bool Inside(const Circle2&) const;
  bool Inside(const Box2&) const;

  bool InsideRange(const Vector2&, const double&) const;
  Vector2 Center() const;
  Vector2 Vertex(const double&) const;

  double Radius() const;

  void Extend(const Vector2&);
  Circle2 Extended(const double&) const;

  bool Intersect(const Box2&) const;
  bool Intersect(const Circle2&) const;
  bool Intersect(const Ray2&, double&, double&) const;
  bool Intersect(const Ray2&) const;
  bool Intersect(const Segment2&) const;

  Vector2 Normal(const Vector2&) const;
  double R(const Vector2&) const;

  Box2 GetBox() const;

  Circle2 Translated(const Vector&) const;
  Circle2 Rotated(const Matrix2&) const;
  Circle2 Transformed(const Frame2&) const;

  void Translate(const Vector2&);
  void Rotate(const Matrix2&);
  void Scale(const double&);

  double Area() const;
  double Area(const Circle2&) const;

  Vector2 RandomInside(Random & = Random::R239) const;
  Vector2 RandomOn(Random & = Random::R239) const;
  Vector2 Vogel(int, int) const;

  static Vector2 RandomUnit(Random & = Random::R239);

  // Stream
  friend std::ostream& operator<<(std::ostream&, const Circle2&);

  QVector<Vector2> Poisson(const double&, int, bool = false, Random & = Random::R239) const;
  QuadricCurve2 QuadricBezierArc(const double&) const;
  CubicCurve2 CubicBezierArc(const double&) const;

  // Drawing
  void Draw(QGraphicsScene&, const QPen & = QPen(), const QBrush & = QBrush()) const;
public:
  static const Circle2 Unit; //!< Unit circle.
  static const Circle2 Infinite; //!< Infinite circle.
protected:
  static constexpr const double Epsilon = 1.e-6; //!< %Epslion for intersection tests
};

//! Center of the circle.
inline Vector2 Circle2::Center() const
{
  return c;
}

//! Radius of the circle.
inline double Circle2::Radius() const
{
  return r;
}

//! Area of the circle.
inline double Circle2::Area() const
{
  return Math::Pi * r * r;
}

/*!
\brief Compute the axis aligned bounding box of a circle.
*/
inline Box2 Circle2::GetBox() const
{
  return Box2(c, r);
}

/*!
\brief Translate a circle.
\param t Translation vector.
*/
inline Circle2 Circle2::Translated(const Vector& t) const
{
  return Circle2(c + t, r);
}

/*!
\brief Rotates a circle.

\param r Rotation matrix.
*/
inline Circle2 Circle2::Rotated(const Matrix2& r) const
{
  return Circle2(r * c, Circle2::r);
}

/*!
\brief Return a circle transformed by a frame.
\param f Transformation.
*/
inline Circle2 Circle2::Transformed(const Frame2& f) const
{
  return Circle2(f.Transform(c), r);
}

/*!
\brief Translate a circle.

\param t Translation vector.
*/
inline void Circle2::Translate(const Vector2& t)
{
  c += t;
}

/*!
\brief Rotates a circle.

\param m Rotation matrix.
*/
inline void Circle2::Rotate(const Matrix2& m)
{
  c = m * c;
}

/*!
\brief Scales a circle.

\param s Scaling factor.
*/
inline void Circle2::Scale(const double& s)
{
  c *= s;
  r *= s;
}

class Disc2 :public Circle2
{
protected:
public:
  //! Empty.
  Disc2() {}
  explicit Disc2(const double&);
  explicit Disc2(const Circle2&);
  explicit Disc2(const Vector2&, const double&);
  explicit Disc2(const Vector2&, const Vector2&, const Vector2&);

  //! Empty.
  ~Disc2() {}

  double R(const Vector2&) const;
  double Signed(const Vector2&) const;
  Vector2 Normal(const Vector2&) const;

  // Stream
  friend std::ostream& operator<<(std::ostream&, const Disc2&);
};



