// Ellipse
#pragma once

class QGraphicsScene;

#include "libs/box.h"
#include "libs/circle.h"

class Ellipsoid
{
protected:
  Vector c = Vector::Null;//!< Center of the ellipsoid.
  Vector axis = Vector::Z; //!< %Axis.
  double a = 1.0, b = 1.0; //!< %Axis length, and radial length (revolution).
public:
  //! Empty.
  Ellipsoid() {}
  explicit Ellipsoid(const Vector&, const Vector&, const double&, const double&);
  explicit Ellipsoid(const Vector&, const Vector&, const double&);

  //! Empty.
  ~Ellipsoid() {}

  Vector Center() const;
  Vector GetAxis() const;
  double A() const;
  double B() const;

  Box GetBox() const;
  Sphere GetSphere() const;

  Ellipsoid Rotated(const Matrix&) const;
  Ellipsoid Translated(const Vector&) const;

  double R(const Vector&) const;

  int Intersect(const Ray&, double&, double&) const;

  double Area() const;
  double Volume() const;
protected:
  double Value(const Vector&) const;
};

//! Gets the center.
inline Vector Ellipsoid::Center() const
{
  return c;
}

//! Gets the axis.
inline Vector Ellipsoid::GetAxis() const
{
  return axis;
}

//! Major radius of the ellipsoid.
inline double Ellipsoid::A() const
{
  return a;
}

//! Minor radius of the ellipsoid.
inline double Ellipsoid::B() const
{
  return b;
}

/*
\brief Compute the surface area.

While the surface area of a general (tri-axial) ellipsoid involves incomplete elliptic integrals of the first and second kind,
the case of an ellipsoid of revolution (or spheroid) may be expressed in terms of elementary functions.

If radial length is greater than axial length, one has an oblate spheroid; otherwize one has a prolate spheroid.
*/
inline double Ellipsoid::Area() const
{
  const double a2 = a * a;
  const double b2 = b * b;

  double s = Math::TwoPi * b2;
  // Oblate
  if (b >= a)
  {
    double e = 1.0 - a2 / b2;
    double se = sqrt(e);
    return  s * (1.0 + a2 * atanh(se) / (e * b2));
  }
  // Prolate
  else
  {
    double e = 1.0 - b2 / a2;
    double se = sqrt(e);
    return s * (1.0 + a * asin(se) / (e * b));
  }
}

/*!
\brief Compute the volume of an ellipsoid.
*/
inline double Ellipsoid::Volume() const
{
  return 4.0 * Math::Pi * a * b * b / 3.0;
}

// Ellipse class
class Ellipse2
{
protected:
  Vector2 c = Vector2::Null;    //!< Center of the ellipse.
  Vector2 u = Vector2::X; //!< Major axis.
  double a = 1.0, b = 1.0; //!< %Axes lengths.
public:
  //! Empty.
  Ellipse2() {}
  explicit Ellipse2(const Vector2&, const double&, const double&, const Vector2 & = Vector2::X);
  explicit Ellipse2(const double&, const double&);

  //! Empty.
  ~Ellipse2() {}

  Vector2 Center() const;
  double C() const;
  double P() const;
  double Length() const;
  Vector2 Axis() const;
  Vector2 Focus(bool) const;
  double A() const;
  double B() const;

  Vector2 Vertex(const double&) const;
  double Curvature(const double&) const;

  bool Inside(const Vector2&) const;

  // Distance
 // Vector2 Normal(const Vector2&) const;
  double R(const Vector2&) const;
  double Signed(const Vector2&) const;

  Ellipse2 Translated(const Vector2&) const;
  Ellipse2 Scaled(const double&) const;
  Ellipse2 Scaled(const Vector2&, bool = true) const;
  Ellipse2 Rotated(const double&) const;
  Ellipse2 Rotated(const Matrix2&) const;
  Ellipse2 Transformed(const Frame2&) const;

  double Eccentricity() const;
  double Focus() const;

  Box2 GetBox() const;
  Circle2 GetCircle() const;

  double Area() const;
  double Perimeter() const;

  // Stream
  friend std::ostream& operator<<(std::ostream&, const Ellipse2&);

  static Ellipse2 Lerp(const double&, const Ellipse2&, const Ellipse2&);

  // Random
  Vector2 RandomInside(Random & = Random::R239) const;

  void Draw(QGraphicsScene&, const QPen&, const QBrush&) const;
public:
  double Value(const Vector2&) const;
  Vector2 Gradient(const Vector2&) const;
protected:
  Matrix2 MatrixForm() const;
};

//! Gets the center.
inline Vector2 Ellipse2::Center() const
{
  return c;
}

/*!
\brief Area of an ellipse.
*/
inline double Ellipse2::Area() const
{
  return Math::Pi * a * b;
}

/*!
\brief Eccentricity.
*/
inline double Ellipse2::Eccentricity() const
{
  if (a >= b)
  {
    return sqrt(a * a - b * b) / a;
  }
  else
  {
    return sqrt(b * b - a * a) / b;
  }
}

/*!
\brief Focus of an ellipse.
*/
inline double Ellipse2::Focus() const
{
  return sqrt(fabs(a * a - b * b));
}

/*!
\brief Compute the bounding circle.
*/
inline Circle2 Ellipse2::GetCircle() const
{
  return Circle2(c, a);
}

/*!
\brief Return the major axis.
*/
inline Vector2 Ellipse2::Axis() const
{
  return u;
}

/*!
\brief Return the major axis length.
*/
inline double Ellipse2::A() const
{
  return a;
}

/*!
\brief Return the minor axis length.
*/
inline double Ellipse2::B() const
{
  return b;
}

