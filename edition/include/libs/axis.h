// Axis
#pragma once

#include "libs/matrix.h"
#include "libs/quadric.h"

class Ray;
class Sphere;
class Circle2;
class Box2;

// Axis class
class Axis
{
protected:
  Vector a = Vector::Null, b = Vector::Z;   //!< End vertexes of the axis.
  Vector axis = Vector::Z;  //!< Normalized axis vector.
  double length = 1.0; //!< Length of the axis.
public:
  //! Empty.
  Axis() {}
  explicit Axis(const Vector&, const Vector&);
  //! Empty.
  ~Axis() {}

  void Rotate(const Matrix&);
  void Translate(const Vector&);
  void Scale(const double&);
  void Scale(const Vector&);

  Quadric Equation(const Ray&) const;
  Vector Vertex(int) const;
  Vector Point(const double&) const;
  Vector GetAxis() const;
  double Length() const;

  Vector Symmetric(const Vector&) const;
  Sphere Symmetric(const Sphere&) const;

  friend std::ostream& operator<<(std::ostream&, const Axis&);

  Vector Normal(const Vector&) const;
  double R(const Vector&) const;
  double R(const Axis&) const;
  double R(const Vector&, double&) const;

  Matrix GetFrame() const;

  static Matrix GetFrame(const Vector&);
  static Vector BoxVector(const Vector&);


  static const Axis X,Y,Z; //!< Axes from origin, identical to Axis(Vector::Null,Vector::X) ...
protected:
  double Radial(const Vector&, Vector&, Vector&) const;
  Vector2 Radial(const Vector&) const;
};

/*!
\brief Return one of the end vertexes of the axis.
*/
inline Vector Axis::Vertex(int i) const
{
  if (i == 0) return a;
  else return b;
}

//! Return the axis length.
inline double Axis::Length() const
{
  return length;
}

/*!
\brief Returns the normalized axis vector.
*/
inline Vector Axis::GetAxis() const
{
  return axis;
}

/*!
\brief Compute the box vector extent of a unit circle with a given axis.

\param axis %Axis (should be unit).

This function is used in particular in Cylinder::GetBox(), Cone::GetBox(), Capsule::GetBox().
*/
inline Vector Axis::BoxVector(const Vector& axis)
{
  return Vector(sqrt(1.0 - axis[0] * axis[0]), sqrt(1.0 - axis[1] * axis[1]), sqrt(1.0 - axis[2] * axis[2]));
}

/*!
\brief Compute a point on the axis.
*/
inline Vector Axis::Point(const double& t) const
{
  return a + t * axis;
}

/*!
\brief Compute the radial coordinates of a point.
\param p Point.
\param av Returned axis vector.
\param rv Returned radial vector.
\return Axial coordinate.
*/
inline double Axis::Radial(const Vector& p, Vector& av, Vector& rv) const
{
  Vector n = p - a;
  double t = (n * axis);

  // Axis
  av = t * axis;

  // Radial
  rv = n - av;

  return t;
}

/*!
\brief Compute the radial coordinates of a point.

\image html radial.png

\param p Point.
\return A vector (x,y) where x is the radial coordinate and y the axial coordinate.
*/
inline Vector2 Axis::Radial(const Vector& p) const
{
  Vector n = p - a;
  double y = (n * axis);

  double x = sqrt(n * n - y * y);

  return Vector2(x, y);
}

// Axis class
class Axis2
{
protected:
  Vector2 a = Vector2(0.0, 0.0), b = Vector2(1.0, 0.0);   //!< End vertexes of the axis.
  Vector2 axis = Vector2::X;  //!< Normalized axis vector.
  double length = 1.0; //!< Length of the axis.
public:
  //! Empty.
  Axis2() {}
  explicit Axis2(const Vector2&, const Vector2&);
  //! Empty.
  ~Axis2() {}

  Vector2 Vertex(int) const;
  Vector2 Point(const double&) const;
  Vector2 GetAxis() const;
  constexpr double Length() const;

  Vector2 Symmetric(const Vector2&) const;
  Circle2 Symmetric(const Circle2&) const;
  Box2 Symmetric(const Box2&) const;

  friend std::ostream& operator<<(std::ostream&, const Axis2&);

protected:
  double Radial(const Vector2&, Vector2&, Vector2&) const;
  Vector2 Radial(const Vector2&) const;
};

/*!
\brief Return one of the end vertexes of the axis.
*/
inline Vector2 Axis2::Vertex(int i) const
{
  if (i == 0) return a;
  else return b;
}

//! Return the axis length.
inline constexpr double Axis2::Length() const
{
  return length;
}

/*!
\brief Returns the normalized axis vector.
*/
inline Vector2 Axis2::GetAxis() const
{
  return axis;
}

/*!
\brief Compute a point on the axis.
*/
inline Vector2 Axis2::Point(const double& t) const
{
  return a + t * axis;
}

/*!
\brief Compute the radial coordinates of a point.
\param p Point.
\param av Returned axis vector.
\param rv Returned radial vector.
\return Axial coordinate.
*/
inline double Axis2::Radial(const Vector2& p, Vector2& av, Vector2& rv) const
{
  Vector2 n = p - a;
  double t = (n * axis);

  // Axis
  av = t * axis;

  // Radial
  rv = n - av;

  return t;
}

/*!
\brief Compute the radial coordinates of a point.
\param p Point.
\return A vector (x,y) where x is the radial coordinate and y the axial coordinate.
*/
inline Vector2 Axis2::Radial(const Vector2& p) const
{
  Vector2 n = p - a;
  double y = (n * axis);

  double x = sqrt(n * n - y * y);

  return Vector2(x, y);
}
