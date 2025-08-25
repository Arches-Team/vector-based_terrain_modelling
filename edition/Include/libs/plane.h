// Plane

#pragma once

#include "libs/matrix.h"
#include "libs/ray.h"

// Forward declare 
class Line;
class Segment;
class Quadric;
class Frame;
class Box;
class Sphere;
class Axis;

// Plane 
class Plane {
protected:
  Vector n = Vector::Z;  //!< %Plane normal.
  double c = 0.0;  //!< %Plane coefficient.
public:
  //! Empty.
  Plane() {}
  explicit Plane(const double&, const double&, const double&, const double&);
  explicit Plane(const Vector&, const double&);
  explicit Plane(const Vector&, const Vector&);
  explicit Plane(const Vector&);

  //! Empty.
  ~Plane() {}
  bool Intersect(const Ray&, double&) const;
  bool Inside(const Vector&) const;

  int Side(const Vector&) const;
  Vector Normal() const;
  Vector Vertex() const;
  double Coeff() const;

  double R(const Vector&) const;
  double Signed(const Vector&) const;

  double Eval(const Vector&) const;

  void Transform(const Matrix4&);
  void Translate(const double&);
  void Rotate(const Matrix&);

  Vector Symmetry(const Vector&) const;
  Box Symmetric(const Box&) const;
  Segment Symmetric(const Segment&) const;
  Sphere Symmetric(const Sphere&) const;

  Vector Reflect(const Vector&) const;
  bool Refract(const Vector&, const double&, Vector&) const;

  Vector Intersect(const Line&) const;
  bool Intersect(const Segment&, double&) const;

  Quadric Equation(const Ray&) const;

  static void Cast(const Vector&, double&, double&);

  Frame GetFrame() const;

  friend std::ostream& operator<<(std::ostream&, const Plane&);

  // Intersection between planes
  static Vector Intersection(const Plane&, const Plane&, const Plane&);
  static QVector<Vector> ConvexPoints(const QVector<Plane>&);
  Axis Intersect(const Plane&) const;

  // Reflection and refraction
  static Vector Reflect(const Vector&, const Vector&);
  static bool Refract(const Vector&, const Vector&, const double&, Vector&);
public:
  static const double epsilon; //!< Small value for checking intersection with a ray.
  static const Plane XY;   //!< Reference horizontal plane.
};

/*!
\brief Creates a plane given normal and distance.

\param n %Plane normal, should be unit.
\param c %Plane constant.
*/
inline Plane::Plane(const Vector& n, const double& c)
{
  Plane::n = n;
  Plane::c = c;
}

/*!
\brief Creates a plane given its analytic equation.

\param a,b,c Coefficients of the plane normal, the plane normal should be unit..
\param d %Plane constant.
*/
inline Plane::Plane(const double& a, const double& b, const double& c, const double& d)
{
  Plane::n = Vector(a, b, c);
  Plane::c = d;
}

/*!
\brief Returns the normal to the plane.
*/
inline Vector Plane::Normal() const
{
  return n;
}

/*!
\brief Returns the coefficient of the plane.
*/
inline double Plane::Coeff() const
{
	return c;
}

/*!
\brief Evaluates the equation of the plane.

This can be interpreted as a signed distance to the plane.

\param p Point.
\sa Plane::R(), Plane::Signed()
*/
inline double Plane::Eval(const Vector& p) const
{
  return n * p - c;
}

/*!
\brief Compute the reflected direction.

\param d Direction.
*/
inline Vector Plane::Reflect(const Vector& d) const
{
  return d - 2.0 * (d * n) * n;
}

/*!
\brief Given an input (normalized) normal, compute a reflected direction.
The code is the same as:
\code
Vector r=d-2.0*(n*d)*n;
\endcode
*/
inline Vector Plane::Reflect(const Vector& d, const Vector& n)
{
  return d - 2.0 * (n * d) * n;
}
