// Rays

#pragma once

#include "libs/evector.h"

// Forward declare Matrix class for ray rotation
class Matrix;
class Frame;
class Frame2;

// Ray class
class Ray
{
protected:
  Vector c = Vector::Null; //!< Origin of the ray.
  Vector n = Vector::Z; //!< Direction.
public:
  //! Empty.
  Ray() {}
  explicit Ray(const Vector&, const Vector&);
  //! Empty.
  ~Ray() {}

  Ray Reflect(const Vector&, const Vector&);

  Vector operator()(const double&) const;

  double Length() const;
  double LengthSquared() const;

  // Ray transformations
  Ray Scaled(const Vector&) const;
  Ray Translated(const Vector&) const;
  Ray Rotated(const Matrix&) const;
  Ray InverseTransform(const Frame&) const;

  // Functions to access Vector class components
  Vector Origin() const;
  Vector Direction() const;

  friend std::ostream& operator<<(std::ostream&, const Ray&);
};

/*!
\brief Creates a ray.

The direction should be unit:
\code
Ray ray(Vector(0.0,0.0,0.0),Normalized(Vector(2.0,-1.0,3.0)));
\endcode
\param p Origin.
\param d Direction (should be unit vector).
*/
inline Ray::Ray(const Vector& p, const Vector& d)
{
  c = p;
  n = d;
}

/*!
\brief Return the origin of the ray.
*/
inline Vector Ray::Origin() const
{
  return c;
}

/*!
\brief Return the direction of the ray.
*/
inline Vector Ray::Direction() const
{
  return n;
}

/*!
\brief Computes the location of a vertex along the ray.
\param t Parameter.
*/
inline Vector Ray::operator()(const double& t) const
{
  return c + t * n;
}

/*!
\brief Gets the length of the direction vector.
*/
inline double Ray::Length() const
{
  return Norm(n);
}

/*!
\brief Gets the squared length of the direction vector.
*/
inline double Ray::LengthSquared() const
{
  return n * n;
}

// Ray class
class Ray2
{
protected:
  Vector2 c = Vector2::Null; //!< Origin of the ray.
  Vector2 n = Vector2::X; //!< Direction.
public:
  //! Empty.
  Ray2() {}
  explicit Ray2(const Vector2&, const Vector2&);
  //! Empty.
  ~Ray2() {}

  Vector2 operator()(const double&) const;

  // Functions to access Vector class components
  Vector2 Origin() const;
  Vector2 Direction() const;

  Ray2 InverseTransform(const Frame2&) const;

  friend std::ostream& operator<<(std::ostream&, const Ray2&);
};

/*!
\brief Creates a ray.

\sa Ray::Ray
\param p Origin.
\param d Direction (should be unit vector).
*/
inline Ray2::Ray2(const Vector2& p, const Vector2& d)
{
  c = p;
  n = d;
}

/*!
\brief Return the origin of the ray.
*/
inline Vector2 Ray2::Origin() const
{
  return c;
}

/*!
\brief Return the direction of the ray.
*/
inline Vector2 Ray2::Direction() const
{
  return n;
}

/*!
\brief Computes the location of a vertex along the ray.
\param t Parameter.
*/
inline Vector2 Ray2::operator()(const double& t) const
{
  return c + t * n;
}
