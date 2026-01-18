// Hemispheres

#pragma once

#include "libs/box.h"

class HemiSphere {
protected:
  Vector c = Vector::Null; //!< Center.
  double r = 1.0; //!< Radius.
  Vector axis = Vector::Z; //!< Axis.
public:
  explicit HemiSphere(const Vector&, const Vector&, const double&);
  explicit HemiSphere(const double&);
  //! Empty.
  ~HemiSphere() {}

  Vector Normal(const Vector&) const;
  double R(const Vector&) const;
  double Signed(const Vector&) const;

  Box GetBox() const;

  // Parameters
  Vector Center() const;
  double Radius() const;
  Vector GetAxis() const;

  Vector Fibonacci(int, int);

  friend std::ostream& operator<<(std::ostream&, const HemiSphere&);
public:
  static Vector FibonacciUnit(int, int);
  static Vector RandomDirection(const Vector&, Random & = Random::R239);
};

//! Gets the center of a hemisphere.
inline Vector HemiSphere::Center() const
{
  return c;
}

//! Gets the radius of a hemisphere.
inline double HemiSphere::Radius() const
{
  return r;
}

//! Gets the axis of the hemisphere.
inline Vector HemiSphere::GetAxis() const
{
  return axis;
}
