// Torus

#pragma once

#include "libs/circle.h"

// Torus class
class Torus :public Circle {
protected:
  double s = 1.0; //!< Small radius of the torus.
public:
  explicit Torus(const Vector&, const Vector&, const double& = 1.0, const double& = 1.0);
  explicit Torus(const Circle&, const double& = 1.0);
  explicit Torus(const double&, const double& = 1.0);

  //! Empty.
  ~Torus() {}

  int Intersect(const Ray&, double*, Vector*) const;
  bool Inside(const Vector&) const;
  Box GetBox() const;

  double R(const Vector&) const;
  double Signed(const Vector&) const;
  Vector Normal(const Vector&) const;

  double Area() const;
  double Volume() const;
  double Small() const;

  void Scale(const double&);

  Vector RandomInside(Random & = Random::R239) const;

  // Stream
  friend std::ostream& operator<<(std::ostream&, const Torus&);
public:
  static constexpr const double epsilon = 1e-06; //!< Epsilon constant for some inside-outside functions.
  static const Torus Unit; //!< %Unit vertical torus, major and minor radii equal to 1.
};

/*!
\brief Return the small radius of the torus.

Use Torus::Radius() to get the major radius.
*/
inline double Torus::Small() const
{
  return s;
}
