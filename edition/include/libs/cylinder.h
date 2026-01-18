// Cylinder
#pragma once


#include "libs/box.h"
#include "libs/axis.h"
#include "libs/circle.h"

// Cylinder class
class Cylinder :public Axis
{
protected:
  double r = 1.0; //!< Radius.
public:
  //! Empty.
  Cylinder() {}
  explicit Cylinder(const Vector&, const Vector&, const double&);
  explicit Cylinder(const double&, const double&, const double& = 1.0);
  explicit Cylinder(const Vector&, const double&, const double& = 1.0);
  explicit Cylinder(const Circle&, const double&, bool = false);
  explicit Cylinder(const Box&, const Axis&);

  //! Empty.
  ~Cylinder() {}

  // Parameters
  using Axis::Vertex;
  double Radius() const;

  // Ray-tracing
  int Intersect(const Ray&, double&, double&) const;
  int Intersect(const Ray&, double&, double&, Vector&, Vector&) const;
  int Intersect(const Ray&, double&, Vector&) const;
  bool Intersect(const Ray&) const;
  bool Inside(const Vector&) const;

  double Volume() const;
  double Area() const;

  void Scale(const double&);
  Cylinder Transformed(const FrameScaled&) const;

  // Mapping
  static Vector2 Cast(const Vector&);
  friend std::ostream& operator<<(std::ostream&, const Cylinder&);
  Vector RandomInside(Random & = Random::R239) const;

  Box GetBox() const;

  Cylinder Rotated(const Matrix&) const;
  Cylinder Translated(const Vector&) const;
  Cylinder Scaled(const double&) const;

  double R(const Vector&) const;
  double R(const Vector&, double&) const;
  double Signed(const Vector&) const;
  Vector Normal(const Vector&) const;
  Vector SignedNormal(const Vector&) const;
public:
  static const double epsilon; //!< Epsilon value for intersection tests.
  static const Cylinder Unit; //!< %Unit vertical cylinder.
};

//! Gets the radius of a cylinder.
inline double Cylinder::Radius() const
{
  return r;
}

//! Computes the volume of a cylinder.
inline double Cylinder::Volume() const
{
  return Math::Pi * r * r * length;
}

//! Computes the surface area of a cylinder.
inline double Cylinder::Area() const
{
  return Math::TwoPi * r * (length + r);
}
