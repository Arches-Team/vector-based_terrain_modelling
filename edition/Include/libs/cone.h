// Cone

#pragma once

#include "libs/box.h"

// Cone class
class Cone :public Axis
{
protected:
  double ra = 1.0, rb = 0.0;   //!< Radius of the cone at the first and second vertexes of the axis.
  double rlength = -1.0; //!< Internal parameter.
  double conelength = Math::Sqrt2; //!< Side length of the cone.
  Vector2 side = Vector2(-1.0, 1.0) / Math::Sqrt2; //!< %Vector representing the side.
public:
  //! Empty.
  Cone() {}
  explicit Cone(const Vector&, const Vector&, const double&, const double& = 0.0);
  explicit Cone(const double&, const double&, const double&);
  //! Empty.
  ~Cone() {}

  // Parameters
  using Axis::Vertex;
  constexpr double Radius(int) const;

  // Ray-tracing
  int Intersect(const Ray&, double&, double&) const;
  int Intersect(const Ray&, double&, double&, Vector&, Vector&) const;
  int Intersect(const Ray&, double&, Vector&) const;
  bool Inside(const Vector&) const;

  double R(const Vector&) const;
  double Signed(const Vector&) const;
  Vector Normal(const Vector&) const;

  void Rotate(const Matrix&);
  void Translate(const Vector&);
  void Scale(const double&);

  double Area() const;
  double Volume() const;

  friend std::ostream& operator<<(std::ostream&, const Cone&);
  Vector RandomInside(Random & = Random::R239) const;

  Box GetBox() const;

  static Vector RandomDirection(const Vector&, const double&, Random & = Random::R239);

public:
  static const double epsilon; //!< Epsilon value.
  static const Cone Unit; //!< %Unit vertical cone.
};

//! Gets either end-radius of the cone.
inline constexpr double Cone::Radius(int i) const
{
  if (i == 0)
    return ra;
  else
    return rb;
}

// Cone class
class NewCone :public Axis
{
protected:
  double ra = 1.0;   //!< Base radius of the cone.
  double conelength = Math::Sqrt2; //!< Side length of the cone.
  Vector2 side = Vector2(-1.0, 1.0) / Math::Sqrt2; //!< %Vector representing the side.
public:
  //! Empty.
  NewCone() {}
  explicit NewCone(const Vector&, const Vector&, const double&);
  explicit NewCone(const double&, const double&);
  //! Empty.
  ~NewCone() {}

  // Parameters
  Vector operator()(int) const;
  constexpr double Radius() const;

  // Ray-tracing
  int Intersect(const Ray&, double&, double&) const;
  int Intersect(const Ray&, double&, double&, Vector&, Vector&) const;
  int Intersect(const Ray&, double&, Vector&) const;
  bool Inside(const Vector&) const;

  double R(const Vector&) const;
  double Signed(const Vector&) const;
  Vector Normal(const Vector&) const;

  void Rotate(const Matrix&);
  void Translate(const Vector&);
  void Scale(const double&);

  double Area() const;
  double Volume() const;

  friend std::ostream& operator<<(std::ostream&, const NewCone&);
  Vector RandomInside(Random & = Random::R239) const;

  Box GetBox() const;

  static Vector RandomDirection(const Vector&, const double&, Random & = Random::R239);

public:
  static const double epsilon; //!< Epsilon value.
  static const NewCone Unit; //!< %Unit vertical cone.
};

//! Gets either end vertex of the cone.
inline Vector NewCone::operator()(int i) const
{
  if (i == 0) { return a; }
  else { return b; }
}

//! Radius of the cone.
inline constexpr double NewCone::Radius() const
{
  return ra;
}
