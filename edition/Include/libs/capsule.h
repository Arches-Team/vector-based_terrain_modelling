// Cylinder Sphere
#pragma once

#include "libs/cylinder.h"

class Capsule :protected Cylinder
{
protected:
public:
  //! Empty.
  Capsule() {}
  explicit Capsule(const Vector&, const Vector&, const double&);
  //! Empty.
  ~Capsule() {}

  // Specific access to inherited members
  using Cylinder::Radius;
  using Cylinder::Vertex;
  using Cylinder::GetAxis;

  // Inside
  bool Inside(const Vector&) const;
  double R(const Vector&) const;
  double Signed(const Vector&) const;

  bool Intersect(const Ray&) const;
  int Intersect(const Ray&, double&, double&) const;
  int Intersect(const Ray&, double&, double&, Vector&, Vector&) const;
  int Intersect(const Ray&, double&, Vector&) const;

  double R(const Capsule&, double&, double&) const;
  bool Intersect(const Capsule&, double&, double&) const;

  Capsule Rotated(const Matrix&) const;
  Capsule Translated(const Vector&) const;
  Capsule Scaled(const double&) const;
  Capsule Transformed(const Frame&) const;

  Box GetBox() const;

  double Volume() const;
  double Area() const;

  friend std::ostream& operator<<(std::ostream&, const Capsule&);
public:
  static const Capsule Unit; //!< Unit vertical capsule.
};

class OrientedBox2;

class Capsule2 :protected Axis2
{
protected:
  double r = 1.0; //!< Radius.
public:
  //! Empty.
  Capsule2() {}
  explicit Capsule2(const Vector2&, const Vector2&, const double&);
  //! Empty.
  ~Capsule2() {}

  // Parameters
  using Axis2::Vertex;

  // Inside
  bool Inside(const Vector2&) const;
  double R(const Vector2&) const;
  double Signed(const Vector2&) const;

  bool Intersect(const Ray2&) const;
  int Intersect(const Ray2&, double&, double&) const;

  Box2 GetBox() const;

  double Area() const;

  friend std::ostream& operator<<(std::ostream&, const Capsule2&);

  void Draw(QGraphicsScene&, const QPen&, const QBrush&) const;
protected:
  OrientedBox2 GetOrientedBox() const;
};
