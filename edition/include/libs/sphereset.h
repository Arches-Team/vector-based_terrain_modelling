// Set of spheres

#pragma once

#include <QtCore/QVector>
#include "libs/sphere.h"

// Sphere class
class SphereSet
{
protected:
  QVector<Sphere> spheres; //!< Set of spheres.
public:
  //! Empty.
  explicit SphereSet() {}
  explicit SphereSet(const Sphere&);
  explicit SphereSet(const QVector<Sphere>&);
  explicit SphereSet(const SphereSet&, const Sphere&);
  //! Empty.
  ~SphereSet() {}

  // Append another sphere set
  void Append(const SphereSet&);

  // Shape intersection
  bool Intersect(const SphereSet&) const;

  bool Inside(const Vector&) const;

  // Is the set empty ?
  bool IsEmpty() const;

  // Distance
  double R(const SphereSet&) const;

  // Box
  Box GetBox() const;

  Vector Normal(const Vector&) const;
  double R(const Vector&) const;
  double Signed(const Vector&) const;

  void Rotate(const Matrix&);
  void Translate(const Vector&);
  void Scale(const double&);
  void Transform(const Frame&);

  // Getters
  int Size(void) const;
  Sphere GetSphere(int) const;

  void RemoveDuplicates(const double&);

public:
  static const double epsilon;
};

