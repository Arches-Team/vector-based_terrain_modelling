// Quaternions
#pragma once

#include "libs/matrix.h"

// Definition of the Quaternion structure
class Quaternion
{
protected:
  double w = 1.0, x = 0.0, y = 0.0, z = 0.0; //!< %Quaternion components.
public:
  //! Empty
  Quaternion() {}
  explicit Quaternion(const double&, const double&, const double&, const double&);
  explicit Quaternion(const Vector&, const double&);
  explicit Quaternion(const Vector&, const Vector&);
  explicit Quaternion(const Matrix&);

  // Conversion between quaternions, matrices, and angle-axes
  Matrix RotationMatrix() const;
  void AngleAxis(double&, Vector&) const;

  double& operator[] (int);
  constexpr double operator[] (int) const;

  // Arithmetic operations
  Quaternion operator+ (const Quaternion&) const;
  Quaternion operator- (const Quaternion&) const;
  Quaternion operator* (const Quaternion&) const;
  Quaternion operator* (const double&) const;
  friend Quaternion operator* (const double&, const Quaternion&);
  Quaternion operator- () const;

  // Unit quaternion functions
  Quaternion Inverse() const;
  Quaternion Conjugate() const;

  friend double Norm(const Quaternion&);
  Quaternion Exp() const;
  Quaternion Log() const;

  // Rotation of a point by a quaternion
  Vector operator* (const Vector&) const;

  // Linear interpolation
  static Quaternion Lerp(const double&, const Quaternion&, const Quaternion&);

  static Quaternion FromAngles(const double&, const double&, const double&);

  friend std::ostream& operator<<(std::ostream&, const Quaternion&);
public:
  static const Quaternion Null; //!< Null quaternion.
  static const Quaternion Identity; //!< Identity.
protected:
  static const double Epsilon; //!< Cutoff value for sine function whenever the angle gets near zero.
};

/*!
\brief Create a quaternion.
\param w Angular parameter.
\param x, y, z %Axis coordinates.
*/
inline Quaternion::Quaternion(const double& w, const double& x, const double& y, const double& z) :x(x), y(y), z(z), w(w)
{
}

/*!
\brief Overloaded.
*/
inline double& Quaternion::operator[] (int i)
{
  if (i == 0) return x;
  else if (i == 1) return y;
  else if (i == 2) return z;
  else return w;
}

/*!
\brief Returns the i-th coefficient of a quaternion.

Terms are x, y, and z if index is in [0,2], and w otherwise.
*/
inline constexpr double Quaternion::operator[] (int i) const
{
  if (i == 0) return x;
  else if (i == 1) return y;
  else if (i == 2) return z;
  else return w;
}

/*!
\brief Computes the conjugatge of a quaternion.
*/
inline Quaternion Quaternion::Conjugate() const
{
  // This is unit length
  return Quaternion(w, -x, -y, -z);
}

/*!
\brief Multiply a quaternion by a double value.
\param c Real.
*/
inline Quaternion Quaternion::operator* (const double& c) const
{
  return Quaternion(c * w, c * x, c * y, c * z);
}

/*!
\brief Left multiply a quaternion by a double value.
\param c Real.
\param q %Quaternion.
*/
inline Quaternion operator* (const double& c, const Quaternion& q)
{
  return Quaternion(c * q.w, c * q.x, c * q.y, c * q.z);
}
