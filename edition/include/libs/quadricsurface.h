// Quadric surfaces

#pragma once

#include "libs/evector.h"
#include "libs/matrix.h"

class QuadricSurface {
protected:
  double c[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; //!< %Quadric coefficients.
public:
  QuadricSurface() {}
  explicit QuadricSurface(const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);
  explicit QuadricSurface(const Matrix&, const double& = 1.0);

  double Value(const Vector2&) const;
  Vector2 Gradient(const Vector2&) const;
  Vector2 Normal(const Vector2&) const;

  QuadricSurface operator/(double) const;
  QuadricSurface operator*(double) const;
  QuadricSurface operator+(const QuadricSurface&) const;
  QuadricSurface operator-(const QuadricSurface&) const;

  QuadricSurface Translated(const Vector2&) const;
  QuadricSurface Scaled(const double&) const;


  double operator()(int, int) const;
  double& operator()(int, int);
protected:
  constexpr const int Index(int, int) const;
protected:
  static constexpr int terms[9] = { 8, 6, 3, 7, 5, 2, 4, 1, 0 }; //! Coefficients for locating the terms of the quadric equation.
};

/*!
\brief Return the index of the coefficient x<sup>a</sup> y<sup>b</sup>
\param a,b Exponents.
*/
inline constexpr const int QuadricSurface::Index(int a, int b) const
{
  return terms[a + 3 * b];
}

/*!
\brief Compute the value.
\param p Point.
*/
inline double QuadricSurface::Value(const Vector2& p) const
{
  const double x = p[0];
  const double y = p[1];

  // Could be: return c[Index(2,2)] * x * x * y * y + c[Index(2,1)] * x * x * y + ...
  return c[0] * x * x * y * y + c[1] * x * x * y + c[2] * x * y * y + c[3] * x * x + c[4] * y * y + c[5] * x * y + c[6] * x + c[7] * y + c[8];
}

/*!
\brief Overloaded operator.
\param a Real.
*/
inline QuadricSurface QuadricSurface::operator/(double a) const
{
  double ia = 1.0 / a;
  return QuadricSurface(c[0] * ia, c[1] * ia, c[2] * ia, c[3] * ia, c[4] * ia, c[5] * ia, c[6] * ia, c[7] * ia, c[8] * ia);
}

/*!
\brief Overloaded operator.
\param a Real.
*/
inline QuadricSurface QuadricSurface::operator*(double a) const
{
  return QuadricSurface(c[0] * a, c[1] * a, c[2] * a, c[3] * a, c[4] * a, c[5] * a, c[6] * a, c[7] * a, c[8] * a);
}

/*!
\brief Overloaded operator.
\param q %Quadric surface.
*/
inline QuadricSurface QuadricSurface::operator+(const QuadricSurface& q) const
{
  return QuadricSurface(c[0] + q.c[0], c[1] + q.c[1], c[2] + q.c[2], c[3] + q.c[3], c[4] + q.c[4], c[5] + q.c[5], c[6] + q.c[6], c[7] + q.c[7], c[8] + q.c[8]);
}

/*!
\brief Overloaded operator.
\param q %Quadric surface.
*/
inline QuadricSurface QuadricSurface::operator-(const QuadricSurface& q) const
{
  return QuadricSurface(c[0] - q.c[0], c[1] - q.c[1], c[2] - q.c[2], c[3] - q.c[3], c[4] - q.c[4], c[5] - q.c[5], c[6] - q.c[6], c[7] - q.c[7], c[8] - q.c[8]);
}

/*!
\brief Return the coefficient.
\param i,j Coefficient indexes.
*/
inline double QuadricSurface::operator()(int i, int j) const
{
  return c[Index(i, j)];
}

/*!
\brief Return the coefficient.
\param i,j Coefficient indexes.
*/
inline double& QuadricSurface::operator()(int i, int j)
{
  return c[Index(i, j)];
}
