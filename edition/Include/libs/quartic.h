// Quartic polynomials

#pragma once

#include "libs/cubic.h"

class Quartic {
protected:
  double c[5] = { 0.0,0.0,0.0,0.0,0.0 }; //!< Array of coefficients.
public:
  Quartic() {} //! Empty.
  explicit Quartic(const double&, const double&, const double&, const double&, const double&);
  Quartic(const Quadric&);
  Quartic(const Cubic&);
  ~Quartic() {} //! Empty.

  constexpr double& operator[] (int);
  constexpr double operator[] (int) const;

  int CheckDegree() const;

  // Unary operators
  Quartic operator+ () const { return *this; }
  Quartic operator- () const;

  Quartic& operator+= (const Quartic&);
  Quartic& operator-= (const Quartic&);
  Quartic& operator*= (const double&);
  Quartic& operator/= (const double&);

  // Binary operators
  friend Quartic operator+ (const Quartic&, const Quartic&);
  friend Quartic operator- (const Quartic&, const Quartic&);

  friend Quartic operator* (const Quartic&, const double&);
  friend Quartic operator* (const double&, const Quartic&);
  friend Quartic operator/ (const Quartic&, const double&);

  friend Quartic operator* (const Quadric&, const Quadric&);

  // Derivatives
  Cubic Prime() const;
  Quadric Second() const;
  Linear Third() const;

  // Lipschitz
  double K(const double&, const double&) const;

  // Composition
  static Quartic Compose(const Quadric&, const Quadric&);

  // Evaluates polynomial
  constexpr double operator()(const double&) const;
  constexpr double Derivative(const double&) const;

  // Solve
  int Solve(double*);
  int Solve(double*, const double&, const double&);

  void Range(double&, double&, const double& = 0.0, const double& = 1.0) const;

  friend std::ostream& operator<<(std::ostream&, const Quartic&);
protected:
  bool Analyze();
public:
  static double epsilon; //!< Epsilon value used to check discriminant terms in the root finding process.
private:
  static double limit; //!< Epsilon value used to limit coefficients.
};

//! Creates a quartic.
inline Quartic::Quartic(const double& a, const double& b, const double& c, const double& d, const double& e)
{
  Quartic::c[4] = a;
  Quartic::c[3] = b;
  Quartic::c[2] = c;
  Quartic::c[1] = d;
  Quartic::c[0] = e;
}

//! Creates a quartic given a quadric. 
inline Quartic::Quartic(const Quadric& p)
{
  c[4] = c[3] = 0.0;
  c[2] = p[2];
  c[1] = p[1];
  c[0] = p[0];
}

//! Creates a quartic given a cubic.
inline Quartic::Quartic(const Cubic& p)
{
  c[4] = 0.0;
  c[3] = p[3];
  c[2] = p[2];
  c[1] = p[1];
  c[0] = p[0];
}

//! Access to the coefficients of the quartic.
inline constexpr double& Quartic::operator[] (int i)
{
  return c[i];
}

//! Overloaded.
inline constexpr double Quartic::operator[] (int i) const
{
  return c[i];
}

/*!
\brief Evaluates the quartic.
\param x Argument value of the function.
*/
inline constexpr double Quartic::operator()(const double& x) const
{
  return c[0] + x * (c[1] + x * (c[2] + x * (c[3] + x * c[4])));
}

/*!
\brief Computes the derivative value at a given point.

If only one evaluation is to be made, this function is more efficient than Prime() and
then evaluating the derivative for a given input value.

\param x Real.

\sa Quartic::Prime
*/
inline constexpr double Quartic::Derivative(const double& x) const
{
  return c[1] + x * (2.0 * c[2] + x * (3.0 * c[3] + x * 4.0 * c[4]));
}

//! Computes the first derivative of a cubic, which is a cubic.
inline Cubic Quartic::Prime() const
{
  return Cubic(4.0 * c[4], 3.0 * c[3], 2.0 * c[2], c[1]);
}

/*!
\brief Computes the second derivative of a quartic, which is a quadric.
*/
inline Quadric Quartic::Second() const
{
  return Quadric(12.0 * c[4], 6.0 * c[3], 2.0 * c[2]);
}
/*!
\brief Computes the third derivative of a quartic, which is a linear.
*/
inline Linear Quartic::Third() const
{
  return Linear(24.0 * c[4], 6.0 * c[3]);
}

//! Overloaded.
inline Quartic operator* (const Quartic& u, const double& e)
{
  return Quartic(e * u[4], e * u[3], e * u[2], e * u[1], e * u[0]);
}

//! Overloaded.
inline Quartic operator+ (const Quartic& u, const Quartic& v)
{
  return Quartic(v[4] + u[4], v[3] + u[3], v[2] + u[2], v[1] + u[1], v[0] + u[0]);
}

//! Overloaded.
inline Quartic operator- (const Quartic& v, const Quartic& u)
{
  return Quartic(v[4] - u[4], v[3] - u[3], v[2] - u[2], v[1] - u[1], v[0] - u[0]);
}

/*!
\brief Multiply two quadrics which generates a quartic.

\param u, v Argument quadrics.
*/
inline Quartic operator*(const Quadric& u, const Quadric& v)
{
  return Quartic(
    u[2] * v[2],
    u[2] * v[1] + u[1] * v[2],
    u[2] * v[0] + v[2] * u[0] + u[1] * v[1],
    u[1] * v[0] + u[0] * v[1],
    u[0] * v[0]);
}

/*!
\brief Overloaded.
*/
inline Quartic operator*(const double& a, const Quartic& p)
{
  return p * a;
}

/*!
\brief Overloaded.
*/
inline Quartic operator/(const Quartic& p, const double& a)
{
  return p * (1.0 / a);
}

/*!
\brief Unary.
*/
inline Quartic Quartic::operator- () const
{
  return Quartic(-c[4], -c[3], -c[2], -c[1], -c[0]);
}
