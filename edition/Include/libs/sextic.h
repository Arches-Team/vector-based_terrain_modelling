// Sextic polynomials

#pragma once

#include "libs/quintic.h"

class Sextic {
protected:
  double c[7] = { 0.0,0.0,0.0,0.0,0.0,0.0,0.0 }; //!< Array of coefficients.
public:
  //! Empty.
  Sextic() {}
  explicit Sextic(const double&, const double&, const double&, const double&, const double&, const double&, const double&);
  Sextic(const Quintic&);
  Sextic(const Quadric&);
  Sextic(const Cubic&);
  Sextic(const Quartic&);
  //! Empty.
  ~Sextic() {}

  constexpr double& operator[] (int);
  constexpr double operator[] (int) const;

  int CheckDegree() const;

  // Unary operators
  Sextic operator+ () const { return *this; }
  Sextic operator- () const;

  Sextic& operator+= (const Sextic&);
  Sextic& operator-= (const Sextic&);
  Sextic& operator*= (const double&);
  Sextic& operator/= (const double&);

  // Binary operators
  friend Sextic operator+ (const Sextic&, const Sextic&);
  friend Sextic operator- (const Sextic&, const Sextic&);

  friend Sextic operator* (const Sextic&, const double&);
  friend Sextic operator* (const double&, const Sextic&);
  friend Sextic operator/ (const Sextic&, const double&);

  friend Sextic operator* (const Quartic&, const Quadric&);
  friend Sextic operator* (const Cubic&, const Cubic&);

  // Derivatives
  Quintic Prime() const;
  Quartic Second() const;
  Cubic Third() const;

  // Lipschitz
  double K(const double&, const double&) const;

  // Evaluates polynomial
  constexpr double operator()(const double&) const;
  constexpr double Derivative(const double&) const;

  // Solve
  int Solve(double*) const;
  int Solve(double*, const double&, const double&) const;
  void Range(double&, double&, const double& = 0.0, const double& = 1.0) const;

  static Sextic Compose(const Quadric&, const Cubic&);
  static Sextic Compose(const Cubic&, const Quadric&);

  friend std::ostream& operator<<(std::ostream&, const Sextic&);
public:
  static double epsilon; //!< Epsilon value used to check discriminant terms in the root finding process.
};

//! Creates a sextic.
inline Sextic::Sextic(const double& a, const double& b, const double& c, const double& d, const double& e, const double& f, const double& g) :c{ g,f,e,d,c,b,a }
{
}

/*!
\brief Creates a sextic given a quadric.

Simply set higher coefficients to 0.
\param p %Quadric.
*/
inline Sextic::Sextic(const Quadric& p) :c{ p[0],p[1],p[2],0.0,0.0,0.0,0.0 }
{
}

/*!
\brief Creates a sextic given a cubic.

Simply set higher coefficients to 0.
\param p %Cubic.
*/
inline Sextic::Sextic(const Cubic& p) :c{ p[0],p[1],p[2],p[3],0.0,0.0,0.0 }
{
}

/*!
\brief Creates a sextic given a quartic.

Simply set highest coefficient to 0.
\param p %Quartic.
*/
inline Sextic::Sextic(const Quartic& p) :c{ p[0],p[1],p[2],p[3],p[4],0.0,0.0 }
{
}

/*!
\brief Creates a sextic given a quintic.

Simply set higher coefficients to 0.
\param p %Quintic.
*/
inline Sextic::Sextic(const Quintic& p) :c{ p[0],p[1],p[2],p[3],p[4],p[5],0.0 }
{
}

//! Access to the coefficients of the sextic.
inline constexpr double& Sextic::operator[] (int i)
{
  return c[i];
}

//! Overloaded.
inline constexpr double Sextic::operator[] (int i) const
{
  return c[i];
}

/*!
\brief Evaluates the sextic.
\param x Argument value of the function.
*/
inline constexpr double Sextic::operator()(const double& x) const
{
  return c[0] + x * (c[1] + x * (c[2] + x * (c[3] + x * (c[4] + x * (c[5] + x * c[6])))));
}

/*!
\brief Computes the derivative of the sextic.

This function is more efficient than using Prime() and
then evaluating the derivative for a given input value.
\param x Real.
*/
inline constexpr double Sextic::Derivative(const double& x) const
{
  return c[1] + x * (2.0 * c[2] + x * (3.0 * c[3] + x * (4.0 * c[4] + x * (5.0 * c[5] + x * 6.0 * c[6]))));
}

/*!
\brief Computes the first derivative of a sextic, which is a quintic.
*/
inline Quintic Sextic::Prime() const
{
  return Quintic(6.0 * c[6], 5.0 * c[5], 4.0 * c[4], 3.0 * c[3], 2.0 * c[2], c[1]);
}

/*!
\brief Computes the second derivative of a sextic, which is a quartic.
*/
inline Quartic Sextic::Second() const
{
  return Quartic(30.0 * c[6], 20.0 * c[5], 12.0 * c[4], 6.0 * c[3], 2.0 * c[2]);
}

/*!
\brief Computes the third derivative of a sextic.
*/
inline Cubic Sextic::Third() const
{
  return Cubic(120.0 * c[6], 60.0 * c[5], 24.0 * c[4], 6.0 * c[3]);
}

//! Overloaded.
inline Sextic operator* (const Sextic& u, const double& e)
{
  return Sextic(e * u[6], e * u[5], e * u[4], e * u[3], e * u[2], e * u[1], e * u[0]);
}

//! Overloaded.
inline Sextic operator+ (const Sextic& u, const Sextic& v)
{
  return Sextic(v[6] + u[6], v[5] + u[5], v[4] + u[4], v[3] + u[3], v[2] + u[2], v[1] + u[1], v[0] + u[0]);
}

//! Overloaded.
inline Sextic operator- (const Sextic& v, const Sextic& u)
{
  return Sextic(v[6] - u[6], v[5] - u[5], v[4] - u[4], v[3] - u[3], v[2] - u[2], v[1] - u[1], v[0] - u[0]);
}

/*!
\brief Overloaded.
*/
inline Sextic operator*(const double& a, const Sextic& p)
{
  return p * a;
}

/*!
\brief Overloaded.
*/
inline Sextic operator/(const Sextic& p, const double& a)
{
  return p * (1.0 / a);
}

/*!
\brief Multiply two cubics which generates a sextic.

\param u, v Arguments.
*/
inline Sextic operator*(const Cubic& u, const Cubic& v)
{
  return Sextic(u[3] * v[3], u[3] * v[2] + u[2] * v[3], u[3] * v[1] + v[2] * u[2] + u[1] * v[2], u[3] * v[0] + u[2] * v[1] + u[1] * v[2] + u[0] * v[3], u[2] * v[0] + u[1] * v[1] + u[0] * v[2], u[1] * v[0] + u[0] * v[1], u[0] * v[0]);
}

/*!
\brief Multiply a quartic by a quadric.

\param u, v Arguments.
*/
inline Sextic operator*(const Quartic& u, const Quadric& v)
{
  return Sextic(u[4] * v[2], u[4] * v[1] + u[3] * v[2], u[4] * v[0] + v[3] * u[1] + u[2] * v[2], u[3] * v[0] + u[2] * v[1] + u[1] * v[2], u[2] * v[0] + u[1] * v[1] + u[0] * v[2], u[1] * v[0] + u[0] * v[1], u[0] * v[0]);
}

/*!
\brief Unary.
*/
inline Sextic Sextic::operator- () const
{
  return Sextic(-c[6], -c[5], -c[4], -c[3], -c[2], -c[1], -c[0]);
}

