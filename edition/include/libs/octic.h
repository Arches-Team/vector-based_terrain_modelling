// Octic polynomials

#pragma once


#include "libs/septic.h"

class Octic {
protected:
  double c[9] = { 0.0 }; //!< Array of coefficients.
public:
  //! Empty.
  Octic() {}
  explicit Octic(const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);

  //! Empty.
  ~Octic() {}

  constexpr double& operator[] (int);
  constexpr double operator[] (int) const;

  int CheckDegree() const;

  // Unary operators
  Octic operator+ () const { return *this; }
  Octic operator- () const;

  Octic& operator+= (const Octic&);
  Octic& operator-= (const Octic&);
  Octic& operator*= (const double&);
  Octic& operator/= (const double&);

  // Binary operators
  friend Octic operator+ (const Octic&, const Octic&);
  friend Octic operator- (const Octic&, const Octic&);

  friend Octic operator* (const Octic&, const double&);
  friend Octic operator* (const double&, const Octic&);
  friend Octic operator/ (const Octic&, const double&);

  friend Octic operator* (const Quartic&, const Quartic&);
  friend Octic operator* (const Sextic&, const Quadric&);

  int Solve(double*) const;

  // Derivatives
  Septic Prime() const;
  Sextic Second() const;
  Quintic Third() const;

  // Evaluates polynomial
  constexpr double operator()(const double&) const;
  constexpr double Derivative(const double&) const;

  friend std::ostream& operator<<(std::ostream&, const Octic&);
public:
  static double epsilon; //!< Epsilon value used to check discriminant terms in the root finding process.
};

//! Creates a Octic polynomial.
inline Octic::Octic(const double& a, const double& b, const double& c, const double& d, const double& e, const double& f, const double& g, const double& h, const double& i) :c{ i,h,g,f,e,d,c,b,a }
{
}

//! Access to the coefficients of the sextic.
inline constexpr double& Octic::operator[] (int i)
{
  return c[i];
}

//! Overloaded.
inline constexpr double Octic::operator[] (int i) const
{
  return c[i];
}

/*!
\brief Evaluates the Octic.
\param x Argument value of the function.
*/
inline constexpr double Octic::operator()(const double& x) const
{
  return c[0] + x * (c[1] + x * (c[2] + x * (c[3] + x * (c[4] + x * (c[5] + x * (c[6] + x * (c[7] + x * c[8])))))));
}

/*!
\brief Computes the derivative of the sextic.

This function is more efficient than using Prime() and
then evaluating the derivative for a given input value.
\param x Real.
*/
inline constexpr double Octic::Derivative(const double& x) const
{
  return c[1] + x * (2.0 * c[2] + x * (3.0 * c[3] + x * (4.0 * c[4] + x * (5.0 * c[5] + x * (6.0 * c[6] + x * (7.0 * c[7] + x * 8.0 * c[8]))))));
}

/*!
\brief Computes the first derivative of a Octic.
*/
inline Septic Octic::Prime() const
{
  return Septic(8.0 * c[8], 7.0 * c[7], 6.0 * c[6], 5.0 * c[5], 4.0 * c[4], 3.0 * c[3], 2.0 * c[2], c[1]);
}

/*!
\brief Computes the second derivative of a Octic.
*/
inline Sextic Octic::Second() const
{
  return Sextic(56.0 * c[8], 42.0 * c[7], 30.0 * c[6], 20.0 * c[5], 12.0 * c[4], 6.0 * c[3], 2.0 * c[2]);
}

/*!
\brief Computes the third derivative of a Octic.
*/
inline Quintic Octic::Third() const
{
  return Quintic(336.0 * c[8], 210.0 * c[7], 120.0 * c[6], 60.0 * c[5], 25.0 * c[4], 6.0 * c[3]);
}


//! Overloaded.
inline Octic operator* (const Octic& u, const double& e)
{
  return Octic(e * u[8], e * u[7], e * u[6], e * u[5], e * u[4], e * u[3], e * u[2], e * u[1], e * u[0]);
}

//! Overloaded.
inline Octic operator+ (const Octic& u, const Octic& v)
{
  return Octic(v[8] + u[8], v[7] + u[7], v[6] + u[6], v[5] + u[5], v[4] + u[4], v[3] + u[3], v[2] + u[2], v[1] + u[1], v[0] + u[0]);
}

//! Overloaded.
inline Octic operator- (const Octic& v, const Octic& u)
{
  return Octic(v[8] - u[8], v[7] - u[7], v[6] - u[6], v[5] - u[5], v[4] - u[4], v[3] - u[3], v[2] - u[2], v[1] - u[1], v[0] - u[0]);
}

/*!
\brief Overloaded.
*/
inline Octic operator*(const double& a, const Octic& p)
{
  return p * a;
}

/*!
\brief Overloaded.
*/
inline Octic operator/(const Octic& p, const double& a)
{
  return p * (1.0 / a);
}

/*!
\brief Unary.
*/
inline Octic Octic::operator- () const
{
  return Octic(-c[8], -c[7], -c[6], -c[5], -c[4], -c[3], -c[2], -c[1], -c[0]);
}
