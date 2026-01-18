// Nonic polynomials

#pragma once

#include "libs/octic.h"

class Nonic {
protected:
  double c[10] = { 0.0 }; //!< Array of coefficients.
public:
  //! Empty.
  Nonic() {}
  explicit Nonic(const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);

  //! Empty.
  ~Nonic() {}

  constexpr double& operator[] (int);
  constexpr double operator[] (int) const;

  int CheckDegree() const;

  // Unary operators
  Nonic operator+ () const { return *this; }
  Nonic operator- () const;

  Nonic& operator+= (const Nonic&);
  Nonic& operator-= (const Nonic&);
  Nonic& operator*= (const double&);
  Nonic& operator/= (const double&);

  // Binary operators
  friend Nonic operator+ (const Nonic&, const Nonic&);
  friend Nonic operator- (const Nonic&, const Nonic&);

  friend Nonic operator* (const Nonic&, const double&);
  friend Nonic operator* (const double&, const Nonic&);
  friend Nonic operator/ (const Nonic&, const double&);

  int Solve(double*) const;

  // Composition
  static Nonic Compose(const Cubic&, const Cubic&);

  // Derivatives
  Octic Prime() const;
  Septic Second() const;
  Sextic Third() const;

  // Evaluates polynomial
  constexpr double operator()(const double&) const;
  constexpr double Derivative(const double&) const;

  static double Smooth(const double&);
  static double SmoothStep(const double&, const double&, const double&);
  static double Gaussian(const double&, const double&, const double&);

  friend std::ostream& operator<<(std::ostream&, const Nonic&);
public:
  static double epsilon; //!< Epsilon value used to check discriminant terms in the root finding process.
};

//! Creates a Nonic polynomial.
inline Nonic::Nonic(const double& a, const double& b, const double& c, const double& d, const double& e, const double& f, const double& g, const double& h, const double& i, const double& j) :c{ j,i,h,g,f,e,d,c,b,a }
{
}

//! Access to the coefficients of the sextic.
inline constexpr double& Nonic::operator[] (int i)
{
  return c[i];
}

//! Overloaded.
inline constexpr double Nonic::operator[] (int i) const
{
  return c[i];
}

/*!
\brief Evaluates the Nonic.
\param x Argument value of the function.
*/
inline constexpr double Nonic::operator()(const double& x) const
{
  return c[0] + x * (c[1] + x * (c[2] + x * (c[3] + x * (c[4] + x * (c[5] + x * (c[6] + x * (c[7] + x * (c[8] + x * c[9]))))))));
}

/*!
\brief Computes the derivative of the sextic.

This function is more efficient than using Prime() and
then evaluating the derivative for a given input value.
\param x Real.
*/
inline constexpr double Nonic::Derivative(const double& x) const
{
  return c[1] + x * (2.0 * c[2] + x * (3.0 * c[3] + x * (4.0 * c[4] + x * (5.0 * c[5] + x * (6.0 * c[6] + x * (7.0 * c[7] + x * (8.0 * c[8] + x * 9.0 * c[9])))))));
}

//! Overloaded.
inline Nonic operator* (const Nonic& u, const double& e)
{
  return Nonic(e * u[9], e * u[8], e * u[7], e * u[6], e * u[5], e * u[4], e * u[3], e * u[2], e * u[1], e * u[0]);
}

//! Overloaded.
inline Nonic operator+ (const Nonic& u, const Nonic& v)
{
  return Nonic(v[9] + u[9], v[8] + u[8], v[7] + u[7], v[6] + u[6], v[5] + u[5], v[4] + u[4], v[3] + u[3], v[2] + u[2], v[1] + u[1], v[0] + u[0]);
}

//! Overloaded.
inline Nonic operator- (const Nonic& v, const Nonic& u)
{
  return Nonic(v[9] - u[9], v[8] - u[8], v[7] - u[7], v[6] - u[6], v[5] - u[5], v[4] - u[4], v[3] - u[3], v[2] - u[2], v[1] - u[1], v[0] - u[0]);
}

/*!
\brief Overloaded.
*/
inline Nonic operator*(const double& a, const Nonic& p)
{
  return p * a;
}

/*!
\brief Overloaded.
*/
inline Nonic operator/(const Nonic& p, const double& a)
{
  return p * (1.0 / a);
}

/*!
\brief Unary.
*/
inline Nonic Nonic::operator- () const
{
  return Nonic(-c[9], -c[8], -c[7], -c[6], -c[5], -c[4], -c[3], -c[2], -c[1], -c[0]);
}

/*!
\brief Compute the value of smooth C<SUP>4</SUP> interpolating function over unit interval.

The Nonic is defined as x<SUP>4</SUP>(-20 x<SUP>3</SUP> + 70 x<SUP>2</SUP> -84 x + 35).
Its first four derivatives at 0.0 and 1.0 are 0.0.

The Lipschitz constant of the smooth nonic over [0,1] is &lambda;=315/128.

\param x Argument in [0,1].
\sa Septic::Smooth(), Quintic::Smooth(), Cubic::Smooth()
*/
inline double Nonic::Smooth(const double& x)
{
  return x * x * x * x * x * (x * (x * (x * (x * 70.0 - 315.0) - 540.0) - 420.0) + 126);
}

/*!
\brief Compute a Nonic smooth step.

The code is slightly more efficient than:
\code
double y=Nonic::Smooth(Linear::Step(x,a,b));
\endcode

\param x Input value.
\param a, b Interval values.
\sa Septic::SmoothStep(), Cubic::SmoothStep(), Quintic::SmoothStep()
*/
inline double Nonic::SmoothStep(const double& x, const double& a, const double& b)
{
  if (x < a)
  {
    return 0.0;
  }
  else if (x > b)
  {
    return 1.0;
  }
  else
  {
    return Nonic::Smooth((x - a) / (b - a));
  }
}

/*!
\brief Compute a compactly supported Gaussian-like pulse.

The function has C<SUP>4</SUP> continuity.
\sa Cubic::Gaussian(), Quintic::Gaussian() and Septic::Gaussian()

\param c Center.
\param r Radius.
\param x Value.
*/
inline double Nonic::Gaussian(const double& c, const double& r, const double& x)
{
  double xc = fabs(x - c);
  if (xc > r) return 0.0;
  xc /= r;
  return Nonic::Smooth(1.0 - xc);
}


/*!
\brief Computes the first derivative of a septic, which is a quintic.
*/
inline Octic Nonic::Prime() const
{
  return Octic(9.0 * c[9], 8.0 * c[8], 7.0 * c[7], 6.0 * c[6], 5.0 * c[5], 4.0 * c[4], 3.0 * c[3], 2.0 * c[2], c[1]);
}

/*!
\brief Computes the second derivative of a Nonic.
*/
inline Septic Nonic::Second() const
{
  return Septic(72.0 * c[9], 56.0 * c[8], 42.0 * c[7], 30.0 * c[6], 20.0 * c[5], 12.0 * c[4], 6.0 * c[3], 2.0 * c[2]);
}

/*!
\brief Computes the third derivative of a Nonic.
*/
inline Sextic Nonic::Third() const
{
  return Sextic(504.0 * c[9], 336.0 * c[8], 210.0 * c[7], 120.0 * c[6], 60.0 * c[5], 25.0 * c[4], 6.0 * c[3]);
}
