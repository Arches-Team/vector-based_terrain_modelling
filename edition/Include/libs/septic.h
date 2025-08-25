// Septic polynomials

#pragma once

#include "libs/sextic.h"

class Septic {
protected:
  double c[8] = { 0.0 }; //!< Array of coefficients.
public:
  //! Empty.
  Septic() {}
  explicit Septic(const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);
  Septic(const Sextic&);
  Septic(const Quintic&);
  Septic(const Quadric&);
  Septic(const Cubic&);
  Septic(const Quartic&);
  //! Empty.
  ~Septic() {}

  constexpr double& operator[] (int);
  constexpr double operator[] (int) const;

  int CheckDegree() const;

  // Unary operators
  Septic operator+ () const { return *this; }
  Septic operator- () const;

  Septic& operator+= (const Septic&);
  Septic& operator-= (const Septic&);
  Septic& operator*= (const double&);
  Septic& operator/= (const double&);

  // Binary operators
  friend Septic operator+ (const Septic&, const Septic&);
  friend Septic operator- (const Septic&, const Septic&);

  friend Septic operator* (const Septic&, const double&);
  friend Septic operator* (const double&, const Septic&);
  friend Septic operator/ (const Septic&, const double&);

  // Derivatives
  Sextic Prime() const;
  Quintic Second() const;
  Quartic Third() const;

  // Lipschitz
  double K(const double&, const double&) const;

  // Evaluates polynomial
  constexpr double operator()(const double&) const;
  constexpr double Derivative(const double&) const;

  // Solve
  int Solve(double*) const;
  int Solve(double*, const double&, const double&) const;

  void Range(double&, double&, const double& = 0.0, const double& = 1.0) const;

  static double Smooth(const double&);
  static double SmoothStep(const double&, const double&, const double&);
  static double Gaussian(const double&, const double&, const double&);

  // Hermite septic
  static Septic Hermite(const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);

  friend std::ostream& operator<<(std::ostream&, const Septic&);
public:
  static double epsilon; //!< Epsilon value used to check discriminant terms in the root finding process.
};

//! Creates a septic polynomial.
inline Septic::Septic(const double& a, const double& b, const double& c, const double& d, const double& e, const double& f, const double& g, const double& h) :c{ h,g,f,e,d,c,b,a }
{
}

/*!
\brief Creates a septic given a quadric.

Simply set higher coefficients to 0.
\param p Quadric.
*/
inline Septic::Septic(const Quadric& p) :c{ p[0],p[1],p[2],0.0,0.0,0.0,0.0,0.0 }
{
}

/*!
\brief Creates a septic given a cubic.

Simply set higher coefficients to 0.
\param p Cubic.
*/
inline Septic::Septic(const Cubic& p) :c{ p[0],p[1],p[2],p[3],0.0,0.0,0.0,0.0 }
{
}

/*!
\brief Creates a septic given a quartic.

\param p Quartic.
*/
inline Septic::Septic(const Quartic& p) :c{ p[0],p[1],p[2],p[3],p[4],0.0,0.0,0.0 }
{
}

/*!
\brief Creates a septic given a quintic.

Simply set higher coefficients to 0.
\param p Quintic.
*/
inline Septic::Septic(const Quintic& p) :c{ p[0],p[1],p[2],p[3],p[4],p[5],0.0,0.0 }
{
}

/*!
\brief Creates a septic given a sextic.

Simply set higher coefficients to 0.
\param p Quintic.
*/
inline Septic::Septic(const Sextic& p) :c{ p[0],p[1],p[2],p[3],p[4],p[5],p[6],0.0 }
{
}

//! Access to the coefficients of the septic.
inline constexpr double& Septic::operator[] (int i)
{
  return c[i];
}

//! Overloaded.
inline constexpr double Septic::operator[] (int i) const
{
  return c[i];
}

/*!
\brief Evaluates the septic.
\param x Argument value of the function.
*/
inline constexpr double Septic::operator()(const double& x) const
{
  return c[0] + x * (c[1] + x * (c[2] + x * (c[3] + x * (c[4] + x * (c[5] + x * (c[6] + x * c[7]))))));
}

/*!
\brief Computes the derivative of the septic.

This function is more efficient than using Prime() and
then evaluating the derivative for a given input value.
\param x Real.
*/
inline constexpr double Septic::Derivative(const double& x) const
{
  return c[1] + x * (2.0 * c[2] + x * (3.0 * c[3] + x * (4.0 * c[4] + x * (5.0 * c[5] + x * (6.0 * c[6] + x * 7.0 * c[7])))));
}

/*!
\brief Computes the first derivative of a septic, which is a quintic.
*/
inline Sextic Septic::Prime() const
{
  return Sextic(7.0 * c[7], 6.0 * c[6], 5.0 * c[5], 4.0 * c[4], 3.0 * c[3], 2.0 * c[2], c[1]);
}

/*!
\brief Computes the second derivative of a septic.
*/
inline Quintic Septic::Second() const
{
  return Quintic(42.0 * c[7], 30.0 * c[6], 20.0 * c[5], 12.0 * c[4], 6.0 * c[3], 2.0 * c[2]);
}

/*!
\brief Computes the third derivative of a septic.
*/
inline Quartic Septic::Third() const
{
  return Quartic(210.0 * c[7], 120.0 * c[6], 60.0 * c[5], 24.0 * c[4], 6.0 * c[3]);
}

//! Overloaded.
inline Septic operator* (const Septic& u, const double& e)
{
  return Septic(e * u[7], e * u[6], e * u[5], e * u[4], e * u[3], e * u[2], e * u[1], e * u[0]);
}

//! Overloaded.
inline Septic operator+ (const Septic& u, const Septic& v)
{
  return Septic(v[7] + u[7], v[6] + u[6], v[5] + u[5], v[4] + u[4], v[3] + u[3], v[2] + u[2], v[1] + u[1], v[0] + u[0]);
}

//! Overloaded.
inline Septic operator- (const Septic& v, const Septic& u)
{
  return Septic(v[7] - u[7], v[6] - u[6], v[5] - u[5], v[4] - u[4], v[3] - u[3], v[2] - u[2], v[1] - u[1], v[0] - u[0]);
}

/*!
\brief Overloaded.
*/
inline Septic operator*(const double& a, const Septic& p)
{
  return p * a;
}

/*!
\brief Overloaded.
*/
inline Septic operator/(const Septic& p, const double& a)
{
  return p * (1.0 / a);
}

/*!
\brief Unary.
*/
inline Septic Septic::operator- () const
{
  return Septic(-c[7], -c[6], -c[5], -c[4], -c[3], -c[2], -c[1], -c[0]);
}

/*!
\brief Compute the value of smooth C<SUP>3</SUP> interpolating function over unit interval.

The septic is defined as x<SUP>4</SUP>(-20 x<SUP>3</SUP> + 70 x<SUP>2</SUP> -84 x + 35).
Its first, second and third derivatives at 0.0 and 1.0 are 0.0.

The Lipschitz constant of the smooth septic over [0,1] is &lambda;=35/16.

\param x Argument in [0,1].
\sa Quintic::Smooth(), Cubic::Smooth()
*/
inline double Septic::Smooth(const double& x)
{
  return x * x * x * x * (x * (x * (x * -20.0 + 70.0) - 84.0) + 35.0);
}

/*!
\brief Compute a septic smooth step.

The code is slightly more efficient than:
\code
double y=Septic::Smooth(Linear::Step(x,a,b));
\endcode

\param x Input value.
\param a, b Interval values.
\sa Cubic::SmoothStep(), Quintic::SmoothStep()
*/
inline double Septic::SmoothStep(const double& x, const double& a, const double& b)
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
    return Septic::Smooth((x - a) / (b - a));
  }
}

/*!
\brief Compute a compactly supported Gaussian-like pulse.

The function has C<SUP>3</SUP> continuity.
\sa Cubic::Gaussian() and Quintic::Gaussian()

\param c Center.
\param r Radius.
\param x Value.
*/
inline double Septic::Gaussian(const double& c, const double& r, const double& x)
{
  double xc = fabs(x - c);
  if (xc > r) return 0.0;
  xc /= r;
  return Septic::Smooth(1.0 - xc);
}

