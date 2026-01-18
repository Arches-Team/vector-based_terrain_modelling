// Quintics

#pragma once

#include "libs/quartic.h"

class Quintic {
protected:
  double c[6] = { 0.0,0.0,0.0,0.0,0.0,0.0 }; //!< Array of coefficients.
public:
  //! Empty.
  Quintic() {}
  explicit Quintic(const double&, const double&, const double&, const double&, const double&, const double&);
  Quintic(const Quadric&);
  Quintic(const Cubic&);
  Quintic(const Quartic&);
  //! Empty.
  ~Quintic() {}

  constexpr double& operator[] (int);
  constexpr double operator[] (int) const;

  int CheckDegree() const;

  // Unary operators
  Quintic operator+ () const { return *this; }
  Quintic operator- () const;

  Quintic& operator+= (const Quintic&);
  Quintic& operator-= (const Quintic&);
  Quintic& operator*= (const double&);
  Quintic& operator/= (const double&);

  // Binary operators
  friend Quintic operator+ (const Quintic&, const Quintic&);
  friend Quintic operator- (const Quintic&, const Quintic&);

  friend Quintic operator* (const Quintic&, const double&);
  friend Quintic operator* (const double&, const Quintic&);
  friend Quintic operator/ (const Quintic&, const double&);

  // Derivatives
  Quartic Prime() const;
  Cubic Second() const;
  Quadric Third() const;

  // Lipschitz
  double K(const double&, const double&) const;

  // Evaluates polynomial
  constexpr double operator()(const double&) const;
  constexpr double Derivative(const double&) const;

  // Hermite quintinc
  static Quintic Hermite(const double&, const double&, const double&, const double&, const double&, const double&);

  // Solve
  int Solve(double*) const;
  int Solve(double*, const double&, const double&) const;

  void Range(double&, double&, const double& = 0.0, const double& = 1.0) const;

  static double Smooth(const double&);
  static double SmoothStep(const double&, const double&, const double&);
  static double Gaussian(const double&, const double&, const double&);
  static double GaussianThick(const double&, const double&, const double&, const double&);

  static Quintic FromRoots(const QVector<double>&);

  friend std::ostream& operator<<(std::ostream&, const Quintic&);
public:
  static double epsilon; //!< \htmlonly\epsilon;\endhtmlonly value used to check discriminant terms in the root finding process.
};

//! Creates a quintic polynomial.
inline Quintic::Quintic(const double& a, const double& b, const double& c, const double& d, const double& e, const double& f)
{
  Quintic::c[5] = a;
  Quintic::c[4] = b;
  Quintic::c[3] = c;
  Quintic::c[2] = d;
  Quintic::c[1] = e;
  Quintic::c[0] = f;
}

/*!
\brief Creates a quintic from a quadric.

Simply set higher coefficients to 0.
\param p Quadric.
*/
inline Quintic::Quintic(const Quadric& p)
{
  c[5] = c[4] = c[3] = 0.0;
  c[2] = p[2];
  c[1] = p[1];
  c[0] = p[0];
}

/*!
\brief Creates a quintic from a cubic.

Simply set higher coefficients to 0.
\param p Cubic.
*/
inline Quintic::Quintic(const Cubic& p)
{
  c[5] = c[4] = 0.0;
  c[3] = p[3];
  c[2] = p[2];
  c[1] = p[1];
  c[0] = p[0];
}

/*!
\brief Creates a quintic from a quartic.

\param p Quartic.
*/
inline Quintic::Quintic(const Quartic& p)
{
  c[5] = 0.0;
  c[4] = p[4];
  c[3] = p[3];
  c[2] = p[2];
  c[1] = p[1];
  c[0] = p[0];
}

//! Access to the coefficients of the quintic.
inline constexpr double& Quintic::operator[] (int i)
{
  return c[i];
}

//! Overloaded.
inline constexpr double Quintic::operator[] (int i) const
{
  return c[i];
}

/*!
\brief Evaluates the quintic.
\param x Argument value of the function.
*/
inline constexpr double Quintic::operator()(const double& x) const
{
  return c[0] + x * (c[1] + x * (c[2] + x * (c[3] + x * (c[4] + x * c[5]))));
}

/*!
\brief Computes the derivative value at a given point.

If only one evaluation is to be made, this function is more efficient than Prime() and
then evaluating the derivative for a given input value.

\param x Real.

\sa Quintic::Prime
*/
inline constexpr double Quintic::Derivative(const double& x) const
{
  return c[1] + x * (2.0 * c[2] + x * (3.0 * c[3] + x * (4.0 * c[4] + x * c[5])));
}

/*!
\brief Computes the first derivative of a quintic, which is a quartic.
*/
inline Quartic Quintic::Prime() const
{
  return Quartic(5.0 * c[5], 4.0 * c[4], 3.0 * c[3], 2.0 * c[2], c[1]);
}

/*!
\brief Computes the second derivative of a quintic, which is a cubic.
*/
inline Cubic Quintic::Second() const
{
  return Cubic(20.0 * c[5], 12.0 * c[4], 6.0 * c[3], 2.0 * c[2]);
}

/*!
\brief Computes the third derivative of a quintic.
*/
inline Quadric Quintic::Third() const
{
  return Quadric(60.0 * c[5], 24.0 * c[4], 6.0 * c[3]);
}

//! Overloaded.
inline Quintic operator* (const Quintic& u, const double& e)
{
  return Quintic(e * u[5], e * u[4], e * u[3], e * u[2], e * u[1], e * u[0]);
}

//! Overloaded.
inline Quintic operator+ (const Quintic& u, const Quintic& v)
{
  return Quintic(v[5] + u[5], v[4] + u[4], v[3] + u[3], v[2] + u[2], v[1] + u[1], v[0] + u[0]);
}

//! Overloaded.
inline Quintic operator- (const Quintic& v, const Quintic& u)
{
  return Quintic(v[5] - u[5], v[4] - u[4], v[3] - u[3], v[2] - u[2], v[1] - u[1], v[0] - u[0]);
}

/*!
\brief Overloaded.
*/
inline Quintic operator*(const double& a, const Quintic& p)
{
  return p * a;
}

/*!
\brief Overloaded.
*/
inline Quintic operator/(const Quintic& p, const double& a)
{
  return p * (1.0 / a);
}

/*!
\brief Unary.
*/
inline Quintic Quintic::operator- () const
{
  return Quintic(-c[5], -c[4], -c[3], -c[2], -c[1], -c[0]);
}

/*!
\brief Compute the value of smooth C<SUP>2</SUP> interpolating function over unit interval.

The quintic is defined as x<SUP>3</SUP>(6 x<SUP>2</SUP>-15 x + 10).
Its first and second derivatives at 0.0 and 1.0 are 0.0.

The Lipschitz constant of the smooth quintic over [0,1] is &lambda;=15/8.

\param x Argument in [0,1].
\sa Septic::Smooth(), Cubic::Smooth()
*/
inline double Quintic::Smooth(const double& x)
{
  return x * x * x * (x * (x * 6.0 - 15.0) + 10.0);
}

/*!
\brief Compute a quintic smooth step.

The code is slightly more efficient than:
\code
double y=Quintic::Smooth(Linear::Step(x,a,b));
\endcode

\param x Input value.
\param a, b Interval values.
\sa Cubic::SmoothStep, Septic::SmoothStep()
*/
inline double Quintic::SmoothStep(const double& x, const double& a, const double& b)
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
    return Quintic::Smooth((x - a) / (b - a));
  }
}

/*!
\brief Compute a compactly supported Gaussian-like pulse.

The function has C<SUP>2</SUP> continuity.
\sa Cubic::Gaussian()

\param c Center.
\param r Radius.
\param x Value.
*/
inline double Quintic::Gaussian(const double& c, const double& r, const double& x)
{
  double xc = fabs(x - c);
  if (xc > r) return 0.0;
  xc /= r;
  return Quintic::Smooth(1.0 - xc);
}

/*!
\brief Compute a compactly supported Gaussian-like pulse with a thick plateau.

\sa Cubic::GaussianThick()
\sa Quintic::Gaussian()

\param c Center.
\param t,r Thickness (plateau) and radius.
\param x Value.
*/
inline double Quintic::GaussianThick(const double& c, const double& t, const double& r, const double& x)
{
  double y = fabs(x - c);
  if (y > r + t) return 0.0;
  y -= t;
  if (y < 0.0) return 1.0;
  y /= r;
  return Quintic::Smooth(1.0 - y);
}
