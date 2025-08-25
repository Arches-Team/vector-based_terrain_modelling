// Quadric polynomials

#pragma once

#include "libs/linear.h"
#include <QtCore/QVector>

class Quadric {
protected:
  double c[3] = { 0.0,0.0,0.0 }; //!< Array of coefficients.
public:
  //! Empty.
  Quadric() {}
  explicit Quadric(const double&, const double&, const double&);
  explicit Quadric(const Linear&);
  //! Empty.
  ~Quadric() {}

  constexpr double& operator[] (int);
  constexpr double operator[] (int) const;

  // Unary operators
  Quadric operator+ () const { return *this; }
  Quadric operator- () const;

  Quadric& operator+= (const Quadric&);
  Quadric& operator-= (const Quadric&);
  Quadric& operator*= (const double&);
  Quadric& operator/= (const double&);

  // Binary operators
  friend Quadric operator+ (const Quadric&, const Quadric&);
  friend Quadric operator- (const Quadric&, const Quadric&);

  friend Quadric operator* (const Quadric&, const double&);
  friend Quadric operator* (const double&, const Quadric&);
  friend Quadric operator/ (const Quadric&, const double&);
  friend Quadric operator* (const Linear&, const Linear&);

  // Evaluates polynomial
  constexpr double operator()(const double&) const;
  constexpr double Derivative(const double&) const;

  // Derivative
  Linear Prime() const;

  // Solve
  int Solve(double&, double&) const;
  int Solve(double*) const;
  int Solve(double*, const double&, const double&) const;

  void Range(double&, double&, const double& = 0.0, const double& = 1.0) const;
  double Maximum(const double& = 0.0, const double& = 1.0) const;
  double Minimum(const double& = 0.0, const double& = 1.0) const;

  // Lipschitz
  double K(const double&, const double&) const;

  static Quadric Bezier(const double&, const double&, const double&);

  static double Smooth(const double&, const double&);
  static double Smooth(const double&);
  static double SmoothCompact(const double&, const double&);
  static double SmoothCompact(const double&, const double&, const double&);

  static double Warp(const double&);

  // Bernstein
  static Quadric Bernstein(int);
  static double Bernstein(int, const double&);

  static Quadric FromRoots(const QVector<double>&);

  friend std::ostream& operator<<(std::ostream&, const Quadric&);
public:
  static double epsilon; //!< Epsilon value used to check b<SUP>2</SUP>-4ac term in the root finding process.
};

//! Creates a quadric.
inline Quadric::Quadric(const double& a, const double& b, const double& c) :c{c,b,a}
{
}

//! Creates a quadric from a Linear polynomial.
inline Quadric::Quadric(const Linear& l) : c{ l[0],l[1],0.0 }
{
}

//! Access class components.
inline constexpr double& Quadric::operator[] (int i)
{
  return c[i];
}

//! Overloaded.
inline constexpr double Quadric::operator[] (int i) const
{
  return c[i];
}

/*!
\brief Compute the derivative of a quadric, which is a linear expression.

\sa Quadric::Derivative
*/
inline Linear Quadric::Prime() const
{
  return Linear(2.0 * c[2], c[1]);
}

//! Evaluates the quadric. 
inline constexpr double Quadric::operator()(const double& x) const
{
  return c[0] + x * (c[1] + x * c[2]);
}

/*!
\brief Computes the derivative value at a given point.

If only one evaluation is to be made, this function is more efficient than Prime() and
then evaluating the derivative for a given input value.

\param x Real.

\sa Quadric::Prime
*/
inline constexpr double Quadric::Derivative(const double& x) const
{
  return 2.0 * c[2] * x + c[1];
}

/*!
\brief Destructive sum of two polynomials.
*/
inline Quadric& Quadric::operator+=(const Quadric& u)
{
  c[0] += u[0];
  c[1] += u[1];
  c[2] += u[2];

  return *this;
}

/*!
\brief Destructive difference of two polynomials.
*/
inline Quadric& Quadric::operator-=(const Quadric& u)
{
  c[0] -= u[0];
  c[1] -= u[1];
  c[2] -= u[2];

  return *this;
}

/*!
\brief Scale a polynomial by a double value.
*/
inline Quadric& Quadric::operator*=(const double& e)
{
  c[0] *= e;
  c[1] *= e;
  c[2] *= e;
  return *this;
}

/*!
\brief Scale a polynomial by a double value.
*/
inline Quadric& Quadric::operator/=(const double& e)
{
  c[0] /= e;
  c[1] /= e;
  c[2] /= e;
  return *this;
}

/*!
\brief Multiply a quadric by a scalar value.
*/
inline Quadric operator*(const Quadric& p, const double& e)
{
  return Quadric(e * p[2], e * p[1], e * p[0]);
}

//! Overloaded.
inline Quadric operator+(const Quadric& u, const Quadric& v)
{
  return Quadric(v[2] + u[2], v[1] + u[1], v[0] + u[0]);
}

//! Overloaded.
inline Quadric operator-(const Quadric& v, const Quadric& u)
{
  return Quadric(v[2] - u[2], v[1] - u[1], v[0] - u[0]);
}

/*!
\brief Multiply two linear polynomials.
*/
inline Quadric operator*(const Linear& u, const Linear& v)
{
  return Quadric(u[1] * v[1], u[0] * v[1] + u[1] * v[0], u[0] * v[0]);
}

//! Overloaded
inline Quadric Quadric::operator- () const
{
  return Quadric(-c[2], -c[1], -c[0]);
}

/*!
\brief Overloaded.
*/
inline Quadric operator*(const double& a, const Quadric& p)
{
  return Quadric(a * p[2], a * p[1], a * p[0]);
}

/*!
\brief Overloaded.
*/
inline Quadric operator/(const Quadric& p, const double& a)
{
  return p * (1.0 / a);
}

/*!
\brief Compute the value of the C<SUP>1</SUP> smooth interpolating function (1-x/r)<SUP>2</SUP>.

\sa Cubic::Smooth(), Quadric::SmoothCompact

\param x Argument value.
\param r Radius.
*/
inline double Quadric::Smooth(const double& x, const double& r)
{
  return Math::Sqr(1.0 - x / r);
}

/*!
\brief Compute the value of the C<SUP>1</SUP> smooth interpolating function (1-x)<SUP>2</SUP>.

\sa Quadric::Smooth(const double&, const double&)

\param x Argument value.
*/
inline double Quadric::Smooth(const double& x)
{
  return Math::Sqr(1.0 - x);
}

/*!
\brief Compact support version of the smooth interpolating function Quadric::Smooth().

Return Quadric::Smooth(x/r) if x<r and 0 otherwise.

\sa Cubic::Smooth(), Quadric::Smooth()

\param x Argument value.
\param r Radius.
*/
inline double Quadric::SmoothCompact(const double& x, const double& r)
{
  return (x > r) ? 0.0 : Math::Sqr(1.0 - x / r);
}

/*!
\brief Compact support version of the smooth interpolating function Quadric::Smooth() with interior and exterior radii.

\sa Cubic::Smooth(), Quadric::Smooth()

\param x Argument value.
\param i, e Interior and exterior radii.
*/
inline double Quadric::SmoothCompact(const double& x, const double& e, const double& i)
{
  if (x > e)
  {
    return 0.0;
  }
  else if (x < i)
  {
    return 1.0;
  }
  else
  {
    return Quadric::Smooth(Math::Sqr((x - i) / (e - i)));
  }
}

/*!
\brief Unit interval warping function.

Remaps the unit interval into the unit interval by expanding the sides and compressing the center, keeping 1/2 mapped to 1/2.

\param x Real in [0,1].
*/
inline double Quadric::Warp(const double& x)
{
  const double y = 2.0 * x;
  if (x < 0.5)
  {
    return 0.5 * Math::Sqr(y);
  }
  else
  {
    return 1.0 - 0.5 * Math::Sqr(2.0 - y);
  }
}

