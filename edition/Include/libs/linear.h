// Linear polynomials

#pragma once

#include "libs/mathematics.h"

class Linear {
protected:
  double c[2] = { 0.0 }; //!< Array of coefficients.
public:
  //! Empty.
  Linear() {}
  explicit Linear(const double&, const double&);
  explicit Linear(const double&);
  //! Empty.
  ~Linear() {}

  constexpr double& operator[] (int);
  constexpr double operator[] (int) const;

  // Unary operators
  Linear operator+ () const { return *this; }
  Linear operator- () const;

  Linear& operator+= (const Linear&);
  Linear& operator-= (const Linear&);
  Linear& operator*= (const double&);
  Linear& operator/= (const double&);

  // Binary operators
  friend Linear operator+ (const Linear&, const Linear&);
  friend Linear operator- (const Linear&, const Linear&);

  friend Linear operator* (const Linear&, const double&);
  friend Linear operator* (const double&, const Linear&);
  friend Linear operator/ (const Linear&, const double&);

  // Evaluates polynomial
  constexpr double operator()(const double&) const;

  // Solve
  int Solve(double&) const;
  int Solve(double&, const double&, const double&) const;
  void Range(double&, double&, const double& = 0.0, const double& = 1.0) const;

  static double Step(const double&, const double&, const double&);
  static double Step(const double&, const double&, const double&, const double&, const double&);
  static double Affine(const double&, const double&, const double&);
  static double Solve(const double&, const double&, const double&, const double&);
  static double SolveAlpha(const double&, const double&, const double&, const double&);

  friend std::ostream& operator<<(std::ostream&, const Linear&);
public:
  static const Linear Id; // Identity.
};

/*!
\brief Creates a linear function.

The function is defined as f(x) = a x + b.
\param a,b Coefficients.
*/
inline Linear::Linear(const double& a, const double& b) :c{ b,a }
{
}

/*!
\brief Creates a constant function.
\param a Constant.
*/
inline Linear::Linear(const double& a) :c{ a,0.0 }
{
}

//! Access class components.
inline constexpr double& Linear::operator[] (int i)
{
  return c[i];
}

//! Overloaded.
inline constexpr double Linear::operator[] (int i) const
{
  return c[i];
}

/*!
\brief Evaluates the linear function.
\param x Argument.
*/
inline constexpr double Linear::operator()(const double& x) const
{
  return c[0] + x * c[1];
}

//! Multiply a polynomial by a scalar value.
inline Linear operator* (const Linear& u, const double& e)
{
  return Linear(e * u[1], e * u[0]);
}

//! Overloaded.
inline Linear operator+ (const Linear& u, const Linear& v)
{
  return Linear(v[1] + u[1], v[0] + u[0]);
}

//! Overloaded.
inline Linear operator- (const Linear& v, const Linear& u)
{
  return Linear(v[1] - u[1], v[0] - u[0]);
}

//! Destructive sum of two linear polynomials.
inline Linear& Linear::operator+=(const Linear& u)
{
  for (int i = 0; i < 2; i++)
  {
    c[i] += u[i];
  }

  return *this;
}

//! Destructive difference of two linear polynomials.
inline Linear& Linear::operator-=(const Linear& u)
{
  for (int i = 0; i < 2; i++)
  {
    c[i] -= u[i];
  }
  return *this;
}

//! Scale a linear polynomial by a double value.
inline Linear& Linear::operator*=(const double& e)
{
  for (int i = 0; i < 2; i++)
  {
    c[i] *= e;
  }
  return *this;
}

//! Scales a linear polynomial by a double value.
inline Linear& Linear::operator/=(const double& e)
{
  for (int i = 0; i < 2; i++)
  {
    c[i] /= e;
  }
  return *this;
}

//! Unary.
inline Linear Linear::operator- () const
{
  return Linear(-c[1], -c[0]);
}

//! Overloaded.
inline Linear operator*(const double& a, const Linear& p)
{
  return p * a;
}

//! Overloaded.
inline Linear operator/(const Linear& p, const double& a)
{
  return p * (1.0 / a);
}

/*!
\brief Create the linear part of a linear step

This function does not test limits.

\param x Value
\param a, b Interval values.
\return Real in unit inverval.
\sa Linear::Step
*/
inline double Linear::Affine(const double& x, const double& a, const double& b)
{
  return (x - a) / (b - a);
}

/*!
\brief Create a linear step.
\param x Value
\param a, b Interval values.
\return Real in unit inverval.
\sa Cubic::SmoothStep, Quintic::SmoothStep
*/
inline double Linear::Step(const double& x, const double& a, const double& b)
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
    return (x - a) / (b - a);
  }
}

/*!
\brief Create a linear step.

This function is the same as coding:
\code
double y=Math::Lerp(ya,yb,Linear::Step(x,a,b));
\endcode
\param x Value
\param a, b Interval values.
\param ya, yb Output interval values.
\sa Cubic::SmoothStep, Quintic::SmoothStep
*/
inline double Linear::Step(const double& x, const double& a, const double& b, const double& ya, const double& yb)
{
  if (x < a)
  {
    return ya;
  }
  else if (x > b)
  {
    return yb;
  }
  else
  {
    return ya + (x - a) * (yb - ya) / (b - a);
  }
}

/*!
\brief Compute the root of a linear function such that f(a)=va and f(b)=vb.

Assuming that the root exists, this function reads easier than the following code, and avoids constructor and parameter computation.

\code
double t;
int n=Linear((vb-va)/(b-a),b*va+a*vb).Solve(t);
\endcode

\sa Vector::Solve(const Vector& a, const Vector& b, const double&, const double&);
*/
inline double Linear::Solve(const double& a, const double& b, const double& va, const double& vb)
{
  return (vb * a - va * b) / (vb - va);
}

