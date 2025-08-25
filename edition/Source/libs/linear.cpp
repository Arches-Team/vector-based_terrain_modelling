// Linear polynomials

#include "libs/linear.h"

/*!
\class Linear linear.h
\brief %Linear polynomials.

When using constructors, the coefficients are given in
descending order:
\code
Linear p(-1.0,1.0); // Linear polynomial -x + 1
\endcode

\ingroup MathGroup

*/

/*!
\page Transfer Transfer functions

Core implements a variety of compact and useful transfer functions. Most of them are grouped in the Linear, Quadric, Cubic classes, although some
specific more complex transfer functions, in general providing higher regularity properties, can be found in other classes such as Quintic, Septic or Nonic polynomials.

\section TransferLinear %Linear functions
<P><I>How do I compute the linear mapping from one interval to another?</I>
<BR>In general, use  double Linear::Step(const double& , const double& , const double& , const double& , const double& ),
for instance for mapping [-3, 7] onto [-1,1] call:
\code
double y = Linear::Step(x,-3.0,7.0,-1.0,1.0);
\endcode
<P><I>How do I compute a gain-like function?</I>
The function double Quadric::Warp(const double&) implements a fast symmetric function based on a quadric;
a higher degree Cubic::Warp(const double&) and lower degree Math::Warp(const double&) implementations are also implemented.
*/

const Linear Linear::Id(1.0, 0.0);

/*!
\brief This function computes the range of values taken by a
Linear over a given interval.
\param a,b Interval.
\param x,y Returned range.
*/
void Linear::Range(double& x, double& y, const double& a, const double& b) const
{
  x = (*this)(a);
  y = (*this)(b);

  Math::Sort(x, y);
}

/*!
\brief Search the roots of a linear equation over a given interval.
\param x Solution.
\param a,b Interval.
*/
int Linear::Solve(double& x, const double& a, const double& b) const
{
  if (Solve(x))
  {
    if ((x > a) && (x < b))
      return 1;
  }

  return 0;
}

/*!
\brief Solve linear equations.

This function store the root (ifit exists) in a double and returns the number of roots.
\param x Solution.
*/
int Linear::Solve(double& x) const
{
  if (c[1] == 0.0)
    return 0;

  x = -c[0] / c[1];
  return 1;
}

/*!
\brief Overloaded output-stream operator.
\param s Stream.
\param p Linear polynomial.
*/
std::ostream& operator<<(std::ostream& s, const Linear& p)
{
  s << "Linear (" << p[1] << ',' << p[0] << ')';
  return s;
}

/*!
\brief Compute the root of a linear function.

The linear function satisfies the constraints: f(a)=va and f(b)=vb;
the function returns an interpolant t such that the root satisties r=(1-t) a+t b.

\sa Linear::Solve(const Vector& a, const Vector& b, const double&, const double&);
*/
double Linear::SolveAlpha(const double& a, const double& b, const double& va, const double& vb)
{
  double t = Solve(a, b, va, vb);
  return (t - a) / (b - a);
}