// Sextic polynomials

#include "libs/sextic.h"
#include "libs/polynomial.h"

/*!
\class Sextic sextic.h
\brief %Sextic polynomials, also known as hextics.

Constructors should provide the coefficients in descending order.
Example of how to code the sextic 4x<SUP>6</SUP>-x+1:
\code
Sextic p(4.0,0.0,0.0,0.0,-1.0,1.0);
\endcode

\ingroup MathGroup
*/

double Sextic::epsilon = 1.0e-10;

/*!
\brief Check the degree of the sextic.
*/
int Sextic::CheckDegree() const
{
  int n = 6;
  while (n > 0)
  {
    if (c[n] != 0.0)
      return n;
    n--;
  }
  return n;
}

/*!
\brief Destructive sum of two sextics.
*/
Sextic& Sextic::operator+= (const Sextic& u)
{
  for (int i = 0; i < 7; i++)
  {
    c[i] += u[i];
  }

  return *this;
}

/*!
\brief Destructive difference of two sextics.
*/
Sextic& Sextic::operator-= (const Sextic& u)
{
  for (int i = 0; i < 7; i++)
  {
    c[i] -= u[i];
  }

  return *this;
}

/*!
\brief Scale a sextic by a double value.
*/
Sextic& Sextic::operator*= (const double& e)
{
  for (int i = 0; i < 7; i++)
  {
    c[i] *= e;
  }
  return *this;
}

/*!
\brief Scale a sextic by a double value.
*/
Sextic& Sextic::operator/= (const double& e)
{
  const double ie = 1.0 / e;
  for (int i = 0; i < 7; i++)
  {
    c[i] *= ie;
  }
  return *this;
}


/*!
\brief Solve the sextic equation over a given interval.

This function stores the
sorted roots in an array and returns the number of roots.

This function calls Polynomial::Solve() if the highest
coefficient is not nul, otherwise it calls Quartic::Solve()
and other lower degree polynomial solvers which are more efficient.
\param roots The array of roots.
*/
int Sextic::Solve(double* roots) const
{
  if (c[6] == 0.0)
  {
    if (c[5] == 0.0)
    {
      if (c[4] == 0.0)
      {
        if (c[3] == 0.0)
        {
          return Quadric(c[2], c[1], c[0]).Solve(roots);
        }
        return Cubic(c[3], c[2], c[1], c[0]).Solve(roots);
      }
      return Quartic(c[4], c[3], c[2], c[1], c[0]).Solve(roots);
    }
    return Polynomial(c[5], c[4], c[3], c[2], c[1], c[0]).Solve(roots);
  }
  return Polynomial(c[6], c[5], c[4], c[3], c[2], c[1], c[0]).Solve(roots);
}

/*!
\brief Search the roots of a sextic equation over a given interval.
\param roots Array for storing the roots.
\param a, b Interval range.
*/
int Sextic::Solve(double* roots, const double& a, const double& b) const
{
  int r = Solve(roots);
  int j = 0;

  for (int i = 0; i < r; i++)
  {
    if ((roots[i] > a) && (roots[i] < b))
    {
      roots[j] = roots[i];
      j++;
    }
  }

  return j;
}

/*!
\brief Compute the range of values taken by a sextic over a given interval.

\param a, b Interval.
\param x, y Returned range.
*/
void Sextic::Range(double& x, double& y, const double& a, const double& b) const
{
  double r[5];

  x = (*this)(a);
  y = (*this)(b);

  Math::Sort(x, y);

  Quintic p = Prime();
  int n = p.Solve(r);

  for (int i = 0; i < n; i++)
  {
    if ((r[i] > a) && (r[i] < b))
    {
      double s = (*this)(r[i]);
      if (s < x)
      {
        x = s;
      }
      else if (s > y)
      {
        y = s;
      }
    }
  }
}

/*!
\brief Overloaded output-stream operator.
\param s Stream.
\param p The sextic.
*/
std::ostream& operator<<(std::ostream& s, const Sextic& p)
{
  s << "Sextic(" << p[6] << ',' << p[5] << ',' << p[4] << ',' << p[3] << ',' << p[2] << ',' << p[1] << ',' << p[0] << ')';
  return s;
}

/*!
\brief Compose the quadric by a cubic.

This function computes p ( q(x) ), where p denotes the current quadric and q the argument cubic.
\param p The quadric.
\param q The cubic.
*/
Sextic Sextic::Compose(const Quadric& p, const Cubic& q)
{
  return Sextic(
    p[2] * q[3] * q[3],
    2.0 * p[2] * q[3] * q[2],
    p[2] * (2.0 * q[3] * q[1] + q[2] * q[2]),
    p[2] * (2.0 * q[3] * q[0] + 2.0 * q[2] * q[1]) + p[1] * q[3],
    p[2] * (2.0 * q[2] * q[0] + q[1] * q[1]) + p[1] * q[2],
    p[2] * (2.0 * q[1] * q[0]) + p[1] * q[1],
    (p[2] * q[0] + p[1]) * q[0] + p[0]);
}

/*!
\brief Compose the cubic by a quadric.

This function computes p ( q(x) ), where p denotes the cubic and q the quadric.
\param p The cubic.
\param q The quadric.
*/
Sextic Sextic::Compose(const Cubic& p, const Quadric& q)
{
  const double q10 = q[1] * q[0];
  const double q21 = q[2] * q[1];
  const double q20 = q[2] * q[0];
  const double q22 = q[2] * q[2];
  const double q11 = q[1] * q[1];
  const double q00 = q[0] * q[0];

  return Sextic(p[3] * (q22 * q[2]),
    p[3] * (3.0 * q22 * q[1]),
    p[2] * q22 + p[3] * (3.0 * q22 * q[0] + 3.0 * q[2] * q11),
    p[2] * 2.0 * q21 + p[3] * (q11 * q[1] + 6.0 * q21 * q[0]),
    p[1] * q[2] + p[3] * (3.0 * q11 * q[0] + 3.0 * q[2] * q00),
    p[1] * q[1] + p[2] * 2.0 * q10 + p[2] * (2.0 * q20 + q11) + p[3] * (3.0 * q[1] * q00),
    p[0] + p[1] * q[0] + p[2] * q00 + p[3] * (q00 * q[0]));
}

/*!
\brief Compute the Lipschitz constant of the sextic.

\param a,b Interval.
*/
double Sextic::K(const double& a, const double& b) const
{
  double x, y;

  Quintic quintic = Prime();
  quintic.Range(x, y, a, b);
  return Math::Max(fabs(x), fabs(y));
}
