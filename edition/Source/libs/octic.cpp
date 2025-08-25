// Octic polynomials

#include "libs/nonic.h"
#include "libs/polynomial.h"

/*!
\class Octic nonic.h
\brief %Octic polynomials.

Constructors should provide the coefficients in descending order.
Example of how to code the octic 19x<SUP>8</SUP>-x+1:
\code
Octic p(19.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0);
\endcode

\ingroup MathGroup
*/

double Octic::epsilon = 1.0e-10;

/*!
\brief Check the degree of the octic.
*/
int Octic::CheckDegree() const
{
  int n = 8;
  while (n > 0)
  {
    if (c[n] != 0.0)
      return n;
    n--;
  }
  return n;
}

/*!
\brief Destructive sum of two octics.
*/
Octic& Octic::operator+= (const Octic& u)
{
  for (int i = 0; i < 9; i++)
  {
    c[i] += u[i];
  }

  return *this;
}

/*!
\brief Destructive difference of two octics.
*/
Octic& Octic::operator-= (const Octic& u)
{
  for (int i = 0; i < 9; i++)
  {
    c[i] -= u[i];
  }

  return *this;
}

/*!
\brief Scale an octic by a double value.
*/
Octic& Octic::operator*= (const double& e)
{
  for (int i = 0; i < 9; i++)
  {
    c[i] *= e;
  }
  return *this;
}

/*!
\brief Scale an octic by a double value.
*/
Octic& Octic::operator/= (const double& e)
{
  const double ie = 1.0 / e;
  for (int i = 0; i < 9; i++)
  {
    c[i] *= ie;
  }
  return *this;
}

/*!
\brief Overloaded output-stream operator.
\param s Stream.
\param p The Octic.
*/
std::ostream& operator<<(std::ostream& s, const Octic& p)
{
  s << "Octic(" << p[8] << ',' << p[7] << ',' << p[6] << ',' << p[5] << ',' << p[4] << ',' << p[3] << ',' << p[2] << ',' << p[1] << ',' << p[0] << ')';
  return s;
}

/*!
\brief Multiply two quartics.

\param u, v Argument quartics.
*/
Octic operator*(const Quartic& u, const Quartic& v)
{
  return Octic(
    u[4] * v[4],
    u[4] * v[3] + u[3] * v[4],
    u[4] * v[2] + u[3] * v[3] + u[2] * v[4],
    u[4] * v[1] + u[3] * v[2] + u[2] * v[3] + u[1] * v[4],
    u[4] * v[0] + u[3] * v[1] + u[2] * v[2] + u[1] * v[3] + u[0] * v[4],
    u[3] * v[0] + u[2] * v[1] + u[1] * v[2] + u[0] * v[3],
    u[2] * v[0] + u[1] * v[1] + u[0] * v[2],
    u[1] * v[0] + u[0] * v[1],
    u[0] * v[0]);
}

/*!
\brief Multiply a sextic by a quadric.

\param s %Sextic.
\param q %Quadric.
*/
Octic operator*(const Sextic& s, const Quadric& q)
{
  return Octic(
    s[6] * q[2],
    s[5] * q[2] + s[6] * q[1],
    s[4] * q[2] + s[5] * q[1] + s[6] * q[0],
    s[3] * q[2] + s[4] * q[1] + s[5] * q[0],
    s[2] * q[2] + s[3] * q[1] + s[4] * q[0],
    s[1] * q[2] + s[2] * q[1] + s[3] * q[0],
    s[0] * q[2] + s[1] * q[1] + s[2] * q[0],
    s[0] * q[1] + s[1] * q[0],
    s[0] * q[0]);
}

/*!
\brief Solve the octic equation over a given interval.

\param roots Roots.
*/
int Octic::Solve(double* roots) const
{
  if (c[8] == 0.0)
  {
    if (c[7] == 0.0)
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
    return Polynomial(c[7], c[6], c[5], c[4], c[3], c[2], c[1], c[0]).Solve(roots);
  }
  return Polynomial(c[8], c[7], c[6], c[5], c[4], c[3], c[2], c[1], c[0]).Solve(roots);
}
