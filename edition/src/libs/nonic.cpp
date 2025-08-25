// Nonic polynomials

#include "libs/nonic.h"
#include "libs/polynomial.h"

/*!
\class Nonic nonic.h
\brief %Nonic polynomials.

Constructors should provide the coefficients in descending order.
Example of how to code the Nonic 19x<SUP>9</SUP>-x+1:
\code
Nonic p(19.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0);
\endcode

\ingroup MathGroup
*/

double Nonic::epsilon = 1.0e-10;

/*!
\brief Check the degree of the Nonic.
*/
int Nonic::CheckDegree() const
{
  int n = 9;
  while (n > 0)
  {
    if (c[n] != 0.0)
      return n;
    n--;
  }
  return n;
}

/*!
\brief Destructive sum of two Nonics.
*/
Nonic& Nonic::operator+= (const Nonic& u)
{
  for (int i = 0; i < 10; i++)
  {
    c[i] += u[i];
  }

  return *this;
}

/*!
\brief Destructive difference of two Nonics.
*/
Nonic& Nonic::operator-= (const Nonic& u)
{
  for (int i = 0; i < 10; i++)
  {
    c[i] -= u[i];
  }

  return *this;
}

/*!
\brief Scale a Nonics by a double value.
*/
Nonic& Nonic::operator*= (const double& e)
{
  for (int i = 0; i < 10; i++)
  {
    c[i] *= e;
  }
  return *this;
}

/*!
\brief Scale a Nonics by a double value.
*/
Nonic& Nonic::operator/= (const double& e)
{
  const double ie = 1.0 / e;
  for (int i = 0; i < 10; i++)
  {
    c[i] *= ie;
  }
  return *this;
}

/*!
\brief Compose the cubic by another one.

This function computes p ( q(x) ), where p denotes the cubic and q the argument.
\param q The cubic.
\param p The argument cubic.
*/
Nonic Nonic::Compose(const Cubic& p, const Cubic& q)
{
  return Nonic(
    p[3] * (q[3] * q[3] * q[3]),
    p[3] * (3.0 * q[3] * q[3] * q[2]),
    p[3] * (3.0 * q[3] * (q[2] * q[2] + q[3] * q[1])),
    p[3] * (q[2] * q[2] * q[2] + 6.0 * q[3] * q[2] * q[1] + 3.0 * q[3] * q[3] * q[0]) + p[2] * (q[3] * q[3]),
    p[3] * (3.0 * (q[2] * q[2] * q[1] + q[3] * q[1] * q[1] + 2.0 * q[3] * q[2] * q[0])) + p[2] * (2.0 * q[3] * q[2]),
    p[3] * (3.0 * (q[2] * q[1] * q[1] + q[2] * q[2] * q[0] + 2.0 * q[3] * q[1] * q[0])) + p[2] * (q[2] * q[2] + 2.0 * q[3] * q[1]),
    p[3] * (q[1] * q[1] * q[1] + 6.0 * q[2] * q[1] * q[0] + 3.0 * q[3] * q[0] * q[0]) + p[2] * (2.0 * q[2] * q[1] + 2.0 * q[3] * q[0]) + p[1] * (q[3]),
    p[3] * (3.0 * q[0] * (q[1] * q[1] + q[2] * q[0])) + p[2] * (q[1] * q[1] + 2.0 * q[2] * q[0]) + p[1] * (q[2]),
    p[3] * (3.0 * q[1] * q[0] * q[0]) + p[2] * (2.0 * q[1] * q[0]) + p[1] * (q[1]),
    p[3] * (q[0] * q[0] * q[0]) + p[2] * (q[0] * q[0]) + p[1] * (q[0]) + p[0]
  );
}



/*!
\brief Overloaded output-stream operator.
\param s Stream.
\param p The Nonic.
*/
std::ostream& operator<<(std::ostream& s, const Nonic& p)
{
  s << "Nonic(" << p[9] << ',' << p[8] << ',' << p[7] << ',' << p[6] << ',' << p[5] << ',' << p[4] << ',' << p[3] << ',' << p[2] << ',' << p[1] << ',' << p[0] << ')';
  return s;
}

/*!
\brief Solve the octic equation over a given interval.

\param roots Roots.
*/
int Nonic::Solve(double* roots) const
{
  if (c[9] == 0.0)
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
  return Polynomial(c[9], c[8], c[7], c[6], c[5], c[4], c[3], c[2], c[1], c[0]).Solve(roots);
}
