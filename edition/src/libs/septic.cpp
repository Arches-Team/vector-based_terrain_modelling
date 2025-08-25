// Septic polynomials

#include "libs/septic.h"
#include "libs/polynomial.h"

/*!
\class Septic septic.h
\brief %Septic (heptic) polynomials.

Constructors should provide the coefficients in descending order.
Example of how to code the septic 12x<SUP>7</SUP>-x+1:
\code
Septic p(12.0,0.0,0.0,0.0,0.0,-1.0,1.0);
\endcode

\ingroup MathGroup
*/

double Septic::epsilon = 1.0e-10;

/*!
\brief Check the degree of the septic.
*/
int Septic::CheckDegree() const
{
  int n = 7;
  while (n > 0)
  {
    if (c[n] != 0.0)
      return n;
    n--;
  }
  return n;
}

/*!
\brief Destructive sum of two septics.
*/
Septic& Septic::operator+= (const Septic& u)
{
  for (int i = 0; i < 8; i++)
  {
    c[i] += u[i];
  }

  return *this;
}

/*!
\brief Destructive difference of two septics.
*/
Septic& Septic::operator-= (const Septic& u)
{
  for (int i = 0; i < 8; i++)
  {
    c[i] -= u[i];
  }

  return *this;
}

/*!
\brief Scale a septics by a double value.
*/
Septic& Septic::operator*= (const double& e)
{
  for (int i = 0; i < 8; i++)
  {
    c[i] *= e;
  }
  return *this;
}

/*!
\brief Scale a septics by a double value.
*/
Septic& Septic::operator/= (const double& e)
{
  const double ie = 1.0 / e;
  for (int i = 0; i < 8; i++)
  {
    c[i] *= ie;
  }
  return *this;
}


/*!
\brief Solve the septics equation over a given interval.

This function stores the
sorted roots in an array and returns the number of roots.

This function calls Polynomial::Solve() if the highest
coefficient is not nul, otherwise it calls Quartic::Solve()
and other lower degree polynomial solvers which are more efficient.
\param roots The array of roots.
*/
int Septic::Solve(double* roots) const
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

/*!
\brief Search the roots of a septic equation over a given interval.
\param roots Array for storing the roots.
\param a, b Interval range.
*/
int Septic::Solve(double* roots, const double& a, const double& b) const
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
\brief Compute the range of values taken by a septic over a given interval.

\param a, b Interval.
\param x, y Returned range.
*/
void Septic::Range(double& x, double& y, const double& a, const double& b) const
{
  double r[6];

  x = (*this)(a);
  y = (*this)(b);

  Math::Sort(x, y);

  Sextic p = Prime();
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
\param p The septic.
*/
std::ostream& operator<<(std::ostream& s, const Septic& p)
{
  s << "Septic(" << p[7] << ',' << p[6] << ',' << p[5] << ',' << p[4] << ',' << p[3] << ',' << p[2] << ',' << p[1] << ',' << p[0] << ')';
  return s;
}

/*!
\brief Creates an Hermite septic polynomial on interval [0,1].

\author LoÃ¯s Paulin

\param a, b Values at t=0 and t=1.
\param da, db <B>D</B>erivatives.
\param sa, sb <B>S</B>econd derivatives.
\param ta, tb <B>T</B>hird derivatives.
*/
Septic Septic::Hermite(const double& a, const double& b, const double& da, const double& db, const double& sa, const double& sb, const double& ta, const double& tb)
{
  return
    Septic(
      (120.0 * a - 120.0 * b + 60.0 * da + 60.0 * db + 12.0 * sa - 12.0 * sb + ta + tb) / 6.0,
      (-420.0 * a + 420.0 * b - 216.0 * da - 204.0 * db - 45.0 * sa + 39.0 * sb - 4.0 * ta - 3.0 * tb) / 6.0,
      (168.0 * a - 168.0 * b + 90.0 * da + 78.0 * db + 20.0 * sa - 14.0 * sb + 2.0 * ta + tb) / 2.0,
      (-210.0 * a + 210.0 * b - 120.0 * da - 90.0 * db - 30.0 * sa + 15.0 * sb - 4.0 * ta - tb) / 6.0,
      ta / 6.0,
      sa / 2.0,
      da,
      a);
}

/*!
\brief Compute the Lipschitz constant of the septic.

\param a,b Interval.
*/
double Septic::K(const double& a, const double& b) const
{
  double x, y;

  Sextic sextic = Prime();
  sextic.Range(x, y, a, b);
  return Math::Max(fabs(x), fabs(y));
}
