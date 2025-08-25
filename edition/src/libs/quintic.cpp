// Quintic polynomials

#include "libs/quintic.h"
#include "libs/polynomial.h"

/*!
\class Quintic quintic.h
\brief %Quintic polynomials.

Closed form expression of roots does not exist for polynomials with a degree greater than 5.

Constructors should provide the coefficients in descending order.
Example of how to code the quintic x<SUP>5</SUP>-x+1:
\code
Quintic p(1.0,0.0,0.0,-1.0,1.0);
\endcode

Quintic implements some inline static member function such as:
\code
inline double Quintic::Smooth(const double& x)
{
return x * x * x * (x * (x * 6 - 15) + 10);
}
\endcode

\ingroup MathGroup
*/

double Quintic::epsilon = 1.0e-10;

/*!
\brief Check the degree of the quintic.
*/
int Quintic::CheckDegree() const
{
  int n = 5;
  while (n > 0)
  {
    if (c[n] != 0.0)
      return n;
    n--;
  }
  return n;
}

/*!
\brief Destructive sum.
*/
Quintic& Quintic::operator+= (const Quintic& u)
{
  for (int i = 0; i < 6; i++)
  {
    c[i] += u[i];
  }

  return *this;
}

/*!
\brief Destructive difference.
*/
Quintic& Quintic::operator-= (const Quintic& u)
{
  for (int i = 0; i < 6; i++)
  {
    c[i] -= u[i];
  }

  return *this;
}

/*!
\brief Scale.
*/
Quintic& Quintic::operator*= (const double& e)
{
  for (int i = 0; i < 6; i++)
  {
    c[i] *= e;
  }
  return *this;
}

/*!
\brief Scale.
*/
Quintic& Quintic::operator/= (const double& e)
{
  const double ie = 1.0 / e;
  for (int i = 0; i < 6; i++)
  {
    c[i] *= ie;
  }
  return *this;
}

/*!
\brief Compute the range of values taken by a quintic over a given interval.

\param a, b Interval.
\param x, y Returned range.
*/
void Quintic::Range(double& x, double& y, const double& a, const double& b) const
{
  double r[4];

  x = (*this)(a);
  y = (*this)(b);

  Math::Sort(x, y);

  Quartic p = Prime();
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
\brief Solve the quintic equation over a given interval.

This function stores the
sorted roots in an array and returns the number of roots.

This function calls Polynomial::Solve() if the highest
coefficient is not nul, otherwise it calls Quartic::Solve()
and other lower degree polynomial solvers which are more efficient.
\param roots The array of roots.
*/
int Quintic::Solve(double* roots) const
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

/*!
\brief Search the roots of a quintic equation over a given interval.
\param roots Array for storing the roots.
\param a, b Interval range.
*/
int Quintic::Solve(double* roots, const double& a, const double& b) const
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
\brief Overloaded output-stream operator.
\param s Stream.
\param p The quintic.
*/
std::ostream& operator<<(std::ostream& s, const Quintic& p)
{
  s << "Quintic(" << p[5] << ',' << p[4] << ',' << p[3] << ',' << p[2] << ',' << p[1] << ',' << p[0] << ')';
  return s;
}

/*!
\brief Creates an Hermite quintic polynomial on interval [0,1].

\param a, b Values at t=0 and t=1.
\param ta, tb Derivatives at t=0 and t=1.
\param na, nb Second derivatives at t=0 and t=1.
*/
Quintic Quintic::Hermite(const double& a, const double& b, const double& ta, const double& tb, const double& na, const double& nb)
{
  return
    Quintic(-6.0, 15.0, -10.0, 0.0, 0.0, 1.0) * a +
    Quintic(6.0, -15.0, 10.0, 0.0, 0.0, 0.0) * b +
    Quintic(-3.0, 8.0, -6.0, 0.0, 1.0, 0.0) * ta +
    Quintic(-3.0, 7.0, -4.0, 0.0, 0.0, 0.0) * tb +
    Quintic(-0.5, 1.5, -1.5, 0.5, 0.0, 0.0) * na +
    Quintic(0.5, -1.0, 0.5, 0.0, 0.0, 0.0) * nb;
}

/*!
\brief Compute the Lipschitz constant of the quintic.

\param a,b Interval.
*/
double Quintic::K(const double& a, const double& b) const
{
  double x, y;

  Quartic quartic = Prime();
  quartic.Range(x, y, a, b);
  return Math::Max(fabs(x), fabs(y));
}

/*!
\brief Create a quintic from a set of roots.
\param r Set of roots.

Example of how to code the quintic (x-1)<SUP>3</SUP>(x+1)(x+2):
\code
Quintic p=Quintic::FromRoots({1.0,1.0,1.0,-1.0,-2.0});
\endcode
*/
Quintic Quintic::FromRoots(const QVector<double>& r)
{
  return Quintic(
    1.0,
    -r[0] - r[1] - r[2] - r[3] - r[4],
    r[0] * r[1] + r[2] * r[3] + (r[0] + r[1]) * (r[2] + r[3]) + (r[0] + r[1] + r[2] + r[3]) * r[4],
    -(r[0] + r[1]) * r[2] * r[3] - r[0] * r[1] * (r[2] + r[3]) - (r[0] * r[1] + r[2] * r[3] + (r[0] + r[1]) * (r[2] + r[3])) * r[4],
    r[0] * r[1] * r[2] * r[3] + ((r[0] + r[1]) * r[2] * r[3] + r[0] * r[1] * (r[2] + r[3])) * r[4],
    -r[0] * r[1] * r[2] * r[3] * r[4]);
}
