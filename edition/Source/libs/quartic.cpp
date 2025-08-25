// Quartic polynomials

#include "libs/quartic.h"

/*!
\class Quartic quartic.h
\brief %Quartic polynomials.

Closed form expression of roots exist for quartics.

Constructors should provide the coefficients in descending order.
Example of how to code the quartic x<SUP>4</SUP>+x<SUP>3</SUP>-x+1:
\code
Quartic p(3.0,1.0,0.0,-1.0,1.0);
\endcode

\ingroup MathGroup
*/

double Quartic::epsilon = 1.0e-10;
double Quartic::limit = 1.0e-16;

/*!
\brief Check the degree of the quartic.
*/
int Quartic::CheckDegree() const
{
  int n = 4;
  while (n > 0)
  {
    if (c[n] != 0.0)
      return n;
    n--;
  }
  return n;
}

/*!
\brief Destructive sum of two quartics.
*/
Quartic& Quartic::operator+= (const Quartic& u)
{
  for (int i = 0; i < 5; i++)
  {
    c[i] += u[i];
  }

  return *this;
}

/*!
\brief Destructive difference.
*/
Quartic& Quartic::operator-= (const Quartic& u)
{
  for (int i = 0; i < 5; i++)
  {
    c[i] -= u[i];
  }

  return *this;
}

/*!
\brief Scale a quartic by a double value.
*/
Quartic& Quartic::operator*= (const double& e)
{
  for (int i = 0; i < 5; i++)
  {
    c[i] *= e;
  }
  return *this;
}

/*!
\brief Scale a quartic by a double value.
*/
Quartic& Quartic::operator/= (const double& e)
{
  const double ie = 1.0 / e;
  for (int i = 0; i < 5; i++)
  {
    c[i] *= ie;
  }
  return *this;
}

/*!
\brief Overloaded output-stream operator.
\param s Stream.
\param p The quartic.
*/
std::ostream& operator<<(std::ostream& s, const Quartic& p)
{
  s << "Quartic(" << p[4] << ',' << p[3] << ',' << p[2] << ',' << p[1] << ',' << p[0] << ')';
  return s;
}

/*!
\brief Compute the range of values taken by a quartic over a given interval.

Compute the roots of the first derivative, and evaluates the cubic
at the roots if they are within the interval bounds.

\param a, b Interval.
\param x, y Returned range.
*/
void Quartic::Range(double& x, double& y, const double& a, const double& b) const
{
  double r[3];

  x = (*this)(a);
  y = (*this)(b);

  Math::Sort(x, y);

  Cubic p = Prime();
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
\brief Search the roots of a quartic equation.

The roots are not sorted.

This function calls other root solving functions,
invoking Cubic::Solve() or Quadric::Solve() as appropriate
if the degree is lower than four (which means that closed
form expressions of the roots exist).

Otherwise, solve the quartic using the method of Francois
Vieta (Circa 1735).

This function checks if the higher coefficients are nul,
possibly calling Quadric::Solve() and Cubic::Solve() if need be.

\param results Array for storing the roots.
*/
int Quartic::Solve(double* results)
{
  if (c[4] == 0.0)
  {
    if (c[3] == 0.0)
    {
      return Quadric(c[2], c[1], c[0]).Solve(results);
    }
    return Cubic(c[3], c[2], c[1], c[0]).Solve(results);
  }
  double roots[3];
  double z, d2;

  // Figure out the size difference between coefficients
  if (Analyze())
  {
    if (fabs(c[4]) < limit)
    {
      if (fabs(c[3]) < limit)
      {
        return Quadric(c[2], c[1], c[0]).Solve(results);
      }
      else
      {
        return Cubic(c[3], c[2], c[1], c[0]).Solve(results);
      }
    }
    else
    {
      return Solve(results);
    }
  }

  // See if the high order term has vanished 
  double c0 = c[4];
  double c1, c2, c3, c4;

  if (fabs(c0) < limit)
  {
    return Cubic(c[3], c[2], c[1], c[0]).Solve(results);
  }
  // See if the constant term has vanished 
  if (fabs(c[0]) < limit)
  {
    results[0] = 0.0;
    return 1 + Cubic(c[4], c[3], c[2], c[1]).Solve(results + 1);
  }
  // Make sure the quartic has a leading coefficient of 1.0 
  if (c0 != 1.0)
  {
    c1 = c[3] / c0;
    c2 = c[2] / c0;
    c3 = c[1] / c0;
    c4 = c[0] / c0;
  }
  else
  {
    c1 = c[3];
    c2 = c[2];
    c3 = c[1];
    c4 = c[0];
  }

  // Compute the cubic resolvant 
  double c12 = c1 * c1;
  double p = -0.375 * c12 + c2;
  double q = 0.125 * c12 * c1 - 0.5 * c1 * c2 + c3;
  double r = -0.01171875 * c12 * c12 + 0.0625 * c12 * c2 - 0.25 * c1 * c3 + c4;

  Cubic cubic(1.0, -0.5 * p, -r, 0.5 * r * p - 0.125 * q * q);

  int i = cubic.Solve(roots);
  if (i > 0)
  {
    z = roots[0];
  }
  else
  {
    return 0;
  }

  double d1 = 2.0 * z - p;

  if (d1 < 0.0)
  {
    if (d1 > -epsilon)
      d1 = 0.0;
    else
      return 0;
  }
  if (d1 < epsilon)
  {
    d2 = z * z - r;
    if (d2 < 0.0)
      return 0;
    d2 = sqrt(d2);
  }
  else
  {
    d1 = sqrt(d1);
    d2 = 0.5 * q / d1;
  }

  // Set up useful values for the quadratic factors 
  double q1 = d1 * d1;
  double q2 = -0.25 * c1;
  i = 0;

  // Solve the first quadratic
  p = q1 - 4.0 * (z - d2);
  if (p == 0)
  {
    results[i++] = -0.5 * d1 - q2;
  }
  else if (p > 0)
  {
    p = sqrt(p);
    results[i++] = -0.5 * (d1 + p) + q2;
    results[i++] = -0.5 * (d1 - p) + q2;
  }
  // Solve the second quadratic 
  p = q1 - 4.0 * (z + d2);
  if (p == 0)
  {
    results[i++] = 0.5 * d1 - q2;
  }
  else if (p > 0)
  {
    p = sqrt(p);
    results[i++] = 0.5 * (d1 + p) + q2;
    results[i++] = 0.5 * (d1 - p) + q2;
  }
  return i;
}

/*!
\brief Search the roots of a quartic equation over a given interval.
\param roots Array for storing the roots.
\param a, b Interval range.
*/
int Quartic::Solve(double* roots, const double& a, const double& b)
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
\brief Check if any coefficient of the quartic is more than 10<SUP>12</SUP> times larger than the smallest.

Small coefficients are set to 0.
*/
bool Quartic::Analyze()
{
  bool flag = false;

  double biggest = fabs(c[0]);
  for (int i = 1; i < 5; i++)
  {
    double t = fabs(c[i]);
    if (t > biggest)
    {
      biggest = t;
    }
  }

  // Everything is zero no sense in doing any more
  if (biggest == 0.0)
    return false;

  for (int j = 0; j < 5; j++)
  {
    if (c[j] != 0.0)
    {
      if (fabs(biggest / c[j]) > 1.0e12)
      {
        c[j] = 0.0;
        flag = true;
      }
    }
  }
  return flag;
}

/*!
\brief Compose the quadric by another one.

This function computes p ( q(x) ), where p denotes the quadric and q the argument.
\param q The quadric.
\param p The argument quadric.
*/
Quartic Quartic::Compose(const Quadric& p, const Quadric& q)
{
  return Quartic(
    p[2] * q[2] * q[2],
    2.0 * p[2] * q[2] * q[1],
    p[2] * (2.0 * q[2] * q[0] + q[1] * q[1]) + p[1] * q[2],
    2.0 * p[2] * q[1] * q[0] + p[1] * q[1],
    (p[2] * q[0] + p[1]) * q[0] + p[0]);
}

/*!
\brief Compute the Lipschitz constant of the quartic.

\param a,b Interval.
*/
double Quartic::K(const double& a, const double& b) const
{
  double x, y;

  Cubic cubic = Prime();
  cubic.Range(x, y, a, b);
  return Math::Max(fabs(x), fabs(y));
}
