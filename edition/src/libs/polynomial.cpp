// Polynomials

#include "libs/polynomial.h"

const double Polynomial::tiny = 1.0e-16;
const double Polynomial::epsilon = 1.0e-10;
const int Polynomial::iterations = 50;

/*!
\class Polynomial polynomial.h
\brief Polynomials of degree lower than or equal to 11.

This class is defined as a fixed size data-structure to speed-up and
facilitate memory management. It is more memory and computationnally
efficient to use Linear, Quadric, Cubic, and Quartic classes whenever possible.

When using constructors, the coefficients are given in descending order,
i.e. start from the highest degree coefficient to the lowest degree.

Example of how to code the quartic polynomial 3 x<SUP>4</SUP>-x+1:
\code
Polynomial p(3.0,0.0,0.0,-1.0,1.0);
\endcode
\sa Septic, Sextic, Quintic, Quartic, Cubic, Quadric, Linear.

\ingroup MathGroup
*/

/*!
\brief Creates a null polynomial.
*/
Polynomial::Polynomial()
{
  c[0] = 0.0;
}

/*!
\brief Set the null polynomial member constant.
*/
const Polynomial Polynomial::Null(0.0);

/*!
\brief Evaluates the polynomial.
\param x Parameter.
*/
double Polynomial::operator()(const double& x) const
{
  if (n < 4)
  {
    if (n == -1)
      return 0.0;
    else if (n == 0)
      return (c[0]);
    else if (n == 1)
      return (c[0] + x * c[1]);
    else if (n == 2)
      return (c[0] + x * (c[1] + x * c[2]));
    else
      return (c[0] + x * (c[1] + x * (c[2] + x * c[3])));
  }
  else
  {
    double y = c[n];

    for (int i = n - 1; i >= 0; i--)
    {
      y *= x;
      y += c[i];
    }
    return y;
  }
}

/*!
\brief Check the degree of the polynomial.

The function decreases the degree if higher coefficients are null.
*/
void Polynomial::Check()
{
  while (n > 0)
  {
    if (c[n] != 0.0)
      break;
    n--;
  }
}

int Polynomial::modp(const Polynomial* u, const Polynomial* v, Polynomial* r)
{
  int k;

  for (int i = 0; i < u->n; i++)
  {
    r[i] = u[i];
  }

  if (v->c[v->n] < 0.0)
  {
    for (k = u->n - v->n - 1; k >= 0; k -= 2)
    {
      r->c[k] = -r->c[k];
    }

    for (k = u->n - v->n; k >= 0; k--)
    {
      for (int j = v->n + k - 1; j >= k; j--)
      {
        r->c[j] = -r->c[j] - r->c[v->n + k] * v->c[j - k];
      }
    }
  }
  else
  {
    for (k = u->n - v->n; k >= 0; k--)
    {
      for (int j = v->n + k - 1; j >= k; j--)
      {
        r->c[j] -= r->c[v->n + k] * v->c[j - k];
      }
    }
  }

  k = v->n - 1;

  while (k >= 0 && fabs(r->c[k]) < Polynomial::tiny)
  {
    r->c[k] = 0.0;
    k--;
  }

  r->n = (k < 0) ? 0 : k;

  return(r->n);
}

/*!
\brief Compute the roots of a polynomial equation using bissection.
\param a,b Interval.
\param val Value.
\param epsilon Precision.
*/
int Polynomial::Bissection(double a, double b, double& val, const double& epsilon) const
{
  double fa = (*this)(a);
  double fb = (*this)(b);

  if (fa * fb > 0.0)
  {
    return 0;
  }

  if (fabs(fa) < tiny)
  {
    val = a;
    return 1;
  }

  if (fabs(fb) < tiny)
  {
    val = b;
    return 1;
  }

  double lfx = fa;

  // Maximum 50 iterations
  for (int i = 0; i < 50; i++)
  {
    double x = (fb * a - fa * b) / (fb - fa);

    double fx = (*this)(x);

    if (fabs(x) > epsilon)
    {
      if (fabs(fx / x) < epsilon)
      {
        val = x;
        return 1;
      }
    }
    else
    {
      if (fabs(fx) < epsilon)
      {
        val = x;
        return 1;
      }
    }

    if (fa < 0)
    {
      if (fx < 0)
      {
        a = x;
        fa = fx;

        if ((lfx * fx) > 0)
        {
          fb *= 0.5;
        }
      }
      else
      {
        b = x;
        fb = fx;

        if ((lfx * fx) > 0)
        {
          fa *= 0.5;
        }
      }
    }
    else
    {
      if (fx < 0)
      {
        b = x;
        fb = fx;

        if ((lfx * fx) > 0)
        {
          fa *= 0.5;
        }
      }
      else
      {
        a = x;
        fa = fx;

        if ((lfx * fx) > 0)
        {
          fb *= 0.5;
        }
      }
    }

    // Check for underflow in the domain
    if (fabs(b - a) < epsilon)
    {
      val = x;
      return 1;
    }

    lfx = fx;
  }

  return 0;
}

/*!
\brief Destructive sum of two polynomials.
*/
Polynomial& Polynomial::operator+= (const Polynomial& u)
{
  if (u.n > n)
  {
    for (int i = 0; i <= n; i++)
    {
      c[i] += u[i];
    }
    for (int i = n + 1; i <= u.n; i++)
    {
      c[i] = u[i];
    }
    n = u.n;
  }
  else
  {
    for (int i = 0; i <= u.n; i++)
    {
      c[i] += u[i];
    }
  }

  return *this;
}
/*!
\brief Destructive difference of two polynomials.
*/
Polynomial& Polynomial::operator-= (const Polynomial& u)
{
  if (u.n > n)
  {
    for (int i = 0; i <= n; i++)
    {
      c[i] -= u[i];
    }
    for (int i = n + 1; i <= u.n; i++)
    {
      c[i] = -u[i];
    }
    n = u.n;
  }
  else
  {
    for (int i = 0; i <= u.n; i++)
    {
      c[i] -= u[i];
    }
  }

  return *this;
}

/*!
\brief Scale a polynomial by a double value. Optimizations are provided
if scalar is 0.0 or 1.0.
*/
Polynomial& Polynomial::operator*= (const double& e)
{
  if (e == 0.0)
  {
    c[0] = 0.0;
    n = 0;
  }
  else if (e != 1.0)
  {
    for (int i = 0; i <= n; i++)
    {
      c[i] *= e;
    }
  }
  return *this;
}

/*!
\brief Scale a polynomial by a double value.
*/
Polynomial& Polynomial::operator/= (const double& e)
{
  for (int i = 0; i <= n; i++)
  {
    c[i] /= e;
  }
  return *this;
}

/*!
\brief Multiply a polynomial by a scalar value. Optimizations are provided
if scalar is 0.0 or 1.0.
*/
Polynomial operator* (const Polynomial& u, const double& e)
{
  // Identity
  if (e == 1.0)
  {
    return u;
  }
  // Null polynomial
  else if (e == 0.0)
  {
    return Polynomial::Null;
  }
  else
  {
    Polynomial r = u;
    r *= e;
    return r;
  }
}

/*!
\brief Overloaded.
*/
Polynomial operator+ (const Polynomial& u, const Polynomial& v)
{
  Polynomial r = u;
  r += v;

  return r;
}

/*!
\brief Overloaded.
*/
Polynomial operator- (const Polynomial& u, const Polynomial& v)
{
  Polynomial r = u;
  r -= v;

  return r;
}

/*!
\brief Overloaded.
*/
Polynomial& Polynomial::operator*=(const Polynomial& u)
{
  Polynomial r = (*this) * u;
  *this = r;
  return *this;
}

/*!
\brief Multiply two polynomials.
*/
Polynomial operator* (const Polynomial& u, const Polynomial& v)
{
  Polynomial r = Polynomial::Null;

  r.n = min(u.n + v.n, Polynomial::MaxDegree);

  for (int i = 0; i <= u.n; i++)
  {
    for (int j = 0; j <= v.n; j++)
    {
      r[i + j] += u[i] * v[j];
    }
  }
  return r;
}

/*!
\brief Overloaded output-stream operator.
\param s Stream.
\param p The polynomial.
*/
std::ostream& operator<<(std::ostream& s, const Polynomial& p)
{
  s << "Polynomial(";
  for (int i = p.n; i > 0; i--)
  {
    s << p.c[i] << ',';
  }
  s << p.c[0] << ')';
  return s;
}

/*!
\brief Reverse the terms of the polynomial.
*/
Polynomial Polynomial::Reversed() const
{
  Polynomial q;
  for (int i = n; i >= 0; i--)
  {
    q[i] = c[n - i];
  }
  return q;
}

/*!
\brief Compose the polynomial by the argument polynomial.
\param v Argument polynomial.
*/
Polynomial Polynomial::Compose(const Polynomial& v) const
{
  // Null or constant polynomial 
  if (n < 1)
  {
    return *this;
  }
  else
  {
    // Linear term
    Polynomial r = c[1] * v;
    r[0] += c[0];

    // Higher order terms
    Polynomial a = v;
    for (int i = 2; i <= n; i++)
    {
      a *= v;
      r += c[i] * a;
    }
    return r;
  }
}

/*!
\brief Compute the derivative of the polynomial.
*/
Polynomial Polynomial::Prime() const
{
  Polynomial r;

  for (int i = 1; i <= n; i++)
  {
    r[i - 1] = double(i) * c[i];
  }
  r.n = n - 1;
  return r;
}

/*!
\brief Unary.
*/
Polynomial Polynomial::operator- () const
{
  Polynomial p;
  p.n = n;
  for (int i = 0; i <= n; i++)
  {
    p[i] = -c[i];
  }
  return p;
}

/*!
\brief Overloaded.
*/
Polynomial operator*(const double& a, const Polynomial& p)
{
  return p * a;
}

/*!
\brief Overloaded.
*/
Polynomial operator/(const Polynomial& p, const double& a)
{
  return p * (1.0 / a);
}

/*!
\brief This function computes the range of values taken by a
polynomial over a given interval.

Computes the roots of the first derivative, and evaluates the polynomial
at the roots if they are within the interval bounds.
\param a,b Interval.
\param x,y Returned range.
*/
void Polynomial::Range(double& x, double& y, const double& a, const double& b) const
{
  double r[16];

  x = (*this)(a);
  y = (*this)(b);
  Math::Sort(x, y);

  // Compute derivative
  Polynomial p = Prime();

  // Find roots
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
      if (s > y)
      {
        y = s;
      }
    }
  }
}

/*!
\brief Compose the quartic by a quadric.

This function computes p ( q(x) ), where p denotes the quartic and q the quadric.
\param p The quartic.
\param q The quadric.
*/
Polynomial Polynomial::Compose(const Quartic& p, const Quadric& q)
{
  const double q10 = q[1] * q[0];
  const double q21 = q[2] * q[1];
  const double q20 = q[2] * q[0];
  const double q22 = q[2] * q[2];
  const double q11 = q[1] * q[1];
  const double q00 = q[0] * q[0];

  return Polynomial(p[1] * (q22 * q22),
    p[1] * (4.0 * q22 * q20),
    p[3] * (q22 * q[2]) + p[1] * (6.0 * q22 * q11 + 4.0 * q22 * q20),
    p[3] * (3.0 * q22 * q[1]) + p[1] * (4.0 * q21 * q11 + 12.0 * q22 * q10),
    p[2] * q22 + p[3] * (3.0 * q22 * q[0] + 3.0 * q[2] * q11) + p[1] * (6.0 * q22 * q00 + q11 * q11 + 12.0 * q20 * q11),
    p[2] * 2.0 * q21 + p[3] * (q11 * q[1] + 6.0 * q21 * q[0]) + p[1] * (4.0 * q11 * q10 + 12.0 * q21 * q00),
    p[1] * q[2] + p[2] * (2.0 * q20 + q11) + p[3] * (3.0 * q11 * q[0] + 3.0 * q[2] * q00) + p[1] * (6.0 * q11 * q00 + 4.0 * q20 * q00),
    p[1] * q[1] + p[2] * 2.0 * q10 + p[3] * (3.0 * q[1] * q00) + p[1] * (4.0 * q10 * q00),
    p[0] + p[1] * q[0] + p[2] * q00 + p[3] * (q00 * q[0]) + p[1] * (q00 * q00));
}