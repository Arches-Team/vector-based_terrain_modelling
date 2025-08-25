// Interval arithmetic

#include "libs/ia.h"

/*!
\class Ia ia.h
\brief Intervals of two reals.

Most interval arithmetics operators are intensive because
of the tests needed to sort the extrema of the resulting
interval. Only addition and subtraction need not use tests and,
therefore, are efficient.

Several functions have been implemented with the same syntax as
for doubles, such as sqrt() and abs().
\ingroup MathGroup

*/

/*!
\brief Infinity interval.
*/
const Ia Ia::Infinity(-1e16, 1e16);

/*!
\brief Empty interval.
*/
const Ia Ia::Empty(1.0, -1.0);

/*!
\brief Null interval is zero-set.
*/
const Ia Ia::Null(0.0, 0.0);

/*!
\brief Unit interval [0,1].
*/
const Ia Ia::Unit(0.0, 1.0);

/*!
\brief Scales a interval by a double.

\param x Interval.
\param a Real.
*/
Ia operator*(const Ia& x, const double& a)
{
  if (a > 0.0)
  {
    return Ia(a * x.a, a * x.b);
  }
  else
  {
    return Ia(a * x.b, a * x.a);
  }
}

/*!
\brief Left multiply, same as scaling.
\param x Interval.
\param a Real.
*/
Ia operator*(const double& a, const Ia& x)
{
  if (a > 0.0)
  {
    return Ia(a * x.a, a * x.b);
  }
  else
  {
    return Ia(a * x.b, a * x.a);
  }
}

/*!
\brief Left divide by a double, same as scaling.
\param x Interval.
\param a Real.
*/
Ia operator/(const Ia& x, const double& a)
{
  if (a > 0.0)
  {
    return Ia(x.a / a, x.b / a);
  }
  else
  {
    return Ia(x.b / a, x.a / a);
  }
}

/*!
\brief Multiplies two intervals.
\param x, y Intervals.
*/
Ia operator*(const Ia& x, const Ia& y)
{
  if (x.a >= 0.0)
  {
    if (y.a >= 0.0)
    {
      return Ia(x.a * y.a, x.b * y.b);
    }
    else if (y.b <= 0.0)
    {
      return Ia(x.b * y.a, x.a * y.b);
    }
    else
    {
      return Ia(x.b * y.a, x.b * y.b);
    }
  }
  else if (x.b <= 0.0)
  {
    if (y.a >= 0.0)
    {
      return Ia(x.a * y.b, x.b * y.a);
    }
    else if (y.b <= 0.0)
    {
      return Ia(x.b * y.b, x.a * y.a);
    }
    else
    {
      return Ia(x.a * y.b, x.a * y.a);
    }
  }
  else
  {
    if (y.a >= 0.0)
    {
      return Ia(x.a * y.b, x.b * y.b);
    }
    else if (y.b <= 0.0)
    {
      return Ia(x.b * y.a, x.a * y.a);
    }
    else
    {
      return Ia(Math::Min(x.a * y.b, x.b * y.a), Math::Max(x.a * y.a, x.b * y.b));
    }
  }
}

/*!
\brief Divides two intervals.
\param x, y Intervals.
*/
Ia operator/(const Ia& x, const Ia& y)
{
  Ia z;
  if (x.a >= 0.0)
  {
    if (y.a >= 0.0)
    {
      z.a = x.a / y.b;

      z.b = x.b / y.a;
    }
    else if (y.b <= 0.0)
    {
      z.a = x.b / y.b;

      z.b = x.a / y.a;
    }
  }
  else if (x.b <= 0.0)
  {
    if (y.a >= 0.0)
    {
      z.a = x.a / y.a;

      z.b = x.b / y.b;
    }
    else if (y.b <= 0.0)
    {
      z.a = x.b / y.a;

      z.b = x.a / y.b;
    }
  }
  else
  {
    if (y.a >= 0.0)
    {
      z.a = x.a / y.a;

      z.b = x.b / y.a;
    }
    else if (y.b <= 0.0)
    {
      z.a = x.b / y.b;

      z.b = x.a / y.b;
    }
  }
  return z;
}

/*!
\brief Compute the squared interval.

This function is more efficient than multiplication.

\param x Interval.
*/
Ia sqr(const Ia& x)
{
  if (x.a >= 0.0)
  {
    return Ia(x.a * x.a, x.b * x.b);
  }
  else if (x.b <= 0.0)
  {
    return Ia(x.b * x.b, x.a * x.a);
  }
  else
  {
    return Ia(0.0, Math::Max(x.a * x.a, x.b * x.b));
  }
}

/*!
\brief Divides a real number by an interval.
\param a Real.
\param x Interval.
*/
Ia operator/(const double& a, const Ia& x)
{
  if (a > 0.0)
  {
    return Ia(a / x.b, a / x.a);
  }
  else
  {
    return Ia(a / x.a, a / x.b);
  }
}

/*!
\brief Computes the square root of an interval.

Interval should have positive bounds.
\param x Interval.
*/
Ia sqrt(const Ia& x)
{
  return Ia(sqrt(x.a), sqrt(x.b));
}

/*!
\brief Computes the absolute value of an interval.
\param x Interval.
*/
Ia abs(const Ia& x)
{
  if (x.a >= 0.0)
  {
    return x;
  }
  else if (x.b <= 0.0)
  {
    return Ia(-x.b, -x.a);
  }
  else
  {
    return Ia(0.0, Math::Max((-x.a), x.b));
  }
}

/*!
\brief Caculate the minimal 2 Math::Pi periodic interval.
*/
Ia Ia::AngleTwoPi() const
{
  double l = b - a;

  double t = floor(a / (Math::TwoPi));

  double x = a - (t * Math::TwoPi);

  return Ia(x, x + l);
}

/*!
\brief Computes the cosine of an interval.
\param x Interval.
*/
Ia cos(const Ia& x)
{
  double a = x.a / Math::Pi;
  double b = x.b / Math::Pi;

  int ua = int(ceil(a));

  if (1 + ua < b)
  {
    return Ia(-1.0, 1.0);
  }

  double ca = cos(x.a);
  double cb = cos(x.b);
  if (ua < b)
  {
    if (ua % 2 == 1)
    {
      return Ia(-1.0, Math::Max(ca, cb));
    }
    else
    {
      return Ia(Math::Min(ca, cb), 1.0);
    }
  }
  else
  {
    return Ia(Math::Min(ca, cb), Math::Max(ca, cb));
  }
}


/*!
\brief Computes the sine of an interval.

\sa cos(const Ia&)

\param x Interval.
*/
Ia sin(const Ia& x)
{
  return -cos(x + Math::Pi / 2);
}

/*!
\brief Computes the maximum of two intervals.
\param x, y Intervals.
*/
Ia Max(const Ia& x, const Ia& y)
{
  return Ia(Math::Max(x.a, y.a), Math::Max(x.b, y.b));
}

/*!
\brief Computes the minimum of two intervals.
\param x,y Intervals.
*/
Ia Min(const Ia& x, const Ia& y)
{
  return Ia(Math::Min(x.a, y.a), Math::Min(x.b, y.b));
}

/*!
\brief Linear interpolation of two intervals.

\param x, y Intervals.
\param t Real.

This code is equivalent but more efficient than an interpolation with overeloaded operators, as coefficients are know to be within unit interval.
\code
const double t=0.25;
Ia z=(1.0-t)*Ia(2.0,3.0)+t*Ia(-1.0,2.0); // Linear interpolation with overloaded operators
\endcode
*/
Ia Ia::Lerp(const Ia& x, const Ia& y, const double& t)
{
  return Ia((1.0 - t) * x.a + t * y.a, (1.0 - t) * x.b + t * y.b);
}

/*!
\brief Computes the power of an interval.
\param x Interval.
\param n Power.
*/
Ia pow(const Ia& x, const double& n)
{
  double a = (x.a == 0.0) ? 0.0 : ((x.a > 0.0) ? pow(x.a, n) : -pow(-x.a, n));
  double b = (x.b == 0.0) ? 0.0 : ((x.b > 0.0) ? pow(x.b, n) : -pow(-x.b, n));
  if (a < b)
    return Ia(a, b);
  else
    return Ia(b, a);
}

/*!
\brief Overloaded output-stream operator.
\param s Stream.
\param x Interval.
*/
std::ostream& operator<<(std::ostream& s, const Ia& x)
{
  s << "Ia(" << x.a << ',' << x.b << ')' << std::endl;
  return s;
}

/*!
\brief Compute the intersection of two intervals.
\param x Other interval.
*/
Ia Ia::Intersection(const Ia& x) const
{
  if ((b < x.a) || (a > x.b))
  {
    return Ia::Empty;
  }
  return Ia(Math::Max(a, x.a), Math::Min(b, x.b));

}

/*!
\brief Check if two intervals intersect.
\param x Other interval.
*/
bool Ia::Intersect(const Ia& x) const
{
  if ((b < x.a) || (a > x.b))
  {
    return false;
  }
  return true;
}

/*!
\brief Compare interval to a real.

Check if the entire interval is superior to the real. This is equivalent to:
\code
Ia i;
double x;
if (i[0]>x) {} // if (i>x) {}
\endcode

\param x Real.
*/
bool Ia::operator>(const double& x) const
{
  return a > x;
}

/*!
\brief Compare interval to a real.

Check if the entire interval is superior to the real. This is equivalent to:
\code
Ia i;
double x;
if (i[1]<x) {} // if (i<x) {}
\endcode

\param x Real.
*/
bool Ia::operator<(const double& x) const
{
  return b < x;
}

/*!
\brief Compare two intervals.

\param y Interval.
*/
bool Ia::operator>(const Ia& y) const
{
  return a > y.b;
}

/*!
\brief Compare two intervals.

\param y Interval.
*/
bool Ia::operator<(const Ia& y) const
{
  return b < y.a;
}
