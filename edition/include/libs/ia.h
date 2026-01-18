// Interval arithmetic

#pragma once

#include "libs/mathematics.h"

class Ia {
private:
  double a = 0.0, b = 0.0; //!< Interval bounds.
public:
  //! Empty 
  Ia() { }
  explicit Ia(const double&);
  explicit Ia(const double&, const double&);
  explicit Ia(const Ia&, const Ia&);

  // Access Ia
  constexpr double& operator[] (int);
  constexpr double operator[] (int) const;

  // Unary operators
  Ia operator+ () const;
  Ia operator- () const;

  friend Ia operator+ (const Ia&, const Ia&);
  friend Ia operator- (const Ia&, const Ia&);

  friend Ia operator- (const Ia&, const double&);
  friend Ia operator+ (const Ia&, const double&);
  friend Ia operator- (const double&, const Ia&);
  friend Ia operator+ (const double&, const Ia&);

  friend Ia operator* (const Ia&, const Ia&);
  friend Ia operator/ (const Ia&, const Ia&);

  Ia& operator+= (const Ia&);
  Ia& operator-= (const Ia&);

  friend Ia operator* (const Ia&, const double&);
  friend Ia operator* (const double&, const Ia&);
  friend Ia operator/ (const Ia&, const double&);
  friend Ia operator/ (const double&, const Ia&);

  friend Ia Min(const Ia&, const Ia&);
  friend Ia Max(const Ia&, const Ia&);
  static Ia Lerp(const Ia&, const Ia&, const double&);

  friend Ia sqrt(const Ia&);
  friend Ia sqr(const Ia&);
  friend Ia abs(const Ia&);
  friend Ia sin(const Ia&);
  friend Ia cos(const Ia&);

  friend Ia pow(const Ia&, const double&);

  Ia AngleTwoPi() const;

  bool operator>(const double&) const;
  bool operator<(const double&) const;

  bool operator<(const Ia&) const;
  bool operator>(const Ia&) const;

  friend bool operator==(const Ia&, const Ia&);
  friend bool operator!=(const Ia&, const Ia&);

  friend std::ostream& operator<<(std::ostream&, const Ia&);

  Ia Intersection(const Ia&) const;
  bool Intersect(const Ia&) const;

  double Length() const;
  double Center() const;
  double Range() const;
public:
  static const Ia Empty;
  static const Ia Null;
  static const Ia Unit;
  static const Ia Infinity;
};

/*!
\brief Creates an interval.
\param x, y Low and high real values.
*/
inline Ia::Ia(const double& x, const double& y) :a(x), b(y)
{
}

/*!
\brief Creates an embedding interval.

It is the tightest embedding interval.

\param x,y Intervals.
*/
inline Ia::Ia(const Ia& x, const Ia& y) :a(Math::Min(x.a, y.a)), b(Math::Max(x.b, y.b))
{
}

/*!
\brief Creates an interval.

The length of the interval is 0.
\param x Real.
*/
inline Ia::Ia(const double& x) :a(x), b(x)
{
}

//! Access interval bounds.
inline constexpr double& Ia::operator[] (int i)
{
  if (i == 0) return a;
  else return b;
}

//! Access interval bounds.
inline constexpr double Ia::operator[] (int i) const
{
  if (i == 0) return a;
  else return b;
}

//! Overloaded.
inline Ia Ia::operator+ () const
{
  return *this;
}

/*!
\brief Compares two intervals.
\param x, y Argument intervals.
*/
inline bool operator==(const Ia& x, const Ia& y)
{
  return ((x.a == y.a) && (x.b == y.b));
}

/*!
\brief Compares two intervals.
\param x, y Intervals.
*/
inline bool operator!=(const Ia& x, const Ia& y)
{
  return ((x.a != y.a) || (x.b != y.b));
}

/*!
\brief Adds two intervals.
\param x, y Intervals.
*/
inline Ia operator+(const Ia& x, const Ia& y)
{
  return Ia(x.a + y.a, x.b + y.b);
}

/*!
\brief Subtracts two intervals.
\param x, y Intervals.
*/
inline Ia operator-(const Ia& x, const Ia& y)
{
  return Ia(x.a - y.b, x.b - y.a);
}

//! Destructive addition.
inline Ia& Ia::operator+=(const Ia& x)
{
  a += x.a;
  b += x.b;
  return *this;
}

//!  Subtracts two intervals.
inline Ia& Ia::operator-=(const Ia& x)
{
  a -= x.b;
  b -= x.a;
  return *this;
}

//! Overloaded.
inline Ia Ia::operator- () const
{
  return Ia(-b, -a);
}

/*!
\brief Overloaded.
\param x Interval.
\param a Real.
*/
inline Ia operator- (const double& a, const Ia& x)
{
  return Ia(a - x.b, a - x.a);
}

/*!
\brief Overloaded.
\param x Interval.
\param a Real.
*/
inline Ia operator+ (const double& a, const Ia& x)
{
  return x + a;
}

/*!
\brief Translates (shift right) an interval by a real value.
\param x Interval.
\param a Real.
*/
inline Ia operator+(const Ia& x, const double& a)
{
  return Ia(x.a + a, x.b + a);
}

/*!
\brief Translates (shift left) an interval by a real value.
\param x Interval.
\param a Real.
*/
inline Ia operator-(const Ia& x, const double& a)
{
  return Ia(x.a - a, x.b - a);
}

/*!
\brief Compute the length of the interval.

Same as:
\code
Ia i;
double l=i.Length(); // double l=i[1]-i[0];
\endcode
*/
inline double Ia::Length() const
{
  return b - a;
}

/*!
\brief Compute the center of the interval.

Same as:
\code
Ia i;
double l=i.Center(); // double l=0.5*(i[1]+i[0]);
\endcode
*/
inline double Ia::Center() const
{
  return 0.5 * (a + b);
}

/*!
\brief Compute the radius r of the interval [-r,r] embedding the interval.

Same as:
\code
Ia i;
double l=i.Range(); // double l=Math::Max(fabs(i[0]),fabs(i[1]));
\endcode
*/
inline double Ia::Range() const
{
  return Math::Max(fabs(a), fabs(b));
}
