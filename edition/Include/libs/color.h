// Color

#pragma once

// Qt Color
#include <QtGui/QColor>
#include <QtCore/QVector>

// Mathematics fundamentals
#include "libs/mathematics.h"
#include "libs/evector.h"


class Color
{
protected:
  double c[4] = { 0.0,0.0,0.0,1.0 }; //!< Array of color components; includes an alpha channel.
public:
  explicit Color(const double& = 0.0);
  explicit Color(unsigned long);
  explicit Color(const double&, const double&, const double&, const double& = 1.0);
  explicit Color(const Vector&, const double& = 1.0);
  explicit Color(int, int, int, int = 255);
  explicit Color(const QColor&);
  explicit Color(const QRgb&);

  //! Empty
  ~Color() {}

  Color Pow(const double& = 2.2) const;

  Color& operator+=(const Color&);

  Color Scale(const Color&) const;

  friend Color operator+(const Color&, const Color&);
  friend Color operator-(const Color&, const Color&);

  friend Color operator*(const Color&, const double&);
  friend Color operator*(const double&, const Color&);
  friend Color operator/(const Color&, const double&);

  Color& operator*= (const double&);

  static Color Min(const Color&, const Color&);
  static Color Max(const Color&, const Color&);

  static Color Lerp(const double&, const Color&, const Color&);
  Color Clamp(const Color & = Color::Black, const Color & = Color::White) const;
  double Luminance() const;

  unsigned long Cast() const;

  constexpr double& operator[] (int);
  constexpr double operator[] (int) const;

  QColor GetQt() const;

  Color Brighten(const double&) const;
  Color Darken(const double&) const;
  Color Saturate(const double&) const;

  // Stream
  friend std::ostream& operator<<(std::ostream&, const Color&);

  static Color Get(const QVector<Color>&, const double& t);

  static QColor LerpQt(const QColor&, const QColor&, const double&, int);

  static Color BiCubic(const double&, const double&, const Color&, const Color&, const Color&, const Color&, const Color&, const Color&, const Color&, const Color&, const Color&, const Color&, const Color&, const Color&, const Color & = Color::Black, const Color & = Color::Black, const Color & = Color::Black, const Color & = Color::Black);
  static Color Bilinear(const Color&, const Color&, const Color&, const Color&, const double&, const double&);

  static const Color Grey(const double&);
  static const Color Black; //!< %Black.
  static const Color White; //!< %White.
  static const Color Red; //!< %Red.
  static const Color Blue; //!< %Blue.

  static const Color Transparent; //!< Transparent.

  static Color Wheel(const double&);
};

/*!
\brief Returns the i<sup>th</sup> channel of the spectrum
\param i nummber of the channel queried (default=last)
*/
inline constexpr double& Color::operator[] (int i)
{
  return c[i];
}

/*!
\brief Returns a copy of the i<sup>th</sup> channel of the spectrum
\param i nummber of the channel queried (default=last)
*/
inline constexpr double Color::operator[] (int i) const
{
  return c[i];
}

//! Destructive scalar multiply.
inline Color& Color::operator*= (const double& a)
{
  c[0] *= a; c[1] *= a; c[2] *= a; c[3] *= a;
  return *this;
}

/*!
\brief Linear interpolation between two colors.

Interpolation is performed in RGB space.

Note that the following lines are equivalent:
\code
Color a,b;
double t;
Color c=(1.0-t)*a+t*b;
Color c=Lerp(t,a,b);
\endcode
*/
inline Color Color::Lerp(const double& t, const Color& a, const Color& b)
{
  return Color((1.0 - t) * a[0] + t * b[0], (1.0 - t) * a[1] + t * b[1], (1.0 - t) * a[2] + t * b[2], (1.0 - t) * a[3] + t * b[3]);
}

/*!
\brief Creates a greyscale color.

Initializes all the components to the given value except the opacity coefficient wich is set to 1.0.
\param v Grey value.
*/
inline Color::Color(const double& v) :Color(v, v, v, 1.0)
{
}

/*!
\brief Creates a color given a compact color representation
\param x Color compacted into an unsigned long.
*/
inline Color::Color(unsigned long x)
{
  c[0] = ((x >> 24) & 255) / 255.0;
  c[1] = ((x >> 16) & 255) / 255.0;
  c[2] = ((x >> 8) & 255) / 255.0;
  c[3] = (x & 255) / 255.0;
}

/*!
\brief Creates a color given each of these components
\param r,g,b Red, green and blue components.
\param a Alpha channel, set to 1.0 (opaque) as default.
*/
inline Color::Color(const double& r, const double& g, const double& b, const double& a)
{
  c[0] = r;
  c[1] = g;
  c[2] = b;
  c[3] = a;
}

/*!
\brief Creates a color.
\param v Vector, coordinates interpreted as red, green and blue components.
\param a Alpha channel, set to 1.0 (opaque) as default.
*/
inline Color::Color(const Vector& v, const double& a) :Color(v[0], v[1], v[2], a)
{
}

/*!
\brief Create a color given integer components
\param r, g, b Red, green and blue components.
\param a Alpha channel, set to 255 (opaque) as default.
*/
inline Color::Color(int r, int g, int b, int a) :Color(r / 255.0, g / 255.0, b / 255.0, a / 255.0)
{
}

/*!
\brief Overloaded operator.
*/
inline Color& Color::operator+=(const Color& ac)
{
  c[0] += ac[0];
  c[1] += ac[1];
  c[2] += ac[2];
  c[3] += ac[3];
  return *this;
}

/*!
\brief Overloaded sum operator
*/
inline Color operator+(const Color& u, const Color& v)
{
  return Color(u[0] + v[0], u[1] + v[1], u[2] + v[2], u[3] + v[3]);
}

/*!
\brief Scale.
\param v %Color.
*/
inline Color Color::Scale(const Color& v) const
{
  return Color(c[0] * v[0], c[1] * v[1], c[2] * v[2], c[3] * v[3]);
}

/*!
\brief Overloaded difference operator
*/
inline Color operator-(const Color& u, const Color& v)
{
  return Color(u[0] - v[0], u[1] - v[1], u[2] - v[2], u[3] - v[3]);
}

/*!
\brief Overloaded product by a scalar operator.
*/
inline Color operator*(const Color& u, const double& a)
{
  return Color(u[0] * a, u[1] * a, u[2] * a, u[3] * a);
}

/*!
\brief Overloaded product by a scalar operator

\param c Color.
\param a Scalar.
*/
inline Color operator*(const double& a, const Color& c)
{
  return c * a;
}

/*!
\brief Overloaded division by a scalar operator.
\param c Color.
\param a Scalar.
*/
inline Color operator/(const Color& c, const double& a)
{
  return c * (1.0 / a);
}

/*!
\brief Clamp a color between two bounds.
\param a, b %Color bounds.
*/
inline Color Color::Clamp(const Color& a, const Color& b) const
{
  return Color(Math::Clamp(c[0], a[0], b[0]), Math::Clamp(c[1], a[1], b[1]), Math::Clamp(c[2], a[2], b[2]), Math::Clamp(c[3], a[3], b[3]));
}

/*!
\brief Create a Qt color from a color.

Example how to write a pixel with a given color in an image:
\code
QImage image;
Color c;
image.setPixel(i, j, c.GetQt().rgb());
\endcode
*/
inline QColor Color::GetQt() const
{
  return QColor(int(255.0 * Math::Clamp(c[0])), int(255.0 * Math::Clamp(c[1])), int(255.0 * Math::Clamp(c[2])), int(255.0 * Math::Clamp(c[3])));
}

/*!
\brief Compute the luminance.
*/
inline double Color::Luminance() const
{
  return 0.212671 * c[0] + 0.715160 * c[1] + 0.072169 * c[2];
}

/*!
\brief Grey.

Same as Color(const double&)

\param g Grey level.
*/
inline const Color Color::Grey(const double& g)
{
  return Color(g);
}

// Hue, saturation, luminance
class Hsl {
protected:
  double h = 0.0, s = 1.0, l = 0.0; //!< Hue, saturation, luminance.
public:
  Hsl() {}
  explicit Hsl(const double&, const double&, const double&);
  explicit Hsl(const Color&);

  void Brighten(const double&);
  Color ToColor() const;
  friend class Color;
protected:
  static double ToRgb(double, double, double);
};

// 
class Xyz {
protected:
  double x = 0.0, y = 0.0, z = 0.0; //!< Color <b>xyz</b> coefficients.
public:
  Xyz() {}
  explicit Xyz(const double&, const double&, const double&);
  explicit Xyz(const Color&);

  Color ToColor() const;
  friend class Lab;
};

class Lab {
protected:
  double l = 0.0, a = 0.0, b = 0.0; //!< Color <b>%Lab</b> coefficients.
public:
  Lab() {}
  explicit Lab(const double&, const double&, const double&);
  explicit Lab(const Xyz&);

  Xyz ToXyz() const;
  friend class Xyz;
};

