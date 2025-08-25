// Color 

#include "libs/color.h"

/*!
\class Color color.h
\brief %Color with alpha channel in RGB space.

Colors are represented using double.

Several constructors exist, two of them may be confusing:
\code
Color c(127,24,216); // Use constructor with integer parameters in [0,255] range.
Color c(0.5,0.1,0.84); // Use double parameters, in [0.0,1.0] range.
\endcode
\ingroup ColorGroup
*/

const Color Color::Black(0, 0, 0);
const Color Color::White(1.0, 1.0, 1.0);
const Color Color::Red(221, 220, 219);
const Color Color::Blue(97, 130, 234);

const Color Color::Transparent(0.0, 0.0, 0.0, 0.0);

/*!
\brief Create a color from Qt representation.
\param c %Color in Qt.
*/
Color::Color(const QColor& c)
{
  Color::c[0] = c.red() / 255.0;
  Color::c[1] = c.green() / 255.0;
  Color::c[2] = c.blue() / 255.0;
  Color::c[3] = c.alpha() / 255.0;
}

/*!
\brief Create a color from Qt representation.
\param c %Color in Qt.
*/
Color::Color(const QRgb& c)
{
  Color::c[0] = qRed(c) / 255.0;
  Color::c[1] = qGreen(c) / 255.0;
  Color::c[2] = qBlue(c) / 255.0;
  Color::c[3] = qAlpha(c) / 255.0;
}

/*!
\brief Returns the minimum value of each component
\param a, b Argument colors.
*/
Color Color::Min(const Color& a, const Color& b)
{
  return Color(Math::Min(a[0], b[0]), Math::Min(a[1], b[1]), Math::Min(a[2], b[2]), Math::Min(a[3], b[3]));
}

/*!
\brief Returns the maximum value of each component
\param a, b Argument colors.
*/
Color Color::Max(const Color& a, const Color& b)
{
  return Color(Math::Max(a[0], b[0]), Math::Max(a[1], b[1]), Math::Max(a[2], b[2]), Math::Max(a[3], b[3]));
}

/*!
\brief Computes the absolute distance between two colors in the system-space.

Alpha values are not used.
\param a, b Argument colors.
*/
double Norm(const Color& a, const Color& b)
{
  return (Math::Abs(a[0] - b[0]) + Math::Abs(a[1] - b[1]) + Math::Abs(a[2] - b[2]));
}

/*!
\brief Apply gamma correction.

\param gamma Exponential coefficient.
*/
Color Color::Pow(const double& gamma) const
{
  return Color(pow(c[0], gamma), pow(c[1], gamma), pow(c[2], gamma), c[3]);
}

/*!
\brief Returns a compact color format of the RGB-alpha system value of the color.
*/
unsigned long Color::Cast() const
{
  unsigned long a = ((unsigned long)(c[3] * 255.0)) << 24;
  unsigned long r = ((unsigned long)(c[0] * 255.0)) << 16;
  unsigned long g = ((unsigned long)(c[1] * 255.0)) << 8;
  unsigned long b = ((unsigned long)(c[2] * 255.0));
  return a | r | g | b;
}

/*!
\brief Compute a color in a discretized array of colors.
\param c Array of colors.
\param t Interpolation factor.
*/
Color Color::Get(const QVector<Color>& c, const double& t)
{
  if (c.size() == 0)
    return Color::Black;

  if (c.size() == 1)
    return c.at(0);

  if (t <= 0.0)
    return c.at(0);

  if (t >= 1.0)
    return c.at(c.size() - 1);

  double u = 1.0 / (c.size() - 1);

  // If t == 1
  int i = 0;
  double tt = t;
  while (tt >= u)
  {
    tt -= u;
    i++;
  }

  return Color::Lerp(tt / u, c.at(i), c.at(i + 1));
}

/*!
\brief Interpolate two Qt colors in the prescribed space.

The color space can be QColor::Rgb, QColor::Hsv, QColor::Cmyk or QColor::Hsl.

\param a, b Colors.
\param t Interpolation parameter.
\param s  Color space.
*/
QColor Color::LerpQt(const QColor& a, const QColor& b, const double& t, int s)
{
  if (s == QColor::Rgb)
  {
    return QColor(
      Math::Lerp(a.red(), b.red(), t),
      Math::Lerp(a.green(), b.green(), t),
      Math::Lerp(a.blue(), b.blue(), t),
      Math::Lerp(a.alpha(), b.alpha(), t));
  }
  else if (s == QColor::Hsl)
  {
    return QColor::fromHsl(
      Math::Lerp(a.hslHue(), b.hslHue(), t),
      Math::Lerp(a.hslSaturation(), b.hslSaturation(), t),
      Math::Lerp(a.lightness(), b.lightness(), t),
      Math::Lerp(a.alpha(), b.alpha(), t));
  }
  else if (s == QColor::Hsv)
  {
    return QColor::fromHsv(
      Math::Lerp(a.hsvHue(), b.hsvHue(), t),
      Math::Lerp(a.hsvSaturation(), b.hsvSaturation(), t),
      Math::Lerp(a.value(), b.value(), t),
      Math::Lerp(a.alpha(), b.alpha(), t));
  }
  else if (s == QColor::Cmyk)
  {
    return QColor::fromCmyk(
      Math::Lerp(a.cyan(), b.cyan(), t),
      Math::Lerp(a.magenta(), b.magenta(), t),
      Math::Lerp(a.yellow(), b.yellow(), t),
      Math::Lerp(a.black(), b.black(), t));
  }
  return QColor("transparent");
}

/*!
\brief Overloaded output-stream operator.
\param c Color.
\param s Stream.
*/
std::ostream& operator<<(std::ostream& s, const Color& c)
{
  s << "Vector(" << c.c[0] << ',' << c.c[1] << ',' << c.c[2] << ',' << c.c[3] << ')';
  return s;
}

/*!
\brief Bi-cubic interpolation between four values, given partial derivatives.

The values are given in trigonometric order.

\sa Math::BiCubic

\param u,v Interpolation coefficients.
\param a00,a10,a11,a01 Interpolated values.
\param u00,u10,u11,u01,v00,v10,v11,v01 Partial derivatives with respect to u and v.
\param x00,x10,x11,x01 Cross derivatives.
*/
Color Color::BiCubic(const double& u, const double& v, const Color& a00, const Color& a10, const Color& a11, const Color& a01, const Color& u00, const Color& u10, const Color& u11, const Color& u01, const Color& v00, const Color& v10, const Color& v11, const Color& v01, const Color& x00, const Color& x10, const Color& x11, const Color& x01)
{
  double u2 = u * u;
  double v2 = v * v;

  double bu0 = u2 * (2.0 * u - 3.0) + 1.0;
  double bu1 = u2 * (-2.0 * u + 3.0);
  double bu2 = (u2 - 2.0 * u + 1.0) * u;
  double bu3 = u2 * (u - 1.0);

  double bv0 = v2 * (2.0 * v - 3.0) + 1.0;
  double bv1 = v2 * (-2.0 * v + 3.0);
  double bv2 = (v2 - 2.0 * v + 1.0) * v;
  double bv3 = v2 * (v - 1.0);

  return bu0 * (a00 * bv0 + a01 * bv1 + v00 * bv2 + v01 * bv3) +
    bu1 * (a10 * bv0 + a11 * bv1 + v10 * bv2 + v11 * bv3) +
    bu2 * (u00 * bv0 + u01 * bv1 + x00 * bv2 + x01 * bv3) +
    bu3 * (u10 * bv0 + u11 * bv1 + x10 * bv2 + x11 * bv3);
}

/*!
\brief Bi-linear interpolation between four colors.

\sa Math::Bilinear

\param a00, a10, a11, a01 Interpolated colors.
\param u,v Interpolation coefficients.
*/
Color Color::Bilinear(const Color& a00, const Color& a10, const Color& a11, const Color& a01, const double& u, const double& v)
{
  return (1 - u) * (1 - v) * a00 + (1 - u) * v * a01 + u * (1 - v) * a10 + u * v * a11;
}

/*
\brief Compute the color wheel in the hue-saturation-value space.
\param a Angle, should be in unit interval.
*/
Color Color::Wheel(const double& a)
{
  return Hsl(a, 1.0, 1.0).ToColor();
}

/*
\brief Brighten the color.
\param c Coefficient added to luminance in Hsl space.
*/
Color Color::Brighten(const double& c) const
{
  Hsl hsl = Hsl(*this);
  hsl.l = Math::Clamp(hsl.l + c);
  return hsl.ToColor();
}

/*
\brief Darken the color.

Convenience function, this is the same as:
\code
  Color c;c=c.Darken(0.1); // Same as c=c.Brighten(-0.1);

\endcode
\sa Brighten

\param c Coefficient added to luminance in Hsl space.
*/
Color Color::Darken(const double& c) const
{
  return Brighten(-c);
}

/*
\brief Saturate the color.
\param c Coefficient added to saturation in Hsl space.
*/
Color Color::Saturate(const double& c) const
{
  Hsl hsl = Hsl(*this);
  hsl.s = Math::Clamp(hsl.s + c);
  return hsl.ToColor();
}