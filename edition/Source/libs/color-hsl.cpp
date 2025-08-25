// Color 

#include "libs/color.h"

/*!
\class Hsl color.h
\brief %Color in HSL space.
\ingroup ColorGroup
*/

/*!
\brief Create a color in HSL space.
\param h Hue.
\param s Saturation.
\param l Lightness.
*/
Hsl::Hsl(const double& h, const double& s, const double& l) :h(h), s(s), l(l)
{
}

/*!
\brief Create a color in HSL space.
\param c Color in RGB space.
*/
Hsl::Hsl(const Color& c)
{
  const double r = c[0];
  const double g = c[1];
  const double b = c[2];

  double min = Math::Min(r, g, b);
  double max = Math::Max(r, g, b);
  double delta = max - min;

  l = (max + min) / 2.0;
  if (delta == 0.0)
  {
    h = s = 0.0;
  }
  else
  {
    if (l < 0.5)
    {
      s = delta / (max + min);
    }
    else
    {
      s = delta / (1 - fabs(2.0 * l - 1.0));
    }

    if (r == max)
    {
      h = (g - b) / delta;
    }
    else if (g == max) {
      h = (b - r) / delta + 2;
    }
    else if (b == max)
    {
      h = (r - g) / delta + 4;
    }
    h = fmod(60.0 * h + 360.0, 360.0);
    h /= 360.0;
  }
}

/*!
\brief Conversion function.
\param a, b Values.
\param h Hue.
*/
double Hsl::ToRgb(double a, double b, double h)
{
  if (h < 0.0)
  {
    h += 1.0;
  }
  if (h > 1.0)
  {
    h -= 1.0;
  }
  if (6.0 * h < 1.0) return a + (b - a) * 6.0 * h;
  if (2.0 * h < 1.0) return b;
  if (3.0 * h < 2.0) return a + (b - a) * (2.0 / 3.0 - h) * 6.0;
  return b;
}

/*!
\brief Convert a color in HSL space to RGB space.
*/
Color Hsl::ToColor() const
{
  if (s == 0.0)
  {
    return Color(l);
  }
  else
  {
    double v = (l < 0.5) ? (l * (1 + s)) : (l + s - (s * l));
    double u = 2 * l - v;

    return Color(ToRgb(u, v, h + 1.0 / 3.0), ToRgb(u, v, h), ToRgb(u, v, h - 1.0 / 3.0));
  }
}
