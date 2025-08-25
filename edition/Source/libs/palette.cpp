// Palette 

#include "libs/palette.h"
#include "libs/linear.h"
#include "libs/cpu.h"

/*!
\class Palette palette.h
\brief %Palette with array of colors.
\ingroup ColorGroup
*/

/*!
\brief Create an empty palette.
*/
Palette::Palette()
{
}

/*!
\brief Create a palette.
\param c Array of colors.
\param a Array of indices.
*/
Palette::Palette(const QVector<Color>& c, const QVector<double>& a) :c(c), a(a)
{
}

/*!
\brief Create a palette.
\param c Set of evenly distributed colors points.
*/
Palette::Palette(const QVector<Color>& c) :c(c)
{
  // Indexes
  for (int i = 0; i < c.size(); i++)
  {
    a.append(Math::Unit(i, c.size()));
  }
}

/*!
\brief Create a palette from an analytic palette.
\param p Analytic palette.
\param n Sampling.
*/
Palette::Palette(const AnalyticPalette& p, int n)
{
  for (int i = 0; i < n; i++)
  {
    c.append(p.GetColor(Math::Unit(i, n)));
  }

  // Indexes
  for (int i = 0; i < c.size(); i++)
  {
    a.append(Math::Unit(i, c.size()));
  }

  // Define as a diverging palette
  if (n == 3)
  {
    type = 1;
  }
}

/*!
\brief Create a palette.

This constructeur allows to create a simple color ramp easily:
\code
  Palette palette({ Color(0.5),Color(1.0) }); // Grey to white
\endcode

\param qtc Set of evenly distributed Qt colors points.
*/
Palette::Palette(const QVector<QColor>& qtc)
{
  for (int i = 0; i < c.size(); i++)
  {
    c.append(Color(qtc.at(i)));
  }

  // Indexes
  for (int i = 0; i < qtc.size(); i++)
  {
    a.append(Math::Unit(i, qtc.size()));
  }
}

/*!
\brief Create a diverging palette.
\param ca, cb, cc End colors and diverging color.
\param ab Diverging entry.
*/
Palette::Palette(const Color& ca, const Color& cb, const Color& cc, const double& ab) :c({ ca,cb,cc }), a({ 0.0,ab,1.0 })
{
  if (ab == 0.5)
  {
    type = 1;
  }
}

/*!
\brief Create a palette from an image.

The function chekcs the larger dimension of the image, and scans it accordingly.

\param image Input palette colors.
*/
Palette::Palette(const QImage& image)
{
  int n = max(image.height(), image.width());
  c.resize(n);

  if (image.height() > image.width())
  {
    for (int j = 0; j < image.height(); j++)
    {
      c.append(Color(image.pixel(0, j)));
    }
  }
  else
  {
    for (int i = 0; i < image.width(); i++)
    {
      c[i] = Color(image.pixel(i, 0));
    }
  }

  // Indexes
  for (int i = 0; i < n; i++)
  {
    a.append(Math::Unit(i, n));
  }
}

/*!
\brief Scale the range of values of the palette so that it lies in the prescribed interval.
\param x,y Interval.
*/
void Palette::ScaleTo(const double& x, const double& y)
{
  // Range
  const double a = Palette::a.at(0);
  const double b = Palette::a.at(Palette::a.size() - 1);

  // Indexes
  for (int i = 0; i < Palette::a.size(); i++)
  {
    Palette::a[i] = Math::Lerp(x, y, Linear::Step(Palette::a.at(i), a, b));
  }
}

/*!
\brief Compute a color in the palette.
\param t Interpolation parameter.
*/
Color Palette::GetColor(double t) const
{
  // Diverging palette with three colors over unit interval
  if (type == 1)
  {
    if (t < 0.5)
    {
      return Color::Lerp(2.0 * t, c.at(0), c.at(1));
    }
    else
    {
      return Color::Lerp(2.0 * t - 1.0, c.at(1), c.at(2));
    }
  }

  // General case
  if (c.size() == 0)
    return Color::White;

  if (c.size() == 1)
    return c.at(0);

  if (t < a.at(0))
    return c.at(0);

  if (t > a.at(a.size() - 1))
    return c.at(c.size() - 1);

  for (int i = 0; i < c.size() - 1; i++)
  {
    if (t < a.at(i + 1))
    {
      double s = Linear::Step(t, a.at(i), a.at(i + 1));
      return Color::Lerp(s, c.at(i), c.at(i + 1));
    }
  }
  return Color::White;
}


const Palette Palette::BrownDesert(QImage(System::GetArchesLib() + QString("/LibTerra/Gradients/browndesert.png")));

const Palette Palette::ShadedRelief(
  { Color(0.627, 0.863, 0.411), Color(1.00, 0.90, 0.45), Color(0.659, 0.607, 0.541), Color(0.95, 0.95, 0.95) },
  { 0.0,150.0,250.0,400.0 });

const Palette Palette::HugoShading(
  { Color(0.45, 0.52, 0.30), Color(0.57, 0.55, 0.38), Color(0.75, 0.71, 0.48), Color(0.75, 0.71, 0.48),
    Color(0.85, 0.80, 0.59),  Color(0.97, 0.89, 0.71), Color(0.98, 0.85, 0.66),Color(0.76, 0.69, 0.54),
    Color(0.73, 0.67, 0.57), Color(0.70, 0.65, 0.61), Color(0.97, 0.97, 0.97) });

const Palette Palette::JoshuaShading(
  { Color(0.54, 0.54, 0.40), Color(0.59, 0.60, 0.44), Color(0.67, 0.64, 0.47), Color(0.76, 0.69, 0.54),
    Color(0.79, 0.74, 0.59), Color(0.84, 0.79, 0.61), Color(0.95, 0.86, 0.70), Color(0.97, 0.97, 0.97) });

const Palette Palette::SimonGaussian(
    { Color(0.627, 0.763, 0.411), Color(1.00, 0.90, 0.45),Color(0.659, 0.607, 0.541), Color(0.95, 0.95, 0.95) });


/*!
\brief Saturate all colors of the palette.
\param x Coefficient.
*/
void Palette::Saturate(const double& x)
{
  for (int i = 0; i < c.size(); i++)
  {
    c[i] = c[i].Saturate(x);
  }
}

/*!
\brief Brighten all colors of the palette.
\param x Coefficient.
*/
void Palette::Brighten(const double& x)
{
  for (int i = 0; i < c.size(); i++)
  {
    c[i] = c[i].Brighten(x);
  }
}

/*!
\brief Reverse the enrties of the palette.
*/
Palette Palette::Reverse() const
{
  QVector<double> ra(a.size());
  QVector<Color> rc(c.size());

  // Indexes
  for (int i = 0; i < a.size(); i++)
  {
    ra[i] = a[a.size() - 1 - i];
  }
  // Colors
  for (int i = 0; i < c.size(); i++)
  {
    rc[i] = c[c.size() - 1 - i];
  }
  return Palette(rc, ra);
}

