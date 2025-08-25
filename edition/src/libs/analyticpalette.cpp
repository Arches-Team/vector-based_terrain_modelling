// Analytic palettes 

#include "libs/palette.h"
#include "libs/array.h"
#include "libs/linear.h"


/*!
\class GenericPalette palette.h

\brief Core class for palettes.
\ingroup ColorGroup
*/

/*!
\brief Compute the color.

\param u %Palette entry.
*/
#pragma warning(push)
#pragma warning(disable: 4100)  
Color GenericPalette::GetColor(double u) const
{
  return Color::Black;
}
#pragma warning(pop)

/*!
\brief Draws an image with the colors of the palette.

The palette is drawn vertically.

\param width, height Width and height of the image.
*/
QImage GenericPalette::Draw(int width, int height) const
{
  QImage image(width, height, QImage::Format_RGB32);

  for (int i = 0; i < height; i++)
  {
    QColor color = GetColor(Math::Unit(i, height)).GetQt();
    for (int j = 0; j < width; j++)
    {
      image.setPixelColor(QPoint(j, i), color);
    }
  }

  return image;
}


/*!
\brief Create an image form a palette.
\param s Size of the image.
\param t Thickness of color ramp.
\param v %Vertical flag, boolean set to true to generate a vertical ramp.
*/
QImage GenericPalette::CreateImage(int s, int t, bool v) const
{
  int w = s;
  int h = t;
  if (v)
  {
    Swap(h, w);
  }

  QImage image(w, h, QImage::Format_ARGB32);

  if (v)
  {
    for (int j = 0; j < h; j++)
    {
      QColor color = GetColor(Math::Unit(j, h)).GetQt();
      for (int i = 0; i < w; i++)
      {

        image.setPixelColor(i, j, color);
      }
    }
  }
  else
  {
    for (int i = 0; i < w; i++)
    {
      QColor color = GetColor(Math::Unit(i, w)).GetQt();
      for (int j = 0; j < h; j++)
      {
        image.setPixelColor(i, j, color);
      }

    }
  }
  return image;
}

void GenericPalette::Draw(QGraphicsScene& scene, const Box2& box) const
{
  const int n = 257;

  Array2 a(box, n, 2);
  for (int i = 0; i < n - 1; i++)
  {
    Color c = GetColor(Math::Unit(i, n));
    QBrush brush(c.GetQt());
    QPen pen(c.GetQt());
    pen.setWidthF(box.Perimeter() / 1000.0);
    Box2(a.Cell(i, 0)).Draw(scene, pen, brush);
  }
  QPen pen;
  pen.setWidthF(box.Perimeter() / 1000.0);
  box.Draw(scene, pen);
}

/*!
\class AnalyticPalette palette.h

\brief Procedurally defined palettes.

Palette entries are as follows:
<B>0</B> Brown-Grey-Green
\image html palette-0.png

<B>1</B> Green-Brown-Grey
\image html palette-1.png

<B>2</B> Matlab Jet
\image html palette-2.png

<B>3</B> Cool-Warm
\image html palette-3.png

<B>4</B> White-Red
\image html palette-4.png

<B>5</B> Blue-Green
\image html palette-5.png

<B>6</B> Blue-Grey-Brown
\image html palette-6.png

<B>7-8</B> Green-Yellow variants
\image html palette-7.png

\image html palette-8.png

\image html palette-9.png

\image html palette-10.png

\image html palette-11.png

\image html palette-12.png

<B>9</B> Brown-Green

<B>10</B> White-Blue


\ingroup ColorGroup
*/

/*!
\brief Create a palette.
\param n Palette identifier.
\param r Reverse flag.
*/
AnalyticPalette::AnalyticPalette(int n, bool r) :n(n), r(r)
{
}


/*!
\brief Compute color for a diverging palette with three colors.
\param a,b,c %Colors.
\param t Interpolation parameter.
*/
inline Color AnalyticPalette::Diverging(const Color& a, const Color& b, const Color& c, double t)
{
  if (t < 0.5)
  {
    return Color::Lerp(2.0 * t, a, c);
  }
  else
  {
    return Color::Lerp(2.0 * t - 1.0, c, b);
  }
}


/*!
\brief Custom divergent palette.

Compute a linear interpolation between brown, grey and green.
\param u Interpolant.
*/
Color AnalyticPalette::BrownGreyGreen(double u)
{
  static const Color Brown = Color(153, 93, 18);
  static const Color Green = Color(12, 112, 104);
  static const Color Grey = Color(244, 244, 244);

  return Diverging(Brown, Green, Grey, u);
}

/*!
\brief Custom divergent palette.

Inspired by Kenneth Moreland's diverging color map.

\image html palette-6.png

Compute a linear interpolation between blue, grey, and brown.
\param u Interpolant.
*/
Color AnalyticPalette::BlueGreyBrown(double u)
{
  static const Color Blue = Color(0.247, 0.522, 0.937);
  static const Color Grey = Color(0.863);
  static const Color Brown = Color(0.67, 0.502, 0.000);

  return Diverging(Blue, Brown, Grey, u);
}

/*!
\brief Custom divergent palette.

\param u Interpolant.
*/
Color AnalyticPalette::GeologyGreenYellow(double u)
{
  static const Color Cool(34.0 / 255.0, 52.0 / 255.0, 63.0 / 255.0);
  static const Color White(108.0 / 255.0, 184.0 / 255.0, 117.0 / 255.0);
  static const Color Warm(236.0 / 255.0, 237.0 / 255.0, 105 / 255.0);

  return Diverging(Cool, Warm, White, u);
}

/*!
\brief Custom divergent palette.

\param u Interpolant.
*/
Color AnalyticPalette::GeologyGreenYellow2(double u)
{
  static const Color White(108.0 / 255.0, 184.0 / 255.0, 117.0 / 255.0);
  static const Color Warm(236.0 / 255.0, 237.0 / 255.0, 105 / 255.0);
  return Color::Lerp(u, White, Warm);
}

/*!
\brief Diverging palette from brown to green.
*/
Color AnalyticPalette::BrownGreen(double u)
{
  static const Color Brown(153, 93, 18);
  static const Color White(244, 244, 244);
  static const Color Green(12, 112, 104);

  return Diverging(Brown, Green, White, u);
}

/*!
\brief Compute the color.

\param u %Palette entry.
*/
Color AnalyticPalette::GetColor(double u) const
{
  u = Reverse(u);

  switch (n)
  {
  case 0:
    return BrownGreyGreen(u);
    break;
  case 1:
    return GreenBrownGrey(u);
    break;
  case 2:
    return MatlabJet(u);
    break;
  case 3:
    return CoolWarm(u);
    break;
  case 4:
    return WhiteRed(u);
    break;
  case 5:
    return BlueGreen(u);
    break;
  case 6:
    return BlueGreyBrown(u);
    break;
  case 7:
    return GeologyGreenYellow(u);
    break;
  case 8:
    return GeologyGreenYellow2(u);
  case 9:
    return BrownGreen(u);
    break;
  case 10:
    return WhiteBlue(u);
    break;
  case 11:
    return WhiteBrown(u);
    break;
  case 12:
    return GreenOrange(u);
    break;
  }
  return Color::Black;
}

/*!
\brief Diverging palette from blue to red.

\param u %Palette entry.
*/
Color AnalyticPalette::CoolWarm(double u)
{
  static const Color Cool(97, 130, 234);
  static const Color White(221, 221, 221);
  static const Color Warm(220, 94, 75);

  if (u < 0.5)
  {
    return Color::Lerp(u / 0.5, Cool, White);
  }
  else
  {
    return Color::Lerp((u - 0.5) / 0.5, White, Warm);
  }
}

/*!
\brief Green to brown coloring.
\param u %Palette entry.
*/
Color AnalyticPalette::GreenBrownGrey(double u)
{
  static const Color Brown(167, 159, 118);
  static const Color Green(77, 97, 73);
  static const Color Tan(126, 128, 84);
  static const Color Grey(202, 195, 191);

  if (u < 0.33)
  {
    return Color::Lerp(u / 0.33, Green, Tan);
  }
  else if (u < 0.67)
  {
    return Color::Lerp((u - 0.33) / 0.33, Tan, Brown);
  }
  else
  {
    return Color::Lerp((u - 0.67) / 0.33, Brown, Grey);
  }
}

/*!
\brief Utility function for AnalyticPalette::MatlabJet
\author Mathieu Gaillard
*/
double MatlabJetBase(double val)
{
  if (val <= 0.125)
  {
    return 0.0;
  }
  else if (val <= 0.375)
  {
    return Linear::Step(val, 0.125, 0.375, 0.0, 1.0);
  }
  else if (val <= 0.625)
  {
    return 1.0;
  }
  else if (val <= 0.875)
  {
    return Linear::Step(val, 0.625, 0.875, 1.0, 0.0);
  }
  else
  {
    return 0.0;
  }
}

/*!
\brief Equivalent of the Jet coloring in Matlab.
\author Mathieu Gaillard
\param u %Palette entry.
*/
Color AnalyticPalette::MatlabJet(double u)
{
  double r = MatlabJetBase(u - 0.25);
  double g = MatlabJetBase(u);
  double b = MatlabJetBase(u + 0.25);

  return Color(r, g, b);
}

/*!
\brief Blue to green smooth palette.
\param u %Palette entry.
*/
Color AnalyticPalette::BlueGreen(double u)
{
  double r = 0.0;
  double g = 0.215 + 0.715 * u;
  double b = 0.825 - 0.625 * u;

  return Color(r, g, b);
}

/*!
\brief White to red.
\param u %Palette entry.
*/
Color AnalyticPalette::WhiteRed(double u)
{
  return Color::Lerp(u, Color::White, Color(1.0, 0.0, 0.0));
}

/*!
\brief White to blue.
\param u %Palette entry.
*/
Color AnalyticPalette::WhiteBlue(double u)
{
  return Color::Lerp(u, Color::White, Color(97, 130, 234));
}

Color AnalyticPalette::WhiteBrown(double u)
{
  return Color::Lerp(u, Color(1.0), Color(158 / 255.0, 137 / 255.0, 114 / 255.0));
}

Color AnalyticPalette::GreenOrange(double u)
{
  static const Color Orange(247 / 255.0, 150 / 255.0, 70 / 255.0);
  static const Color Green(146 / 255.0, 208 / 255.0, 80 / 255.0);

  if (u < 0.5)
    return Color::Lerp(u / 0.5, Orange, Color::White);
  else
    return Color::Lerp((u - 0.5) / 0.5, Color::White, Green);
}

LookupPalette::LookupPalette(const QVector<Color>& c) :c(c)
{
}

Color LookupPalette::GetColor(int i) const
{
  if ((i < 0) || (i >= c.size()))
  {
    return Color::Black;
  }
  return c.at(i);
}

