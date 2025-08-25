// Palettes

#pragma once

#include "libs/color.h"
#include "libs/box.h"

class QImage;
class QGraphicsScene;

class LookupPalette {
protected:
  QVector<Color> c; //!< Set of colors.
public:
  LookupPalette(const QVector<Color>&);
  virtual Color GetColor(int) const;
};

class GenericPalette {
public:
  virtual Color GetColor(double) const;
  // Draw
  virtual QImage Draw(int = 32, int = 1024) const;
  virtual QImage CreateImage(int, int = 16, bool = false) const;

  virtual void Draw(QGraphicsScene&, const Box2&) const;
};


class AnalyticPalette : public GenericPalette {
protected:
  int n = 0; //!< Palette identifier.
  bool r = false; //!< Reverse flag.
public:
  explicit AnalyticPalette(int = 0, bool = false);
  virtual Color GetColor(double) const;

public:
  static Color BrownGreyGreen(double);
  static Color GreenBrownGrey(double);
  static Color CoolWarm(double);
  static Color MatlabJet(double);
  static Color BlueGreen(double);
  static Color WhiteRed(double);
  static Color BlueGreyBrown(double);
  static Color GeologyGreenYellow(double);
  static Color GeologyGreenYellow2(double);
  static Color BrownGreen(double);
  static Color WhiteBlue(double);
  static Color WhiteBrown(double);
  static Color GreenOrange(double);

  AnalyticPalette Reversed() const;
private:
  static Color Diverging(const Color&, const Color&, const Color&, double);
  double Reverse(double) const;
};

/*!
\brief Return the reversed palette.
*/
inline AnalyticPalette AnalyticPalette::Reversed() const
{
  return AnalyticPalette(n, !r);
}

inline double AnalyticPalette::Reverse(double u) const
{
  return r == false ? u : 1.0 - u;
}


class QImage;

class Palette : public GenericPalette
{
protected:
  QVector<Color> c;  //!< %Array of colors.
  QVector<double> a; //!< Anchors.
  int type = 0;      //!< Type, used to speed-up queries.
public:
  Palette();
  Palette(const AnalyticPalette&, int = 16);
  Palette(const QVector<QColor>&);
  Palette(const QVector<Color>&);
  Palette(const QVector<Color>&, const QVector<double>&);
  Palette(const Color&, const Color&, const Color&, const double& = 0.5);
  Palette(const QImage&);

  Palette Reverse() const;
  void ScaleTo(const double& = 0.0, const double& = 1.0);

  Color GetColor(double) const;
  void Saturate(const double&);
  void Brighten(const double&);
public:
  static const Palette BrownDesert;
  static const Palette SimonGaussian;
  static const Palette ShadedRelief;
  static const Palette HugoShading;
  static const Palette JoshuaShading;
};


