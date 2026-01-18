#pragma once

#include "libs/heightfield.h"

class HeightFieldShader :public HeightField
{
protected:
public:
  explicit HeightFieldShader(const HeightField&);
  //! Empty
  ~HeightFieldShader() {}

  QImage Relief(bool = false) const;
  QImage ReliefLighting() const;
  QImage Mitsuba() const;

  QImage BrownShading() const;
  QImage BrownShading(double, double) const;

  QImage ShadedRelief() const;
  QImage ImageShadeRadianceScaling(const Color & = Color::White, const Color & = Color(142, 92, 31)) const;
  ScalarField2 RadianceScaling(const double&, const double&) const;
  QImage IsoLines(const QColor & = QColor(80, 60, 60), const QColor & = QColor(40, 20, 20), const double& = 100.0, int = 4) const;
  QImage IsoLinesAA(const QColor & = QColor(80, 60, 60), const QColor & = QColor(40, 20, 20), const double& = 100.0, int = 4) const;

  QImage CartographicGreen() const;
  QImage CartographicGreenYellow() const;

  QImage RiverWidth(const double&, int) const;

  QImage LaplacianShading(const Palette&, double = 1.6, double = 0.3) const;

};


