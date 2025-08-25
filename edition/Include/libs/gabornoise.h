// Noise

#pragma once

#include "libs/scalarfield.h"

// Gabor noise
class GaborNoise :public AnalyticScalarField2
{
protected:
  double K;  //!< Amplitude.
  double a;  //!< Sigma.
  double F0; //!< Frequency.
  double r;  //!< Kernel radius.
  double d;  //!< Impulse density.
  unsigned int o; //!< Random offset.
  double c0, s0;  //!< Cosine and sine of pulsation.
public:
  GaborNoise(const double&, const double&, const double&, const double&, const double&, unsigned int);
  double Value(const Vector2&) const;
  double Cell(int, int, const double&, const double&) const;
  double Variance() const;
protected:
  double Gabor(const double&, const double&) const;
};


