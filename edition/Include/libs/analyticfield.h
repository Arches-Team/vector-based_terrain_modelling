#pragma once

#include "libs/heightfield.h"
#include "libs/turbulence.h"

class AnalyticHeightField :public AnalyticScalarField2 {
protected:
public:
  //! Empty.
  AnalyticHeightField() {}
  virtual double Height(const Vector2&) const;
  virtual Vector Normal(const Vector2&) const;
  HeightField CreateHeightField(const Box2&, int, int) const;

  // Slope
  double Slope(const Vector2&) const;
  double AverageSlope(const Vector2&, int = 8) const;
  double Slope(const Vector2&, const Vector2&) const;

  QVector<double> Cross(const Vector2&, const Vector2&, int) const;

public:
  static double Ridge(const double&, const double& = 0.0);
protected:
  static double Epsilon; //!< Epsilon value used for gradient computation.
};

/*!
\brief Ridge function.
\param z Elevation.
\param r Ridge elevation.
*/
inline double AnalyticHeightField::Ridge(const double& z, const double& r)
{
  if (z < r) return z;
  else return 2.0 * r - z;
}

/*!
\brief Compute elevation at a given point.
\param p Point.
*/
inline double AnalyticHeightField::Height(const Vector2& p) const
{
  return Value(p);
}

#include "libs/noise.h"

class AnalyticFieldDune : public AnalyticHeightField
{
protected:
  static Noise2 noise; //! %Noise function.
protected:
  double e;//!< Base elevation.
  double a; //!< Amplitude.
  double w; //!< Average wavelength between dunes.
public:
  AnalyticFieldDune(const double&, const double&, const double&);
  ~AnalyticFieldDune();
  double Value(const Vector2&) const;
  int Memory() const;
};


class AnalyticMusgrave :public AnalyticHeightField
{
protected:
  static const Noise2 noise; //! Base noise function, the same for all Musgrave fractal terrains.
protected:
  double z0;//!< Base elevation.
  double amp; //!< Amplitude.
  double wave; //!< Wavelength.
  double H; //! Hurst exponent.
  double lacunarity;
  int octaves; //!< Number of octaves.
  double gain;
  Vector2 tt;
public:
  explicit AnalyticMusgrave(const double&, const double&, const double&, const double&, int, const double&, const double&, const Vector2 & = Vector2::Null);

  double Value(const Vector2&) const;
public:
  static double Ridge(double);
};

/*!
\brief Musgrave ridge function.

Should be applied to a noise function value.
Ridges are created by taking the absolute value of the input signal, and squaring the result to amplify the ridged effect.

\param y Signal, should be a noise value in [-1,1].
*/
inline double AnalyticMusgrave::Ridge(double y)
{
  return Math::Sqr(1.0 - fabs(y));
}


class AnalyticRidged :public AnalyticHeightField, protected SimplexTurbulence2
{
protected:
public:
  explicit AnalyticRidged(const double&, const double&, const double&, const double&, const double&, int, const Vector2 & = Vector::Null);
  double Value(const Vector2&) const;
};

