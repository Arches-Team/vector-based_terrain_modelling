// Noise

#pragma once

#include "libs/noise.h"
#include "libs/matrix.h"
#include "libs/voronoise.h"

class NoiseTurbulence :protected Noise {
protected:
  double v0 = 0.0; //!< Base value.
  double a0 = 1.0; //!< Base amplitude.
  double l0 = 1.0; //!< Wavelength.
  double alpha = 0.5; //!< Amplitude attenuation coefficient.
  double lambda = 0.5;  //!< Wavelength amplification coefficient.
  int octaves = 8; //!< Number of octaves.
  Vector t = Vector::Null; //!< Translation.
public:
  explicit NoiseTurbulence(const double& = 0.0, const double& = 1.0, const double& = 1.0, const double& = 0.5, const double& = 0.5, int = 5, const Vector & = Vector::Null);
  ~NoiseTurbulence();
  double Value(const Vector&) const;
  Vector AtVector(const Vector&) const;

  double Omega() const;
  double Lambda() const;
  int Octaves() const;

  double K() const;
  double Maximum() const;

  friend std::ostream& operator<<(std::ostream&, const NoiseTurbulence&);

protected:
  static const Matrix R3; //!< %Random incremental rotation.
};


class NoiseTurbulence2 :public Noise2 {
protected:
  double v0 = 0.0; //!< Base value.
  double a0 = 1.0; //!< Base amplitude.
  double l0 = 1.0; //!< Wavelength.
  double alpha = 0.5; //!< Amplitude attenuation coefficient.
  double lambda = 0.5;  //!< Wavelength amplification coefficient.
  int octaves = 8; //!< Number of octaves.
  Vector t = Vector::Null; //!< Translation.
public:
  explicit NoiseTurbulence2(const double& = 0.0, const double& = 1.0, const double& = 1.0, const double& = 0.5, const double& = 0.5, int = 5, const Vector & = Vector::Null);
  ~NoiseTurbulence2();
  double Value(const Vector2&) const;

  double Omega() const;
  double Lambda() const;
  int Octaves() const;

  double K() const;
  double Maximum() const;

  friend std::ostream& operator<<(std::ostream&, const NoiseTurbulence2&);

  using AnalyticScalarField2::Sample;
protected:
  static const Matrix2 R2; //!< %Random incremental rotation in the plane.
};

class SimplexTurbulence2 :protected SimplexNoise2
{
protected:
  double alpha = 0.5; //!< Amplitude attenuation coefficient.
  double lambda = 0.5;  //!< Wavelength amplification coefficient.
  int octaves = 8; //!< Number of octaves.
  Vector2 t = Vector2::Null; //!< Translation.
  double v0 = 0.0; //!< Base value.
  double a0 = 1.0; //!< Base amplitude.
  double l0 = 1.0; //!< Wavelength.
public:
  explicit SimplexTurbulence2(const double&, const double&, const double&, const double& = 0.5, const double& = 0.5, int = 8, const Vector2 & = Vector::Null);
  //explicit SimplexTurbulence2(const double&, const double&);
  explicit SimplexTurbulence2();

  //! Empty.
  ~SimplexTurbulence2() {}

  double Value(const Vector2&) const;

  double GetAmplitude() const;
  double GetWavelength() const;

  double GetAlpha() const;
  double GetLambda() const;
  double GetOctaves() const;

  double K() const;
  double Maximum() const;

  friend std::ostream& operator<<(std::ostream&, const SimplexTurbulence2&);

  using AnalyticScalarField2::Sample;
protected:
  static const Matrix2 R2; //!< %Random incremental rotation in the plane.
protected:
  double UnitAt(const Vector2&) const;
};


class CellularTurbulence2 :protected CellularNoise2
{
protected:
  double alpha = 0.5; //!< Amplitude attenuation coefficient.
  double lambda = 0.5;  //!< Wavelength amplification coefficient.
  int octaves = 4; //!< Number of octaves.
  Vector2 t = Vector2::Null; //!< Translation.
  double v0 = 0.0; //!< Base value.
  double a0 = 1.0; //!< Base amplitude.
  double l0 = 1.0; //!< Wavelength.
public:
  explicit CellularTurbulence2(const double&, const double&, const double&, const double& = 0.5, const double& = 0.5, int = 8, const Vector2 & = Vector2::Null);
  explicit CellularTurbulence2();

  //! Empty.
  ~CellularTurbulence2() {}

  double Value(const Vector2&) const;
  //using AnalyticScalarField2::Sample;
protected:
  double UnitAt(const Vector2&) const;
protected:
  static const Matrix2 R2; //!< %Random incremental rotation in the plane.
};

class SimplexTurbulence :protected SimplexNoise
{
protected:
  double alpha = 0.5; //!< Amplitude attenuation coefficient.
  double lambda = 0.5;  //!< Wavelength amplification coefficient.
  int octaves = 8; //!< Number of octaves.
  Vector t = Vector::Null; //!< Translation.
  double v0 = 0.0; //!< Base value.
  double a0 = 1.0; //!< Base amplitude.
  double l0 = 1.0; //!< Wavelength.
public:
  explicit SimplexTurbulence(const double&, const double&, const double&, const double& = 0.5, const double& = 0.5, int = 8, const Vector & = Vector::Null);
  explicit SimplexTurbulence();

  //! Empty.
  ~SimplexTurbulence() {}

  double Value(const Vector&) const;

  double GetAmplitude() const;
  double GetWavelength() const;

  double GetAlpha() const;
  double GetLambda() const;
  double GetOctaves() const;

  double K() const;
  double Maximum() const;

  friend std::ostream& operator<<(std::ostream&, const SimplexTurbulence&);

protected:
  static const Matrix R3; //!< %Random incremental rotation.
protected:
  double UnitAt(const Vector&) const;
};


