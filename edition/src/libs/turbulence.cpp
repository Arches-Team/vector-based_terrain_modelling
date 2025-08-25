
#include "libs/turbulence.h"

/*!
\class NoiseTurbulence turbulence.h
\brief %Turbulence computed from a modified noise function.

A turbulence field is characterized by the number of
octaves, the attenuation of the amplitude of the noise summed at every
step, and by the amplification of the wavelength.
*/

const Matrix NoiseTurbulence::R3 = Matrix::Rotation(Vector(4.0 * Math::Pi / 19.0, 3.0 * Math::Pi / 19.0, 2.0 * Math::Pi / 19.0));

/*!
\brief Creates a turbulence field.
\param v, a, l Base value, amplitude and wavelength.
\param alpha, lambda Amplitude and wavelength attenuation (should be less than 1.0 in general).
\param octaves Number of octaves.
\param t Offset vector.
*/
NoiseTurbulence::NoiseTurbulence(const double& v, const double& a, const double& l, const double& alpha, const double& lambda, int octaves, const Vector& t) :v0(v), a0(a), l0(l), alpha(alpha), lambda(lambda), octaves(octaves), t(t)
{
}

/*!
\brief Destroys a turbulence instance.
*/
NoiseTurbulence::~NoiseTurbulence()
{
}

/*!
\brief Computes the turbulence value at a given point.

Sums a series of scaled noise functions.
\param p Point.
*/
double NoiseTurbulence::Value(const Vector& p) const
{
  Vector q = p - t;

  // First wavelength
  double tu = v0;
  double l = l0;
  double a = a0;

  // Next waves
  for (int i = 0; i < octaves; i++)
  {
    tu += a * Noise::Value(q / l);
    l *= lambda;
    a *= alpha;
    q = R3 * q;
  }

  return tu;
}

/*!
\brief Computes a turbulence displacement vector at a given point.

Sums a series of scaled noise displacement functions.
\param p Point.
*/
Vector NoiseTurbulence::AtVector(const Vector& p) const
{
  Vector q = p - t;
  // First wave
  Vector tu = Vector::Null;

  double l = l0;
  double a = a0;
  // Next waves
  for (int i = 0; i < octaves; i++)
  {
    tu += a * Noise::AtVector(q / l);
    l *= lambda;
    a *= alpha;
    q = R3 * q;
  }
  return tu;
}

/*!
\brief Computes the maximum of the turbulence.

Sums the maximum of the terms of the series of scaled noise.
*/
double NoiseTurbulence::Maximum() const
{
  return v0 + a0 * Math::Geometric(alpha, octaves);
}

/*!
\brief Return the global Lipschitz constant.
\sa Noise::K()
*/
double NoiseTurbulence::K() const
{
  return (a0 * Noise::K() / l0) * Math::Geometric(alpha / lambda, octaves);
}

/*!
\brief Return the amplitude dampening coefficient.
*/
double NoiseTurbulence::Omega() const
{
  return alpha;
}

/*!
\brief Return the wavelength dampening coefficient.
*/
double NoiseTurbulence::Lambda() const
{
  return lambda;
}

/*!
\brief Number of octaves.
*/
int NoiseTurbulence::Octaves() const
{
  return octaves;
}

/*!
\brief Overloaded.
\param s Stream.
\param noise The noise.
*/
std::ostream& operator<<(std::ostream& s, const NoiseTurbulence& noise)
{
  s << "NoiseTurbulence(" << noise.v0 << ',' << noise.a0 << "," << noise.l0 << ',' << noise.alpha << ',' << noise.lambda << ',' << noise.octaves << ')';
  return s;
}

/*!
\class SimplexTurbulence turbulence.h
\brief Fractal Brownian motion from a simplex noise function.

A turbulence field is characterized by the number of
octaves, the attenuation of the amplitude of the noise summed at every
step, and by the amplification of the wavelength.
*/

const Matrix SimplexTurbulence::R3 = Matrix::Rotation(Vector(4.0 * Math::Pi / 19.0, 3.0 * Math::Pi / 19.0, 2.0 * Math::Pi / 19.0));

/*!
\brief Create a scaled turbulence.
\param v, a, l  Base value, amplitude and wavelength.
\param alpha Amplitude attenuation factor (~0.5, should be less than 1.0).
\param lambda Wavelength attenuation factor (~0.5, should be less than 1.0).
\param octaves Number of octaves.
\param t Translation offset.
*/
SimplexTurbulence::SimplexTurbulence(const double& v, const double& a, const double& l, const double& alpha, const double& lambda, int octaves, const Vector& t) :v0(v), a0(a), l0(l), alpha(alpha), lambda(lambda), octaves(octaves), t(t)
{
}

/*!
\brief Unit turbulence.
*/
SimplexTurbulence::SimplexTurbulence() :SimplexTurbulence(0.0, 1.0, 1.0)
{
}

/*!
\brief Computes the turbulence value at a given point.

Sums a series of scaled simplex noise functions.
\param p Point.
*/
double SimplexTurbulence::Value(const Vector& p) const
{
  return v0 + a0 * UnitAt(p / l0);
}

double SimplexTurbulence::GetAmplitude() const
{
  return a0;
}

double SimplexTurbulence::GetWavelength() const
{
  return l0;
}

/*!
\brief Return the maximum of the turbulence.
*/
double SimplexTurbulence::Maximum() const
{
  std::cout << v0 << std::endl;
  std::cout << (1.0 - pow(alpha, octaves)) << std::endl;
  std::cout << a0 << std::endl;
  return v0 + a0 * (1.0 - pow(alpha, octaves)) / (1.0 - alpha);
}

/*!
\brief Return the global Lipschitz constant.
\sa Noise::K()
*/
double SimplexTurbulence::K() const
{
  return  SimplexNoise::K() * Math::Geometric(alpha / lambda, octaves) * (a0 / l0);
}

/*!
\brief Return the amplitude dampening coefficient.
*/
double SimplexTurbulence::GetAlpha() const
{
  return alpha;
}

/*!
\brief Return the wavelength dampening coefficient.
*/
double SimplexTurbulence::GetLambda() const
{
  return lambda;
}

/*!
\brief Number of octaves.
*/
double SimplexTurbulence::GetOctaves() const
{
  return octaves;
}


/*!
\brief Computes the turbulence value at a given point.

Sums a series of scaled noise functions.
\param p Point.
*/
double SimplexTurbulence::UnitAt(const Vector& p) const
{
  Vector q = p - t;

  // First wavelength
  double tu = SimplexNoise::Value(q);
  double l = lambda;
  double a = alpha;

  // Next waves
  for (int i = 1; i < octaves; i++)
  {
    q = R3 * q;
    tu += a * SimplexNoise::Value(q / l);
    l *= lambda;
    a *= alpha;
  }

  return tu;
}



/*!
\brief Overloaded.
\param s Stream.
\param noise The noise.
*/
std::ostream& operator<<(std::ostream& s, const SimplexTurbulence& noise)
{
  s << "SimplexTurbulence(" << noise.v0 << ',' << noise.a0 << "," << noise.l0 << ',' << noise.alpha << ',' << noise.lambda << ',' << noise.octaves << ')';
  return s;
}
