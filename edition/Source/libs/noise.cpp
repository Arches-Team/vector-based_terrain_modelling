// Noise

#include "libs/noise.h"

#include "libs/quintic.h"

/*
\brief Hashing function.
\param a,b Integer values
*/
inline short Noise::Hash1d(int a, int b) const
{
  return hashTable[(int)(a) ^ ((b) & 0xfffL)];
}

inline short Noise::Hash2d(int a, int b) const
{
  return hashTable[(int)(hashTable[(int)((a) & 0xfffL)] ^ ((b) & 0xfffL))];
}

inline short Noise::Hash3d(int a, int b, int c) const
{
  return hashTable[(int)(hashTable[(int)(hashTable[(int)((a) & 0xfffL)] ^ ((b) & 0xfffL))] ^ ((c) & 0xfffL))];
}

inline double Noise::IncrSum(int m, double s, double x, double y, double z) const
{
  return s * (rtable[m] * 0.5 + rtable[m + 1] * x + rtable[m + 2] * y + rtable[m + 3] * z);
}

inline double Noise::IncrSumP(double* mp, double s, double x, double y, double z) const
{
  return s * (mp[0] * 0.5 + mp[1] * x + mp[2] * y + mp[3] * z);
}

inline double Noise::GradientVertex(double* mp, double x, double y, double z) const
{
  return mp[0] * 0.5 + mp[1] * x + mp[2] * y + mp[3] * z;
}

inline double Noise::GradientVertex(double* mp, double x, double y) const
{
  return mp[0] * 0.5 + mp[1] * x + mp[2] * y;
}

inline Vector Noise::GradientVertexGrad(double* mp) const
{
  return Vector(mp[1], mp[2], mp[3]);
}
inline double Noise::ValueVertex(double* mp) const
{
  return mp[0] * 0.5;
}

/*!
\class Noise noise.h
\brief Perlin noise.

This class is used as a core class to create various types
of fractal noise, such as sums of scaled noises.

\sa SimplexNoise
*/

/*!
\brief Create a noise function.
\param quintic Flag defining whether quintic interpolation will be used, noise uses cubic interpolation if set to false.
*/
Noise::Noise(bool quintic)
{
  Noise::quintic = quintic;
}

/*!
\brief Compute noise at a given point.

\param p Point.
*/
double Noise::Value(const Vector& p) const
{
  Vector u = p;
  int ix = Math::Integer(u[0]);
  int iy = Math::Integer(u[1]);
  int iz = Math::Integer(u[2]);

  int jx = ix + 1;
  int jy = iy + 1;
  int jz = iz + 1;

  const double x = u[0] - ix;
  const double y = u[1] - iy;
  const double z = u[2] - iz;

  // Grid vertices
  int ixiy_hash = Hash2d(ix, iy);
  int jxiy_hash = Hash2d(jx, iy);
  int ixjy_hash = Hash2d(ix, jy);
  int jxjy_hash = Hash2d(jx, jy);

  int n000 = (int)Hash1d(ixiy_hash, iz) & 0xFF;
  int n100 = (int)Hash1d(jxiy_hash, iz) & 0xFF;
  int n010 = (int)Hash1d(ixjy_hash, iz) & 0xFF;
  int n110 = (int)Hash1d(jxjy_hash, iz) & 0xFF;

  int n001 = (int)Hash1d(ixiy_hash, jz) & 0xFF;
  int n101 = (int)Hash1d(jxiy_hash, jz) & 0xFF;
  int n011 = (int)Hash1d(ixjy_hash, jz) & 0xFF;
  int n111 = (int)Hash1d(jxjy_hash, jz) & 0xFF;

  // Gradients and values at grid vertices
  double g000 = GradientVertex(&rtable[n000], x, y, z);
  double g100 = GradientVertex(&rtable[n100], x - 1.0, y, z);
  double g010 = GradientVertex(&rtable[n010], x, y - 1.0, z);
  double g110 = GradientVertex(&rtable[n110], x - 1.0, y - 1.0, z);

  double g001 = GradientVertex(&rtable[n001], x, y, z - 1.0);
  double g101 = GradientVertex(&rtable[n101], x - 1.0, y, z - 1.0);
  double g011 = GradientVertex(&rtable[n011], x, y - 1.0, z - 1.0);
  double g111 = GradientVertex(&rtable[n111], x - 1.0, y - 1.0, z - 1.0);

  // Interpolation coefficients
  double sx, sy, sz;
  if (quintic)
  {
    sx = Quintic::Smooth(x);
    sy = Quintic::Smooth(y);
    sz = Quintic::Smooth(z);
  }
  else
  {
    sx = Cubic::Smooth(x);
    sy = Cubic::Smooth(y);
    sz = Cubic::Smooth(z);
  }

  // The complement values
  double tx = 1.0 - sx;
  double ty = 1.0 - sy;
  double tz = 1.0 - sz;

  // Factorize some interpolation terms
  double txty = tx * ty;
  double sxty = sx * ty;
  double txsy = tx * sy;
  double sxsy = sx * sy;

  double sum = (txty * tz) * g000 + (sxty * tz) * g100 + (txsy * tz) * g010 + (sxsy * tz) * g110 + (txty * sz) * g001 + (sxty * sz) * g101 + (txsy * sz) * g011 + (sxsy * sz) * g111;

  // Range at -0.5, 0.5
  sum = sum + 0.5;

  return sum;
}

/*!
\brief Compute the gradient of the noise function at a given point.

\param p Point.
*/
Vector Noise::Gradient(const Vector& p) const
{
  int ix = Math::Integer(p[0]);
  int iy = Math::Integer(p[1]);
  int iz = Math::Integer(p[2]);

  int jx = ix + 1;
  int jy = iy + 1;
  int jz = iz + 1;

  double x = p[0] - ix;
  double y = p[1] - iy;
  double z = p[2] - iz;

  // Grid vertices
  int ixiy_hash = Hash2d(ix, iy);
  int jxiy_hash = Hash2d(jx, iy);
  int ixjy_hash = Hash2d(ix, jy);
  int jxjy_hash = Hash2d(jx, jy);

  int n000 = (int)Hash1d(ixiy_hash, iz) & 0xFF;
  int n100 = (int)Hash1d(jxiy_hash, iz) & 0xFF;
  int n010 = (int)Hash1d(ixjy_hash, iz) & 0xFF;
  int n110 = (int)Hash1d(jxjy_hash, iz) & 0xFF;

  int n001 = (int)Hash1d(ixiy_hash, jz) & 0xFF;
  int n101 = (int)Hash1d(jxiy_hash, jz) & 0xFF;
  int n011 = (int)Hash1d(ixjy_hash, jz) & 0xFF;
  int n111 = (int)Hash1d(jxjy_hash, jz) & 0xFF;

  // Gradients at grid vertices
  Vector g000 = GradientVertexGrad(&rtable[n000]);
  Vector g100 = GradientVertexGrad(&rtable[n100]);
  Vector g010 = GradientVertexGrad(&rtable[n010]);
  Vector g110 = GradientVertexGrad(&rtable[n110]);

  Vector g001 = GradientVertexGrad(&rtable[n001]);
  Vector g101 = GradientVertexGrad(&rtable[n101]);
  Vector g011 = GradientVertexGrad(&rtable[n011]);
  Vector g111 = GradientVertexGrad(&rtable[n111]);

  // Values at grid vertices
  double v000 = ValueVertex(&rtable[n000]) + g000 * Vector(x, y, z);
  double v100 = ValueVertex(&rtable[n100]) + g100 * Vector(x - 1.0, y, z);
  double v010 = ValueVertex(&rtable[n010]) + g010 * Vector(x, y - 1.0, z);
  double v110 = ValueVertex(&rtable[n110]) + g110 * Vector(x - 1.0, y - 1.0, z);

  double v001 = ValueVertex(&rtable[n001]) + g001 * Vector(x, y, z - 1.0);
  double v101 = ValueVertex(&rtable[n101]) + g101 * Vector(x - 1.0, y, z - 1.0);
  double v011 = ValueVertex(&rtable[n011]) + g011 * Vector(x, y - 1.0, z - 1.0);
  double v111 = ValueVertex(&rtable[n111]) + g111 * Vector(x - 1.0, y - 1.0, z - 1.0);

  // Interpolation coefficients
  double sx, sy, sz;
  double sxp, syp, szp;
  if (quintic)
  {
    sx = Quintic::Smooth(p[0] - ix);
    sy = Quintic::Smooth(p[1] - iy);
    sz = Quintic::Smooth(p[2] - iz);
    sxp = 30.0 * x * x * Math::Sqr(1.0 - x);
    syp = 30.0 * y * y * Math::Sqr(1.0 - y);
    szp = 30.0 * z * z * Math::Sqr(1.0 - z);
  }
  else
  {
    sx = Cubic::Smooth(p[0] - ix);
    sy = Cubic::Smooth(p[1] - iy);
    sz = Cubic::Smooth(p[2] - iz);
    sxp = 6.0 * x * (1.0 - x);
    syp = 6.0 * y * (1.0 - y);
    szp = 6.0 * z * (1.0 - z);
  }

  // The complement values
  double tx = 1.0 - sx;
  double ty = 1.0 - sy;
  double tz = 1.0 - sz;

  double nx = g000[0]
    + sxp * (v100 - v000)
    + sx * (g100[0] - g000[0])
    + sy * (g010[0] - g000[0])
    + sz * (g001[0] - g000[0])
    + sxp * sy * (v110 - v010 - v100 + v000)
    + sx * sy * (g110[0] - g010[0] - g100[0] + g000[0])
    + sxp * sz * (v101 - v001 - v100 + v000)
    + sx * sz * (g101[0] - g001[0] - g100[0] - g000[0])
    + sy * sz * (g011[0] - g001[0] - g010[0] + g000[0])
    + sxp * sy * sz * (v111 - v011 - v101 + v001 - v110 + v010 + v100 - v000)
    + sx * sy * sz * (g111[0] - g011[0] - g101[0] + g001[0] - g110[0] + g010[0] + g100[0] - g000[0]);

  double ny = g000[1]
    + sx * (g100[1] - g000[1])
    + syp * (v010 - v000)
    + sy * (g010[1] - g000[1])
    + sz * (g001[1] - g000[1])
    + sx * syp * (v110 - v010 - v100 + v000)
    + sx * sy * (g110[1] - g010[1] - g100[1] + g000[1])
    + sx * sz * (g101[1] - g001[1] - g100[1] + g000[1])
    + syp * sz * (v011 - v001 - v010 + v000)
    + sy * sz * (g011[1] - g001[1] - g010[1] + g000[1])
    + sx * syp * sz * (v111 - v011 - v101 + v001 - v110 + v010 + v100 - v000)
    + sx * sy * sz * (g111[1] - g011[1] - g101[1] + g001[1] - g110[1] + g010[1] + g100[1] - g000[1]);

  double nz = g000[2]
    + sx * (g100[2] - g000[2])
    + sy * (g010[2] - g000[2])
    + szp * (v001 - v000)
    + sz * (g001[2] - g000[2])
    + sx * sy * (g110[2] - g010[2] - g100[2] + g000[2])
    + sx * szp * (v101 - v001 - v100 + v000)
    + sx * sz * (g101[2] - g001[2] - g100[2] + g000[2])
    + sy * szp * (v011 - v001 - v010 + v000)
    + sy * sz * (g011[2] - g001[2] - g010[2] + g000[2])
    + sx * sy * szp * (v111 - v011 - v101 + v001 - v110 + v010 + v100 - v000)
    + sx * sy * sz * (g111[2] - g011[2] - g101[2] + g001[2] - g110[2] + g010[2] + g100[2] - g000[2]);

  // Factorize some interpolation terms
  double txty = tx * ty;
  double sxty = sx * ty;
  double txsy = tx * sy;
  double sxsy = sx * sy;

  double sum = (txty * tz) * v000 + (sxty * tz) * v100 + (txsy * tz) * v010 + (sxsy * tz) * v110 + (txty * sz) * v001 + (sxty * sz) * v101 + (txsy * sz) * v011 + (sxsy * sz) * v111;

  // Range at -0.5, 0.5
  sum = sum + 0.5;

  return Vector(nx, ny, nz);
}

/*!
\brief Compute the gradient of the noise function at a given point.

Note that it is much faster to call this function rather than calling Noise::Gradient() and NoiseAt() consecutively.

\param p Point.
\param n Returned gradient.
\return The noise value.
*/
double Noise::AtGradient(const Vector& p, Vector& n) const
{
  int ix = Math::Integer(p[0]);
  int iy = Math::Integer(p[1]);
  int iz = Math::Integer(p[2]);

  int jx = ix + 1;
  int jy = iy + 1;
  int jz = iz + 1;

  double x = p[0] - ix;
  double y = p[1] - iy;
  double z = p[2] - iz;

  // Grid vertices
  int ixiy_hash = Hash2d(ix, iy);
  int jxiy_hash = Hash2d(jx, iy);
  int ixjy_hash = Hash2d(ix, jy);
  int jxjy_hash = Hash2d(jx, jy);

  int n000 = (int)Hash1d(ixiy_hash, iz) & 0xFF;
  int n100 = (int)Hash1d(jxiy_hash, iz) & 0xFF;
  int n010 = (int)Hash1d(ixjy_hash, iz) & 0xFF;
  int n110 = (int)Hash1d(jxjy_hash, iz) & 0xFF;

  int n001 = (int)Hash1d(ixiy_hash, jz) & 0xFF;
  int n101 = (int)Hash1d(jxiy_hash, jz) & 0xFF;
  int n011 = (int)Hash1d(ixjy_hash, jz) & 0xFF;
  int n111 = (int)Hash1d(jxjy_hash, jz) & 0xFF;

  // Gradients at grid vertices
  Vector g000 = GradientVertexGrad(&rtable[n000]);
  Vector g100 = GradientVertexGrad(&rtable[n100]);
  Vector g010 = GradientVertexGrad(&rtable[n010]);
  Vector g110 = GradientVertexGrad(&rtable[n110]);

  Vector g001 = GradientVertexGrad(&rtable[n001]);
  Vector g101 = GradientVertexGrad(&rtable[n101]);
  Vector g011 = GradientVertexGrad(&rtable[n011]);
  Vector g111 = GradientVertexGrad(&rtable[n111]);

  // Values at grid vertices
  double v000 = ValueVertex(&rtable[n000]) + g000 * Vector(x, y, z);
  double v100 = ValueVertex(&rtable[n100]) + g100 * Vector(x - 1.0, y, z);
  double v010 = ValueVertex(&rtable[n010]) + g010 * Vector(x, y - 1.0, z);
  double v110 = ValueVertex(&rtable[n110]) + g110 * Vector(x - 1.0, y - 1.0, z);

  double v001 = ValueVertex(&rtable[n001]) + g001 * Vector(x, y, z - 1.0);
  double v101 = ValueVertex(&rtable[n101]) + g101 * Vector(x - 1.0, y, z - 1.0);
  double v011 = ValueVertex(&rtable[n011]) + g011 * Vector(x, y - 1.0, z - 1.0);
  double v111 = ValueVertex(&rtable[n111]) + g111 * Vector(x - 1.0, y - 1.0, z - 1.0);

  // Interpolation coefficients and primes
  double sx, sy, sz;
  double sxp, syp, szp;
  if (quintic)
  {
    sx = Quintic::Smooth(p[0] - ix);
    sy = Quintic::Smooth(p[1] - iy);
    sz = Quintic::Smooth(p[2] - iz);
    sxp = 30.0 * x * x * Math::Sqr(1.0 - x);
    syp = 30.0 * y * y * Math::Sqr(1.0 - y);
    szp = 30.0 * z * z * Math::Sqr(1.0 - z);
  }
  else
  {
    sx = Cubic::Smooth(p[0] - ix);
    sy = Cubic::Smooth(p[1] - iy);
    sz = Cubic::Smooth(p[2] - iz);
    sxp = 6.0 * x * (1.0 - x);
    syp = 6.0 * y * (1.0 - y);
    szp = 6.0 * z * (1.0 - z);
  }

  // The complement values
  double tx = 1.0 - sx;
  double ty = 1.0 - sy;
  double tz = 1.0 - sz;

  double nx = g000[0]
    + sxp * (v100 - v000)
    + sx * (g100[0] - g000[0])
    + sy * (g010[0] - g000[0])
    + sz * (g001[0] - g000[0])
    + sxp * sy * (v110 - v010 - v100 + v000)
    + sx * sy * (g110[0] - g010[0] - g100[0] + g000[0])
    + sxp * sz * (v101 - v001 - v100 + v000)
    + sx * sz * (g101[0] - g001[0] - g100[0] - g000[0])
    + sy * sz * (g011[0] - g001[0] - g010[0] + g000[0])
    + sxp * sy * sz * (v111 - v011 - v101 + v001 - v110 + v010 + v100 - v000)
    + sx * sy * sz * (g111[0] - g011[0] - g101[0] + g001[0] - g110[0] + g010[0] + g100[0] - g000[0]);

  double ny = g000[1]
    + sx * (g100[1] - g000[1])
    + syp * (v010 - v000)
    + sy * (g010[1] - g000[1])
    + sz * (g001[1] - g000[1])
    + sx * syp * (v110 - v010 - v100 + v000)
    + sx * sy * (g110[1] - g010[1] - g100[1] + g000[1])
    + sx * sz * (g101[1] - g001[1] - g100[1] + g000[1])
    + syp * sz * (v011 - v001 - v010 + v000)
    + sy * sz * (g011[1] - g001[1] - g010[1] + g000[1])
    + sx * syp * sz * (v111 - v011 - v101 + v001 - v110 + v010 + v100 - v000)
    + sx * sy * sz * (g111[1] - g011[1] - g101[1] + g001[1] - g110[1] + g010[1] + g100[1] - g000[1]);

  double nz = g000[2]
    + sx * (g100[2] - g000[2])
    + sy * (g010[2] - g000[2])
    + szp * (v001 - v000)
    + sz * (g001[2] - g000[2])
    + sx * sy * (g110[2] - g010[2] - g100[2] + g000[2])
    + sx * szp * (v101 - v001 - v100 + v000)
    + sx * sz * (g101[2] - g001[2] - g100[2] + g000[2])
    + sy * szp * (v011 - v001 - v010 + v000)
    + sy * sz * (g011[2] - g001[2] - g010[2] + g000[2])
    + sx * sy * szp * (v111 - v011 - v101 + v001 - v110 + v010 + v100 - v000)
    + sx * sy * sz * (g111[2] - g011[2] - g101[2] + g001[2] - g110[2] + g010[2] + g100[2] - g000[2]);

  // Factorize some interpolation terms
  double txty = tx * ty;
  double sxty = sx * ty;
  double txsy = tx * sy;
  double sxsy = sx * sy;

  double v = (txty * tz) * v000 + (sxty * tz) * v100 + (txsy * tz) * v010 + (sxsy * tz) * v110 + (txty * sz) * v001 + (sxty * sz) * v101 + (txsy * sz) * v011 + (sxsy * sz) * v111;

  // Range at -0.5, 0.5
  v += 0.5;
  n = Vector(nx, ny, nz);
  return v;
}

/*!
\brief Return the global Lipschitz constant of the noise.
*/
double Noise::K() const
{
  // Did not compute the Lipschitz constant, should be bounded by 2
  return 2.0;
}

/*!
\brief Vector version of Perlin-style noise function.
\param p Point
*/
Vector Noise::AtVector(const Vector& p) const
{
  Vector result;

  int ix = Math::Integer(p[0]);
  int iy = Math::Integer(p[1]);
  int iz = Math::Integer(p[2]);
  int jx = ix + 1;
  int jy = iy + 1;
  int jz = iz + 1;

  double sx = Cubic::Smooth(p[0] - ix);
  double sy = Cubic::Smooth(p[1] - iy);
  double sz = Cubic::Smooth(p[2] - iz);

  // The complement values
  double tx = 1.0 - sx;
  double ty = 1.0 - sy;
  double tz = 1.0 - sz;

  // Interpolate
  double txty = tx * ty;
  double sxty = sx * ty;
  double txsy = tx * sy;
  double sxsy = sx * sy;
  int ixiy_hash = Hash2d(ix, iy);
  int jxiy_hash = Hash2d(jx, iy);
  int ixjy_hash = Hash2d(ix, jy);
  int jxjy_hash = Hash2d(jx, jy);

  double* mp = &rtable[int(Hash1d(ixiy_hash, iz)) & 0xFF];
  double px = p[0] - ix;
  double py = p[1] - iy;
  double pz = p[2] - iz;
  double s = txty * tz;
  result[0] = IncrSumP(mp, s, px, py, pz);
  mp += 4;
  result[1] = IncrSumP(mp, s, px, py, pz);
  mp += 4;
  result[2] = IncrSumP(mp, s, px, py, pz);

  mp = &rtable[int(Hash1d(jxiy_hash, iz)) & 0xFF];
  px = p[0] - jx;
  s = sxty * tz;
  result[0] += IncrSumP(mp, s, px, py, pz);
  mp += 4;
  result[1] += IncrSumP(mp, s, px, py, pz);
  mp += 4;
  result[2] += IncrSumP(mp, s, px, py, pz);

  mp = &rtable[int(Hash1d(jxjy_hash, iz)) & 0xFF];
  py = p[1] - jy;
  s = sxsy * tz;
  result[0] += IncrSumP(mp, s, px, py, pz);
  mp += 4;
  result[1] += IncrSumP(mp, s, px, py, pz);
  mp += 4;
  result[2] += IncrSumP(mp, s, px, py, pz);

  mp = &rtable[int(Hash1d(ixjy_hash, iz)) & 0xFF];
  px = p[0] - ix;
  s = txsy * tz;
  result[0] += IncrSumP(mp, s, px, py, pz);
  mp += 4;
  result[1] += IncrSumP(mp, s, px, py, pz);
  mp += 4;
  result[2] += IncrSumP(mp, s, px, py, pz);

  mp = &rtable[int(Hash1d(ixjy_hash, jz)) & 0xFF];
  pz = p[2] - jz;
  s = txsy * sz;
  result[0] += IncrSumP(mp, s, px, py, pz);
  mp += 4;
  result[1] += IncrSumP(mp, s, px, py, pz);
  mp += 4;
  result[2] += IncrSumP(mp, s, px, py, pz);

  mp = &rtable[int(Hash1d(jxjy_hash, jz)) & 0xFF];
  px = p[0] - jx;
  s = sxsy * sz;
  result[0] += IncrSumP(mp, s, px, py, pz);
  mp += 4;
  result[1] += IncrSumP(mp, s, px, py, pz);
  mp += 4;
  result[2] += IncrSumP(mp, s, px, py, pz);

  mp = &rtable[int(Hash1d(jxiy_hash, jz)) & 0xFF];
  py = p[1] - iy;
  s = sxty * sz;
  result[0] += IncrSumP(mp, s, px, py, pz);
  mp += 4;
  result[1] += IncrSumP(mp, s, px, py, pz);
  mp += 4;
  result[2] += IncrSumP(mp, s, px, py, pz);

  mp = &rtable[int(Hash1d(ixiy_hash, jz)) & 0xFF];
  px = p[0] - ix;
  s = txty * sz;
  result[0] += IncrSumP(mp, s, px, py, pz);
  mp += 4;
  result[1] += IncrSumP(mp, s, px, py, pz);
  mp += 4;
  result[2] += IncrSumP(mp, s, px, py, pz);

  return result;
}

