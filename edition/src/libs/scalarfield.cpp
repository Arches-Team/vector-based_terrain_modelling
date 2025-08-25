// Fields

#include "libs/scalarfield.h"
#include "libs/cubic.h"
#include <immintrin.h>
#include "libs/cpu.h"
#include <omp.h>
#include "libs/mesh.h"

/*!
\class ScalarField scalarfield.h
\brief A base three-dimensional field of real values.

\ingroup StructureGroup
*/

/*!
\brief Create the field structure.
\param a The array.
\param v Default value of field.
\param s Sign, by default, set sign convention to true, i.e., negative values inside and positive outside.
*/
ScalarField::ScalarField(const Array& a, const double& v, bool s) :Array(a), sign(s)
{
  field.fill(v, nx * ny * nz);
}

/*!
\brief Create the field structure.
\param box The box.
\param x,y,z Size of the array.
\param v Default value of field.
\param s Sign, by default, set sign convention to true, i.e., negative values inside and positive outside.
*/
ScalarField::ScalarField(const Box& box, int x, int y, int z, const double& v, bool s) :Array(box, x, y, z), sign(s)
{
  field.fill(v, nx * ny * nz);
}

/*!
\brief Create the field structure.
\param box The box.
\param x,y,z Size of the array.
\param v Array of scalar values.
\param s Sign, by default, set sign convention to true, i.e., negative values inside and positive outside.
*/
ScalarField::ScalarField(const Box& box, int x, int y, int z, const QVector<double>& v, bool s) :Array(box, x, y, z), sign(s)
{
  field = v;
}

/*!
\brief Compute the gradient at a given sample.

\param i,j,k Integer coordinates of the sample.
*/
Vector ScalarField::Gradient(int i, int j, int k) const
{
  Vector n;

  // Gradient along x axis
  if (i == 0)
  {
    n[0] = (at(i + 1, j, k) - at(i, j, k)) * inversecelldiagonal[0];
  }
  else if (i == nx - 1)
  {
    n[0] = (at(i, j, k) - at(i - 1, j, k)) * inversecelldiagonal[0];
  }
  else
  {
    n[0] = (at(i + 1, j, k) - at(i - 1, j, k)) * 0.5 * inversecelldiagonal[0];
  }

  // Gradient along y axis
  if (j == 0)
  {
    n[1] = (at(i, j + 1, k) - at(i, j, k)) * inversecelldiagonal[1];
  }
  else if (j == ny - 1)
  {
    n[1] = (at(i, j, k) - at(i, j - 1, k)) * inversecelldiagonal[1];
  }
  else
  {
    n[1] = (at(i, j + 1, k) - at(i, j - 1, k)) * 0.5 * inversecelldiagonal[1];
  }

  // Gradient along z axis
  if (k == 0)
  {
    n[2] = (at(i, j, k + 1) - at(i, j, k)) * inversecelldiagonal[2];
  }
  else if (j == ny - 1)
  {
    n[2] = (at(i, j, k) - at(i, j, k - 1)) * inversecelldiagonal[2];
  }
  else
  {
    n[2] = (at(i, j, k + 1) - at(i, j, k - 1)) * 0.5 * inversecelldiagonal[2];
  }

  return n;
}

/*!
\brief Compute the normal.

\param p Point (should be on the surface).
*/
Vector ScalarField::Normal(const Vector& p) const
{
  Vector normal = Normalized(Gradient(p));
  if (sign)
  {
    return normal;
  }
  else
  {
    return -normal;
  }
}

/*!
\brief Access to array of data.
\param i,j,k Integer coordinates.
*/
double ScalarField::Value(int i, int j, int k) const
{
  return field.at(VertexIndex(i, j, k));

}

/*!
\brief Return the intensity at a point.

The function returns -1.0 if point sample is outside
the bounding box. Otherwise, the intensity is computed
by blending the intensity at the vertices of the voxel
grid with cubic polynomials.

\param p Point.
*/
double ScalarField::BicubicValue(const Vector& p) const
{
  // Compute the cell index and local coordinates
  Vector u;
  int i, j, k;
  CellInteger(p, i, j, k, u[0], u[1], u[2]);

  if (!InsideCellIndex(i, j, k))
    return -1.0;

  int c = VertexIndex(i, j, k);

  // Store bernstein coefficients
  double x[4];
  double y[4];
  double z[4];

  Blend(u[0], x);
  Blend(u[1], y);
  Blend(u[2], z);

  // Offset to lower diagonal vertex
  c -= nx * ny + nx + 1;

  // Sum
  double r = 0.0;

  for (i = 0; i < 4; i++)
  {
    for (j = 0; j < 4; j++)
    {
      for (k = 0; k < 4; k++)
      {
        if ((c > 0) && (c < nx * ny * nz))
        {
          r += field[c] * x[i] * y[j] * z[k];
        }
        c += nx * ny;
      }
      c += nx - 4 * nx * ny;
    }
    c += 1 - 4 * nx;
  }
  return r;
}

/*!
\brief Compute the blending coefficients using cubic polynomials.

After Allen Van Gelder and Jane Wilhels, <I>Topological considerations in iso-surface generation</I>,
<B>ACM Transactions on Graphics</B>, <B>13</B>(4):21-28, 1994.
*/
void ScalarField::Blend(const double& u, double b[4]) const
{
  double a = u - 1.0;
  b[0] = -0.5 * u * a * a;
  b[1] = 0.5 * (2.0 + u * u * (3.0 * u - 5.0));
  b[2] = 0.5 * u * (-3.0 * u * u + 4.0 * u + 1.0);
  b[3] = 0.5 * u * u * a;
}

/*!
\brief Compute the Lipschitz constant of the field.
*/
double ScalarField::K() const
{
  double a = 0.0;

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      for (int k = 0; k < nz; k++)
      {
        a = Math::Max(a, Norm(Gradient(i, j, k)));
      }
    }
  }
  return a;
}

/*!
\brief Set all values in the field that are lower than epsilon to true zero.
\param e Epsilon coefficient.
*/
void ScalarField::CutEpsilon(const double& e)
{
  // Process all elements
  for (int i = 0; i < nx * ny * nz; i++)
  {
    if (field.at(i) < e)
    {
      field[i] = 0.0;
    }
  }
}

/*!
\brief Get the integral of the scalar field.

Simply compute the sum of all the elements and make the product by the area of small cell element.

\sa Sum()
*/
double ScalarField::Integral() const
{
  double s = Sum();

  // Scale by the size of a cell element
  s *= CellVolume();

  return s;
}


/*!
\brief Compute the sum of the elements in the scalar field.

\avx

\sa Integral(), Average()
*/
double ScalarField::Sum() const
{
  const int size = nx * ny * nz;

  double s = 0.0;

#ifdef _MSC_VER
  if (System::Avx() && field.size() > 8)
  {
    const double* p = field.data();

    __m256d avxa = _mm256_load_pd(p);
    __m256d avxb = _mm256_load_pd(p + 4);

    const int offset = size % 8;
    for (int i = 1; i < size / 8; i++)
    {
      __m256d avxi = _mm256_load_pd(p + i * 8);
      avxa = _mm256_add_pd(avxa, avxi);
      __m256d avxj = _mm256_load_pd(p + i * 8 + 4);
      avxb = _mm256_add_pd(avxb, avxj);
    }
    avxa = _mm256_add_pd(avxa, avxb);
    double e[4];
    _mm256_store_pd(e, avxa);
    s = e[0] + e[1] + e[2] + e[3];

    // Could try this instead of the previous four lines
    // avxa = _mm256_hadd_pd(avxa, avxa);
    // s= ((double*)&avxa)[0] + ((double*)&avxa)[2];

    for (int i = size - offset; i < size; i++)
    {
      s += field.at(i);
    }
  }
  else
#endif
  {
    // Process all elements

    for (int i = 0; i < size; i++)
    {
      s += field.at(i);
    }
  }
  return s;
}

/*!
\brief Get the range of the field.
\param a,b Returned minimum and maximum.
*/
void ScalarField::GetRange(double& a, double& b) const
{
  a = field.at(0);
  b = a;

  for (int i = 1; i < field.size(); i++)
  {
    double x = field.at(i);
    if (x < a)
    {
      a = x;
    }
    else if (x > b)
    {
      b = x;
    }
  }
}

/*!
\brief Get the field value with world coordinate system.

This function computes a tri-linear interpolation of the values.

\sa Math::Trilinear
\param p Point position (should be strictly within bounding box of the domain).
*/
double ScalarField::Value(const Vector& p) const
{
  double u, v, w;
  int i, j, k;
  CellInteger(p, i, j, k, u, v, w);

  // Test position
  if (!InsideCellIndex(i, j, k))
  {
    return 0.0;
  }

  return Math::Trilinear(at(i, j, k), at(i + 1, j, k), at(i + 1, j + 1, k), at(i, j + 1, k), at(i, j, k + 1), at(i + 1, j, k + 1), at(i + 1, j + 1, k + 1), at(i, j + 1, k + 1), u, v, w);
}

/*!
\brief Sum two scalar fields.

Scalar fields should have the same size.

\param s Argument scalar field.
*/
ScalarField& ScalarField::operator+=(const ScalarField& s)
{
  for (int i = 0; i < field.size(); i++)
  {
    field[i] += s.at(i);
  }
  return *this;
}

/*!
\brief Scale the values of a scalar field.

\param s Scaling factor.
*/
ScalarField& ScalarField::operator*=(const double& s)
{
  for (int i = 0; i < field.size(); i++)
  {
    field[i] *= s;
  }
  return *this;
}

/*!
\brief Translate the domain of the scalar field.

\param t Translation vector.
*/
void ScalarField::Translate(const Vector& t)
{
  Box::Translate(t);
}

/*!
\brief Scale the domain of the scalar field.

\param s Scaling factor.
*/
void ScalarField::Scale(const Vector& s)
{
  Array::Scale(s);
}

/*!
\brief Compute the size of the scalarfield.
*/
int ScalarField::Memory() const
{
  return Array::Memory() + sizeof(double) * int(field.size());
}


int ScalarField::edgeTable[256] = {
  0, 273, 545, 816, 1042, 1283, 1587, 1826, 2082, 2355, 2563, 2834, 3120, 3361, 3601, 3840,
  324, 85, 869, 628, 1366, 1095, 1911, 1638, 2406, 2167, 2887, 2646, 3444, 3173, 3925, 3652,
  644, 917, 165, 436, 1686, 1927, 1207, 1446, 2726, 2999, 2183, 2454, 3764, 4005, 3221, 3460,
  960, 721, 481, 240, 2002, 1731, 1523, 1250, 3042, 2803, 2499, 2258, 4080, 3809, 3537, 3264,
  1096, 1369, 1641, 1912, 90, 331, 635, 874, 3178, 3451, 3659, 3930, 2168, 2409, 2649, 2888,
  1292, 1053, 1837, 1596, 286, 15, 831, 558, 3374, 3135, 3855, 3614, 2364, 2093, 2845, 2572,
  1740, 2013, 1261, 1532, 734, 975, 255, 494, 3822, 4095, 3279, 3550, 2812, 3053, 2269, 2508,
  1928, 1689, 1449, 1208, 922, 651, 443, 170, 4010, 3771, 3467, 3226, 3000, 2729, 2457, 2184,
  2184, 2457, 2729, 3000, 3226, 3467, 3771, 4010, 170, 443, 651, 922, 1208, 1449, 1689, 1928,
  2508, 2269, 3053, 2812, 3550, 3279, 4095, 3822, 494, 255, 975, 734, 1532, 1261, 2013, 1740,
  2572, 2845, 2093, 2364, 3614, 3855, 3135, 3374, 558, 831, 15, 286, 1596, 1837, 1053, 1292,
  2888, 2649, 2409, 2168, 3930, 3659, 3451, 3178, 874, 635, 331, 90, 1912, 1641, 1369, 1096,
  3264, 3537, 3809, 4080, 2258, 2499, 2803, 3042, 1250, 1523, 1731, 2002, 240, 481, 721, 960,
  3460, 3221, 4005, 3764, 2454, 2183, 2999, 2726, 1446, 1207, 1927, 1686, 436, 165, 917, 644,
  3652, 3925, 3173, 3444, 2646, 2887, 2167, 2406, 1638, 1911, 1095, 1366, 628, 869, 85, 324,
  3840, 3601, 3361, 3120, 2834, 2563, 2355, 2082, 1826, 1587, 1283, 1042, 816, 545, 273, 0
};

int ScalarField::TriangleTable[256][16] = {
  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 8, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 5, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 5, 8, 4, 9, 8, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 4, 10, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 10, 1, 8, 10, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 5, 9, 0, 1, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 5, 10, 1, 5, 9, 10, 9, 8, 10, -1, -1, -1, -1, -1, -1, -1 },
  { 5, 1, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 8, 4, 5, 1, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 1, 11, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 8, 4, 1, 11, 8, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1 },
  { 4, 11, 5, 10, 11, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 11, 5, 0, 8, 11, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1 },
  { 4, 9, 0, 4, 10, 9, 10, 11, 9, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 8, 11, 11, 8, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 2, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 2, 4, 0, 6, 4, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 5, 9, 8, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 2, 5, 9, 2, 6, 5, 6, 4, 5, -1, -1, -1, -1, -1, -1, -1 },
  { 8, 2, 6, 4, 10, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 10, 2, 6, 10, 1, 2, 1, 0, 2, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 0, 5, 8, 2, 6, 1, 4, 10, -1, -1, -1, -1, -1, -1, -1 },
  { 2, 6, 10, 9, 2, 10, 9, 10, 1, 9, 1, 5, -1, -1, -1, -1 },
  { 5, 1, 11, 8, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 4, 2, 6, 4, 0, 2, 5, 1, 11, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 1, 11, 9, 0, 1, 8, 2, 6, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 11, 9, 1, 9, 6, 1, 6, 4, 6, 9, 2, -1, -1, -1, -1 },
  { 4, 11, 5, 4, 10, 11, 6, 8, 2, -1, -1, -1, -1, -1, -1, -1 },
  { 5, 10, 11, 5, 2, 10, 5, 0, 2, 6, 10, 2, -1, -1, -1, -1 },
  { 2, 6, 8, 9, 0, 10, 9, 10, 11, 10, 0, 4, -1, -1, -1, -1 },
  { 2, 6, 10, 2, 10, 9, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 7, 2, 0, 8, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 7, 2, 5, 7, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 8, 7, 2, 8, 4, 7, 4, 5, 7, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 7, 2, 1, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 10, 1, 0, 8, 10, 2, 9, 7, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 7, 2, 0, 5, 7, 1, 4, 10, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 5, 7, 1, 7, 8, 1, 8, 10, 2, 8, 7, -1, -1, -1, -1 },
  { 5, 1, 11, 9, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 4, 0, 8, 5, 1, 11, 2, 9, 7, -1, -1, -1, -1, -1, -1, -1 },
  { 7, 1, 11, 7, 2, 1, 2, 0, 1, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 11, 7, 4, 1, 7, 4, 7, 2, 4, 2, 8, -1, -1, -1, -1 },
  { 11, 4, 10, 11, 5, 4, 9, 7, 2, -1, -1, -1, -1, -1, -1, -1 },
  { 2, 9, 7, 0, 8, 5, 8, 11, 5, 8, 10, 11, -1, -1, -1, -1 },
  { 7, 2, 0, 7, 0, 10, 7, 10, 11, 10, 0, 4, -1, -1, -1, -1 },
  { 7, 2, 8, 7, 8, 11, 11, 8, 10, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 6, 8, 7, 6, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 4, 0, 9, 7, 4, 7, 6, 4, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 6, 8, 0, 5, 6, 5, 7, 6, -1, -1, -1, -1, -1, -1, -1 },
  { 5, 7, 4, 4, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 6, 9, 7, 6, 8, 9, 4, 10, 1, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 7, 6, 9, 6, 1, 9, 1, 0, 1, 6, 10, -1, -1, -1, -1 },
  { 1, 4, 10, 0, 5, 8, 5, 6, 8, 5, 7, 6, -1, -1, -1, -1 },
  { 10, 1, 5, 10, 5, 6, 6, 5, 7, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 6, 8, 9, 7, 6, 11, 5, 1, -1, -1, -1, -1, -1, -1, -1 },
  { 11, 5, 1, 9, 7, 0, 7, 4, 0, 7, 6, 4, -1, -1, -1, -1 },
  { 8, 0, 1, 8, 1, 7, 8, 7, 6, 11, 7, 1, -1, -1, -1, -1 },
  { 1, 11, 7, 1, 7, 4, 4, 7, 6, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 7, 8, 8, 7, 6, 11, 5, 4, 11, 4, 10, -1, -1, -1, -1 },
  { 7, 6, 0, 7, 0, 9, 6, 10, 0, 5, 0, 11, 10, 11, 0, -1 },
  { 10, 11, 0, 10, 0, 4, 11, 7, 0, 8, 0, 6, 7, 6, 0, -1 },
  { 10, 11, 7, 6, 10, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 6, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 4, 0, 8, 10, 6, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 5, 9, 10, 6, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 8, 5, 9, 8, 4, 5, 10, 6, 3, -1, -1, -1, -1, -1, -1, -1 },
  { 6, 1, 4, 3, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 6, 0, 8, 6, 3, 0, 3, 1, 0, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 6, 3, 1, 4, 6, 0, 5, 9, -1, -1, -1, -1, -1, -1, -1 },
  { 5, 3, 1, 5, 8, 3, 5, 9, 8, 8, 6, 3, -1, -1, -1, -1 },
  { 11, 5, 1, 3, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 5, 1, 11, 4, 0, 8, 3, 10, 6, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 9, 0, 1, 11, 9, 3, 10, 6, -1, -1, -1, -1, -1, -1, -1 },
  { 3, 10, 6, 1, 11, 4, 11, 8, 4, 11, 9, 8, -1, -1, -1, -1 },
  { 11, 6, 3, 11, 5, 6, 5, 4, 6, -1, -1, -1, -1, -1, -1, -1 },
  { 11, 6, 3, 5, 6, 11, 5, 8, 6, 5, 0, 8, -1, -1, -1, -1 },
  { 0, 4, 6, 0, 6, 11, 0, 11, 9, 3, 11, 6, -1, -1, -1, -1 },
  { 6, 3, 11, 6, 11, 8, 8, 11, 9, -1, -1, -1, -1, -1, -1, -1 },
  { 3, 8, 2, 10, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 4, 3, 10, 4, 0, 3, 0, 2, 3, -1, -1, -1, -1, -1, -1, -1 },
  { 8, 3, 10, 8, 2, 3, 9, 0, 5, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 2, 3, 9, 3, 4, 9, 4, 5, 10, 4, 3, -1, -1, -1, -1 },
  { 8, 1, 4, 8, 2, 1, 2, 3, 1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 2, 1, 2, 3, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 5, 9, 0, 1, 4, 2, 1, 2, 3, 2, 4, 8, -1, -1, -1, -1 },
  { 5, 9, 2, 5, 2, 1, 1, 2, 3, -1, -1, -1, -1, -1, -1, -1 },
  { 3, 8, 2, 3, 10, 8, 1, 11, 5, -1, -1, -1, -1, -1, -1, -1 },
  { 5, 1, 11, 4, 0, 10, 0, 3, 10, 0, 2, 3, -1, -1, -1, -1 },
  { 2, 10, 8, 2, 3, 10, 0, 1, 9, 1, 11, 9, -1, -1, -1, -1 },
  { 11, 9, 4, 11, 4, 1, 9, 2, 4, 10, 4, 3, 2, 3, 4, -1 },
  { 8, 5, 4, 8, 3, 5, 8, 2, 3, 3, 11, 5, -1, -1, -1, -1 },
  { 11, 5, 0, 11, 0, 3, 3, 0, 2, -1, -1, -1, -1, -1, -1, -1 },
  { 2, 3, 4, 2, 4, 8, 3, 11, 4, 0, 4, 9, 11, 9, 4, -1 },
  { 11, 9, 2, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 2, 9, 7, 6, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 8, 4, 2, 9, 7, 10, 6, 3, -1, -1, -1, -1, -1, -1, -1 },
  { 7, 0, 5, 7, 2, 0, 6, 3, 10, -1, -1, -1, -1, -1, -1, -1 },
  { 10, 6, 3, 8, 4, 2, 4, 7, 2, 4, 5, 7, -1, -1, -1, -1 },
  { 6, 1, 4, 6, 3, 1, 7, 2, 9, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 7, 2, 0, 8, 3, 0, 3, 1, 3, 8, 6, -1, -1, -1, -1 },
  { 4, 3, 1, 4, 6, 3, 5, 7, 0, 7, 2, 0, -1, -1, -1, -1 },
  { 3, 1, 8, 3, 8, 6, 1, 5, 8, 2, 8, 7, 5, 7, 8, -1 },
  { 9, 7, 2, 11, 5, 1, 6, 3, 10, -1, -1, -1, -1, -1, -1, -1 },
  { 3, 10, 6, 5, 1, 11, 0, 8, 4, 2, 9, 7, -1, -1, -1, -1 },
  { 6, 3, 10, 7, 2, 11, 2, 1, 11, 2, 0, 1, -1, -1, -1, -1 },
  { 4, 2, 8, 4, 7, 2, 4, 1, 7, 11, 7, 1, 10, 6, 3, -1 },
  { 9, 7, 2, 11, 5, 3, 5, 6, 3, 5, 4, 6, -1, -1, -1, -1 },
  { 5, 3, 11, 5, 6, 3, 5, 0, 6, 8, 6, 0, 9, 7, 2, -1 },
  { 2, 0, 11, 2, 11, 7, 0, 4, 11, 3, 11, 6, 4, 6, 11, -1 },
  { 6, 3, 11, 6, 11, 8, 7, 2, 11, 2, 8, 11, -1, -1, -1, -1 },
  { 3, 9, 7, 3, 10, 9, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1 },
  { 4, 3, 10, 0, 3, 4, 0, 7, 3, 0, 9, 7, -1, -1, -1, -1 },
  { 0, 10, 8, 0, 7, 10, 0, 5, 7, 7, 3, 10, -1, -1, -1, -1 },
  { 3, 10, 4, 3, 4, 7, 7, 4, 5, -1, -1, -1, -1, -1, -1, -1 },
  { 7, 8, 9, 7, 1, 8, 7, 3, 1, 4, 8, 1, -1, -1, -1, -1 },
  { 9, 7, 3, 9, 3, 0, 0, 3, 1, -1, -1, -1, -1, -1, -1, -1 },
  { 5, 7, 8, 5, 8, 0, 7, 3, 8, 4, 8, 1, 3, 1, 8, -1 },
  { 5, 7, 3, 1, 5, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 5, 1, 11, 9, 7, 10, 9, 10, 8, 10, 7, 3, -1, -1, -1, -1 },
  { 0, 10, 4, 0, 3, 10, 0, 9, 3, 7, 3, 9, 5, 1, 11, -1 },
  { 10, 8, 7, 10, 7, 3, 8, 0, 7, 11, 7, 1, 0, 1, 7, -1 },
  { 3, 10, 4, 3, 4, 7, 1, 11, 4, 11, 7, 4, -1, -1, -1, -1 },
  { 5, 4, 3, 5, 3, 11, 4, 8, 3, 7, 3, 9, 8, 9, 3, -1 },
  { 11, 5, 0, 11, 0, 3, 9, 7, 0, 7, 3, 0, -1, -1, -1, -1 },
  { 0, 4, 8, 7, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 11, 7, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 11, 3, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 8, 4, 7, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 0, 5, 7, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 5, 8, 4, 5, 9, 8, 7, 11, 3, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 4, 10, 11, 3, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 10, 0, 8, 10, 1, 0, 11, 3, 7, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 5, 9, 1, 4, 10, 7, 11, 3, -1, -1, -1, -1, -1, -1, -1 },
  { 7, 11, 3, 5, 9, 1, 9, 10, 1, 9, 8, 10, -1, -1, -1, -1 },
  { 5, 3, 7, 1, 3, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 5, 3, 7, 5, 1, 3, 4, 0, 8, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 3, 7, 9, 0, 3, 0, 1, 3, -1, -1, -1, -1, -1, -1, -1 },
  { 7, 9, 8, 7, 8, 1, 7, 1, 3, 4, 1, 8, -1, -1, -1, -1 },
  { 3, 4, 10, 3, 7, 4, 7, 5, 4, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 8, 10, 0, 10, 7, 0, 7, 5, 7, 10, 3, -1, -1, -1, -1 },
  { 4, 10, 3, 0, 4, 3, 0, 3, 7, 0, 7, 9, -1, -1, -1, -1 },
  { 3, 7, 9, 3, 9, 10, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1 },
  { 7, 11, 3, 2, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 2, 4, 0, 2, 6, 4, 3, 7, 11, -1, -1, -1, -1, -1, -1, -1 },
  { 5, 9, 0, 7, 11, 3, 8, 2, 6, -1, -1, -1, -1, -1, -1, -1 },
  { 11, 3, 7, 5, 9, 6, 5, 6, 4, 6, 9, 2, -1, -1, -1, -1 },
  { 4, 10, 1, 6, 8, 2, 11, 3, 7, -1, -1, -1, -1, -1, -1, -1 },
  { 7, 11, 3, 2, 6, 1, 2, 1, 0, 1, 6, 10, -1, -1, -1, -1 },
  { 0, 5, 9, 2, 6, 8, 1, 4, 10, 7, 11, 3, -1, -1, -1, -1 },
  { 9, 1, 5, 9, 10, 1, 9, 2, 10, 6, 10, 2, 7, 11, 3, -1 },
  { 3, 5, 1, 3, 7, 5, 2, 6, 8, -1, -1, -1, -1, -1, -1, -1 },
  { 5, 1, 7, 7, 1, 3, 4, 0, 2, 4, 2, 6, -1, -1, -1, -1 },
  { 8, 2, 6, 9, 0, 7, 0, 3, 7, 0, 1, 3, -1, -1, -1, -1 },
  { 6, 4, 9, 6, 9, 2, 4, 1, 9, 7, 9, 3, 1, 3, 9, -1 },
  { 8, 2, 6, 4, 10, 7, 4, 7, 5, 7, 10, 3, -1, -1, -1, -1 },
  { 7, 5, 10, 7, 10, 3, 5, 0, 10, 6, 10, 2, 0, 2, 10, -1 },
  { 0, 7, 9, 0, 3, 7, 0, 4, 3, 10, 3, 4, 8, 2, 6, -1 },
  { 3, 7, 9, 3, 9, 10, 2, 6, 9, 6, 10, 9, -1, -1, -1, -1 },
  { 11, 2, 9, 3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 2, 11, 3, 2, 9, 11, 0, 8, 4, -1, -1, -1, -1, -1, -1, -1 },
  { 11, 0, 5, 11, 3, 0, 3, 2, 0, -1, -1, -1, -1, -1, -1, -1 },
  { 8, 4, 5, 8, 5, 3, 8, 3, 2, 3, 5, 11, -1, -1, -1, -1 },
  { 11, 2, 9, 11, 3, 2, 10, 1, 4, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 8, 1, 1, 8, 10, 2, 9, 11, 2, 11, 3, -1, -1, -1, -1 },
  { 4, 10, 1, 0, 5, 3, 0, 3, 2, 3, 5, 11, -1, -1, -1, -1 },
  { 3, 2, 5, 3, 5, 11, 2, 8, 5, 1, 5, 10, 8, 10, 5, -1 },
  { 5, 2, 9, 5, 1, 2, 1, 3, 2, -1, -1, -1, -1, -1, -1, -1 },
  { 4, 0, 8, 5, 1, 9, 1, 2, 9, 1, 3, 2, -1, -1, -1, -1 },
  { 0, 1, 2, 2, 1, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 8, 4, 1, 8, 1, 2, 2, 1, 3, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 3, 2, 9, 4, 3, 9, 5, 4, 10, 3, 4, -1, -1, -1, -1 },
  { 8, 10, 5, 8, 5, 0, 10, 3, 5, 9, 5, 2, 3, 2, 5, -1 },
  { 4, 10, 3, 4, 3, 0, 0, 3, 2, -1, -1, -1, -1, -1, -1, -1 },
  { 3, 2, 8, 10, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 6, 11, 3, 6, 8, 11, 8, 9, 11, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 6, 4, 0, 11, 6, 0, 9, 11, 3, 6, 11, -1, -1, -1, -1 },
  { 11, 3, 6, 5, 11, 6, 5, 6, 8, 5, 8, 0, -1, -1, -1, -1 },
  { 11, 3, 6, 11, 6, 5, 5, 6, 4, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 4, 10, 11, 3, 8, 11, 8, 9, 8, 3, 6, -1, -1, -1, -1 },
  { 1, 0, 6, 1, 6, 10, 0, 9, 6, 3, 6, 11, 9, 11, 6, -1 },
  { 5, 8, 0, 5, 6, 8, 5, 11, 6, 3, 6, 11, 1, 4, 10, -1 },
  { 10, 1, 5, 10, 5, 6, 11, 3, 5, 3, 6, 5, -1, -1, -1, -1 },
  { 5, 1, 3, 5, 3, 8, 5, 8, 9, 8, 3, 6, -1, -1, -1, -1 },
  { 1, 3, 9, 1, 9, 5, 3, 6, 9, 0, 9, 4, 6, 4, 9, -1 },
  { 6, 8, 0, 6, 0, 3, 3, 0, 1, -1, -1, -1, -1, -1, -1, -1 },
  { 6, 4, 1, 3, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 8, 9, 3, 8, 3, 6, 9, 5, 3, 10, 3, 4, 5, 4, 3, -1 },
  { 0, 9, 5, 10, 3, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 6, 8, 0, 6, 0, 3, 4, 10, 0, 10, 3, 0, -1, -1, -1, -1 },
  { 6, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 10, 7, 11, 6, 7, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 10, 7, 11, 10, 6, 7, 8, 4, 0, -1, -1, -1, -1, -1, -1, -1 },
  { 7, 10, 6, 7, 11, 10, 5, 9, 0, -1, -1, -1, -1, -1, -1, -1 },
  { 11, 6, 7, 11, 10, 6, 9, 8, 5, 8, 4, 5, -1, -1, -1, -1 },
  { 1, 7, 11, 1, 4, 7, 4, 6, 7, -1, -1, -1, -1, -1, -1, -1 },
  { 8, 1, 0, 8, 7, 1, 8, 6, 7, 11, 1, 7, -1, -1, -1, -1 },
  { 9, 0, 5, 7, 11, 4, 7, 4, 6, 4, 11, 1, -1, -1, -1, -1 },
  { 9, 8, 1, 9, 1, 5, 8, 6, 1, 11, 1, 7, 6, 7, 1, -1 },
  { 10, 5, 1, 10, 6, 5, 6, 7, 5, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 8, 4, 5, 1, 6, 5, 6, 7, 6, 1, 10, -1, -1, -1, -1 },
  { 9, 6, 7, 9, 1, 6, 9, 0, 1, 1, 10, 6, -1, -1, -1, -1 },
  { 6, 7, 1, 6, 1, 10, 7, 9, 1, 4, 1, 8, 9, 8, 1, -1 },
  { 5, 4, 7, 4, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 8, 6, 0, 6, 5, 5, 6, 7, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 0, 4, 9, 4, 7, 7, 4, 6, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 8, 6, 7, 9, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 7, 8, 2, 7, 11, 8, 11, 10, 8, -1, -1, -1, -1, -1, -1, -1 },
  { 7, 0, 2, 7, 10, 0, 7, 11, 10, 10, 4, 0, -1, -1, -1, -1 },
  { 0, 5, 9, 8, 2, 11, 8, 11, 10, 11, 2, 7, -1, -1, -1, -1 },
  { 11, 10, 2, 11, 2, 7, 10, 4, 2, 9, 2, 5, 4, 5, 2, -1 },
  { 1, 7, 11, 4, 7, 1, 4, 2, 7, 4, 8, 2, -1, -1, -1, -1 },
  { 7, 11, 1, 7, 1, 2, 2, 1, 0, -1, -1, -1, -1, -1, -1, -1 },
  { 4, 11, 1, 4, 7, 11, 4, 8, 7, 2, 7, 8, 0, 5, 9, -1 },
  { 7, 11, 1, 7, 1, 2, 5, 9, 1, 9, 2, 1, -1, -1, -1, -1 },
  { 1, 7, 5, 1, 8, 7, 1, 10, 8, 2, 7, 8, -1, -1, -1, -1 },
  { 0, 2, 10, 0, 10, 4, 2, 7, 10, 1, 10, 5, 7, 5, 10, -1 },
  { 0, 1, 7, 0, 7, 9, 1, 10, 7, 2, 7, 8, 10, 8, 7, -1 },
  { 9, 2, 7, 1, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 8, 2, 7, 8, 7, 4, 4, 7, 5, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 2, 7, 5, 0, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 8, 2, 7, 8, 7, 4, 9, 0, 7, 0, 4, 7, -1, -1, -1, -1 },
  { 9, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 2, 10, 6, 2, 9, 10, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 8, 4, 2, 9, 6, 9, 10, 6, 9, 11, 10, -1, -1, -1, -1 },
  { 5, 11, 10, 5, 10, 2, 5, 2, 0, 6, 2, 10, -1, -1, -1, -1 },
  { 4, 5, 2, 4, 2, 8, 5, 11, 2, 6, 2, 10, 11, 10, 2, -1 },
  { 1, 9, 11, 1, 6, 9, 1, 4, 6, 6, 2, 9, -1, -1, -1, -1 },
  { 9, 11, 6, 9, 6, 2, 11, 1, 6, 8, 6, 0, 1, 0, 6, -1 },
  { 4, 6, 11, 4, 11, 1, 6, 2, 11, 5, 11, 0, 2, 0, 11, -1 },
  { 5, 11, 1, 8, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 2, 10, 6, 9, 10, 2, 9, 1, 10, 9, 5, 1, -1, -1, -1, -1 },
  { 9, 6, 2, 9, 10, 6, 9, 5, 10, 1, 10, 5, 0, 8, 4, -1 },
  { 10, 6, 2, 10, 2, 1, 1, 2, 0, -1, -1, -1, -1, -1, -1, -1 },
  { 10, 6, 2, 10, 2, 1, 8, 4, 2, 4, 1, 2, -1, -1, -1, -1 },
  { 2, 9, 5, 2, 5, 6, 6, 5, 4, -1, -1, -1, -1, -1, -1, -1 },
  { 2, 9, 5, 2, 5, 6, 0, 8, 5, 8, 6, 5, -1, -1, -1, -1 },
  { 2, 0, 4, 6, 2, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 2, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 11, 8, 11, 10, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 4, 0, 9, 4, 9, 10, 10, 9, 11, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 5, 11, 0, 11, 8, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1 },
  { 4, 5, 11, 10, 4, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 4, 8, 1, 8, 11, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 11, 1, 0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 4, 8, 1, 8, 11, 0, 5, 8, 5, 11, 8, -1, -1, -1, -1 },
  { 5, 11, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 5, 1, 10, 5, 10, 9, 9, 10, 8, -1, -1, -1, -1, -1, -1, -1 },
  { 4, 0, 9, 4, 9, 10, 5, 1, 9, 1, 10, 9, -1, -1, -1, -1 },
  { 0, 1, 10, 8, 0, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 4, 1, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 5, 4, 8, 9, 5, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 9, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 4, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 } };

/*!
\brief Compute the configuration for the voxel corners of the cell.
\param cx, cy, cz Cell coordinates.
*/
int ScalarField::DualCellCode(int cx, int cy, int cz) const
{
  // Determine for each cube corner if it is outside or inside
  int code = 0;
  if (at(cx, cy, cz) >= 0.0)
    code |= 1;
  if (at(cx + 1, cy, cz) >= 0.0)
    code |= 2;
  if (at(cx, cy + 1, cz) >= 0.0)
    code |= 4;
  if (at(cx + 1, cy + 1, cz) >= 0.0)
    code |= 8;
  if (at(cx, cy, cz + 1) >= 0.0)
    code |= 16;
  if (at(cx + 1, cy, cz + 1) >= 0.0)
    code |= 32;
  if (at(cx, cy + 1, cz + 1) >= 0.0)
    code |= 64;
  if (at(cx + 1, cy + 1, cz + 1) >= 0.0)
    code |= 128;
  return code;
}

/*!
\brief Compute the gradient.
\param p Point.
*/
Vector ScalarField::Gradient(const Vector& p) const
{
  static const double Epsilon = 1e-4;

  double x = Value(Vector(p[0] + Epsilon, p[1], p[2])) - Value(Vector(p[0] - Epsilon, p[1], p[2]));
  double y = Value(Vector(p[0], p[1] + Epsilon, p[2])) - Value(Vector(p[0], p[1] - Epsilon, p[2]));
  double z = Value(Vector(p[0], p[1], p[2] + Epsilon)) - Value(Vector(p[0], p[1], p[2] - Epsilon));
  return Vector(x, y, z) * (0.5 / Epsilon);
}

/*!
\brief Compute the Hessian symmetric matrix of the field function.

\param p Point.
\return Hessian matrix.
*/
Matrix ScalarField::Hessian(const Vector& p) const
{
  double a[3][3][3];
  Matrix H;

  for (int i = -1; i < 2; i++)
  {
    for (int j = -1; j < 2; j++)
    {
      for (int k = -1; k < 2; k++)
      {
        a[1 + i][1 + j][1 + k] = Value(p + Vector(i, j, k) * Epsilon);
      }
    }
  }
  H[0] = (a[2][1][1] - 2.0 * a[1][1][1] + a[0][1][1]) / (Epsilon * Epsilon);
  H[1] = (a[2][2][1] - a[2][0][1] - a[0][2][1] + a[0][0][1]) / (2.0 * Epsilon * Epsilon);
  H[2] = (a[2][1][2] - a[2][1][0] - a[0][1][2] + a[0][1][0]) / (2.0 * Epsilon * Epsilon);

  H[3] = H[1];
  H[4] = (a[1][2][1] - 2.0 * a[1][1][1] + a[1][0][1]) / (Epsilon * Epsilon);
  H[5] = (a[1][2][2] - a[1][2][0] - a[1][0][2] + a[1][0][0]) / (2.0 * Epsilon * Epsilon);

  H[6] = H[2];
  H[7] = H[5];
  H[8] = (a[1][1][2] - 2.0 * a[1][1][1] + a[1][1][0]) / (Epsilon * Epsilon);

  return H;
}

/*!
\brief Compute the Hessian symmetric matrix of the field.

\param i,j,k Integer coordinates of the point.
\return Hessian matrix.
*/
Matrix ScalarField::Hessian(int i, int j, int k) const
{
  double a[3][3][3];
  Matrix H;

  for (int ai = -1; ai < 2; ai++)
  {
    for (int aj = -1; aj < 2; aj++)
    {
      for (int ak = -1; ak < 2; ak++)
      {
        a[1 + i][1 + j][1 + k] = Value(i + ai, j + aj, k + ak);
      }
    }
  }
  H[0] = (a[2][1][1] - 2.0 * a[1][1][1] + a[0][1][1]) / (Epsilon * Epsilon);
  H[1] = (a[2][2][1] - a[2][0][1] - a[0][2][1] + a[0][0][1]) / (2.0 * Epsilon * Epsilon);
  H[2] = (a[2][1][2] - a[2][1][0] - a[0][1][2] + a[0][1][0]) / (2.0 * Epsilon * Epsilon);

  H[3] = H[1];
  H[4] = (a[1][2][1] - 2.0 * a[1][1][1] + a[1][0][1]) / (Epsilon * Epsilon);
  H[5] = (a[1][2][2] - a[1][2][0] - a[1][0][2] + a[1][0][0]) / (2.0 * Epsilon * Epsilon);

  H[6] = H[2];
  H[7] = H[5];
  H[8] = (a[1][1][2] - 2.0 * a[1][1][1] + a[1][1][0]) / (Epsilon * Epsilon);

  return H;
}

static const int EDGE0 = 1;
static const int EDGE1 = 1 << 1;
static const int EDGE2 = 1 << 2;
static const int EDGE3 = 1 << 3;
static const int EDGE4 = 1 << 4;
static const int EDGE5 = 1 << 5;
static const int EDGE6 = 1 << 6;
static const int EDGE7 = 1 << 7;
static const int EDGE8 = 1 << 8;
static const int EDGE9 = 1 << 9;
static const int EDGE10 = 1 << 10;
static const int EDGE11 = 1 << 11;

#include <unordered_map>

static std::unordered_map<int64_t, int> pointToIndex; //!< Hash map for shared vertex index computations

static int64_t Codec(int id, int p)
{
  return int64_t(id) | (int64_t(p) << 32u);
}

int const dualvertex[256][4] = { //!< Dual Marching Cubes table encoding the edge vertices for the 256 marching cubes cases.
{0, 0, 0, 0}, // 0
{EDGE0 | EDGE3 | EDGE8, 0, 0, 0}, // 1
{EDGE0 | EDGE1 | EDGE9, 0, 0, 0}, // 2
{EDGE1 | EDGE3 | EDGE8 | EDGE9, 0, 0, 0}, // 3
{EDGE4 | EDGE7 | EDGE8, 0, 0, 0}, // 4
{EDGE0 | EDGE3 | EDGE4 | EDGE7, 0, 0, 0}, // 5
{EDGE0 | EDGE1 | EDGE9, EDGE4 | EDGE7 | EDGE8, 0, 0}, // 6
{EDGE1 | EDGE3 | EDGE4 | EDGE7 | EDGE9, 0, 0, 0}, // 7
{EDGE4 | EDGE5 | EDGE9, 0, 0, 0}, // 8
{EDGE0 | EDGE3 | EDGE8, EDGE4 | EDGE5 | EDGE9, 0, 0}, // 9
{EDGE0 | EDGE1 | EDGE4 | EDGE5, 0, 0, 0}, // 10
{EDGE1 | EDGE3 | EDGE4 | EDGE5 | EDGE8, 0, 0, 0}, // 11
{EDGE5 | EDGE7 | EDGE8 | EDGE9, 0, 0, 0}, // 12
{EDGE0 | EDGE3 | EDGE5 | EDGE7 | EDGE9, 0, 0, 0}, // 13
{EDGE0 | EDGE1 | EDGE5 | EDGE7 | EDGE8, 0, 0, 0}, // 14
{EDGE1 | EDGE3 | EDGE5 | EDGE7, 0, 0, 0}, // 15
{EDGE2 | EDGE3 | EDGE11, 0, 0, 0}, // 16
{EDGE0 | EDGE2 | EDGE8 | EDGE11, 0, 0, 0}, // 17
{EDGE0 | EDGE1 | EDGE9, EDGE2 | EDGE3 | EDGE11, 0, 0}, // 18
{EDGE1 | EDGE2 | EDGE8 | EDGE9 | EDGE11, 0, 0, 0}, // 19
{EDGE4 | EDGE7 | EDGE8, EDGE2 | EDGE3 | EDGE11, 0, 0}, // 20
{EDGE0 | EDGE2 | EDGE4 | EDGE7 | EDGE11, 0, 0, 0}, // 21
{EDGE0 | EDGE1 | EDGE9, EDGE4 | EDGE7 | EDGE8, EDGE2 | EDGE3 | EDGE11, 0}, // 22
{EDGE1 | EDGE2 | EDGE4 | EDGE7 | EDGE9 | EDGE11, 0, 0, 0}, // 23
{EDGE4 | EDGE5 | EDGE9, EDGE2 | EDGE3 | EDGE11, 0, 0}, // 24
{EDGE0 | EDGE2 | EDGE8 | EDGE11, EDGE4 | EDGE5 | EDGE9, 0, 0}, // 25
{EDGE0 | EDGE1 | EDGE4 | EDGE5, EDGE2 | EDGE3 | EDGE11, 0, 0}, // 26
{EDGE1 | EDGE2 | EDGE4 | EDGE5 | EDGE8 | EDGE11, 0, 0, 0}, // 27
{EDGE5 | EDGE7 | EDGE8 | EDGE9, EDGE2 | EDGE3 | EDGE11, 0, 0}, // 28
{EDGE0 | EDGE2 | EDGE5 | EDGE7 | EDGE9 | EDGE11, 0, 0, 0}, // 29
{EDGE0 | EDGE1 | EDGE5 | EDGE7 | EDGE8, EDGE2 | EDGE3 | EDGE11, 0, 0}, // 30
{EDGE1 | EDGE2 | EDGE5 | EDGE7 | EDGE11, 0, 0, 0}, // 31
{EDGE1 | EDGE2 | EDGE10, 0, 0, 0}, // 32
{EDGE0 | EDGE3 | EDGE8, EDGE1 | EDGE2 | EDGE10, 0, 0}, // 33
{EDGE0 | EDGE2 | EDGE9 | EDGE10, 0, 0, 0}, // 34
{EDGE2 | EDGE3 | EDGE8 | EDGE9 | EDGE10, 0, 0, 0}, // 35
{EDGE4 | EDGE7 | EDGE8, EDGE1 | EDGE2 | EDGE10, 0, 0}, // 36
{EDGE0 | EDGE3 | EDGE4 | EDGE7, EDGE1 | EDGE2 | EDGE10, 0, 0}, // 37
{EDGE0 | EDGE2 | EDGE9 | EDGE10, EDGE4 | EDGE7 | EDGE8, 0, 0}, // 38
{EDGE2 | EDGE3 | EDGE4 | EDGE7 | EDGE9 | EDGE10, 0, 0, 0}, // 39
{EDGE4 | EDGE5 | EDGE9, EDGE1 | EDGE2 | EDGE10, 0, 0}, // 40
{EDGE0 | EDGE3 | EDGE8, EDGE4 | EDGE5 | EDGE9, EDGE1 | EDGE2 | EDGE10, 0}, // 41
{EDGE0 | EDGE2 | EDGE4 | EDGE5 | EDGE10, 0, 0, 0}, // 42
{EDGE2 | EDGE3 | EDGE4 | EDGE5 | EDGE8 | EDGE10, 0, 0, 0}, // 43
{EDGE5 | EDGE7 | EDGE8 | EDGE9, EDGE1 | EDGE2 | EDGE10, 0, 0}, // 44
{EDGE0 | EDGE3 | EDGE5 | EDGE7 | EDGE9, EDGE1 | EDGE2 | EDGE10, 0, 0}, // 45
{EDGE0 | EDGE2 | EDGE5 | EDGE7 | EDGE8 | EDGE10, 0, 0, 0}, // 46
{EDGE2 | EDGE3 | EDGE5 | EDGE7 | EDGE10, 0, 0, 0}, // 47
{EDGE1 | EDGE3 | EDGE10 | EDGE11, 0, 0, 0}, // 48
{EDGE0 | EDGE1 | EDGE8 | EDGE10 | EDGE11, 0, 0, 0}, // 49
{EDGE0 | EDGE3 | EDGE9 | EDGE10 | EDGE11, 0, 0, 0}, // 50
{EDGE8 | EDGE9 | EDGE10 | EDGE11, 0, 0, 0}, // 51
{EDGE4 | EDGE7 | EDGE8, EDGE1 | EDGE3 | EDGE10 | EDGE11, 0, 0}, // 52
{EDGE0 | EDGE1 | EDGE4 | EDGE7 | EDGE10 | EDGE11, 0, 0, 0}, // 53
{EDGE0 | EDGE3 | EDGE9 | EDGE10 | EDGE11, EDGE4 | EDGE7 | EDGE8, 0, 0}, // 54
{EDGE4 | EDGE7 | EDGE9 | EDGE10 | EDGE11, 0, 0, 0}, // 55
{EDGE4 | EDGE5 | EDGE9, EDGE1 | EDGE3 | EDGE10 | EDGE11, 0, 0}, // 56
{EDGE0 | EDGE1 | EDGE8 | EDGE10 | EDGE11, EDGE4 | EDGE5 | EDGE9, 0, 0}, // 57
{EDGE0 | EDGE3 | EDGE4 | EDGE5 | EDGE10 | EDGE11, 0, 0, 0}, // 58
{EDGE4 | EDGE5 | EDGE8 | EDGE10 | EDGE11, 0, 0, 0}, // 59
{EDGE5 | EDGE7 | EDGE8 | EDGE9, EDGE1 | EDGE3 | EDGE10 | EDGE11, 0, 0}, // 60
{EDGE0 | EDGE1 | EDGE5 | EDGE7 | EDGE9 | EDGE10 | EDGE11, 0, 0, 0}, // 61
{EDGE0 | EDGE3 | EDGE5 | EDGE7 | EDGE8 | EDGE10 | EDGE11, 0, 0, 0}, // 62
{EDGE5 | EDGE7 | EDGE10 | EDGE11, 0, 0, 0}, // 63
{EDGE6 | EDGE7 | EDGE11, 0, 0, 0}, // 64
{EDGE0 | EDGE3 | EDGE8, EDGE6 | EDGE7 | EDGE11, 0, 0}, // 65
{EDGE0 | EDGE1 | EDGE9, EDGE6 | EDGE7 | EDGE11, 0, 0}, // 66
{EDGE1 | EDGE3 | EDGE8 | EDGE9, EDGE6 | EDGE7 | EDGE11, 0, 0}, // 67
{EDGE4 | EDGE6 | EDGE8 | EDGE11, 0, 0, 0}, // 68
{EDGE0 | EDGE3 | EDGE4 | EDGE6 | EDGE11, 0, 0, 0}, // 69
{EDGE0 | EDGE1 | EDGE9, EDGE4 | EDGE6 | EDGE8 | EDGE11, 0, 0}, // 70
{EDGE1 | EDGE3 | EDGE4 | EDGE6 | EDGE9 | EDGE11, 0, 0, 0}, // 71
{EDGE4 | EDGE5 | EDGE9, EDGE6 | EDGE7 | EDGE11, 0, 0}, // 72
{EDGE0 | EDGE3 | EDGE8, EDGE4 | EDGE5 | EDGE9, EDGE6 | EDGE7 | EDGE11, 0}, // 73
{EDGE0 | EDGE1 | EDGE4 | EDGE5, EDGE6 | EDGE7 | EDGE11, 0, 0}, // 74
{EDGE1 | EDGE3 | EDGE4 | EDGE5 | EDGE8, EDGE6 | EDGE7 | EDGE11, 0, 0}, // 75
{EDGE5 | EDGE6 | EDGE8 | EDGE9 | EDGE11, 0, 0, 0}, // 76
{EDGE0 | EDGE3 | EDGE5 | EDGE6 | EDGE9 | EDGE11, 0, 0, 0}, // 77
{EDGE0 | EDGE1 | EDGE5 | EDGE6 | EDGE8 | EDGE11, 0, 0, 0}, // 78
{EDGE1 | EDGE3 | EDGE5 | EDGE6 | EDGE11, 0, 0, 0}, // 79
{EDGE2 | EDGE3 | EDGE6 | EDGE7, 0, 0, 0}, // 80
{EDGE0 | EDGE2 | EDGE6 | EDGE7 | EDGE8, 0, 0, 0}, // 81
{EDGE0 | EDGE1 | EDGE9, EDGE2 | EDGE3 | EDGE6 | EDGE7, 0, 0}, // 82
{EDGE1 | EDGE2 | EDGE6 | EDGE7 | EDGE8 | EDGE9, 0, 0, 0}, // 83
{EDGE2 | EDGE3 | EDGE4 | EDGE6 | EDGE8, 0, 0, 0}, // 84
{EDGE0 | EDGE2 | EDGE4 | EDGE6, 0, 0, 0}, // 85
{EDGE0 | EDGE1 | EDGE9, EDGE2 | EDGE3 | EDGE4 | EDGE6 | EDGE8, 0, 0}, // 86
{EDGE1 | EDGE2 | EDGE4 | EDGE6 | EDGE9, 0, 0, 0}, // 87
{EDGE4 | EDGE5 | EDGE9, EDGE2 | EDGE3 | EDGE6 | EDGE7, 0, 0}, // 88
{EDGE0 | EDGE2 | EDGE6 | EDGE7 | EDGE8, EDGE4 | EDGE5 | EDGE9, 0, 0}, // 89
{EDGE0 | EDGE1 | EDGE4 | EDGE5, EDGE2 | EDGE3 | EDGE6 | EDGE7, 0, 0}, // 90
{EDGE1 | EDGE2 | EDGE4 | EDGE5 | EDGE6 | EDGE7 | EDGE8, 0, 0, 0}, // 91
{EDGE2 | EDGE3 | EDGE5 | EDGE6 | EDGE8 | EDGE9, 0, 0, 0}, // 92
{EDGE0 | EDGE2 | EDGE5 | EDGE6 | EDGE9, 0, 0, 0}, // 93
{EDGE0 | EDGE1 | EDGE2 | EDGE3 | EDGE5 | EDGE6 | EDGE8, 0, 0, 0}, // 94
{EDGE1 | EDGE2 | EDGE5 | EDGE6, 0, 0, 0}, // 95
{EDGE1 | EDGE2 | EDGE10, EDGE6 | EDGE7 | EDGE11, 0, 0}, // 96
{EDGE0 | EDGE3 | EDGE8, EDGE1 | EDGE2 | EDGE10, EDGE6 | EDGE7 | EDGE11, 0}, // 97
{EDGE0 | EDGE2 | EDGE9 | EDGE10, EDGE6 | EDGE7 | EDGE11, 0, 0}, // 98
{EDGE2 | EDGE3 | EDGE8 | EDGE9 | EDGE10, EDGE6 | EDGE7 | EDGE11, 0, 0}, // 99
{EDGE4 | EDGE6 | EDGE8 | EDGE11, EDGE1 | EDGE2 | EDGE10, 0, 0}, // 100
{EDGE0 | EDGE3 | EDGE4 | EDGE6 | EDGE11, EDGE1 | EDGE2 | EDGE10, 0, 0}, // 101
{EDGE0 | EDGE2 | EDGE9 | EDGE10, EDGE4 | EDGE6 | EDGE8 | EDGE11, 0, 0}, // 102
{EDGE2 | EDGE3 | EDGE4 | EDGE6 | EDGE9 | EDGE10 | EDGE11, 0, 0, 0}, // 103
{EDGE4 | EDGE5 | EDGE9, EDGE1 | EDGE2 | EDGE10, EDGE6 | EDGE7 | EDGE11, 0}, // 104
{EDGE0 | EDGE3 | EDGE8, EDGE4 | EDGE5 | EDGE9, EDGE1 | EDGE2 | EDGE10, EDGE6 | EDGE7 | EDGE11}, // 105
{EDGE0 | EDGE2 | EDGE4 | EDGE5 | EDGE10, EDGE6 | EDGE7 | EDGE11, 0, 0}, // 106
{EDGE2 | EDGE3 | EDGE4 | EDGE5 | EDGE8 | EDGE10, EDGE6 | EDGE7 | EDGE11, 0, 0}, // 107
{EDGE5 | EDGE6 | EDGE8 | EDGE9 | EDGE11, EDGE1 | EDGE2 | EDGE10, 0, 0}, // 108
{EDGE0 | EDGE3 | EDGE5 | EDGE6 | EDGE9 | EDGE11, EDGE1 | EDGE2 | EDGE10, 0, 0}, // 109
{EDGE0 | EDGE2 | EDGE5 | EDGE6 | EDGE8 | EDGE10 | EDGE11, 0, 0, 0}, // 110
{EDGE2 | EDGE3 | EDGE5 | EDGE6 | EDGE10 | EDGE11, 0, 0, 0}, // 111
{EDGE1 | EDGE3 | EDGE6 | EDGE7 | EDGE10, 0, 0, 0}, // 112
{EDGE0 | EDGE1 | EDGE6 | EDGE7 | EDGE8 | EDGE10, 0, 0, 0}, // 113
{EDGE0 | EDGE3 | EDGE6 | EDGE7 | EDGE9 | EDGE10, 0, 0, 0}, // 114
{EDGE6 | EDGE7 | EDGE8 | EDGE9 | EDGE10, 0, 0, 0}, // 115
{EDGE1 | EDGE3 | EDGE4 | EDGE6 | EDGE8 | EDGE10, 0, 0, 0}, // 116
{EDGE0 | EDGE1 | EDGE4 | EDGE6 | EDGE10, 0, 0, 0}, // 117
{EDGE0 | EDGE3 | EDGE4 | EDGE6 | EDGE8 | EDGE9 | EDGE10, 0, 0, 0}, // 118
{EDGE4 | EDGE6 | EDGE9 | EDGE10, 0, 0, 0}, // 119
{EDGE4 | EDGE5 | EDGE9, EDGE1 | EDGE3 | EDGE6 | EDGE7 | EDGE10, 0, 0}, // 120
{EDGE0 | EDGE1 | EDGE6 | EDGE7 | EDGE8 | EDGE10, EDGE4 | EDGE5 | EDGE9, 0, 0}, // 121
{EDGE0 | EDGE3 | EDGE4 | EDGE5 | EDGE6 | EDGE7 | EDGE10, 0, 0, 0}, // 122
{EDGE4 | EDGE5 | EDGE6 | EDGE7 | EDGE8 | EDGE10, 0, 0, 0}, // 123
{EDGE1 | EDGE3 | EDGE5 | EDGE6 | EDGE8 | EDGE9 | EDGE10, 0, 0, 0}, // 124
{EDGE0 | EDGE1 | EDGE5 | EDGE6 | EDGE9 | EDGE10, 0, 0, 0}, // 125
{EDGE0 | EDGE3 | EDGE8, EDGE5 | EDGE6 | EDGE10, 0, 0}, // 126
{EDGE5 | EDGE6 | EDGE10, 0, 0, 0}, // 127
{EDGE5 | EDGE6 | EDGE10, 0, 0, 0}, // 128
{EDGE0 | EDGE3 | EDGE8, EDGE5 | EDGE6 | EDGE10, 0, 0}, // 129
{EDGE0 | EDGE1 | EDGE9, EDGE5 | EDGE6 | EDGE10, 0, 0}, // 130
{EDGE1 | EDGE3 | EDGE8 | EDGE9, EDGE5 | EDGE6 | EDGE10, 0, 0}, // 131
{EDGE4 | EDGE7 | EDGE8, EDGE5 | EDGE6 | EDGE10, 0, 0}, // 132
{EDGE0 | EDGE3 | EDGE4 | EDGE7, EDGE5 | EDGE6 | EDGE10, 0, 0}, // 133
{EDGE0 | EDGE1 | EDGE9, EDGE4 | EDGE7 | EDGE8, EDGE5 | EDGE6 | EDGE10, 0}, // 134
{EDGE1 | EDGE3 | EDGE4 | EDGE7 | EDGE9, EDGE5 | EDGE6 | EDGE10, 0, 0}, // 135
{EDGE4 | EDGE6 | EDGE9 | EDGE10, 0, 0, 0}, // 136
{EDGE0 | EDGE3 | EDGE8, EDGE4 | EDGE6 | EDGE9 | EDGE10, 0, 0}, // 137
{EDGE0 | EDGE1 | EDGE4 | EDGE6 | EDGE10, 0, 0, 0}, // 138
{EDGE1 | EDGE3 | EDGE4 | EDGE6 | EDGE8 | EDGE10, 0, 0, 0}, // 139
{EDGE6 | EDGE7 | EDGE8 | EDGE9 | EDGE10, 0, 0, 0}, // 140
{EDGE0 | EDGE3 | EDGE6 | EDGE7 | EDGE9 | EDGE10, 0, 0, 0}, // 141
{EDGE0 | EDGE1 | EDGE6 | EDGE7 | EDGE8 | EDGE10, 0, 0, 0}, // 142
{EDGE1 | EDGE3 | EDGE6 | EDGE7 | EDGE10, 0, 0, 0}, // 143
{EDGE2 | EDGE3 | EDGE11, EDGE5 | EDGE6 | EDGE10, 0, 0}, // 144
{EDGE0 | EDGE2 | EDGE8 | EDGE11, EDGE5 | EDGE6 | EDGE10, 0, 0}, // 145
{EDGE0 | EDGE1 | EDGE9, EDGE2 | EDGE3 | EDGE11, EDGE5 | EDGE6 | EDGE10, 0}, // 146
{EDGE1 | EDGE2 | EDGE8 | EDGE9 | EDGE11, EDGE5 | EDGE6 | EDGE10, 0, 0}, // 147
{EDGE4 | EDGE7 | EDGE8, EDGE2 | EDGE3 | EDGE11, EDGE5 | EDGE6 | EDGE10, 0}, // 148
{EDGE0 | EDGE2 | EDGE4 | EDGE7 | EDGE11, EDGE5 | EDGE6 | EDGE10, 0, 0}, // 149
{EDGE0 | EDGE1 | EDGE9, EDGE4 | EDGE7 | EDGE8, EDGE2 | EDGE3 | EDGE11, EDGE5 | EDGE6 | EDGE10}, // 150
{EDGE1 | EDGE2 | EDGE4 | EDGE7 | EDGE9 | EDGE11, EDGE5 | EDGE6 | EDGE10, 0, 0}, // 151
{EDGE4 | EDGE6 | EDGE9 | EDGE10, EDGE2 | EDGE3 | EDGE11, 0, 0}, // 152
{EDGE0 | EDGE2 | EDGE8 | EDGE11, EDGE4 | EDGE6 | EDGE9 | EDGE10, 0, 0}, // 153
{EDGE0 | EDGE1 | EDGE4 | EDGE6 | EDGE10, EDGE2 | EDGE3 | EDGE11, 0, 0}, // 154
{EDGE1 | EDGE2 | EDGE4 | EDGE6 | EDGE8 | EDGE10 | EDGE11, 0, 0, 0}, // 155
{EDGE6 | EDGE7 | EDGE8 | EDGE9 | EDGE10, EDGE2 | EDGE3 | EDGE11, 0, 0}, // 156
{EDGE0 | EDGE2 | EDGE6 | EDGE7 | EDGE9 | EDGE10 | EDGE11, 0, 0, 0}, // 157
{EDGE0 | EDGE1 | EDGE6 | EDGE7 | EDGE8 | EDGE10, EDGE2 | EDGE3 | EDGE11, 0, 0}, // 158
{EDGE1 | EDGE2 | EDGE6 | EDGE7 | EDGE10 | EDGE11, 0, 0, 0}, // 159
{EDGE1 | EDGE2 | EDGE5 | EDGE6, 0, 0, 0}, // 160
{EDGE0 | EDGE3 | EDGE8, EDGE1 | EDGE2 | EDGE5 | EDGE6, 0, 0}, // 161
{EDGE0 | EDGE2 | EDGE5 | EDGE6 | EDGE9, 0, 0, 0}, // 162
{EDGE2 | EDGE3 | EDGE5 | EDGE6 | EDGE8 | EDGE9, 0, 0, 0}, // 163
{EDGE4 | EDGE7 | EDGE8, EDGE1 | EDGE2 | EDGE5 | EDGE6, 0, 0}, // 164
{EDGE0 | EDGE3 | EDGE4 | EDGE7, EDGE1 | EDGE2 | EDGE5 | EDGE6, 0, 0}, // 165
{EDGE0 | EDGE2 | EDGE5 | EDGE6 | EDGE9, EDGE4 | EDGE7 | EDGE8, 0, 0}, // 166
{EDGE2 | EDGE3 | EDGE4 | EDGE5 | EDGE6 | EDGE7 | EDGE9, 0, 0, 0}, // 167
{EDGE1 | EDGE2 | EDGE4 | EDGE6 | EDGE9, 0, 0, 0}, // 168
{EDGE0 | EDGE3 | EDGE8, EDGE1 | EDGE2 | EDGE4 | EDGE6 | EDGE9, 0, 0}, // 169
{EDGE0 | EDGE2 | EDGE4 | EDGE6, 0, 0, 0}, // 170
{EDGE2 | EDGE3 | EDGE4 | EDGE6 | EDGE8, 0, 0, 0}, // 171
{EDGE1 | EDGE2 | EDGE6 | EDGE7 | EDGE8 | EDGE9, 0, 0, 0}, // 172
{EDGE0 | EDGE1 | EDGE2 | EDGE3 | EDGE6 | EDGE7 | EDGE9, 0, 0, 0}, // 173
{EDGE0 | EDGE2 | EDGE6 | EDGE7 | EDGE8, 0, 0, 0}, // 174
{EDGE2 | EDGE3 | EDGE6 | EDGE7, 0, 0, 0}, // 175
{EDGE1 | EDGE3 | EDGE5 | EDGE6 | EDGE11, 0, 0, 0}, // 176
{EDGE0 | EDGE1 | EDGE5 | EDGE6 | EDGE8 | EDGE11, 0, 0, 0}, // 177
{EDGE0 | EDGE3 | EDGE5 | EDGE6 | EDGE9 | EDGE11, 0, 0, 0}, // 178
{EDGE5 | EDGE6 | EDGE8 | EDGE9 | EDGE11, 0, 0, 0}, // 179
{EDGE4 | EDGE7 | EDGE8, EDGE1 | EDGE3 | EDGE5 | EDGE6 | EDGE11, 0, 0}, // 180
{EDGE0 | EDGE1 | EDGE4 | EDGE5 | EDGE6 | EDGE7 | EDGE11, 0, 0, 0}, // 181
{EDGE0 | EDGE3 | EDGE5 | EDGE6 | EDGE9 | EDGE11, EDGE4 | EDGE7 | EDGE8, 0, 0}, // 182
{EDGE4 | EDGE5 | EDGE6 | EDGE7 | EDGE9 | EDGE11, 0, 0, 0}, // 183
{EDGE1 | EDGE3 | EDGE4 | EDGE6 | EDGE9 | EDGE11, 0, 0, 0}, // 184
{EDGE0 | EDGE1 | EDGE4 | EDGE6 | EDGE8 | EDGE9 | EDGE11, 0, 0, 0}, // 185
{EDGE0 | EDGE3 | EDGE4 | EDGE6 | EDGE11, 0, 0, 0}, // 186
{EDGE4 | EDGE6 | EDGE8 | EDGE11, 0, 0, 0}, // 187
{EDGE1 | EDGE3 | EDGE6 | EDGE7 | EDGE8 | EDGE9 | EDGE11, 0, 0, 0}, // 188
{EDGE0 | EDGE1 | EDGE9, EDGE6 | EDGE7 | EDGE11, 0, 0}, // 189
{EDGE0 | EDGE3 | EDGE6 | EDGE7 | EDGE8 | EDGE11, 0, 0, 0}, // 190
{EDGE6 | EDGE7 | EDGE11, 0, 0, 0}, // 191
{EDGE5 | EDGE7 | EDGE10 | EDGE11, 0, 0, 0}, // 192
{EDGE0 | EDGE3 | EDGE8, EDGE5 | EDGE7 | EDGE10 | EDGE11, 0, 0}, // 193
{EDGE0 | EDGE1 | EDGE9, EDGE5 | EDGE7 | EDGE10 | EDGE11, 0, 0}, // 194
{EDGE1 | EDGE3 | EDGE8 | EDGE9, EDGE5 | EDGE7 | EDGE10 | EDGE11, 0, 0}, // 195
{EDGE4 | EDGE5 | EDGE8 | EDGE10 | EDGE11, 0, 0, 0}, // 196
{EDGE0 | EDGE3 | EDGE4 | EDGE5 | EDGE10 | EDGE11, 0, 0, 0}, // 197
{EDGE0 | EDGE1 | EDGE9, EDGE4 | EDGE5 | EDGE8 | EDGE10 | EDGE11, 0, 0}, // 198
{EDGE1 | EDGE3 | EDGE4 | EDGE5 | EDGE9 | EDGE10 | EDGE11, 0, 0, 0}, // 199
{EDGE4 | EDGE7 | EDGE9 | EDGE10 | EDGE11, 0, 0, 0}, // 200
{EDGE0 | EDGE3 | EDGE8, EDGE4 | EDGE7 | EDGE9 | EDGE10 | EDGE11, 0, 0}, // 201
{EDGE0 | EDGE1 | EDGE4 | EDGE7 | EDGE10 | EDGE11, 0, 0, 0}, // 202
{EDGE1 | EDGE3 | EDGE4 | EDGE7 | EDGE8 | EDGE10 | EDGE11, 0, 0, 0}, // 203
{EDGE8 | EDGE9 | EDGE10 | EDGE11, 0, 0, 0}, // 204
{EDGE0 | EDGE3 | EDGE9 | EDGE10 | EDGE11, 0, 0, 0}, // 205
{EDGE0 | EDGE1 | EDGE8 | EDGE10 | EDGE11, 0, 0, 0}, // 206
{EDGE1 | EDGE3 | EDGE10 | EDGE11, 0, 0, 0}, // 207
{EDGE2 | EDGE3 | EDGE5 | EDGE7 | EDGE10, 0, 0, 0}, // 208
{EDGE0 | EDGE2 | EDGE5 | EDGE7 | EDGE8 | EDGE10, 0, 0, 0}, // 209
{EDGE0 | EDGE1 | EDGE9, EDGE2 | EDGE3 | EDGE5 | EDGE7 | EDGE10, 0, 0}, // 210
{EDGE1 | EDGE2 | EDGE5 | EDGE7 | EDGE8 | EDGE9 | EDGE10, 0, 0, 0}, // 211
{EDGE2 | EDGE3 | EDGE4 | EDGE5 | EDGE8 | EDGE10, 0, 0, 0}, // 212
{EDGE0 | EDGE2 | EDGE4 | EDGE5 | EDGE10, 0, 0, 0}, // 213
{EDGE0 | EDGE1 | EDGE9, EDGE2 | EDGE3 | EDGE4 | EDGE5 | EDGE8 | EDGE10, 0, 0}, // 214
{EDGE1 | EDGE2 | EDGE4 | EDGE5 | EDGE9 | EDGE10, 0, 0, 0}, // 215
{EDGE2 | EDGE3 | EDGE4 | EDGE7 | EDGE9 | EDGE10, 0, 0, 0}, // 216
{EDGE0 | EDGE2 | EDGE4 | EDGE7 | EDGE8 | EDGE9 | EDGE10, 0, 0, 0}, // 217
{EDGE0 | EDGE1 | EDGE2 | EDGE3 | EDGE4 | EDGE7 | EDGE10, 0, 0, 0}, // 218
{EDGE4 | EDGE7 | EDGE8, EDGE1 | EDGE2 | EDGE10, 0, 0}, // 219
{EDGE2 | EDGE3 | EDGE8 | EDGE9 | EDGE10, 0, 0, 0}, // 220
{EDGE0 | EDGE2 | EDGE9 | EDGE10, 0, 0, 0}, // 221
{EDGE0 | EDGE1 | EDGE2 | EDGE3 | EDGE8 | EDGE10, 0, 0, 0}, // 222
{EDGE1 | EDGE2 | EDGE10, 0, 0, 0}, // 223
{EDGE1 | EDGE2 | EDGE5 | EDGE7 | EDGE11, 0, 0, 0}, // 224
{EDGE0 | EDGE3 | EDGE8, EDGE1 | EDGE2 | EDGE5 | EDGE7 | EDGE11, 0, 0}, // 225
{EDGE0 | EDGE2 | EDGE5 | EDGE7 | EDGE9 | EDGE11, 0, 0, 0}, // 226
{EDGE2 | EDGE3 | EDGE5 | EDGE7 | EDGE8 | EDGE9 | EDGE11, 0, 0, 0}, // 227
{EDGE1 | EDGE2 | EDGE4 | EDGE5 | EDGE8 | EDGE11, 0, 0, 0}, // 228
{EDGE0 | EDGE1 | EDGE2 | EDGE3 | EDGE4 | EDGE5 | EDGE11, 0, 0, 0}, // 229
{EDGE0 | EDGE2 | EDGE4 | EDGE5 | EDGE8 | EDGE9 | EDGE11, 0, 0, 0}, // 230
{EDGE4 | EDGE5 | EDGE9, EDGE2 | EDGE3 | EDGE11, 0, 0}, // 231
{EDGE1 | EDGE2 | EDGE4 | EDGE7 | EDGE9 | EDGE11, 0, 0, 0}, // 232
{EDGE0 | EDGE3 | EDGE8, EDGE1 | EDGE2 | EDGE4 | EDGE7 | EDGE9 | EDGE11, 0, 0}, // 233
{EDGE0 | EDGE2 | EDGE4 | EDGE7 | EDGE11, 0, 0, 0}, // 234
{EDGE2 | EDGE3 | EDGE4 | EDGE7 | EDGE8 | EDGE11, 0, 0, 0}, // 235
{EDGE1 | EDGE2 | EDGE8 | EDGE9 | EDGE11, 0, 0, 0}, // 236
{EDGE0 | EDGE1 | EDGE2 | EDGE3 | EDGE9 | EDGE11, 0, 0, 0}, // 237
{EDGE0 | EDGE2 | EDGE8 | EDGE11, 0, 0, 0}, // 238
{EDGE2 | EDGE3 | EDGE11, 0, 0, 0}, // 239
{EDGE1 | EDGE3 | EDGE5 | EDGE7, 0, 0, 0}, // 240
{EDGE0 | EDGE1 | EDGE5 | EDGE7 | EDGE8, 0, 0, 0}, // 241
{EDGE0 | EDGE3 | EDGE5 | EDGE7 | EDGE9, 0, 0, 0}, // 242
{EDGE5 | EDGE7 | EDGE8 | EDGE9, 0, 0, 0}, // 243
{EDGE1 | EDGE3 | EDGE4 | EDGE5 | EDGE8, 0, 0, 0}, // 244
{EDGE0 | EDGE1 | EDGE4 | EDGE5, 0, 0, 0}, // 245
{EDGE0 | EDGE3 | EDGE4 | EDGE5 | EDGE8 | EDGE9, 0, 0, 0}, // 246
{EDGE4 | EDGE5 | EDGE9, 0, 0, 0}, // 247
{EDGE1 | EDGE3 | EDGE4 | EDGE7 | EDGE9, 0, 0, 0}, // 248
{EDGE0 | EDGE1 | EDGE4 | EDGE7 | EDGE8 | EDGE9, 0, 0, 0}, // 249
{EDGE0 | EDGE3 | EDGE4 | EDGE7, 0, 0, 0}, // 250
{EDGE4 | EDGE7 | EDGE8, 0, 0, 0}, // 251
{EDGE1 | EDGE3 | EDGE8 | EDGE9, 0, 0, 0}, // 252
{EDGE0 | EDGE1 | EDGE9, 0, 0, 0}, // 253
{EDGE0 | EDGE3 | EDGE8, 0, 0, 0}, // 254
{0, 0, 0, 0} // 255
};

// Encodes the ambiguous face of cube configurations, which can cause non-manifold meshes.
// The first bit of each value actually encodes a positive or negative
// direction while the second and third bit define the axis.

const unsigned char manifold[256] = { //!< Ambiguous face of cube configurations that cause non-manifold meshes.
255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
255,255,255,255,255,255,255,255,255,255,255,255,255,1,0,255,
255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
255,255,255,255,255,255,255,255,255,255,255,3,255,255,2,255,
255,255,255,255,255,255,255,5,255,255,255,255,255,255,5,5,
255,255,255,255,255,255,4,255,255,255,3,3,1,1,255,255,
255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
255,255,255,255,255,255,255,255,255,255,255,5,255,5,255,5,
255,255,255,255,255,255,255,3,255,255,255,255,255,2,255,255,
255,255,255,255,255,3,255,3,255,4,255,255,0,255,0,255,
255,255,255,255,255,255,255,1,255,255,255,0,255,255,255,255,
255,255,255,1,255,255,255,1,255,4,2,255,255,255,2,255,
255,255,255,0,255,2,4,255,255,255,255,0,255,2,255,255,
255,255,255,255,255,255,4,255,255,4,255,255,255,255,255,255
};

/*!
\brief Get the 12-bit dual point code mask.
*/
int getDualPointCode(const ScalarField& sf, int cx, int cy, int cz, int edge)
{
  int cubeCode = sf.DualCellCode(cx, cy, cz);

  // If a problematic C16 or C19 configuration shares the ambiguous face 
  // with another C16 or C19 configuration we simply invert the cube code
  // before looking up dual points. 

  const unsigned char direction = manifold[(unsigned char)(cubeCode)];
  // C16 or C19 configuration
  if (direction != 255)
  {
    // We have to check the neighboring cube, which shares the ambiguous
    // face. For this we decode the direction. This could also be done
    // with another lookup table.
    // copy current cube coordinates into an array.
    int neighborCoords[] = { cx,cy,cz };
    // get the dimension of the non-zero coordinate axis
    unsigned int component = direction >> 1;
    // get the sign of the direction
    int delta = (direction & 1) == 1 ? 1 : -1;
    // modify the correspong cube coordinate
    neighborCoords[component] += delta;
    // have we left the volume in this direction?

    int dims = sf.GetSizeX();
    if (component == 1) dims = sf.GetSizeY();
    if (component == 2) dims = sf.GetSizeZ();
    if (neighborCoords[component] >= 0 && neighborCoords[component] < (dims - 1))
    {
      // get the cube configuration of the relevant neighbor
      int neighborCubeCode = sf.DualCellCode(neighborCoords[0], neighborCoords[1], neighborCoords[2]);
      // Look up the neighbor configuration ambiguous face direction.
      // If the direction is valid we have a C16 or C19 neighbor.
      // As C16 and C19 have exactly one ambiguous face this face is
      // guaranteed to be shared for the pair.
      if (manifold[(unsigned char)(neighborCubeCode)] != 255) {
        // replace the cube configuration with its inverse.
        cubeCode ^= 0xff;
      }
    }
  }

  for (int i = 0; i < 4; ++i)
  {
    if (dualvertex[cubeCode][i] & edge)
    {
      return dualvertex[cubeCode][i];
    }
  }
  return 0;
}

/*!
\brief Compute the dual point.

The dual point is defined as the mean of the face vertices belonging to the marching cubes face.
\param c Configuration.
\param cx,cy,cz Coordinates.
\param v Returned vertex.
*/
void ScalarField::DualVertex(int cx, int cy, int cz, int c, Vector& v) const
{
  v = Vector::Null;
  int points = 0;

  // sum edge intersection vertices using the point code
  if (c & EDGE0) {
    v += Vector::Solve(Vertex(cx, cy, cz), Vertex(cx + 1, cy, cz), at(cx, cy, cz), at(cx + 1, cy, cz));
    points++;
  }

  if (c & EDGE1) {
    v += Vector::Solve(Vertex(cx + 1, cy, cz), Vertex(cx + 1, cy, cz + 1), at(cx + 1, cy, cz), at(cx + 1, cy, cz + 1));
    points++;
  }

  if (c & EDGE2) {
    v += Vector::Solve(Vertex(cx, cy, cz + 1), Vertex(cx + 1, cy, cz + 1), at(cx, cy, cz + 1), at(cx + 1, cy, cz + 1));
    points++;
  }

  if (c & EDGE3) {
    v += Vector::Solve(Vertex(cx, cy, cz), Vertex(cx, cy, cz + 1), at(cx, cy, cz), at(cx, cy, cz + 1));
    points++;
  }

  if (c & EDGE4) {
    v += Vector::Solve(Vertex(cx, cy + 1, cz), Vertex(cx + 1, cy + 1, cz), at(cx, cy + 1, cz), at(cx + 1, cy + 1, cz));
    points++;
  }

  if (c & EDGE5) {
    v += Vector::Solve(Vertex(cx + 1, cy + 1, cz), Vertex(cx + 1, cy + 1, cz + 1), at(cx + 1, cy + 1, cz), at(cx + 1, cy + 1, cz + 1));
    points++;
  }

  if (c & EDGE6) {
    v += Vector::Solve(Vertex(cx, cy + 1, cz + 1), Vertex(cx + 1, cy + 1, cz + 1), at(cx, cy + 1, cz + 1), at(cx + 1, cy + 1, cz + 1));
    points++;
  }

  if (c & EDGE7) {
    v += Vector::Solve(Vertex(cx, cy + 1, cz), Vertex(cx, cy + 1, cz + 1), at(cx, cy + 1, cz), at(cx, cy + 1, cz + 1));
    points++;
  }

  if (c & EDGE8) {
    v += Vector::Solve(Vertex(cx, cy, cz), Vertex(cx, cy + 1, cz), at(cx, cy, cz), at(cx, cy + 1, cz));
    points++;
  }

  if (c & EDGE9) {
    v += Vector::Solve(Vertex(cx + 1, cy, cz), Vertex(cx + 1, cy + 1, cz), at(cx + 1, cy, cz), at(cx + 1, cy + 1, cz));
    points++;
  }

  if (c & EDGE10) {
    v += Vector::Solve(Vertex(cx + 1, cy, cz + 1), Vertex(cx + 1, cy + 1, cz + 1), at(cx + 1, cy, cz + 1), at(cx + 1, cy + 1, cz + 1));
    points++;
  }

  if (c & EDGE11) {
    v += Vector::Solve(Vertex(cx, cy, cz + 1), Vertex(cx, cy + 1, cz + 1), at(cx, cy, cz + 1), at(cx, cy + 1, cz + 1));
    points++;
  }

  v /= double(points);
}

/*!
\brief Get the shared index of a dual point.
*/
int getSharedDualPointIndex(const ScalarField& sf, int cx, int cy, int cz, int edge, std::vector<Vector>& vertices)
{
  // create a key for the dual point from its linearized cell ID and point code
  int pc = getDualPointCode(sf, cx, cy, cz, edge);
  int64_t kkey = Codec(sf.VertexIndex(cx, cy, cz), pc);

  // have we already computed the dual point?
  auto iterator = pointToIndex.find(kkey);
  if (iterator != pointToIndex.end())
  {
    // just return the dual point index
    return iterator->second;
  }
  else
  {
    // create new Vector and Vector id
    int newVertexId = int(vertices.size());
    vertices.emplace_back();
    sf.DualVertex(cx, cy, cz, pc, vertices.back());
    // insert vertex ID into map and also return it
    pointToIndex[kkey] = newVertexId;
    return newVertexId;
  }
}

static void DoQ(std::vector<int>& quads, int i0, int i1, int i2, int i3)
{
  quads.emplace_back(i0);
  quads.emplace_back(i1);
  quads.emplace_back(i2);
  quads.emplace_back(i3);
}

//  Coordinate system
//
//       y
//       |
//       |
//       |
//       0-----x
//      /
//     /
//    z
//

//         Cell Edges
//  
//       o--------4----------o
//      /|                  /|
//     7 |                 5 |
//    /  |                /  |
//   o--------6----------o   |
//   |   8               |   9
//   |   |               |   |
//   |   |               |   |
//   11  |               10  |
//   |   o--------0------|---o
//   |  /                |  /
//   | 3                 | 1
//   |/                  |/
//   o--------2----------o
//

/*!
\brief Polygonization using Dual Marching Cubes.

This function implements a modified version of the dual marching cubes presented in:
G. Nielson. Dual Marching Cubes. <I>Proceedings of the Conference on Visualization</I>, 489-496, 2004.

A vertex in standard marching cubes is shared by 4 faces (as it lies on an edge, thus shared
by four cubes), thus the dual mesh is made from quadrangles. Manifold
mesh is produced by using the modifications described in R. Wenger. <I> Isosurfaces: Geometry, Topology, and Algorithms</I>.

\param g %Returned mesh.
*/
void ScalarField::Dual(Mesh& g) const
{
  std::vector<Vector> vertices;
  std::vector<int> quads;

  pointToIndex.clear();

  // iterate voxels
  for (int z = 0; z < nz - 1; z++)
  {
    for (int y = 0; y < ny - 1; y++)
    {
      for (int x = 0; x < nx - 1; x++)
      {
        // X edge
        if (z > 0 && y > 0)
        {
          bool const entering = at(x, y, z) < 0.0 && at(x + 1, y, z) >= 0.0;
          bool const exiting = at(x, y, z) >= 0.0 && at(x + 1, y, z) < 0.0;
          if (entering || exiting)
          {
            // generate 
            int i0 = getSharedDualPointIndex(*this, x, y, z, EDGE0, vertices);
            int i1 = getSharedDualPointIndex(*this, x, y, z - 1, EDGE2, vertices);
            int i2 = getSharedDualPointIndex(*this, x, y - 1, z - 1, EDGE6, vertices);
            int i3 = getSharedDualPointIndex(*this, x, y - 1, z, EDGE4, vertices);

            if (entering)
            {
              DoQ(quads, i0, i1, i2, i3);
            }
            else
            {
              DoQ(quads, i0, i3, i2, i1);
            }
          }
        }

        // Y edge
        if (z > 0 && x > 0)
        {
          bool const entering = at(x, y, z) < 0.0 && at(x, y + 1, z) >= 0.0;
          bool const exiting = at(x, y, z) >= 0.0 && at(x, y + 1, z) < 0.0;
          if (entering || exiting)
          {
            // generate quad
            int i0 = getSharedDualPointIndex(*this, x, y, z, EDGE8, vertices);
            int i1 = getSharedDualPointIndex(*this, x, y, z - 1, EDGE11, vertices);
            int i2 = getSharedDualPointIndex(*this, x - 1, y, z - 1, EDGE10, vertices);
            int i3 = getSharedDualPointIndex(*this, x - 1, y, z, EDGE9, vertices);

            if (exiting)
            {
              DoQ(quads, i0, i1, i2, i3);
            }
            else
            {
              DoQ(quads, i0, i3, i2, i1);
            }
          }
        }

        // Z edge
        if (x > 0 && y > 0)
        {
          bool const entering = at(x, y, z) < 0.0 && at(x, y, z + 1) >= 0.0;
          bool const exiting = at(x, y, z) >= 0.0 && at(x, y, z + 1) < 0.0;
          if (entering || exiting)
          {
            // generate quad
            int i0 = getSharedDualPointIndex(*this, x, y, z, EDGE3, vertices);
            int i1 = getSharedDualPointIndex(*this, x - 1, y, z, EDGE1, vertices);
            int i2 = getSharedDualPointIndex(*this, x - 1, y - 1, z, EDGE5, vertices);
            int i3 = getSharedDualPointIndex(*this, x, y - 1, z, EDGE7, vertices);

            if (exiting)
            {
              DoQ(quads, i0, i1, i2, i3);
            }
            else
            {
              DoQ(quads, i0, i3, i2, i1);
            }
          }
        }
      }
    }
  }
  QVector<Vector> v;
  QVector<int> t;

  for (int k = 0; k < vertices.size(); k++)
  {
    v.append(vertices[k]);
  }
  for (int k = 0; k < quads.size(); k += 4)
  {
    t.append(quads[k]);
    t.append(quads[k + 1]);
    t.append(quads[k + 2]);

    t.append(quads[k]);
    t.append(quads[k + 2]);
    t.append(quads[k + 3]);
  }

  // Normals
  QVector<Vector> n(v.size());

  for (int i = 0; i < v.size(); i++)
  {
    n[i] = Normal(v[i]);
  }
  // Mesh
  g = Mesh(v, n, t, t);
}

/*!
\brief Polygonization.
\param g %Returned mesh.
*/
void ScalarField::Polygonize(Mesh& g) const
{
  QVector<Vector> vertex;
  QVector<Vector> normal;
  QVector<int> varray;
  QVector<int> narray;
  vertex.reserve(200000);
  normal.reserve(200000);
  varray.reserve(200000);
  narray.reserve(200000);

  int nv = 0;

  Box clipped = GetBox();

  const int size = nx * ny;

  // Intensities
  double* a = new double[size];
  double* b = new double[size];

  // Vertex
  Vector2* u = new Vector2[size];

  // Edges
  int* eax = new int[size];
  int* eay = new int[size];
  int* ebx = new int[size];
  int* eby = new int[size];
  int* ez = new int[size];

  double za = clipped[0][2];

  // Compute field inside lower Oxy plane
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      u[i * ny + j] = clipped[0] + Vector2(i * celldiagonal[0], j * celldiagonal[1]);
      a[i * ny + j] = at(i, j, 0);
    }
  }

  // Compute straddling edges inside lower Oxy plane
  for (int i = 0; i < nx - 1; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      // We need a xor b, which can be implemented a == !b 
      if (!((a[i * ny + j] < 0.0) == !(a[(i + 1) * ny + j] >= 0.0)))
      {
        vertex.append(Vector::Solve(u[i * ny + j].ToVector(za), u[(i + 1) * ny + j].ToVector(za), a[i * ny + j], a[(i + 1) * ny + j]));
        normal.append(Normal(vertex.last()));
        eax[i * ny + j] = nv;
        nv++;
      }
    }
  }
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny - 1; j++)
    {
      if (!((a[i * ny + j] < 0.0) == !(a[i * ny + (j + 1)] >= 0.0)))
      {
        vertex.append(Vector::Solve(u[i * ny + j].ToVector(za), u[i * ny + (j + 1)].ToVector(za), a[i * ny + j], a[i * ny + (j + 1)]));
        normal.append(Normal(vertex.last()));
        eay[i * ny + j] = nv;
        nv++;
      }
    }
  }

  // Array for edge vertices
  int e[12];

  // For all layers
  for (int k = 0; k < nz - 1; k++)
  {
    double zb = za + celldiagonal[2];
    for (int i = 0; i < nx; i++)
    {
      for (int j = 0; j < ny; j++)
        b[i * ny + j] = at(i, j, k + 1);
    }

    // Compute straddling edges inside lower Oxy plane
    for (int i = 0; i < nx - 1; i++)
    {
      for (int j = 0; j < ny; j++)
      {
        if (!((b[i * ny + j] < 0.0) == !(b[(i + 1) * ny + j] >= 0.0)))
        {
          vertex.append(Vector::Solve(u[i * ny + j].ToVector(zb), u[(i + 1) * ny + j].ToVector(zb), b[i * ny + j], b[(i + 1) * ny + j]));
          normal.append(Normal(vertex.last()));
          ebx[i * ny + j] = nv;
          nv++;
        }
      }
    }

    for (int i = 0; i < nx; i++)
    {
      for (int j = 0; j < ny - 1; j++)
      {
        if (!((b[i * ny + j] < 0.0) == !(b[i * ny + (j + 1)] >= 0.0)))
        {
          vertex.append(Vector::Solve(u[i * ny + j].ToVector(zb), u[i * ny + (j + 1)].ToVector(zb), b[i * ny + j], b[i * ny + (j + 1)]));
          normal.append(Normal(vertex.last()));
          eby[i * ny + j] = nv;
          nv++;
        }
      }
    }

    // Create vertical straddling edges
    for (int i = 0; i < nx; i++)
    {
      for (int j = 0; j < ny; j++)
      {
        if (!((a[i * ny + j] < 0.0) == !(b[i * ny + j] >= 0.0)))
        {
          vertex.append(Vector::Solve(u[i * ny + j].ToVector(za), u[i * ny + j].ToVector(zb), a[i * ny + j], b[i * ny + j]));
          normal.append(Normal(vertex.last()));
          ez[i * ny + j] = nv;
          nv++;
        }
      }
    }
    // Create mesh
    for (int i = 0; i < nx - 1; i++)
    {
      for (int j = 0; j < ny - 1; j++)
      {
        int cubeindex = 0;
        if (a[i * ny + j] < 0.0)       cubeindex |= 1;
        if (a[(i + 1) * ny + j] < 0.0)   cubeindex |= 2;
        if (a[i * ny + j + 1] < 0.0)     cubeindex |= 4;
        if (a[(i + 1) * ny + j + 1] < 0.0) cubeindex |= 8;
        if (b[i * ny + j] < 0.0)       cubeindex |= 16;
        if (b[(i + 1) * ny + j] < 0.0)   cubeindex |= 32;
        if (b[i * ny + j + 1] < 0.0)     cubeindex |= 64;
        if (b[(i + 1) * ny + j + 1] < 0.0) cubeindex |= 128;

        // Cube is straddling the surface
        if ((cubeindex != 255) && (cubeindex != 0))
        {
          e[0] = eax[i * ny + j];
          e[1] = eax[i * ny + (j + 1)];
          e[2] = ebx[i * ny + j];
          e[3] = ebx[i * ny + (j + 1)];
          e[4] = eay[i * ny + j];
          e[5] = eay[(i + 1) * ny + j];
          e[6] = eby[i * ny + j];
          e[7] = eby[(i + 1) * ny + j];
          e[8] = ez[i * ny + j];
          e[9] = ez[(i + 1) * ny + j];
          e[10] = ez[i * ny + (j + 1)];
          e[11] = ez[(i + 1) * ny + (j + 1)];

          for (int h = 0; TriangleTable[cubeindex][h] != -1; h += 3)
          {
            int v0 = e[TriangleTable[cubeindex][h + 0]];
            int v1 = e[TriangleTable[cubeindex][h + 1]];
            int v2 = e[TriangleTable[cubeindex][h + 2]];
            varray.append(v0); varray.append(v1); varray.append(v2);
            narray.append(v0); narray.append(v1); narray.append(v2);
          }
        }
      }
    }

    za += celldiagonal[2];
    Math::Swap(a, b);
    Math::Swap(eax, ebx);
    Math::Swap(eay, eby);
  }

  delete[]a;
  delete[]b;
  delete[]u;

  delete[]eax;
  delete[]eay;
  delete[]ebx;
  delete[]eby;
  delete[]ez;

  g = Mesh(vertex, normal, varray, narray);
}

#include <QtCore/QFile>

/*!
\brief Save as data file.

\param name Filename.
*/
void ScalarField::Save(const QString& name) const
{
  QFile file(name);
  if (!file.open(QFile::WriteOnly))
  {
    std::cout << "Could not open file for writing";
    return;
  }
  QDataStream out(&file);

  Array::OutStream(out);

  for (int i = 0; i < field.size(); i++)
  {
    double v = field.at(i);
    v = (1.0 + Math::Clamp(v, -1.0, 1.0)) * 0.5;
    unsigned char c = (unsigned char)(v * 255.0);
    out << c;
  }
  file.close();
}

/*!
\brief Read data file.
\param name Filename.
*/
void ScalarField::Read(const QString& name)
{
  QFile file(name);
  if (!file.open(QFile::ReadOnly))
  {
    std::cout << "Could not open file for reading";
    return;
  }
  QDataStream in(&file);

  Array::InStream(in);

  for (int i = 0; i < field.size(); i++)
  {
    unsigned char c;
    in >> c;
    field[i] = 2.0 * (c / 255.0) - 1.0;
  }
  file.close();
}



/*
\brief Compute the minimum curvature.

See Schmidt et al 2003, Comparison of polynomial models for land surface curvature calculation,
for a complete reference of the formulations of the different curvatures.

Max curvature is also expressed as H - sqrt(H^2 - K) with H being mean curvature and K gaussian curvature.
See Table 2.2 in Florinsky's book "Digital Terrain Analysis in Soil Science and Geology"

\author Oscar Argudo
*/
ScalarField2 ScalarField2::CurvatureMin() const
{
  ScalarField2 curv(Box2(a, b), nx, ny, 0.0);

  ScalarField2 cGauss = CurvatureGaussian();
  ScalarField2 cMean = CurvatureMean();

  for (int i = 1; i < nx - 1; i++)
  {
    for (int j = 1; j < ny - 1; j++)
    {
      double k = cGauss.at(i, j);
      double h = cMean.at(i, j);
      curv(i, j) = h - sqrt(h * h - k);
    }
  }
  return curv;
}

/*
\brief Compute the maximum curvature.

See Schmidt et al 2003, Comparison of polynomial models for land surface curvature calculation,
for a complete reference of the formulations of the different curvatures.

Max curvature is also expressed as H + sqrt(H^2 - K) with H being mean curvature and K gaussian curvature.
See Table 2.2 in Florinsky's book "Digital Terrain Analysis in Soil Science and Geology"

\author Oscar Argudo
*/
ScalarField2 ScalarField2::CurvatureMax() const
{
  ScalarField2 curv(Box2(a, b), nx, ny, 0.0);

  ScalarField2 cGauss = CurvatureGaussian();
  ScalarField2 cMean = CurvatureMean();

  for (int i = 1; i < nx - 1; i++)
  {
    for (int j = 1; j < ny - 1; j++)
    {
      double k = cGauss.at(i, j);
      double h = cMean.at(i, j);
      curv(i, j) = h + sqrt(h * h - k);
    }
  }
  return curv;
}
