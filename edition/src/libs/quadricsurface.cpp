// Quadric surfaces

#include "libs/quadricsurface.h"

/*!
\class QuadricSurface quadricsurface.h
\brief %Quadric surfaces.

\ingroup MathGroup
*/

/*!
\brief Create a quadric surface.

z=a22  x<SUP>2</SUP>y<sup>2</sup> + a21 x<sup>2</sup>y +  a21 xy<sup>2</sup> + a20 x<sup>2</sup> + a02 y<sup>2</sup> +  a11 xy +  a10 x + a01 y + a00.
*/
QuadricSurface::QuadricSurface(const double& a22, const double& a21, const double& a12, const double& a20, const double& a02, const double& a11, const double& a10, const double& a01, const double& a00) : c{ a22, a21, a12, a20, a02, a11, a10, a01, a00 }
{
}

/*!
\brief Create a quadric surface from elevations.
\param m Matrix of elevations.
\param e Distance.
*/
QuadricSurface::QuadricSurface(const Matrix& m, const double& e)
{
  c[0] = ((m(0, 2) + m(2, 2) + m(0, 0) + m(2, 0)) / 4.0 - (m(1, 2) + m(0, 1) + m(2, 1) + m(1, 0)) / 2.0 + m(1, 1)) / (e * e * e * e);
  c[1] = ((m(0, 2) + m(2, 2) - m(0, 0) - m(2, 0)) / 4.0 - (m(1, 2) - m(1, 0)) / 2.0) / (e * e * e);
  c[2] = ((-m(0, 2) + m(2, 2) - m(0, 0) + m(2, 0)) / 4.0 + (m(0, 1) - m(2, 1)) / 2.0) / (e * e * e);
  c[3] = ((m(0, 1) + m(2, 1)) / 2.0 - m(1, 1)) / (e * e);
  c[4] = ((m(1, 2) + m(1, 0)) / 2.0 - m(1, 1)) / (e * e);
  c[5] = (-m(0, 2) + m(2, 2) + m(0, 0) - m(2, 0)) / (4.0 * e * e);
  c[6] = (-m(0, 1) + m(2, 1)) / (2.0 * e);
  c[7] = (m(1, 2) - m(1, 0)) / (2.0 * e);
  c[8] = m(1, 1);
}

/*!
\brief Compute gradient.
\param p Point.
*/
Vector2 QuadricSurface::Gradient(const Vector2& p) const
{
  const double x = p[0];
  const double y = p[1];

  return Vector2(
    2.0 * c[0] * x * y * y + 2.0 * c[1] * x * y + c[2] * y * y + 2.0 * c[3] * x + c[5] * y + c[6],
    2.0 * c[0] * x * x * y + c[1] * x * x + 2.0 * c[2] * x * y + 2.0 * c[4] * y + c[5] * x + c[7]);
}

/*!
\brief Compute normal.
\param p Point.
*/
Vector2 QuadricSurface::Normal(const Vector2& p) const
{
  Vector2 g = Gradient(p);

  return Normalized(Vector(-g[0], -g[1], 1.0));
}

/*!
\brief Translate a quadric surface.
\param t Translation vector.
*/
QuadricSurface QuadricSurface::Translated(const Vector2& t) const
{
  const double& x = t[0];
  const double& y = t[1];

  const double xx = x * x;
  const double yy = y * y;
  const double xy = x * y;

  return QuadricSurface(
    c[0],
    +c[1] + 2.0 * c[0] * y,
    +c[2] + 2.0 * c[0] * x,
    +c[3] + c[1] * y + c[0] * yy,
    +c[4] + c[0] * xx + c[2] * x,
    +c[5] + 2.0 * c[2] * y + 2.0 * c[1] * x + 4.0 * c[0] * x,
    +c[6] + y * c[5] + 2.0 * c[3] * x + c[2] * yy + 2.0 * c[1] * xy + 2.0 * c[0] * x * yy,
    +c[7] + x * c[5] + 2.0 * c[4] * y + 2.0 * c[2] * xy + c[1] * xx + 2.0 * c[0] * xx * y,
    +c[8] + c[7] * y + c[6] * x + xy * c[5] + c[3] * xx + c[2] * x * yy + c[1] * xx * y + c[0] * xx * yy + c[4] * yy
  );
}

/*!
\brief Scale a quadric surface.
\param s Scaling coefficient.
*/
QuadricSurface QuadricSurface::Scaled(const double& s) const
{
  double is = 1.0 / s;
  return QuadricSurface(c[0] * is * is * is * is, c[1] * is * is * is, c[2] * is * is * is, c[3] * is * is, c[4] * is * is, c[5] * is * is, c[6] * is, c[7] * is, c[8]);
}