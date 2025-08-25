// Fields

#include "libs/scalarfield.h"

/*
\brief Compute the profile curvature, i.e., the curvature along flow direction

\sa CurvatureContour

\author Oscar Argudo
*/
ScalarField2 ScalarField2::CurvatureProfile() const
{
  ScalarField2 curv(Box2(a, b), nx, ny, 0.0);

  const double dx = celldiagonal[0];
  const double dy = celldiagonal[1];

  for (int i = 1; i < nx - 1; i++)
  {
    for (int j = 1; j < ny - 1; j++)
    {
      double z_x = (at(i + 1, j) - at(i - 1, j)) / (2.0 * dx);
      double z_y = (at(i, j + 1) - at(i, j - 1)) / (2.0 * dy);
      double z_xx = (at(i + 1, j) + at(i - 1, j) - 2.0 * at(i, j)) / (dx * dx);
      double z_yy = (at(i, j + 1) + at(i, j - 1) - 2.0 * at(i, j)) / (dy * dy);
      double z_xy = (at(i + 1, j + 1) + at(i - 1, j - 1) - at(i + 1, j - 1) - at(i - 1, j + 1)) / (4.0 * dx * dy);
      double p = z_x * z_x + z_y * z_y;
      double q = p + 1.0;

      if (p > 1e-16)
        curv(i, j) = (z_xx * z_x * z_x + 2 * z_xy * z_x * z_y + z_yy * z_y * z_y) / (p * Math::Sqrt32(q));
    }
  }
  return curv;
}

/*
\brief Compute the Contour Curvature, i.e. curvature along the XY plane

See [Schmidt et al 2003, Comparison of polynomial models for land surface curvature calculation]
for a complete reference of the formulations of the different curvatures.

\author Oscar Argudo
*/
ScalarField2 ScalarField2::CurvatureContour() const
{
  ScalarField2 curv(Box2(a, b), nx, ny, 0.0);

  const double dx = celldiagonal[0];
  const double dy = celldiagonal[1];

  for (int i = 1; i < nx - 1; i++)
  {
    for (int j = 1; j < ny - 1; j++)
    {
      double z_x = (at(i + 1, j) - at(i - 1, j)) / (2.0 * dx);
      double z_y = (at(i, j + 1) - at(i, j - 1)) / (2.0 * dy);
      double z_xx = (at(i + 1, j) + at(i - 1, j) - 2.0 * at(i, j)) / (dx * dx);
      double z_yy = (at(i, j + 1) + at(i, j - 1) - 2.0 * at(i, j)) / (dy * dy);
      double z_xy = (at(i + 1, j + 1) + at(i - 1, j - 1) - at(i + 1, j - 1) - at(i - 1, j + 1)) / (4.0 * dx * dy);
      double p = z_x * z_x + z_y * z_y;

      if (p > 1e-16)
        curv(i, j) = (z_xx * z_y * z_y - 2.0 * z_xy * z_x * z_y + z_yy * z_x * z_x) / Math::Sqrt32(p);
    }
  }
  return curv;
}

/*
\brief Compute the Tangential Curvature, i.e. curvature on a plane tangential to the surface.
Similar but not exactly the same as the Contour Curvature.

\sa CurvatureContour

\author Oscar Argudo
*/
ScalarField2 ScalarField2::CurvatureTangential() const
{
  ScalarField2 curv(Box2(a, b), nx, ny, 0.0);

  const double dx = celldiagonal[0];
  const double dy = celldiagonal[1];

  for (int i = 1; i < nx - 1; i++)
  {
    for (int j = 1; j < ny - 1; j++)
    {
      double z_x = (at(i + 1, j) - at(i - 1, j)) / (2.0 * dx);
      double z_y = (at(i, j + 1) - at(i, j - 1)) / (2.0 * dy);
      double z_xx = (at(i + 1, j) + at(i - 1, j) - 2.0 * at(i, j)) / (dx * dx);
      double z_yy = (at(i, j + 1) + at(i, j - 1) - 2.0 * at(i, j)) / (dy * dy);
      double z_xy = (at(i + 1, j + 1) + at(i - 1, j - 1) - at(i + 1, j - 1) - at(i - 1, j + 1)) / (4.0 * dx * dy);
      double p = z_x * z_x + z_y * z_y;
      double q = p + 1.0;

      if (p > 1e-16)
        curv(i, j) = (z_xx * z_y * z_y - 2.0 * z_xy * z_x * z_y + z_yy * z_x * z_x) / (p * sqrt(q));
    }
  }
  return curv;
}

/*
\brief Compute the Gaussian Curvature

\sa CurvatureContour

\author Oscar Argudo
*/
ScalarField2 ScalarField2::CurvatureGaussian() const
{
  ScalarField2 curv(Box2(a, b), nx, ny, 0.0);

  const double dx = celldiagonal[0];
  const double dy = celldiagonal[1];

  for (int i = 1; i < nx - 1; i++)
  {
    for (int j = 1; j < ny - 1; j++)
    {
      double z_x = (at(i + 1, j) - at(i - 1, j)) / (2.0 * dx);
      double z_y = (at(i, j + 1) - at(i, j - 1)) / (2.0 * dy);
      double z_xx = (at(i + 1, j) + at(i - 1, j) - 2.0 * at(i, j)) / (dx * dx);
      double z_yy = (at(i, j + 1) + at(i, j - 1) - 2.0 * at(i, j)) / (dy * dy);
      double z_xy = (at(i + 1, j + 1) + at(i - 1, j - 1) - at(i + 1, j - 1) - at(i - 1, j + 1)) / (4.0 * dx * dy);
      double p = z_x * z_x + z_y * z_y;
      double q = p + 1.0;

      curv(i, j) = (z_xx * z_yy - z_xy * z_xy) / (q * q);
    }
  }
  return curv;
}

/*
\brief Compute the Mean Curvature

\sa CurvatureContour

\author Oscar Argudo
*/
ScalarField2 ScalarField2::CurvatureMean() const
{
  ScalarField2 curv(Box2(a, b), nx, ny, 0.0);

  const double dx = celldiagonal[0];
  const double dy = celldiagonal[1];

  for (int i = 1; i < nx - 1; i++)
  {
    for (int j = 1; j < ny - 1; j++)
    {
      double z_x = (at(i + 1, j) - at(i - 1, j)) / (2.0 * dx);
      double z_y = (at(i, j + 1) - at(i, j - 1)) / (2.0 * dy);
      double z_xx = (at(i + 1, j) + at(i - 1, j) - 2.0 * at(i, j)) / (dx * dx);
      double z_yy = (at(i, j + 1) + at(i, j - 1) - 2.0 * at(i, j)) / (dy * dy);
      double z_xy = (at(i + 1, j + 1) + at(i - 1, j - 1) - at(i + 1, j - 1) - at(i - 1, j + 1)) / (4.0 * dx * dy);
      double p = z_x * z_x + z_y * z_y;
      double q = p + 1.0;

      curv(i, j) = ((1 + z_y * z_y) * z_xx - 2 * z_xy * z_x * z_y + (1 + z_x * z_x) * z_yy) / (2.0 * Math::Sqrt32(q));
    }
  }
  return curv;
}
