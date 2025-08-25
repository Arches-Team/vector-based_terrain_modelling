#include "libs/analyticfield.h"

/*!
\class AnalyticHeightField analyticfield.h
\brief A base two-dimensional analytic field of real values.
*/

double AnalyticHeightField::Epsilon = 1e-6; //!< Epsilon value used for all fields.

/*!
\brief Compute the normal to the field.
\param p Point.
*/
Vector AnalyticHeightField::Normal(const Vector2& p) const
{
  Vector2 g = Gradient(p);
  return Normalized(Vector(-g[0], -g[1], 1.0));
}

/*!
\brief Create a HeightField given an input domain.
\param box The box.
\param nx, ny Resolution.
*/
HeightField AnalyticHeightField::CreateHeightField(const Box2& box, int nx, int ny) const
{
  HeightField field(box, nx, ny);
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      Vector2 q = field.ArrayVertex(i, j);
      field(i, j) = Height(q);
    }
  }
  return field;
}

/*!
\brief Compute the slope at a given point.

This is the maximum slope, which is the norm of the gradient.

\sa AnalyticHeightField::Gradient()
\param p The point.
*/
double AnalyticHeightField::Slope(const Vector2& p) const
{
  return Norm(Gradient(p));
}

/*!
\brief Compute the slope in a given direction and at a given point.

\param p The point.
\param d Direction (should be normalized).
*/
double AnalyticHeightField::Slope(const Vector2& p, const Vector2& d) const
{
  return (Height(p + Epsilon * d) - Height(p - Epsilon * d)) / (2.0 * Epsilon);
}

/*!
\brief Compute the average slope at a given point.

\param p The point.
\param n Number of sampling directions.
*/
double AnalyticHeightField::AverageSlope(const Vector2& p, int n) const
{
  double s = 0.0;
  for (int i = 0; i < n; i++)
  {
    double t = (i * Math::Pi) / n;

    // Direction computed with rotation 
    Vector2 u = Vector2(cos(t), sin(t));

    // Slope
    s += fabs(Slope(p, u));

    // Slope in orthogonal direction
    s += fabs(Slope(p, u.Orthogonal()));
  }
  return s;
}

/*!
\brief Compute the elevation along a segment.

\param a,b %Segment.
\param n Number of samples.
*/
QVector<double> AnalyticHeightField::Cross(const Vector2& a, const Vector2& b, int n) const
{
  QVector<double> cut;

  cut.reserve(n);

  cut.append(Height(a));

  for (int k = 1; k < n - 1; k++)
  {
    Vector2 p = Vector2::Lerp(a, b, Math::Unit(k, n));
    cut.append(Height(p));
  }

  cut.append(Height(b));

  return cut;
}
