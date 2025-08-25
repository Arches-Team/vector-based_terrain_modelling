// Smooth ellipse

#include "libs/smooth.h"
#include "libs/cubic.h"

/*!
\class SmoothEllipse2 smooth.h
\brief An ellipse with a smooth falloff, resembles a compact scaled Gaussian.
\ingroup PlanarGroup
*/

/*!
\brief Create a smooth ellipse.
\param e %Ellipse.
*/

SmoothEllipse2::SmoothEllipse2(const Ellipse2& e) :Ellipse2(e)
{
}


/*!
\brief Compute the intensity.
\param p Point.
*/
double SmoothEllipse2::Value(const Vector2& p) const
{
  Vector2 q = p - c;
  q = Vector2(u * q, u.Orthogonal() * q);

  double x = q[0];
  double y = q[1];

  // Distance to center such that d = Min(a, b) on the border of the ellipse
  // double d = Math::Min(a, b) * sqrt(x * x / (a * a) + y * y / (b * b));
  // Instead, use normalized distance
  return Cubic::SmoothCompact(x * x / (a * a) + y * y / (b * b), 1.0);
}
