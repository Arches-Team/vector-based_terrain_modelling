// Set of cubic curves

#include "libs/curveset.h"

/*!
\class CubicCurveSet curveset.h
\brief Piecewise cubic curves.

\ingroup KernelGroup
*/

/*!
\brief Creates an empty piecewise quadric curve.
*/
CubicCurveSet::CubicCurveSet()
{
}

/*!
\brief Creates the piecewise cubic curve given a set.
\param control Set of cubic curves.
*/
CubicCurveSet::CubicCurveSet(const QVector<CubicCurve>& control) :curve(control)
{
  lengths.resize(curve.size());

  length = 0.0;
  for (int i = 0; i < curve.size(); i++)
  {
    lengths[i] = curve[i].S(1.0, 100);
    length += lengths[i];
  }
}

/*!
\brief Creates the piecewise cubic curve given a set of control points with the Cattmul-Rom construction.
\param control Set of points.
\param ta Starting tangent vector
\param tb Ending tangent vector
*/
CubicCurveSet::CubicCurveSet(const QVector<Vector>& control, const Vector& ta, const Vector& tb)
{
  if (control.size() >= 3)
  {
    for (int i = 0; i < control.size() - 1; i++)
    {
      Vector t1;
      Vector t2;
      if (i == 0)
      {
        t1 = ta;
      }
      else
      {
        t1 = 0.5 * (control[i + 1] - control[i - 1]);
      }
      if (i == control.size() - 2)
      {
        t2 = tb;
      }
      else
      {
        t2 = 0.5 * (control[i + 2] - control[i]);
      }
      curve << CubicCurve::Hermite(control[i], control[i + 1], t1, t2);
    }
    lengths.clear();
    lengths.resize(curve.size());
    length = 0.0;
    for (int i = 0; i < curve.size(); i++)
    {
      lengths[i] = curve[i].S(1.0, 100);
      length += lengths[i];
    }
  }
}

/*!
\brief Compute the bounding box of the curve.
*/
Box CubicCurveSet::GetBox() const
{
  Box box = curve.at(0).GetBox();
  for (int i = 1; i < curve.size(); i++)
  {
    box = Box(box, curve.at(i).GetBox());
  }
  return box;
}

/*!
\brief Compute the Frenet vectors at a given point on the curve.

The columns of the returned matrix contain the tangent, normal and binormal vectors.

\param u Parameter.
*/
Matrix CubicCurveSet::GetMatrix(const double& u) const
{
  int iu = 0;
  if (u < 0.0)
  {
    // iu=0;
  }
  else if (u >= curve.size())
  {
    iu = curve.size() - 1;
  }
  else
  {
    iu = int(u);
  }
  return curve.at(iu).GetMatrix(u - iu);
}

/*!
\brief Compute the Frenet frame at a given point on the curve.

\param u Parameter.
*/
Frame CubicCurveSet::GetFrame(const double& u) const
{
  int iu = 0;
  if (u < 0.0)
  {
    // iu=0;
  }
  else if (u >= curve.size() - 1)
  {
    iu = curve.size() - 1;
  }
  else
  {
    iu = int(u);
  }
  return curve.at(iu).GetFrame(u - iu);
}

/*!
\brief Generates a discretization of the curve with a linear curvilign absisca parameterization.

\param s Stepping distance.
*/
QVector<Vector> CubicCurveSet::GetDiscretisation(const double& s) const
{
  // Total length
  double l = 0.0;

  // Curve identifier
  int k = 0;

  // Set of points
  QVector<Vector> p;

  // Start with the very first vertex of the piecewise curve
  p.append(curve[0](0));
  Vector last = curve[0](0);

  while (true)
  {
    l += s;
    double u;

    k = U(l, u);

    if (k != -1)
    {
      if (last != curve[k](u))
      {
        p.append(curve[k](u));
        last = curve[k](u);
      }
    }
    else
    {
      break;
    }
  }

  // Append last point of the curve
  if (last != curve[curve.size() - 1](1.0))
  {
    p.append(curve[curve.size() - 1](1.0));
  }
  return p;
}

/*!
\brief Generates a discretization of the curve with a linear curvilign absisca parameterization.

\param s Stepping distance.
\param tangents Tangents to the curve at sample points.
*/
QVector<Vector> CubicCurveSet::GetDiscretisation(const double& s, QVector<Vector>& tangents) const
{
  // Total length
  double l = 0.0;

  // Curve identifier
  int k = 0;

  // Set of points
  QVector<Vector> p;
  Vector last;

  // Start with the very first vertex of the piecewise curve
  p.append(curve[0](0));
  tangents.append(curve[0].Tangent(0.0001));
  last = curve[0](0);

  while (true)
  {
    l += s; // Position suivante sur la trajectoire

    double u;
    k = U(l, u);

    if (k != -1)
    {
      if (last != curve[k](u))
      {
        p.append(curve[k](u));  // Condition d'arret
        tangents.append(curve[k].Tangent(u));
        last = curve[k](u);
      }
    }
    else
    {
      break;
    }
  }

  // Append last point of the curve
  if (last != curve[curve.size() - 1](1.0))
  {
    p.append(curve[curve.size() - 1](1.0));  // Condition d'arret
    tangents.append(curve[curve.size() - 1].Tangent(0.9999));
  }
  return p;
}

/*!
\brief Compute the parameter of the curve corresponding to the input length.
\param s Input length.
\param u Parameter of the i-th curve.
\return Identifier of the curve in the piecewise definition.
*/
int CubicCurveSet::U(const double& s, double& u) const
{
  if ((s > length) || (s < 0.0))
  {
    u = 0.0;
    return -1;
  }

  // Cumulative length
  double l = 0.0;

  for (int i = 0; i < curve.size(); i++)
  {
    // We found the right interval
    if (s < l + lengths[i])
    {
      // Find the right location
      u = curve[i].U(s - l, 256);
      return i;
    }
    // Increment the length
    l += lengths[i];
  }

  // Extreme case : we reached the end of the curve
  u = 1.0;
  return curve.size() - 1;
}

/*!
\brief Computes the squared distance between a point and a set of cubic curves.

\sa CubicCurve::R()

\param p Point.
\param u Parameter defining the coordinate of the projection of the argument vertex onto the curve.
\param k Index of the curve for which the minimum distance was found.
*/
double CubicCurveSet::R(const Vector& p, double& u, int& k) const
{
  if (curve.size() == 0)
    return 0.0;

  k = 0;
  double a = curve.at(0).R(p, u);
  for (int i = 1; i < curve.size(); i++)
  {
    double v;
    double b = curve.at(i).R(p, v);

    if (b < a)
    {
      a = b;
      u = v;
      k = i;
    }
  }
  return a;
}

/*!
\brief Approximate a piecewise cubic curve by a piecewise quadric curve.

There are twice as many quadric curves as cubic curves.

\sa CubicCurve::Approximate
\param gamma Subdision parameter, in [0,1].
*/
QuadricCurveSet CubicCurveSet::Approximate(double gamma) const
{
  QVector <QuadricCurve> q;

  for (int i = 0; i < curve.size(); i++)
  {
    QuadricCurve a, b;
    curve.at(i).Approximate(a, b, gamma);
    q.append(a);
    q.append(b);
  }
  return QuadricCurveSet(q);
}

