// Set of quadric curves

#include "libs/curveset.h"

/*!
\class QuadricCurveSet curveset.h
\brief Piecewise quadric curves.

The array of quadric curves is stored using a QVector container so as to simplify memory management.
\ingroup KernelGroup
*/

/*!
\brief Creates an empty piecewise quadric curve.
*/
QuadricCurveSet::QuadricCurveSet()
{
}

/*!
\brief Compute the bounding box of the curve.
*/
Box QuadricCurveSet::GetBox() const
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
Matrix QuadricCurveSet::GetMatrix(const double& u) const
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
Frame QuadricCurveSet::GetFrame(const double& u) const
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
\brief Creates the piecewise quadric curve given a set of control points.
\param p Control vertices.
*/
QuadricCurveSet::QuadricCurveSet(const QVector<Vector>& p)
{
  int n = p.size();

  if (n > 2)
  {
    QVector<Vector> knot;
    knot.clear();

    // Starting vertices
    knot.append(p[0]);
    knot.append(p[1]);

    // Mid vertices
    for (int i = 1; i < n - 2; i++)
    {
      knot.append(0.5 * (p[i] + p[i + 1]));
      knot.append(p[i + 1]);
    }

    // End vertex
    knot.append(p[n - 1]);

    curve.clear();
    curve.resize(n - 2);

    for (int i = 0; i < n - 2; i++)
    {
      curve[i] = QuadricCurve::Bezier(knot[2 * i], knot[2 * i + 1], knot[2 * i + 2]);
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
\brief Create a set of quadric curves.
\param s Set of quadric curves.
*/
QuadricCurveSet::QuadricCurveSet(const QVector<QuadricCurve>& s)
{
  curve = s;
  lengths.clear();
  lengths.resize(curve.size());
  length = 0.0;
  for (int i = 0; i < curve.size(); i++)
  {
    lengths[i] = curve[i].S(1.0, 100);
    length += lengths[i];
  }
}

/*!
\brief Generates a discretization of the curve with a linear curvilign absisca parameterization.
\param step Stepping distance.
*/
QVector<Vector> QuadricCurveSet::GetDiscretisation(const double& step) const
{
  double s = 0.0; // Curvilign abscisca
  int id = 0; // Identifiant de la courbe courante
  double u; // position dans une sous courbe

  // Set of points
  QVector<Vector> p;
  Vector last;

  // Start with the very first vertex of the piecewise curve
  p.append(curve[0](0));
  last = curve[0](0);

  while (true)
  {
    s += step;                              // Position suivante sur la trajectoire

    id = U(s, u);

    if (id != -1)
    {
      if (last != curve[id](u))
      {
        p.append(curve[id](u));  // Condition d'arret
        last = curve[id](u);
      }
    }
    else
    {
      break;
    }
  }

  // Append last point of the curve
  if (last != curve[curve.size() - 1](1.0))
    p.append(curve[curve.size() - 1](1.0));  // Condition d'arret
  return p;
}

/*!
\brief Compute the parameter of the curve corresponding to the input length.
\param s Input length.
\param u Parameter of the i-th curve.
\return Identifier of the curve in the piecewise definition.
*/
int QuadricCurveSet::U(const double& s, double& u) const
{
  if ((s > length) || (s < 0.0))
  {
    u = 0.0;
    return -1;
  }

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

  // Extreme case : we reached the end of the curve (in fact we should almost never get there)
  u = 1.0;
  return curve.size() - 1;
}

/*!
\brief Appends a piecewise cubic curve to an existing piecewise cubic curve.
\param c Argument piecewise cubic curve.
*/
QuadricCurveSet& QuadricCurveSet::operator+=(const QuadricCurveSet& c)
{
  curve += c.curve;
  lengths += c.lengths;
  length += c.length;

  return *this;
}

/*!
\brief This function returns the squared distance between
a point and a set of quadric curves.

\sa QuadicCurve::R()

\param p Point.
\param u Parameter defining the coordinate of the projection of the argument vertex onto the curve.
\param k Index of the curve for which the minimum distance was found.
*/
double QuadricCurveSet::R(const Vector& p, double& u, int& k) const
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
\brief Compute the squared distance between a point and a set of quadric curves.

\sa QuadricCurveSet::R(const Vector&,double&,int&)

\param p Point.
*/
double QuadricCurveSet::R(const Vector& p) const
{
  double u;
  int k;
  return R(p, u, k);
}

/*!
\brief Compute the distance between a point and a set of quadric curves.

\sa QuadricCurveSet::R(const Vector&,double&,int&)

\param p Point.
*/
double QuadricCurveSet::Signed(const Vector& p) const
{
  double u;
  int k;
  return sqrt(R(p, u, k));
}
