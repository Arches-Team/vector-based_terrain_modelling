// Curves

#include "libs/curveset.h"

/*!
\class CubicCurve2Set curveset.h
\brief Piecewise cubic curves in the plane.

The length of every curve CubicCurve2Set::lengths and the total length CubicCurve2Set::length are computed
to speed up computations (internal optimization).
\ingroup PlanarGroup
*/

/*!
\brief Creates an empty piecewise quadric curve.
*/
CubicCurve2Set::CubicCurve2Set()
{
}

/*!
\brief Creates the piecewise quadric curve given a set of CubicCurve2.
\param control Set of cubic curves.
*/
CubicCurve2Set::CubicCurve2Set(const QVector<CubicCurve2>& control)
{
  curve.clear();
  curve = control;

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
\brief Creates the piecewise cubic curve given a set of control points with the Cattmul-Rom construction.
\param control Set of 2D points.
\param ta Starting tangent vector.
\param tb Ending tangent vector.
*/
CubicCurve2Set::CubicCurve2Set(const QVector<Vector2>& control, const Vector2& ta, const Vector2& tb)
{
  if (control.size() >= 3)
  {
    for (int i = 0; i < control.size() - 1; i++)
    {
      Vector2 t1;
      Vector2 t2;
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
      curve << CubicCurve2::Hermite(control[i], control[i + 1], t1, t2);
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
\brief Creates a piecewise cubic curve in the plane from a similar curve.
\param c Curve.
*/
CubicCurve2Set::CubicCurve2Set(const CubicCurveSet& c)
{
  curve.clear();

  for (int i = 0; i < c.Size(); i++)
  {
    curve.append(CubicCurve2(c(i)));
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

/*!
\brief Compute the bounding box of the curve.
*/
Box2 CubicCurve2Set::GetBox() const
{
  Box2 box = curve.at(0).GetBox();
  for (int i = 1; i < curve.size(); i++)
  {
    box = Box2(box, curve.at(i).GetBox());
  }
  return box;
}

/*!
\brief Generates a discretization of the curve with a linear curvilign absisca parameterization.

\param step Stepping distance.
\param tangents Tangents to the curve at sample points.
*/
QVector<Vector2> CubicCurve2Set::GetDiscretisation(const double& step, QVector<Vector2>& tangents) const
{
  double niveau = 0.0;
  int idCurve = 0;

  // Set of points
  QVector<Vector2> p;
  Vector2 last;

  // Start with the very first vertex of the piecewise curve
  p.append(curve[0](0));
  tangents.append(curve[0].Tangent(0.0001));
  last = curve[0](0);

  while (true)
  {
    niveau += step; // Position suivante sur la trajectoire

    double u;
    idCurve = U(niveau, u);

    if (idCurve != -1)
    {
      if (last != curve[idCurve](u))
      {
        p.append(curve[idCurve](u));  // Condition d'arret
        tangents.append(curve[idCurve].Tangent(u));
        last = curve[idCurve](u);
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
\brief Generates a discretization of the piecewise cubic curve with a linear curvilign absisca parameterization.

\param step Stepping distance.
*/
QVector<Vector2> CubicCurve2Set::GetDiscretisation(const double& step) const
{
  double niveau = 0.0;
  int idCurve = 0;

  // Set of points
  QVector<Vector2> p;
  Vector2 last;

  // Start with the very first vertex of the piecewise curve
  p.append(curve[0](0));
  last = curve[0](0);

  //int parcours = 0;                                 // Condition d'arret
  while (true)
  {
    niveau += step;                              // Position suivante sur la trajectoire

    double u;
    idCurve = U(niveau, u);

    if (idCurve != -1)
    {
      if (last != curve[idCurve](u))
      {
        p.append(curve[idCurve](u));  // Condition d'arret
        last = curve[idCurve](u);
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
  }
  return p;
}

/*!
\brief Compute the parameter of the curve corresponding to the input length.
\param s Input length.
\param u Parameter of the i-th curve.
\return Identifier of the curve in the piecewise definition.
*/
int CubicCurve2Set::U(const double& s, double& u) const
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
\brief Compute the squared distance between a point and a set of cubic curves.

\sa CubicCurve::R()

\param p Point.
\param u Parameter defining the coordinate of the projection of the argument vertex onto the curve.
\param k Index of the curve for which the minimum distance was found.
*/
double CubicCurve2Set::R(const Vector2& p, double& u, int& k) const
{
  if (curve.size() == 0)
    return 0.0;

  k = 0;
  double a = curve.at(0).R(p, u);
  for (int i = 1; i < curve.size(); i++)
  {
    double t;
    double b = curve.at(i).R(p, t);

    if (b < a)
    {
      a = b;
      u = t;
      k = i;
    }
  }
  return a;
}

/*!
\brief Compute the signed distance between a point and a set of cubic curves forming a closed path.

\sa CubicCurve2Set::R()
\param p Point.
*/
double CubicCurve2Set::Signed(const Vector2& p) const
{
  // Curve distance
  double u;
  int k;
  double d = sqrt(R(p, u, k));

  int nbIntersect = 0;

  // Check if p is inside with the number if intersections with the curve
  for (int i = 0; i < curve.size(); i++)
  {
    Ray2 ray = Ray2(p, Vector2(1.0, 0.0));
    Cubic c = ray.Direction()[0] * curve.at(i)[1] - ray.Direction()[1] * curve.at(i)[0];
    c[0] += ray.Direction()[1] * ray.Origin()[0] - ray.Direction()[0] * ray.Origin()[1];
    double t[3];
    nbIntersect += Intersect(ray, i, t);
  }

  bool inside = (nbIntersect % 2 != 0);

  if (inside) d = -d;

  return d;
}

/*!
\brief Compute the intersection between the set of curves and a ray.

\param ray The ray.
\param i Index.
\param t Intersection depths.
*/
int CubicCurve2Set::Intersect(const Ray2& ray, int i, double t[3]) const
{
  Cubic c = ray.Direction()[0] * curve.at(i)[1] - ray.Direction()[1] * curve.at(i)[0];
  c[0] += ray.Direction()[1] * ray.Origin()[0] - ray.Direction()[0] * ray.Origin()[1];

  double u[4];

  int n = c.Solve(u, 0.0, 1.0);

  if (n == 0)
  {
    return 0;
  }

  // Compute intersection depths
  if (fabs(ray.Direction()[0]) > fabs(ray.Direction()[1]))
  {
    for (int j = 0; j < n; j++)
    {
      t[j] = (curve.at(i)[0](u[j]) - ray.Origin()[0]) / ray.Direction()[0];
    }
  }
  else
  {
    for (int j = 0; j < n; j++)
    {
      t[j] = (curve.at(i)[1](u[j]) - ray.Origin()[1]) / ray.Direction()[1];
    }
  }

  if (n == 1)
  {
    if (t[0] < 0.0)
    {
      return 0;
    }
    else
    {
      return 1;
    }
  }

  if (n == 2)
  {
    if (t[0] < 0.0)
    {
      if (t[1] < 0.0)
      {
        return 0;
      }
      else
      {
        t[0] = t[1];
        return 1;
      }
    }
    else
    {
      if (t[1] < 0.0)
      {
        return 1;
      }
      else
      {
        return 2;
      }
    }
  }

  if (t[0] < 0.0)
  {
    if (t[1] < 0.0)
    {
      if (t[2] < 0.0)
      {
        return 0;
      }
      else
      {
        t[0] = t[2];
        return 1;
      }
    }
    else
    {
      if (t[2] < 0.0)
      {
        t[0] = t[1];
        return 1;
      }
      else
      {
        t[0] = t[1];
        t[1] = t[2];
        return 2;
      }
    }
  }
  else
  {
    if (t[1] < 0.0)
    {
      if (t[2] < 0.0)
      {
        return 1;
      }
      else
      {
        t[1] = t[2];
        return 2;
      }
    }
    else
    {
      if (t[2] < 0.0)
      {
        return 2;
      }
      else
      {
        return 3;
      }
    }
  }
}


/*!
\brief Compute the inverse mapping for a given input point.
\param p Point.
\param u, v Inverse Coordinates, u will denote the position along the curve and v the signed distance to the curve.
\param k Index of the curve for which the minimum distance was found.
*/
void CubicCurve2Set::UV(const Vector2& p, double& u, double& v, int& k) const
{
  if (curve.size() == 0)
    return;

  // Process first part of the curve
  k = 0;
  curve.at(0).UV(p, u, v);

  // Process other parts
  for (int i = 1; i < curve.size(); i++)
  {
    double tu, tv;
    curve.at(i).UV(p, tu, tv);

    if (tv < v)
    {
      v = tv;
      u = tu;
      k = i;
    }
  }
}

/*!
\brief Draw the curve.

The brush is not used.
\param scene Graphics scene.
\param pen The pen.
\param brush The brush.
*/
void CubicCurve2Set::Draw(QGraphicsScene& scene, const QPen& pen, const QBrush& brush) const
{
  QPen joining = pen;
  joining.setCapStyle(Qt::RoundCap); // Could have been FlatCap, however this yields discontinuities
  joining.setJoinStyle(Qt::RoundJoin);

  for (int i = 0; i < curve.size(); i++)
  {
    curve.at(i).Draw(scene, joining, brush);
  }
}

/*!
\brief Check if a point lies inside the piecewise cubic curve closed contour.

The curve should be closed.

\param p Point.
*/
bool CubicCurve2Set::Inside(const Vector2& p) const
{
  Ray2 ray(p, Vector2::X);
  int c = 0;

  for (int i = 0; i < curve.size(); i++)
  {
    double u[3];
    c += curve.at(i).Intersect(ray, u);
  }

  return (c % 2 == 0) ? false : true;
}

/*!
\brief Compute the sinuosity.

Sinusoity is defined as the ratio between the length of the curve and the distance between end points.
*/
double CubicCurve2Set::Sinuosity() const
{
  return length / Norm(curve.first()(0.0) - curve.last()(1.0));
}


/*!
\brief Approximate a cubic curve by two quadric curves.

N. Truong, C. Yuksel and L. Seiler. Quadratic Approximation of Cubic Curves. <I> Proceedings of the ACM on Computer Graphics and Interactive Techniques.</I> 2020.
\param gamma Subdision parameter, in [0,1].
\param a,b Returned quadric curves.
*/
void CubicCurve2::Approximate(QuadricCurve2& a, QuadricCurve2& b, const double& gamma) const
{
  Vector2 c[4];
  for (int i = 0; i < 4; i++)
  {
    c[i] = BezierControl(i);
  }

  // Control points
  Vector2 q[5];

  q[0] = c[0];
  q[4] = c[3];
  q[1] = c[0] + 1.5 * gamma * (c[1] - c[0]);
  q[3] = c[3] - 1.5 * (1.0 - gamma) * (c[3] - c[2]);
  q[2] = (1.0 - gamma) * q[1] + gamma * q[3];

  a = QuadricCurve2::Bezier(q[0], q[1], q[2]);
  b = QuadricCurve2::Bezier(q[2], q[3], q[4]);
}
