// Curves

#include "libs/curveset.h"

/*!
\class QuadricCurve2Set curveset.h
\brief Piecewise quadric curves.

This class stores an array of quadric curves.
\ingroup PlanarGroup
*/

/*!
\brief Creates an empty piecewise quadric curve.
*/
QuadricCurve2Set::QuadricCurve2Set()
{
}

/*!
\brief Creates the piecewise quadric curve given a set of control points.
\param p Control vertices.
*/
QuadricCurve2Set::QuadricCurve2Set(const QVector<Vector2>& p)
{
  int n = p.size();

  if (n > 2)
  {
    QVector<Vector2> knot;
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
      curve[i] = QuadricCurve2::Bezier(knot[2 * i], knot[2 * i + 1], knot[2 * i + 2]);
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
  else if (n == 1)
  {
    curve.append(QuadricCurve2(Quadric(p[0][0], 0.0, 0.0), Quadric(p[0][1], 0.0, 0.0)));

  }
  else if (n == 2)
  {
    Vector2 d = p[1] - p[0];
    curve.append(QuadricCurve2(Quadric(p[0][0], d[0], 0.0), Quadric(p[0][1], d[1], 0.0)));
  }
}

/*!
\brief Create a set of quadric curves.

Could use QuadricCurve2Set::LengthIntegral() to compute the length of every part.
\param s Set of quadric curves.
*/
QuadricCurve2Set::QuadricCurve2Set(const QVector<QuadricCurve2>& s)
{
  curve = s;
  lengths.clear();
  lengths.resize(curve.size());
  length = 0.0;
  for (int i = 0; i < curve.size(); i++)
  {
    lengths[i] = curve[i].S(1.0, 100);
    //lengths[i] = curve[i].LengthIntegral();
    length += lengths[i];
  }
}

/*!
\brief Creates a piecewise cubic curve in the plane from a similar curve.
\param c Curve.
*/
QuadricCurve2Set::QuadricCurve2Set(const QuadricCurveSet& c)
{
  curve.clear();

  for (int i = 0; i < c.Size(); i++)
  {
    curve.append(QuadricCurve2(c(i)));
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
\brief Destroys a piecewise quadric curve.
*/
QuadricCurve2Set::~QuadricCurve2Set()
{
}

/*!
\brief Compute the bounding box of the curve.
*/
Box2 QuadricCurve2Set::GetBox() const
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
\param step Step.
*/
QVector<Vector2> QuadricCurve2Set::GetDiscretisation(const double& step) const
{
  double niveau = 0.0; // Position sur la trajectoire 
  int idCurve = 0; // Identifiant de la courbe courante

  // Set of points
  QVector<Vector2> p;
  Vector2 last;

  // Start with the very first vertex of the piecewise curve
  p.append(curve[0](0));
  last = curve[0](0);

  while (true)
  {
    double u;
    niveau += step;                              // Position suivante sur la trajectoire

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
    p.append(curve[curve.size() - 1](1.0));  // Condition d'arret
  return p;
}

/*!
\brief Compute the parameter of the curve corresponding to the input length.
\param s Input length.
\param u Parameter of the i-th curve.
\return Identifier of the curve in the piecewise definition.
*/
int QuadricCurve2Set::U(const double& s, double& u) const
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
QuadricCurve2Set& QuadricCurve2Set::operator+=(const QuadricCurve2Set& c)
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
double QuadricCurve2Set::R(const Vector2& p, double& u, int& k) const
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
\brief Compute the inverse mapping for a given input point.
\param p Point.
\param u, v Inverse Coordinates, u will denote the position along the curve and v the signed distance to the curve.
\param k Index defining the sub-curve element.
*/
double QuadricCurve2Set::UV(const Vector2& p, double& u, double& v, int& k) const
{
  if (curve.size() == 0)
    return 0.0;

  // Process first part of the curve
  k = 0;
  double r2 = curve.at(0).UV(p, u, v);

  // Process other parts
  for (int i = 1; i < curve.size(); i++)
  {
    double tu, tv;
    double tr2 = curve.at(i).UV(p, tu, tv);

    if (tr2 < r2)
    {
      v = tv;
      u = tu;
      k = i;
      r2 = tr2;
    }
  }
  return r2;
}

/*!
\brief Translate the curve.
\param t Translation vector.
*/
void QuadricCurve2Set::Translate(const Vector2& t)
{
  for (int i = 0; i < curve.size(); i++)
  {
    curve[i].Translate(t);
  }
}

/*!
\brief Draw the curve.

The brush is not used.
\param scene Graphics scene.
\param pen The pen.
\param brush The brush.
*/
void QuadricCurve2Set::Draw(QGraphicsScene& scene, const QPen& pen, const QBrush& brush) const
{
  for (int i = 0; i < curve.size(); i++)
  {
    curve.at(i).Draw(scene, pen, brush);
  }
}

/*!
\brief Check if a point lies inside the piecewise quadric curve closed contour.

The curve should be closed.

\sa Intersections

\param p Point.
*/
bool QuadricCurve2Set::Inside(const Vector2& p) const
{
  int c = Intersections(Ray2(p, Vector2::X));

  return (c % 2 == 0) ? false : true;
}

/*!
\brief Compute the first intersection between a ray and the curve.

\param ray The ray.
\param t Returned intersection depth.
*/
bool QuadricCurve2Set::Intersect(const Ray2& ray, double& t) const
{
  t = Math::Infinity;
  bool r = false;
  for (int i = 0; i < curve.size(); i++)
  {
    double u[2];

    // Intersection
    int c = curve.at(i).Intersect(ray, u);

    for (int j = 0; j < c; j++)
    {
      if ((u[j] < t) && (u[j] >= 0.0))
      {
        r = true;
        t = u[j];
      }
    }
  }

  return r;
}

/*!
\brief Compute the number of intersections between a ray and a piecewize quadric curve.

\param ray The ray.
*/
int QuadricCurve2Set::Intersections(const Ray2& ray) const
{
  int c = 0;

  for (int i = 0; i < curve.size(); i++)
  {
    double u[2];
    c += curve.at(i).Intersect(ray, u);
  }
  return c;
}

/*!
\brief Compute the sinuosity.

Sinusoity is defined as the ratio between the length of the curve and the distance between end points.
*/
double QuadricCurve2Set::Sinuosity() const
{
  return length / Norm(curve.first()(0.0) - curve.last()(1.0));
}