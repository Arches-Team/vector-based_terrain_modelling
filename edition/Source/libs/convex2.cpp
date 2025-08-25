// Convex

#include "libs/convex.h"

/*!
\class Convex2 convex.h
\brief %Convex polygons in the plane.

Drawing member function is inherited from Polygon2::Draw() const.
\ingroup PlanarGroup
*/

/*!
\brief Create a convex polygon.

Note that convexity is not checked. See Convex2::Hull to
compute the convex hull of a set of points.

\param q Array of points defining the polygon (should be convex).
*/
Convex2::Convex2(const QVector<Vector>& q) :Polygon2(q)
{
}

/*!
\brief Create a convex polygon.

\param q Array of points defining the polygon (should be convex).
*/
Convex2::Convex2(const QVector<Vector2>& q) :Polygon2(q)
{
}

/*!
\brief Create a convex polygon from an array of points.

See Convex2::Hull to compute the convex hull of a set of points.

\param p Array of points.
\param n Size.
*/
Convex2::Convex2(Vector2* p, int n) :Polygon2(p, n)
{
}

/*!
\brief Create a convex polygon from three points.

\param a,b,c Points (should be in trigonometric ordrer).
*/
Convex2::Convex2(const Vector2& a, const Vector2& b, const Vector2& c) :Polygon2(Triangle2(a, b, c))
{
}

/*!
\brief Create a convex from a box.
\param b The box.
*/
Convex2::Convex2(const Box2& b) :Polygon2(b)
{
}

/*!
\brief Create a convex from a hexagon.
\param h The hexagon.
*/
Convex2::Convex2(const Hexagon2& h) :Polygon2(h)
{
}

/*!
\brief Cut the convex polygon by a line.

This function updates the convex polygon and returns true if intersection occured.
\param line The line.
*/
bool Convex2::Cut(const Line2& line)
{
  // LEFT OR MIDDLE = true / RIGHT = false
  unsigned int nb_true = 0;
  unsigned int nb_false = 0;
  QVector<bool> status;
  QVector<Vector2> q2;
  for (int i = 0; i < q.size(); i++)
  {
    Vector2 a = q[i];
    if (line.IsLeftOrOn(a))
    {
      nb_true++;
      status.append(true);
    }
    else
    {
      nb_false++;
      status.append(false);
    }
  }

  if (nb_true == 0)
  {
    q = q2;
    return false;
  }

  else if (nb_false == 0)
  {
    return true;
  }

  for (int i = 0; i < q.size(); i++)
  {
    Vector2 a = q[i];
    Vector2 b = q[(i + 1) % q.size()];

    if (status[i] == true && status[(i + 1) % q.size()] == true)
      q2.append(a);
    else if (status[i] == true && status[(i + 1) % q.size()] == false)
    {
      q2.append(a);
      Vector2 pp;
      line.Intersection(Segment2(a, b), pp);
      q2.append(pp);
    }
    else if (status[i] == false && status[(i + 1) % q.size()] == true)
    {
      Vector2 pp;
      line.Intersection(Segment2(a, b), pp);
      q2.append(pp);
    }
    // else (false -> false) : nothing
  }
  q = q2;

  return true;
}

/*!
\brief Check the intersection between the convex and a segment, if intersection occurs, add the points to the convex.

This function adds the points to the list of vertices of the convex polygon.
\param s The segment.
*/
bool Convex2::AddIntersection(const Segment2& s)
{
  bool e = false;

  for (int i = 0; i < q.size(); i++)
  {
    Vector2 a = q[i];
    Vector2 b = q[(i + 1) % q.size()];

    Vector2 pi;
    if (s.Intersection(Segment2(a, b), pi))
    {
      if (!(NormInfinity(pi - a) < 0.00001) && !(NormInfinity(pi - b) < 0.00001))
      {
        q.insert(q.begin() + i + 1, pi);
        e = true;
      }
    }
  }

  return e;
}

/*!
\brief Create one Voronoi cell.
\param q Set of points.
\param p Point.
\param box %Box embedding the set of points.
*/
Convex2 Convex2::VoronoiCell(const QVector<Vector2>& q, const Box2& box, const Vector2& p)
{
  // Initialize convex hull
  Convex2 c = Convex2(box);

  for (int i = 0; i < q.size(); i++)
  {
    const Vector2& pi = q.at(i);
    if (Vector2::Equal(p, pi, 0.00001))
      continue;

    Vector2 m = (pi + p) / 2.0;
    Vector2 d = Vector2(pi - p).Orthogonal();

    c.Cut(Line2(m, m + d));
  }
  return c;
}

/*!
\brief Create the Voronoi cells from a given set of points.
\param q Set of points.
\param box %Box embedding the set of points.
*/
QVector<Convex2> Convex2::VoronoiCells(const QVector<Vector2>& q, const Box2& box)
{
  // Cells
  QVector<Convex2> c;
  c.reserve(q.size());

  for (int i = 0; i < q.size(); i++)
  {
    c.append(VoronoiCell(q, box, q[i]));
  }
  return c;
}

/*!
\brief Test if a convex is inside another one.

\param convex The convex polygon.
*/
bool Convex2::Inside(const Convex2& convex) const
{
  for (int i = 0; i < convex.Size(); i++)
  {
    if (!Polygon2::Inside(convex.Vertex(i)))
    {
      return false;
    }
  }

  return true;
}
/*!
\brief Create the convex hull of a set of points.
\param p Set of points.
*/
Convex2 Convex2::Hull(QVector<Vector2> p)
{
  // Escape case
  if (p.size() == 0)
    return Convex2();

  int n = p.size();
  int k = 0;
  QVector<Vector2> hull(2 * n);

  // Sort points lexicographically
  std::sort(p.begin(), p.end(), [](const Vector2& a, const Vector2& b) { return a[0] < b[0];  });

  // Build lower hull
  for (int i = 0; i < n; i++)
  {
    // while (k >= 2 && !Clockwise(hull[k - 2], hull[k - 1], p[i])) k--;
    while (k >= 2 && ((hull[k - 1] - hull[k - 2]) / (p[i] - hull[k - 2]) <= 0.0))
    {
      k--;
    }
    hull[k++] = p[i];
  }

  // Build upper hull
  for (int i = n - 2, t = k + 1; i >= 0; i--)
  {
    // while (k >= t && !Clockwise(hull[k - 2], hull[k - 1], p[i])) 
    while (k >= t && ((hull[k - 1] - hull[k - 2]) / (p[i] - hull[k - 2]) <= 0.0))
    {
      k--;
    }
    hull[k++] = p[i];
  }

  hull.resize(k);

  return Convex2(hull);
}

/*!
\brief Compute the squared distance between a point and a convex polygon.
\param p Point.
*/
double Convex2::R(const Vector2& p) const
{
  double r = Segment2(q.at(q.size() - 1), q.at(0)).R(p);

  for (int i = 0; i < q.size() - 1; i++)
  {
    const Vector2& a = q.at(i);
    const Vector2& b = q.at(i + 1);

    double t = Segment2(a, b).R(p);
    if (t < r)
    {
      r = t;
    }
  }
  return r;
}

/*!
\brief Compute the normal vector between a point and a convex polygon.

Returns null vector if the point is inside the polygon.

\param p Point.
*/
Vector2 Convex2::Normal(const Vector2& p) const
{
  Vector2 a = q.at(q.size() - 1);
  Vector2 pa = p - a;
  Vector2 ea = q.at(q.size() - 1) - q.at(q.size() - 2);
  double paea = pa * ea;

  for (int i = 0; i < q.size(); i++)
  {
    Vector2 b = q.at(i);
    Vector2 pb = p - b;
    Vector2 eb = b - a;

    // Orthogonal vector points inward
    Vector2 ebo = eb.Orthogonal();

    double d = ebo * pa;

    // Not on the right side
    if (d > 0.0)
      continue;

    double paeb = pa * eb;
    double pbeb = pb * eb;
    //
    if (paeb < 0.0)
    {
      // Vertex
      if (paea > 0.0)
      {
        return pa;
      }
    }
    else
    {
      // Edge
      if (pbeb < 0.0)
      {
        double r = pa * pa - (paeb * paeb) / (eb * eb);
        return ebo * r / Norm(ebo);
      }
    }

    ea = eb;
    pa = pb;
    a = b;
    paea = pbeb;
  }

  // Face
  return Vector2::Null;
}

/*!
\brief Straightforward implementation of the Minkowski sum of two convex.
\param c %Convex.
*/
Convex2 Convex2::Minkowski(const Convex2& c) const
{
  QVector<Vector2> v;
  for (int i = 0; i < q.size(); i++)
  {
    for (int j = 0; j < c.q.size(); j++)
    {
      v.append(q.at(i) + c.q.at(j));
    }
  }
  return Hull(v);
}
