// Octogon

#include <QtWidgets/QGraphicsScene>

#include "libs/octogon.h"

/*!
\class Octogon2 octogon.h
\brief Regular Octogons (in the plane).

\ingroup PlanarGroup
*/

const Vector2 Octogon2::vertex[8] = {
  Vector2(1.0+Math::Sqrt2, 1.0) / std::sqrt(4+2*Math::Sqrt2),
  Vector2(1.0,1.0+Math::Sqrt2) / std::sqrt(4+2*Math::Sqrt2),
  Vector2(-1.0,1.0+Math::Sqrt2) / std::sqrt(4+2*Math::Sqrt2),
  Vector2(-1.0-Math::Sqrt2, 1.0) / std::sqrt(4+2*Math::Sqrt2),
  Vector2(-1.0-Math::Sqrt2, -1.0) / std::sqrt(4+2*Math::Sqrt2),
  Vector2(-1.0,-1.0-Math::Sqrt2) / std::sqrt(4+2*Math::Sqrt2),
  Vector2(+1.0,-1.0-Math::Sqrt2) / std::sqrt(4+2*Math::Sqrt2),
  Vector2(1.0+Math::Sqrt2, -1.0) / std::sqrt(4+2*Math::Sqrt2),
};

const Vector2 Octogon2::edge[8] = {
  vertex[1] - vertex[0],
  vertex[2] - vertex[1],
  vertex[3] - vertex[2],
  vertex[4] - vertex[3],
  vertex[5] - vertex[4],
  vertex[6] - vertex[5],
  vertex[7] - vertex[6],
  vertex[0] - vertex[7]
};

const Vector2 Octogon2::normal[8] = {
  Normalized(edge[0].Orthogonal()),
  Normalized(edge[1].Orthogonal()),
  Normalized(edge[2].Orthogonal()),
  Normalized(edge[3].Orthogonal()),
  Normalized(edge[4].Orthogonal()),
  Normalized(edge[5].Orthogonal()),
  Normalized(edge[6].Orthogonal()),
  Normalized(edge[7].Orthogonal())
};

/*!
\brief Create an octogon.
\param c Center.
\param r Radius.
*/
Octogon2::Octogon2(const Vector2& c, const double& r) :c(c), r(r)
{
}

/*!
\brief Create an octogon.
\param r Radius.
*/
Octogon2::Octogon2(const double& r) :Octogon2(Vector2::Null, r)
{
}

/*!
\brief Area of the octogon.
*/
double Octogon2::Area() const
{
  return 2.0 * sqrt(2.0) * r * r;
}

/*!
\brief Perimeter of the octogon.
*/
double Octogon2::Perimeter() const
{
  return 8.0 * Edge();
}

/*!
\brief Translate an octogon.
\param t Translation vector.
*/
void Octogon2::Translate(const Vector2& t)
{
  c += t;
}

/*!
\brief Scale an octogon.
\param s Scaling factor.
*/
void Octogon2::Scale(const double& s)
{
  c *= s;
  r *= s;
}

/*!
\brief Computes the bounding box of an octogon.
*/
Box2 Octogon2::GetBox() const
{
  return Box2(c, r);
}

/*!
\brief Overloaded.
\param s Stream.
\param octogon The octogon.
*/
std::ostream& operator<<(std::ostream& s, const Octogon2& octogon)
{
  s << "Octogon2(" << octogon.c << ',' << octogon.r << ")";
  return s;
}

/*!
\brief Draws an octogon.
\param scene Graphics scene.
\param pen The pen.
\param brush The brush.
*/
void Octogon2::Draw(QGraphicsScene& scene, const QPen& pen, const QBrush& brush) const
{
  QPolygonF polygon;
  for (int i = 0; i < 8; i++)
  {
    Vector2 p = Vertex(i);
    polygon.append(QPointF(p[0], p[1]));
  }
  scene.addPolygon(polygon, pen, brush);
}

/*!
\brief Test if a point is inside.
\param p Point.
*/
bool Octogon2::Inside(const Vector2& p) const
{
  // Take advantage of vertical and horizontal symmetry 
  Vector2 cp = Abs(p - c) / r;

  // Diagonal symmetry
  if (cp[1] > cp[0])
  {
    Math::Swap(cp[1], cp[0]);
  }

  // Plane test
  if ((cp - vertex[0]) * normal[0] < 0.0) return false;
  if ((cp - vertex[7]) * normal[7] < 0.0) return false;

  return true;
}

/*!
\brief Compute the squared Euclidean distance to the octogon.
\param p Point.
*/
double Octogon2::R(const Vector2& p) const
{
  // Take advantage of vertical and horizontal symmetry 
  Vector2 q = Abs(p - c);

  // Diagonal symmetry
  if (q[1] > q[0])
  {
    Math::Swap(q[1], q[0]);
  }

  double d = q * normal[0] - r;

  // Inside
  if (d < 0.0) return 0.0;

  Vector2 a = vertex[0] * r;
  Vector2 aq = q - a;

  Vector2 e = edge[0];
  double l = aq * e;
  if (l < 0.0)
  {
    return SquaredNorm(aq);
  }
  else if (l > r)
  {
    Vector2 b = vertex[1] * r;
    Vector2 bq = q - b;
    return SquaredNorm(bq);
  }
  else
  {
    return d * d;
  }
}

/*!
\brief Compute the signed Euclidean distance to the octogon.
\param p Point.
*/
double Octogon2::Signed(const Vector2& p) const
{
  // Take advantage of vertical and horizontal symmetry 
  Vector2 q = Abs(p - c);

  // Diagonal symmetry
  if (q[1] > q[0])
  {
    Math::Swap(q[1], q[0]);
  }

  double d = q * normal[0] - r;

  // Inside
  if (d < 0.0) return d;

  Vector2 a = vertex[0] * r;
  Vector2 aq = q - a;


  Vector2 e = edge[0];
  double l = aq * e;
  if (l < 0.0)
  {
    return Norm(aq);
  }
  else if (l > r)
  {
    Vector2 b = vertex[1] * r;
    Vector2 bq = q - b;
    return Norm(bq);
  }
  else
  {
    return d;
  }
}

/*!
\brief Computes the distance vector between an octogon and a point.
\param p Point.
*/
Vector2 Octogon2::Normal(const Vector2& p) const
{
  // Take advantage of vertical and horizontal symmetry 
  Vector2 q = p - c;
  bool sx = (q[0] < 0.0) ? true : false;
  bool sy = (q[1] < 0.0) ? true : false;

  q = Abs(p - c);

  // Diagonal symmetry
  bool sxy = false;
  if (q[1] > q[0])
  {
    Math::Swap(q[1], q[0]);
    sxy = true;
  }

  double d = q * normal[0] - r;

  // Inside
  if (d < 0.0) return Vector2::Null;

  Vector2 a = vertex[0] * r;
  Vector2 aq = q - a;

  Vector2 n;

  Vector2 e = edge[0];
  double l = aq * e;
  if (l < 0.0)
  {
    n = aq;
  }
  else if (l > r)
  {
    Vector2 b = vertex[1] * r;
    Vector2 bq = q - b;
    n = bq;
  }
  else
  {
    n = d * normal[0];
  }

  // Symmetries, performed in inverse order
  if (sxy == true)
  {
    Math::Swap(n[1], n[0]);
  }
  if (sy == true)
  {
    n[0] = -n[1];
  }
  if (sx == true)
  {
    n[0] = -n[0];
  }

  return n;

}

/*!
\brief Generate a random vector inside the box.
\param random %Random number generator.
*/
Vector2 Octogon2::RandomInside(Random& random) const
{
  Vector2 p;
  while (true)
  {
    p = c + r * Vector2(random.Uniform(-1.0, 1.0), random.Uniform(-1.0, 1.0));
    if (Inside(p)) break;
  }
  return p;
}

/*!
\brief Compute a Poisson disc distribution inside the irregular octogon.

This function uses a simple dart throwing algorithm.

\param r Radius of the discs.
\param n Number of candidate points.
\param s Set of points.
\param a Add set of points flag: set to true it the set points should be added to the sampling, set to false to define constraints (typically for borders).
\param random %Random number generator.
*/
QVector<Vector2> Octogon2::Poisson(const double& r, int n, const QVector<Vector2>& s, bool a, Random& random) const
{
  QVector<Vector2> p = s;

  // Collision radius
  double c = 4.0 * r * r;

  // Create instances
  for (int i = 0; i < n; i++)
  {
    Vector2 t = RandomInside(random);
    bool hit = false;
    for (int j = 0; j < p.size(); j++)
    {
      if (SquaredNorm(t - p.at(j)) < c)
      {
        hit = true;
        break;
      }
    }
    if (hit == false)
      p.append(t);
  }

  // In this case, we need to keep only the generated samples
  if (a == false)
  {
    p = p.mid(s.size());
  }
  return p;
}

/*!
\brief Compute a Poisson disc distribution inside the irregular octogon.

This function uses a simple dart throwing algorithm.

\param r Radius of the discs.
\param n Number of candidate points.
\param random %Random number generator.
*/
QVector<Vector2> Octogon2::Poisson(const double& r, int n, Random& random) const
{
  QVector<Vector2> p;

  // Collision radius
  double c = 4.0 * r * r;

  // Create instances
  for (int i = 0; i < n; i++)
  {
    Vector2 t = RandomInside(random);
    bool hit = false;
    for (int j = 0; j < p.size(); j++)
    {
      if (SquaredNorm(t - p.at(j)) < c)
      {
        hit = true;
        break;
      }
    }
    if (hit == false)
      p.append(t);
  }

  return p;
}
