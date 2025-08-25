// Hashing 

#include <QtWidgets/QGraphicsScene>

#include "libs/hexagon.h"
#include "libs/circle.h"

/*!
\class Hexagon2 hexagon.h
\brief %Flat topped hexagon in the plane.

\ingroup PlanarGroup
*/

const double Hexagon2::Alpha = sqrt(3.0) / 2.0;

const Vector2 Hexagon2::vertex[6] = {
  Vector2(1.0, 0.0),
  Vector2(0.5, Alpha),
  Vector2(-0.5, Alpha),
  Vector2(-1.0, 0.0),
  Vector2(-0.5,-Alpha),
  Vector2(0.5,-Alpha) };

const Vector2 Hexagon2::normal[6] = {
  Vector2(Alpha, 0.5),
  Vector2(0.0, 1.0),
  Vector2(-Alpha, 0.5),
  Vector2(-Alpha, -0.5),
  Vector2(0.0,-1.0),
  Vector2(Alpha, -0.5) };

const Vector2 Hexagon2::edge[6] = {
  Vector2(-0.5, Alpha),
  Vector2(-1.0, 0.0),
  Vector2(-0.5, -Alpha),
  Vector2(0.5, -Alpha),
  Vector2(1.0, 0.0),
  Vector2(0.5, Alpha) };

/*!
\brief Create an hexagon.
\param c Center.
\param r Radius.
*/
Hexagon2::Hexagon2(const Vector2& c, const double& r) :c(c), r(r)
{
}

/*!
\brief Create an hexagon.
\param r Radius.
*/
Hexagon2::Hexagon2(const double& r) :Hexagon2(Vector2::Null, r)
{
}

/*!
\brief Area of the hexagon.
*/
double Hexagon2::Area() const
{
  return 3.0 * Alpha * r * r;
}

/*!
\brief Perimeter of the hexagon.
*/
double Hexagon2::Perimeter() const
{
  return 6.0 * r;
}

/*!
\brief Translate a hexagon.
\param t Translation vector.
*/
void Hexagon2::Translate(const Vector2& t)
{
  c += t;
}

/*!
\brief Scale a hexagon.
\param s Scaling factor.
*/
void Hexagon2::Scale(const double& s)
{
  c *= s;
  r *= s;
}

/*!
\brief Compute the bounding box of a hexagon.
*/
Box2 Hexagon2::GetBox() const
{
  return Box2(Vector2(-r, -Alpha * r), Vector2(r, Alpha * r)).Translated(c);
}

/*!
\brief Overloaded.
\param s Stream.
\param hex The hexagon.
*/
std::ostream& operator<<(std::ostream& s, const Hexagon2& hex)
{
  s << "Hexagon2(" << hex.c << ',' << hex.r << ")";
  return s;
}

/*!
\brief Draws an hexagon.
\param scene Graphics scene.
\param pen The pen.
\param brush The brush.
*/
void Hexagon2::Draw(QGraphicsScene& scene, const QPen& pen, const QBrush& brush) const
{
  QPolygonF polygon;
  for (int i = 0; i < 6; i++)
  {
    Vector2 p = Vertex(i);
    polygon.append(QPointF(p[0], p[1]));
  }
  scene.addPolygon(polygon, pen, brush);
}

/*!
\brief Compute the sector given an input direction.
\param d Direction.
*/
int Hexagon2::Sector(const Vector2& d)
{
  bool a = d[1] >= 0.0; // d * Vector2(0.0,1.0);
  bool b = d * Vector2(-Math::Sqrt3 / 2.0, 0.5) >= 0.0;
  bool c = d * Vector2(-Math::Sqrt3 / 2.0, -0.5) >= 0.0;

  // Line of a = 0Â°
  if (a)
  {
    // Line a = 60Â°
    if (b)
    {
      if (c)
      {
        return 2;
      }
      else
      {
        return 1;
      }
    }
    else
    {
      return 0;
    }
  }
  else
  {
    if (b)
    {
      return 3;
    }
    else
    {
      // Third sector
      if (c)
      {
        return 4;
      }
      // Second
      else
      {
        return 5;
      }
    }
  }
}

/*!
\brief Test if a point is inside.
\param p Point.
*/
bool Hexagon2::Inside(const Vector2& p) const
{
  Vector2 cp = p - c;
  // Take advantage of horizontal symmetry and test only three slabs
  cp[1] = fabs(cp[1]);

  for (int i = 0; i < 3; i++)
  {
    double d = cp * normal[i] - Alpha*r;
    if (d > 0.0) return false;
  }
  return true;
}

/*!
\brief Check if two hexagons intersect.
\param hexagon %Hexagon.
*/
bool Hexagon2::Intersect(const Hexagon2& hexagon) const
{
  Vector2 a = hexagon.c - c;

  // Take advantage of symmetry and test only three slabs
  for (int i = 0; i < 3; i++)
  {
    if (fabs(a * normal[i]) > hexagon.r + r)
    {
      return false;
    }
  }
  return true;
}

/*!
\brief Compute the squared Euclidean distance to the hexagon.
\param p Point.
*/
double Hexagon2::R(const Vector2& p) const
{
  Vector2 pc = p - c;
  int s = Sector(pc);

  Vector2 a = vertex[s] * r;

  Vector2 ap = pc - a;
  Vector2 n = normal[s];

  // Distance to line
  double d = ap * n;

  // Inside
  if (d < 0.0)
  {
    return 0.0;
  }

  Vector2 e = edge[s];
  double l = ap * e;
  if (l < 0.0)
  {
    return SquaredNorm(ap);
  }
  else if (l > r)
  {
    Vector2 b = vertex[(s + 1) % 6] * r;
    Vector2 bp = pc - b;
    return SquaredNorm(bp);
  }
  else
  {
    return d;
  }
}

/*!
\brief Compute the signed Euclidean distance to the hexagon.
\param p Point.
*/
double Hexagon2::Signed(const Vector2& p) const
{
  Vector2 pc = p - c;
  int s = Sector(pc);

  Vector2 a = vertex[s] * r;
  Vector2 ap = pc - a;
  Vector2 n = normal[s];

  // Distance to line
  double d = ap * n;

  // Inside
  if (d < 0.0)
  {
    return d;
  }

  Vector2 e = edge[s];

  double l = ap * e;
  if (l <= 0.0)
  {
    return Norm(ap);
  }
  else if (l >= r)
  {
    Vector2 b = vertex[(s + 1) % 6] * r;
    Vector2 bp = pc - b;
    return Norm(bp);
  }
  else
  {
    return d;
  }
}

/*!
\brief Computes the distance vector between a hexagon and a point.
\param p Point.
*/
Vector2 Hexagon2::Normal(const Vector2& p) const
{
  Vector2 pc = p - c;
  int s = Sector(pc);

  Vector2 a = vertex[s] * r;

  Vector2 ap = pc - a;
  Vector2 n = normal[s];

  // Distance to line
  double d = ap * n;

  // Inside
  if (d < 0.0)
  {
    return Vector2::Null;
  }

  Vector2 e = edge[s];
  double l = ap * e;
  if (l < 0.0)
  {
    return ap;
  }
  else if (l > r)
  {
    Vector2 b = vertex[(s + 1) % 6] * r;
    Vector2 bp = pc - b;
    return bp;
  }
  else
  {
    return d * normal[s];
  }
}

/*!
\brief Check the intersection with a circle.

Accurate intersection test using the distance from the center of the circle to the hexagon.

\param circle The circle.
*/
bool Hexagon2::Intersect(const Circle2& circle) const
{
  return R(circle.Center()) > circle.Radius();
}

/*!
\brief Generate a random vector inside the hexagon.

\sa Triangle::RandomInside

\param random %Random number generator.
*/
Vector2 Hexagon2::RandomInside(Random& random) const
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
\brief Compute a Poisson disc distribution inside the hexagon.

This function uses a simple dart throwing algorithm.

\param r Radius of the discs.
\param n Number of candidate points.
\param random %Random number generator.
*/
QVector<Vector2> Hexagon2::Poisson(const double& r, int n, Random& random) const
{
  QVector<Vector2> p;

  // Collision radius
  double r24 = 4.0 * r * r;

  // Create instances
  for (int i = 0; i < n; i++)
  {
    Vector2 t = RandomInside(random);
    bool hit = false;
    for (int j = 0; j < p.size(); j++)
    {
      if (SquaredNorm(t - p.at(j)) < r24)
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