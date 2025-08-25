// Pentagon

#include <QtWidgets/QGraphicsScene>

#include "libs/pentagon.h"

/*!
\class Pentagon2 pentagon.h
\brief %Flat topped pentagon in the plane.

\image html pentagon.png

\ingroup PlanarGroup
*/

const double Pentagon2::Alpha = 1.0 / (2.0 * sqrt(5.0 - 2.0 * sqrt(5.0)));

const Vector2 Pentagon2::vertex[5] = {
  Vector2(sqrt(10.0 + 2.0 * sqrt(5.0)) / 4.0,(sqrt(5.0) - 1) / 4.0),
  Vector2(0.0, 1.0),
  Vector2(-sqrt(10.0 + 2.0 * sqrt(5.0)) / 4.0,(sqrt(5.0) - 1) / 4.0),
  Vector2(-sqrt(10.0 - 2.0 * sqrt(5.0)) / 4.0,-(sqrt(5.0) + 1) / 4.0),
  Vector2(sqrt(10.0 - 2.0 * sqrt(5.0)) / 4.0,-(sqrt(5.0) + 1) / 4.0)
};

const Vector2 Pentagon2::edge[5] = {
  Normalized(vertex[1] - vertex[0]),
  Normalized(vertex[2] - vertex[1]),
  Normalized(vertex[3] - vertex[2]),
  Normalized(vertex[4] - vertex[3]),
  Normalized(vertex[0] - vertex[4])
};

const Vector2 Pentagon2::normal[5] = {
  Normalized(edge[0]).Orthogonal(),
  Normalized(edge[1]).Orthogonal(),
  Normalized(edge[2]).Orthogonal(),
  Normalized(edge[3]).Orthogonal(),
  Normalized(edge[4]).Orthogonal()
};

/*!
\brief Create a pentagon.
\param c Center.
\param r Radius.
*/
Pentagon2::Pentagon2(const Vector2& c, const double& r) :c(c), r(r)
{
}

/*!
\brief Create an pentagon.
\param r Radius.
*/
Pentagon2::Pentagon2(const double& r) :Pentagon2(Vector2::Null, r)
{
}

/*!
\brief Area of the pentagon.
*/
double Pentagon2::Area() const
{
  return 5.0 * r * r * sqrt((5.0 + sqrt(5.0)) / 2.0) / 4.0;
}

/*!
\brief Perimeter of the pentagon.
*/
double Pentagon2::Perimeter() const
{
  return 5.0 * Edge();
}

/*!
\brief Translate a pentagon.
\param t Translation vector.
*/
void Pentagon2::Translate(const Vector2& t)
{
  c += t;
}

/*!
\brief Scale a pentagon.
\param s Scaling factor.
*/
void Pentagon2::Scale(const double& s)
{
  c *= s;
  r *= s;
}

/*!
\brief Computes the bounding box of a pentagon.
*/
Box2 Pentagon2::GetBox() const
{
  static const double a = (Math::Sqrt5 + 1.0) / 4.0;
  static const double b = sqrt(10.0 + 2.0 * sqrt(5.0)) / 4.0;
  return Box2(Vector2(-b * r, -a * r), Vector2(b * r, r)).Translated(c);
}

/*!
\brief Overloaded.
\param s Stream.
\param pentagon The pentagon.
*/
std::ostream& operator<<(std::ostream& s, const Pentagon2& pentagon)
{
  s << "Pentagon2(" << pentagon.c << ',' << pentagon.r << ")";
  return s;
}

/*!
\brief Draws an pentagon.
\param scene Graphics scene.
\param pen The pen.
\param brush The brush.
*/
void Pentagon2::Draw(QGraphicsScene& scene, const QPen& pen, const QBrush& brush) const
{
  QPolygonF polygon;
  for (int i = 0; i < 5; i++)
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
int Pentagon2::Sector(const Vector2& d)
{
  bool a = d[0] >= 0.0;

  // Vertical symmetry
  Vector2 s(fabs(d[0]), d[1]);

  bool b = s * edge[0] >= 0.0;
  bool c = s * edge[4] >= 0.0;

  // Line of a = 0Â°
  if (a)
  {
    if (b)
    {
      return 1;
    }
    else
    {
      if (c)
      {
        return 0;
      }
      else
      {
        return 4;
      }
    }
  }
  else
  {
    if (b)
    {
      return 1;
    }
    else
    {
      if (c)
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
\brief Test if a point is inside.
\param p Point.
*/
bool Pentagon2::Inside(const Vector2& p) const
{
  Vector2 cp = p - c;

  // Take advantage of symmetry and test only three planes
  cp[0] = fabs(cp[0]);

  // Three planes
  if (cp * normal[0] - r > 0) return false;
  if (cp * normal[4] - r > 0) return false;
  if (cp * normal[3] - r > 0) return false;

  return true;
}

/*!
\brief Compute the squared Euclidean distance to the pentagon.
\param p Point.
*/
double Pentagon2::R(const Vector2& p) const
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
    Vector2 b = vertex[(s + 1) % 5] * r;
    Vector2 bp = pc - b;
    return SquaredNorm(bp);
  }
  else
  {
    return d;
  }
}

/*!
\brief Compute the signed Euclidean distance to the pentagon.
\param p Point.
*/
double Pentagon2::Signed(const Vector2& p) const
{
  /*
  // signed distance to a regular pentagon
  float sdPentagon(in vec2 p, in float r)
  {
    const vec3 k = vec3(0.809016994, 0.587785252, 0.726542528); // pi/5: cos, sin, tan
   // p.y = -p.y;
    p.x = abs(p.x);
    p -= 2.0 * min(dot(vec2(-k.x, k.y), p), 0.0) * vec2(-k.x, k.y);
    p -= 2.0 * min(dot(vec2(k.x, k.y), p), 0.0) * vec2(k.x, k.y);
    p -= vec2(clamp(p.x, -r * k.z, r * k.z), r);
    return length(p) * sign(p.y);
  ////1) p.y = -p.y; // This is not a reflection (folding) in the sense of the others, but just an inversion around the x axis. For convenience and aesthetic orientation, I guess.

//2) p.x = abs(p.x); // Folding about the y axis, so we treat the left side as a mirror reflection of the right side.

//3) p -= 2.0*min(dot(vec2(-k.x,k.y),p),0.0)*vec2(-k.x,k.y); // Fold about the line through the origin whose angle is (3/2 + 4/5)pi.
//This is done by projecting p along a vector perpendicular to that line, [-k.x, k.y] (whose angle is 4/5 pi); and if the projection has a negative length, subtracting double the projection.

//4) p -= 2.0*min(dot(vec2(k.x,k.y),p),0.0)*vec2(k.x,k.y); // Fold about the line through the origin whose angle is (3/2 - 4/5)pi in a similar way.

//5) Finally, compute signed distance to (the closest point on) a horizontal line segment y = r, whose length is 2 r tan pi/5.

  */

  const Vector k(0.809016994, 0.587785252, 0.726542528); // pi/5: cos, sin, tan
  // p.y = -p.y;
  Vector2 q = p;
  q[0] = abs(q[0]);
  q -= 2.0 * Math::Min(Vector2(-k[0], k[1]) * p, 0.0) * Vector2(-k[0], k[1]);
  q -= 2.0 * Math::Min(Vector2(k[0], k[1]) * p, 0.0) * Vector2(k[0], k[1]);
  q -= Vector2(Math::Clamp(p[0], -r * k[2], r * k[2]), r);
  return Norm(p) * Math::Sign(p[1]);


  /*Vector2 pc = p - c;
  int s = Sector(pc);
  int t = (s + 5 - 1) % 5;

  Vector2 a = vertex[s] * r;

  Vector2 ap = pc - a;

  double eu = ap * edge[s];
  double ev = ap * edge[t];

  // Edge
  if (eu > 0.0)
  {
    double d = ap * normal[s];
    if (d > 0.0)
    {
      return d;
    }
    else
    {
      return Math::Max(d, ap * normal[t]);
    }
  }
  else
  {
    // Edge
    if (ev < 0.0)
    {
      double d = ap * normal[t];
      if (d > 0.0)
      {
        return d;
      }
      else
      {
        return Math::Max(d, ap * normal[s]);
      }
    }
    // Vertex
    else
    {
      return Norm(ap);
    }
  }
  */
}

/*!
\brief Computes the distance vector between a pentagon and a point.
\param p Point.
*/
Vector2 Pentagon2::Normal(const Vector2& p) const
{
  Vector2 pc = p - c;
  int s = Sector(pc);
  int t = (s + 5 - 1) % 5;

  Vector2 a = vertex[s] * r;

  Vector2 ap = pc - a;

  double eu = ap * edge[s];
  double ev = ap * edge[t];

  // Edge
  if (eu > 0.0)
  {
    double d = ap * normal[s];
    if (d > 0.0)
    {
      return d * normal[s];
    }
    else
    {
      return Vector2::Null;
    }
  }
  else
  {
    // Edge
    if (ev < 0.0)
    {
      double d = ap * normal[t];
      if (d > 0.0)
      {
        return d * normal[t];
      }
      else
      {
        return Vector2::Null;
      }
    }
    // Vertex
    else
    {
      return ap;
    }
  }
}

const double Pentagon2::epsilon = 1.0e-7;

/*!
\brief Compute the intersection with a ray.

\param ray The ray.
\param ta, tb Intersection depths
*/
bool Pentagon2::Intersect(const Ray2& ray, double& ta, double& tb) const
{
  ta = -Math::Infinity;
  tb = Math::Infinity;

  // Half space intersection
  for (int i = 0; i < 5; i++)
  {
    Vector2 pi = Vertex(i);
    Vector2 ni = Edge(i).Orthogonal();

    double e = ray.Direction() * ni;

    // Lines are parallel
    if (fabs(e) < epsilon) continue;

    double t = (ray.Origin() - pi) * ni / e;

    if (e < 0.0)
    {
      ta = Math::Max(ta, t);
    }
    else
    {
      ta = Math::Min(tb, t);
    }
  }
  return (ta < tb);
}


/*!
\brief Generate a random vector inside the pentagon.

\sa Triangle::RandomInside

\param random %Random number generator.
*/
Vector2 Pentagon2::RandomInside(Random& random) const
{
  Vector2 p;
  while (true)
  {
    p = c + r * Vector2(random.Uniform(), random.Uniform());
    if (Inside(p)) break;
  }
  return p;
}


/*!
\brief Compute a Poisson disc distribution inside the pentagon.

This function uses a simple dart throwing algorithm.

\param r Radius of the discs.
\param n Number of candidate points.
\param random %Random number generator.
*/
QVector<Vector2> Pentagon2::Poisson(const double& r, int n, Random& random) const
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