// Polygon

#include "libs/polygon.h"
#include "libs/particle.h"
#include "libs/quadrangle.h"

// Used for drawing
#include <QtGui/QPolygonF>
#include <QtWidgets/QGraphicsScene>

/*!
\class Polygon2 polygon.h
\brief Polygons in the plane.

The class provides range-based for loops.

Note that polygons can be non-convex or even self intersecting.
\ingroup PlanarGroup
*/

/*!
\brief Create a triangle.
\param t %Triangle.
*/
Polygon2::Polygon2(const Triangle2& t) :q({ t[0], t[1], t[2] })
{
}

/*!
\brief Create a triangle.
\param a,b,c Points.
*/
Polygon2::Polygon2(const Vector2& a, const Vector2& b, const Vector2& c) :q({ a,b,c })
{
}

/*!
\brief Create a quadrangle.
\param a,b,c,d Points.
*/
Polygon2::Polygon2(const Vector2& a, const Vector2& b, const Vector2& c, const Vector2& d) :q({ a,b,c,d })
{
}

/*!
\brief Create a polygon.
\param p Array of points.
*/
Polygon2::Polygon2(const QVector<Vector>& p)
{
  q.resize(p.size());

  for (int i = 0; i < p.size(); i++)
  {
    q[i] = p.at(i);
  }
}

/*!
\brief Create a planar polygon for a polygon.
\param p Array of points.
*/
Polygon2::Polygon2(const Polygonal& p)
{
  q.resize(p.Size());

  for (int i = 0; i < p.Size(); i++)
  {
    q[i] = Vector2(p.Vertex(i));
  }
}

/*!
\brief Create a polygon approximation of an ellipse.
\param e %Ellipse.
\param n Discretization.
*/
Polygon2::Polygon2(const Ellipse2& e, int n)
{
  q.resize(n);

  for (int i = 0; i < n; i++)
  {
    double u = Math::Angle(i, n);

    q[i] = e.Vertex(u);
  }
}

/*!
\brief Create a polygon.
\param p Array of points.
*/
Polygon2::Polygon2(const QVector<Vector2>& p) :q(p)
{
}

/*!
\brief Create a polygon from a set of vertexes and a subset of indexes.
\param p Array of points.
\param indexes Set of indexes.
*/
Polygon2::Polygon2(const QVector<Vector2>& p, const QVector<int>& indexes)
{
  int n = indexes.size();
  q.resize(n);

  for (int i = 0; i < n; i++)
  {
    q[i] = p.at(indexes.at(i));
  }
}

/*!
\brief Create a polygon given an array of points.

Note that the points should be provided in trigonometric order.
\param a Array of vertices.
\param n Size.
*/
Polygon2::Polygon2(Vector2* a, int n)
{
  q.reserve(n);
  for (int i = 0; i < n; i++)
  {
    q.append(a[i]);
  }
}

/*!
\brief Create a polygon from a box.

\sa Polygon2(const Hexagon2&)
\param box The box.
*/
Polygon2::Polygon2(const Box2& box) :q({ box.Vertex(0), box.Vertex(1), box.Vertex(3), box.Vertex(2) })
{
}

#include "libs/hexagon.h"
#include "libs/pentagon.h"

/*!
\brief Create a pentagon.

\sa Polygon2(const Box2&)
\param pentagon %Pentagon.
*/
Polygon2::Polygon2(const Pentagon2& pentagon)
{
  q.reserve(5);
  for (int i = 0; i < 5; i++)
  {
    q.append(pentagon.Vertex(i));
  }
}

/*!
\brief Create a hexagon.

\sa Polygon2(const Box2&)
\param hexagon %Hexagon.
*/
Polygon2::Polygon2(const Hexagon2& hexagon)
{
  q.reserve(6);
  for (int i = 0; i < 6; i++)
  {
    q.append(hexagon.Vertex(i));
  }
}

/*!
\brief Create a quadrangle.

\sa Polygon2(const Box2&)
\param quadrangle %Quadrangle.
*/
Polygon2::Polygon2(const Quadrangle2& quadrangle) :q({ quadrangle.Vertex(0) ,quadrangle.Vertex(1) ,quadrangle.Vertex(2) ,quadrangle.Vertex(3) })
{
}

/*!
\brief Compute the bounding box of the polygon.

The polygon should have some vertices, otherwise the result is undefined.
*/
Box2 Polygon2::GetBox() const
{
  // Escape if empty
  if (q.size() == 0)
  {
    return Box2::Null;
  }

  return Box2(q);
}

/*!
\brief Compute the perimeter of the polygon.
*/
double Polygon2::Length() const
{
  if (q.size() == 0)
  {
    return 0.0;
  }

  const int n = q.size() - 1;
  double length = Norm(q.at(0) - q.at(n));

  for (int i = 0; i < n; i++)
  {
    length += Norm(q.at(i + 1) - q.at(i));
  }

  return length;
}

/*!
\brief Check whether the polygon is convex.
*/
bool Polygon2::IsConvex() const
{
  // Chose clockwise or counter-clockwise
  double side = 0;
  int i = 0;

  Vector2 a, b, c;
  do
  {
    a = q.at(i);
    b = q.at((i + 1) % q.size());
    c = q.at((i + 2) % q.size());
    side = WhichSide(c, a, b);
    i++;
  } while (i < q.size() && side == 0);

  bool bo = Line2(a, b).IsLeftOrOn(c, 0.000001);

  // Check if it is the same
  for (int j = i + 1; j < q.size(); j++)
  {
    a = q[j];
    b = q[(j + 1) % q.size()];
    c = q[(j + 2) % q.size()];

    if (bo != Line2(a, b).IsLeftOrOn(c, 0.000001))
      return false;
  }
  return true;
}

/*!
\brief Return the position of the point on the polygon at a given length from the starting point.

\param length Distance on the perimeter of the polygon. It should be less than the
perimeter of the polygon, still a modulo operation is performed inside.
*/
Vector2 Polygon2::PointAtLength(const double& length) const
{
  double l = Length();

  // Modulus
  double r = fmod(length, l);

  for (int i = 0; i < q.size(); i++)
  {
    const Vector2& a = q.at(i);
    const Vector2& b = q.at((i + 1) % q.size());
    double d = Norm(b - a);

    // In the edge
    if (r < d)
    {
      double t = r / d;
      return a * (1 - t) + b * t;
    }
    else
    {
      r -= d;
    }
  }

  // Should never get there
  return q.at(0);
}

/*!
\brief Return the normal at the position of the point on the polygon.
\param length Perimeter length from the starting point.
*/
Vector2 Polygon2::NormalAtLength(const double& length) const
{
  double l = Length();

  // Modulus
  double r = fmod(length, l);

  for (int i = 0; i < q.size(); i++)
  {
    const Vector2& a = q.at(i);
    const Vector2& b = q.at((i + 1) % q.size());
    double d = Norm(b - a);

    // In the edge
    if (r < d)
    {
      return -Normalized((b - a).Orthogonal());
    }
    else
    {
      r -= d;
    }
  }
  // Should never get there
  return Vector2::Null;
}

/*!
\brief Overloaded stream operator.
\param s Stream.
\param p Polygon.
*/
std::ostream& operator<<(std::ostream& s, const Polygon2& p)
{
  s << "Polygon2(";
  for (int i = 0; i < p.Size() - 1; i++)
  {
    s << p.Vertex(i) << ',';
  }
  s << p.Vertex(p.Size() - 1) << ')';
  return s;
}

/*!
\brief Check if a point is inside or outside of the polygon.
\param p Point.
*/
bool Polygon2::Inside(const Vector2& p) const
{
  int n = 0;

  for (int i = 0; i < q.size(); i++)
  {
    const Vector2& a = q.at(i);
    const Vector2& b = q.at((i + 1) % q.size());

    if (p[1] > Math::Min(a[1], b[1]))
    {
      if (p[1] <= Math::Max(a[1], b[1]))
      {
        if (p[0] <= Math::Max(a[0], b[0]))
        {
          if (a[1] != b[1])
          {
            double t = (p[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1]) + a[0];
            if (a[0] == b[0] || p[0] <= t)
            {
              n++;
            }
          }
        }
      }
    }
  }

  if (n % 2 == 0)
    return false;
  else
    return true;
}

/*!
\brief Compute the distance from a polygon to a line.
\param line %Line.
*/
double Polygon2::R(const Line2& line) const
{
  // Orthogonal vector to line, does not compute any square root
  Vector2 n = line.Orthogonal();

  // Get a normal with first vertex of the polygon on the backside
  if (n * (q.at(0) - line.Vertex(0)) > 0.0)
  {
    n = -n;
  }

  int k = 0;
  // Find the vertex with the maximum distance
  for (int i = 1; i < q.size(); i++)
  {
    if (n * (q.at(i) - q.at(k)) > 0.0)
    {
      k = i;
    }
  }

  // Maximum on the positive side: intersection
  if (n * (q.at(k) - line.Vertex(0)) > 0.0)
  {
    return 0.0;
  }

  // Square roots are postponed here, only if distance was not null
  return line.R(q.at(k));
}

/*!
\brief Compute the squared distance between a point and a polygon.
\param p Point.
*/
double Polygon2::R(const Vector2& p) const
{
  if (Inside(p))
  {
    return 0.0;
  }
  return RC(p);
}

/*!
\brief Compute the signed distance between a point and a polygon.
\param p Point.
*/
double Polygon2::Signed(const Vector2& p) const
{
  double r = sqrt(RC(p));
  if (Inside(p))
  {
    return -r;
  }
  return r;
}

/*!
\brief Check the intersection between the boundaty of the polygon and a circle.

Compute the signed distance between the center and the polygon.

\sa Polygon2::Signed(const Vector2&) const

\param c The circle.
*/
bool Polygon2::Intersect(const Circle2& c) const
{
  const double d = Signed(c.Center());
  const double r = c.Radius();

  if (fabs(d) < r)
  {
    return true;
  }

  return false;
}

/*!
\brief Check the position of a circle againts the polygon.

\return 0 if the circle intersects the boundary of the polygon, -1 if strictly inside, +1 if strictly outside.

\sa Polygon2::Signed(const Vector2&) const

\param c The circle.
*/
int Polygon2::Where(const Circle2& c) const
{
  const double d = Signed(c.Center());
  const double r = c.Radius();

  if (d < -r)
  {
    return -1;
  }
  else if (d > r)
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

/*!
\brief Compute the squared distance between a point and the contour of the polygon.
\param p Point.
*/
double Polygon2::RC(const Vector2& p) const
{
  double r = Segment2(q.at(q.size() - 1), q.at(0)).R(p);

  for (int i = 0; i < q.size() - 1; i++)
  {
    double t = Segment2(q.at(i), q.at((i + 1))).R(p);
    if (t < r)
    {
      r = t;
    }
  }
  return r;
}

/*!
\brief Test if a point is inside the union of a set of polygons.

Polygons should be non-intersecting. Polygons inside others will be considered as holes.

\param s Set of polygons.
\param p Point.
*/
bool Polygon2::Inside(const QVector<Polygon2>& s, const Vector2& p)
{
  int n = 0;

  for (int j = 0; j < s.size(); j++)
  {
    if (s.at(j).Inside(p))
    {
      n++;
    }
  }

  if (n % 2 == 0)
    return false;
  else
    return true;
}

/*!
\brief Test if a segment intersects the polygon.
\param s %Segment.
*/
bool Polygon2::IntersectSegment(const Segment2& s) const
{
  for (int i = 0; i < q.size(); i++)
  {
    Segment2 ab(q.at(i), q.at((i + 1) % q.size()));

    if (s.Intersect(ab))
    {
      return true;
    }

  }
  return false;
}

/*!
\brief Compute the area of the polygon.
*/
double Polygon2::Area() const
{
  // Area
  double a = 0.0;

  for (int i = 0; i < q.size(); i++)
  {
    const Vector2& pa = q.at(i);
    const Vector2& pb = q.at((i + 1) % q.size());

    a += pa[0] * pb[1] - pb[0] * pa[1];
  }

  a = Math::Abs(a);

  return 0.5 * a;
}

/*!
\brief Compute the centroid of the polygon.

Note that the centroid is not the same as the barycenter.
\sa Barycenter
*/
Vector2 Polygon2::Centroid() const
{
  double area = 0.0;
  double x = 0.0;
  double y = 0.0;

  for (int i = 0; i < q.size(); i++)
  {
    const Vector2& pa = q.at(i);
    const Vector2& pb = q.at((i + 1) % q.size());

    double a = pa[0] * pb[1] - pb[0] * pa[1];
    area += a;
    x += (pa[0] + pb[0]) * a;
    y += (pa[1] + pb[1]) * a;
  }
  area = Math::Abs(0.5 * area);

  return Vector2(x, y) / (6.0 * area);
}

/*!
\brief Compute the barycenter of the polygon.

Note that the barycenter is not the same as the centroid.
\sa Centroid
*/
Vector2 Polygon2::Center() const
{
  Vector2 g = Vector2::Null;

  for (int i = 0; i < q.size(); i++)
  {
    g += q.at(i);
  }

  return g / q.size();
}

inline int A(int k, int i, int j, int n, int p)
{
  if (i == 0) return 0;
  else
  {
    if (j < i)
    {
      return 1 + (k % p) * (n * (n + 1) / 2) + i * (i - 1) / 2 + j;
    }
    else
    {
      return 1 + ((k + 1) % p) * (n * (n + 1) / 2) + i * (i - 1) / 2 + 0;
    }
  }
}
/*!
\brief Create an n-adic subdivision of a polygon.

Note that the polygon should be star-shaped.

This function first computes the barycenter of the polygon and uses it to create triangles which are subdivided.
\param n Subdivision level.
\param vertex Array of points.
\param index Array of indexes defining the triangles.
*/
void Polygon2::Subdivide(int n, QVector<Vector2>& vertex, QVector<int>& index) const
{
  int s = q.size();
  Vector2 g = Center();

  // Starting point.
  vertex.append(g);

  // Append vertices
  for (int k = 0; k < q.size(); k++)
  {
    Vector2 e01 = (q[k] - g) / n;
    Vector2 e12 = (q[(k + 1) % q.size()] - q[k]) / n;

    for (int j = 1; j <= n; j++)
    {
      for (int i = 0; i < j; i++)
      {
        vertex.append((g + e01 * j + e12 * i));
      }
    }
  }
  // Indexes
  for (int k = 0; k < q.size(); k++)
  {
    // Counter clockwise oriented triangle connecting the center : i=1
    index.append(0);
    index.append(A(k, 1, 0, n, s));
    index.append(A(k + 1, 1, 0, n, s));

    // Other counter clockwise oriented triangles
    for (int i = 2; i <= n; i++)
    {
      for (int j = 0; j < i; j++)
      {
        index.append(A(k, i, j, n, s));
        index.append(A(k, i, j + 1, n, s));
        index.append(A(k, i - 1, j, n, s));
      }
    }

    // Other clockwise oriented triangles
    for (int i = 2; i <= n; i++)
    {
      for (int j = 1; j < i; j++)
      {
        index.append(A(k, i, j, n, s));
        index.append(A(k, i - 1, j, n, s));
        index.append(A(k, i - 1, j - 1, n, s));
      }
    }
  }
}

/*!
\brief Translate the polygon by a given vector.
\param t %Vector.
*/
void Polygon2::Translate(const Vector2& t)
{
  for (int i = 0; i < q.size(); i++)
  {
    q[i] += t;
  }
}

/*!
\brief Translate the polygon by a given vector.
\param t %Vector.
*/
Polygon2 Polygon2::Translated(const Vector2& t) const
{
  QVector<Vector2> tq(q.size());

  for (int i = 0; i < q.size(); i++)
  {
    tq[i] = q.at(i) + t;
  }

  return Polygon2(tq);
}

/*!
\brief Scale the polygon by a given factor.
\param s Scaling factor.
*/
void Polygon2::Scale(const double& s)
{
  for (int i = 0; i < q.size(); i++)
  {
    q[i] *= s;
  }
}

/*!
\brief Rotate the polygon.
\param r Rotation matrix.
*/
void Polygon2::Rotate(const Matrix2& r)
{
  for (int i = 0; i < q.size(); i++)
  {
    q[i] = r * q[i];
  }
}


/*!
\brief Convert the rectangle into a Qt polygon.
*/
QPolygonF Polygon2::GetQt() const
{
  QVector<QPointF> qp;
  for (int i = 0; i < q.size(); i++)
  {
    // Compute bits to reorder points
    Vector2 qi = q.at(i);
    qp.append(QPointF(qi[0], qi[1]));
  }

  return QPolygonF(qp);
}

/*!
\brief Draw a polygon.
\param scene Graphics scene.
\param pen The pen.
\param brush The brush.
*/
void Polygon2::Draw(QGraphicsScene& scene, const QPen& pen, const QBrush& brush) const
{
  QPolygonF polygon = GetQt();
  scene.addPolygon(polygon, pen, brush);
}

/*!
\brief Compute a simple polygon triangulation using a ear clipping method.
\return Indices triplet which define triangles.
*/
QVector<int> Polygon2::EarClip() const
{
  if (Size() < 3)
    return QVector<int>();

  int size = q.size();
  QVector<int> vertices = QVector<int>(size);
  QVector<double> angles = QVector<double>(size);
  QVector<int> indices;

  for (int i = 0; i < size; i++)
    vertices[i] = i;

  int current = 0;
  for (int i = 0; i < size; i++)
    angles[i] = Normalized(q[vertices[(i + 2) % size]] - q[vertices[(i + 1) % size]]).Angle(Normalized(q[vertices[i]] - q[vertices[(i + 1) % size]]));

  while (size > 3 && current < size)
  {
    if (angles[current] < Math::Pi)
    {
      Triangle2 triangle(q[vertices[current]], q[vertices[(current + 1) % size]], q[vertices[(current + 2) % size]]);
      bool is_ear = true;

      for (int j = current + 3; j < size + current; j++)
      {
        if (angles[(j - 1) % size] > Math::Pi && triangle.Inside(q[vertices[j % size]]))
        {
          is_ear = false;
          break;
        }
      }

      if (is_ear)
      {
        indices.append(vertices[current]);
        indices.append(vertices[(current + 1) % size]);
        indices.append(vertices[(current + 2) % size]);

        angles[(current + size - 1) % size] = Normalized(q[vertices[(current + 2) % size]] - q[vertices[current]]).Angle(Normalized(q[vertices[(current + size - 1) % size]] - q[vertices[current]]));
        angles[(current) % size] = Normalized(q[vertices[(current + 3) % size]] - q[vertices[(current + 2) % size]]).Angle(Normalized(q[vertices[current]] - q[vertices[(current + 2) % size]]));
        angles.remove((current + 1) % size);
        vertices.remove((current + 1) % size);
        size--;
        current %= size;
      }
      else
        current = (current + 1) % size;
    }
    else
      current = (current + 1) % size;
  }

  indices.append(vertices[0]);
  indices.append(vertices[1]);
  indices.append(vertices[2]);

  return indices;
}

/*!
\brief Expand the polygon by a given distance.
\param e Expansion factor.
*/
void Polygon2::Expand(const double& e)
{
  for (int i = 0; i < Size(); i++)
  {
    Vector2 d = Normalized((Vector2(q[(i + Size() - 1) % Size()]) - Vector2(q[(i + 1) % Size()])).Orthogonal());
    q[i] += d * e;
  }
}

/*!
\brief Compute the squared distance between a point and a set of polygons.
\param p Point.
\param s Set of polygons.
*/
double Polygon2::R(const QVector<Polygon2>& s, const Vector2& p)
{
  double r = s.at(0).R(p);
  for (int i = 1; i < s.size(); i++)
  {
    r = Math::Min(r, s.at(i).R((p)));
  }
  return r;
}

/*!
\brief Resample the edges of the polygon.

Generate points separated by a given input distance on the edges.

\param x Distance
*/
Polygon2 Polygon2::Resampled(const double& x) const
{
  Polygon2 polygon;

  for (int i = 0; i < q.size(); i++)
  {
    // Previous vertex
    int j = (i + q.size() - 1) % q.size();

    // Edge
    Segment2 e(q.at(j), q.at(i));

    double dist = e.Length();

    // Number of points that will be on the edge
    int n = int(dist / x) - 1;

    for (int k = 1; k <= n; k++)
    {
      double t = double(k) / double(n + 1);
      Vector2 p = e.VertexAt(t);
      polygon.Append(p);
    }

    polygon.Append(q.at(i));
  }
  return Polygon2(polygon);
}


/*!
\brief Compute a Poisson sphere distribution inside a polygon.

This function uses a simple O(n<SUP>3</SUP>) dart throwing algorithm.

\sa Box2::Poisson

\param ra Radius of the sphere.
\param n Number of candidate points.
\param border True to add vertexes of the polygon into the set (if possible).
\param random %Random number generator.
*/
QVector<Vector2> Polygon2::Poisson(const double& ra, int n, bool border, Random& random) const
{
  Box2 box = GetBox();

  ParticleSet2 p(ra);

  if (border)
  {
    // Vertexes
    for (int i = 0; i < q.size(); i++)
    {
      if (!p.Intersect(Disc2(q.at(i), ra)))
      {
        p.Append(q.at(i));
      }
    }
  }

  // Sample inside and collision detection
  for (int i = 0; i < n; i++)
  {
    Vector2 t = RandomInside(random);
    if (!p.Intersect(Disc2(t, ra)))
    {
      p.Append(t);
    }
  }
  return p.GetCenters();
}
/*!
\brief Generate a random vector inside the box.
\param random %Random number generator.
*/
Vector2 Polygon2::RandomInside(Random& random) const
{
  Vector2 p;

  Box2 box = GetBox();
  while (true)
  {
    p = box.RandomInside(random);
    if (Inside(p)) break;
  }
  return p;
}

/*!
\brief Compute the Hausdorff distance between two polygons.
\param polygon The second polygon.
\param directed Directed distance flag, set to true to compute the directed Hausdorff distance, false to compute the (symmetrized) maximum of directed Hausdorff distances.

Note that the maximum of directed Hausdorff distances is simply implemented as:
\code
Polygon2 a,b;
double d=a.Hausdorff(b,true)+b.Hausdorff(a,true);
\endcode
*/
double Polygon2::Hausdorff(const Polygon2& polygon, bool directed) const
{
  if (directed == true)
  {
    double d = 0.0;
    for (const Vector2& a : (*this))
    {
      double m = Math::Infinity;
      for (const Vector2& b : polygon)
      {
        // Postpone square root
        m = Math::Min(m, SquaredNorm(b - a));
      }
      d = Math::Max(d, m);
    }
    return sqrt(d);
  }
  else
  {
    return Hausdorff(polygon, true) + polygon.Hausdorff(*this, true);
  }
}
