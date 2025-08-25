// Polygon

#include "libs/polygon.h"
#include "libs/matrix.h"

/*!
\class Polygonal polygon.h
\brief Polygons.

The name of this class was chosen to avoid conflicts with a pre-defined Windows open-GL class.

Note that polygons can be non-convex.

\ingroup ExtendedKernelGroup
*/

/*!
\brief Create a polygon.
\param p Array of points.
*/
Polygonal::Polygonal(const QVector<Vector>& p) :vertex(p)
{
}

/*!
\brief Create a planar polygon.
\param p Array of points.
*/
Polygonal::Polygonal(const QVector<Vector2>& p)
{
  vertex.reserve(p.size());

  for (int i = 0; i < p.size(); i++)
  {
    vertex[i] = p.at(i).ToVector();
  }
}

/*!
\brief Create a polygon from a set of segments.

The polygon is created from the first vertices of all segments.
\param e Array of segments.
*/
Polygonal::Polygonal(const QVector<Segment>& e)
{
  for (int i = 0; i < e.size(); i++)
  {
    vertex.append(e[i].Vertex(0));
  }
}

/*!
\brief Create a polygon given an array of points.

Note that the points should be provided in trigonometric order.
\param a Array of vertices.
\param n Size.
*/
Polygonal::Polygonal(Vector* a, int n)
{
  vertex.reserve(n);

  for (int i = 0; i < n; i++)
  {
    vertex.append(a[i]);
  }
}

/*!
\brief Create a polygon from a box.
\param box The box.
*/
Polygonal::Polygonal(const Box2& box)
{
  vertex.append(Vector(box[0][0], box[0][1], 0.0));
  vertex.append(Vector(box[1][0], box[0][1], 0.0));
  vertex.append(Vector(box[1][0], box[1][1], 0.0));
  vertex.append(Vector(box[0][0], box[1][1], 0.0));
}

/*!
\brief Compute the normal from the first three vertexes.
*/
Vector Polygonal::Normal() const
{
  // Face normal
  return Normalized((vertex.at(1) - vertex.at(0)) / (vertex.at(2) - vertex.at(1)));
}

/*!
\brief Compute the bounding box.

The polygon should have some vertices, otherwise the result is undefined.
*/
Box Polygonal::GetBox() const
{
  // Escape if empty
  if (vertex.size() == 0)
  {
    return Box::Null;
  }

  Vector a = vertex.at(0);
  Vector b = vertex.at(0);

  for (int i = 1; i < vertex.size(); i++)
  {
    a = Vector::Min(a, vertex.at(i));
    b = Vector::Max(b, vertex.at(i));
  }

  return Box(a, b);
}

/*!
\brief Compute the perimeter of the polygon.
*/
double Polygonal::Length() const
{
  double length = 0.0;

  for (int i = 0; i < vertex.size(); i++)
  {
    length += Norm(vertex.at((i + 1) % vertex.size()) - vertex.at(i));
  }

  return length;
}

/*!
\brief Compute the perimeter of the projection of the polygon on the horizontal plane.

This is equivalent to projecting the polygon on the plane and computing its perimeter:
\code
Polygonal p;
double length = Polygon2(p).Length();
\endcode
*/
double Polygonal::PlanarLength() const
{
  double length = 0.0;

  for (int i = 0; i < vertex.size(); i++)
  {
    length += Norm(Vector2(vertex.at((i + 1) % vertex.size())) - Vector2(vertex.at(i)));
  }

  return length;
}

/*!
\brief Return the position of the point on the polygon at a given length from the starting point.

\param length Distance on the perimeter of the polygon. Note that it should be less than the perimeter of the polygon, still a modulo operation is performed inside.
*/
Vector Polygonal::PointAtLength(const double& length) const
{
  double r = length;
  double l = Length();

  r = fmod(r, l);

  for (int i = 0; i < vertex.size(); i++)
  {
    const Vector& a = vertex.at(i);
    const Vector& b = vertex.at((i + 1) % vertex.size());
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
  return vertex.at(0);
}

/*!
\brief Return the normal at the position of the point on the polygon.
\param length Given length from the starting point.
*/
Vector Polygonal::NormalAtLength(const double& length) const
{
  double r = length;
  double l = Length();

  r = fmod(r, l);

  for (int i = 0; i < vertex.size(); i++)
  {
    const Vector& a = vertex.at(i);
    const Vector& b = vertex.at((i + 1) % vertex.size());
    double d = Norm(b - a);

    // In the edge
    if (r < d)
    {
      Matrix mt_rot = Matrix::RotationZ(-Math::HalfPi);
      Vector vt_rot = mt_rot * (b - a);
      return Normalized(vt_rot);
    }
    else
    {
      r -= d;
    }
  }
  // Should never get there
  return Vector(-1, 0, 0);
}

/*!
\brief Overloaded stream operator.
\param s Stream.
\param p Polygon.
*/
std::ostream& operator<<(std::ostream& s, const Polygonal& p)
{
  s << "Polygonal(";
  for (int i = 0; i < p.Size() - 1; i++)
  {
    s << p.Vertex(i) << ',';
  }
  s << p.Vertex(p.Size() - 1) << ')';
  return s;
}

/*!
\brief Compute the area of the polygon.
*/
double Polygonal::Area() const
{
  double a = 0.0;

  const Vector& po = vertex.at(0);

  for (int i = 0; i < vertex.size(); i++)
  {
    const Vector& pa = vertex.at(i);
    const Vector& pb = vertex.at((i + 1) % vertex.size());

    a += Norm((pa - po) / (pb - po));
  }

  a = fabs(a);

  return 0.5 * a;
}

/*!
\brief Compute the centroid of the polygon.

Note that the centroid is not the same as the barycenter.
\sa Barycenter
*/
Vector Polygonal::Centroid() const
{
  double area = 0.0;
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;

  const Vector& po = vertex.at(0);

  for (int i = 0; i < vertex.size(); i++)
  {
    const Vector& pa = vertex.at(i);
    const Vector& pb = vertex.at((i + 1) % vertex.size());

    double a = Norm((pa - po) / (pb - po));
    area += a;
    x += (pa[0] + pb[0]) * a;
    y += (pa[1] + pb[1]) * a;
    z += (pa[2] + pb[2]) * a;
  }
  area = fabs(0.5 * area);

  return Vector(x, y, z) / (6.0 * area);
}

/*!
\brief Compute the barycenter of the polygon.

Note that the barycenter is not the same as the centroid.
\sa Centroid
*/
Vector Polygonal::Center() const
{
  Vector g(0.0);

  for (int i = 0; i < vertex.size(); i++)
  {
    g += vertex.at(i);
  }

  return g / vertex.size();
}

static inline int A(int k, int i, int j, int n, int p)
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
void Polygonal::Subdivide(int n, QVector<Vector>& vertex, QVector<int>& index) const
{
  int s = vertex.size();
  Vector2 g = Center();

  // Starting point.
  vertex.append(g.ToVector());

  // Append vertices
  for (int k = 0; k < vertex.size(); k++)
  {
    Vector2 e01 = (vertex.at(k) - g) / n;
    Vector2 e12 = (vertex.at((k + 1) % vertex.size()) - vertex.at(k)) / n;

    for (int j = 1; j <= n; j++)
    {
      for (int i = 0; i < j; i++)
      {
        vertex.append((g + e01 * j + e12 * i).ToVector());
      }
    }
  }
  // Indexes
  for (int k = 0; k < vertex.size(); k++)
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
void Polygonal::Translate(const Vector& t)
{
  for (int i = 0; i < vertex.size(); i++)
  {
    vertex[i] += t;
  }
}

/*!
\brief Scale the polygon by a given factor.
\param s Scaling factor.
*/
void Polygonal::Scale(const double& s)
{
  for (int i = 0; i < vertex.size(); i++)
  {
    vertex[i] *= s;
  }
}

/*!
\brief Rotate the polygon.
\param r Rotation matrix.
*/
void Polygonal::Rotate(const Matrix& r)
{
  for (int i = 0; i < vertex.size(); i++)
  {
    vertex[i] = r * vertex[i];
  }
}

/*!
\brief Computes the vector distance between the polygon and a point in space.

It analyses the location of the point againt the Voronoi diagram of the polygon,
and projects it to derive the right vector distance.

\param p Point.
*/
Vector Polygonal::Normal(const Vector& p) const
{
  // Normal
  Vector up = Normal();

  const int n = vertex.size();

  // Previous edge
  Vector edgej = vertex[0] - vertex[n - 1];

  Vector pi;

  for (int i = 0; i < n; i++)
  {
    // Edge before edge i
    Vector edgei = vertex[(i + 1) % n] - vertex[i];

    pi = p - vertex.at(i);

    // Scalar product to tell side
    double s = pi * edgei;

    if (s <= 0.0)
    {
      if (pi * edgej >= 0.0)
      {
        // Vertex i
        return pi;
      }
    }
    else
    {
      Vector pj = p - vertex.at((i + 1) % n);

      if (pj * edgei < 0.0)
      {
        Vector orthogonal = edgei / up;
        if (pi * orthogonal >= 0.0)
        {
          // Edge 
          return pi - s * edgei;
        }
      }
    }
    edgej = edgei;
  }

  // Face
  double t = up * pi;
  return up * t;
}

/*!
\brief Computes the squered distance between the polygon and a point.

\param p Point.
*/
double Polygonal::R(const Vector& p) const
{
  return SquaredNorm(Normal(p));
}