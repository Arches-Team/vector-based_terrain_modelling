// Box

#include <QtWidgets/QGraphicsScene>

#include "libs/box.h"

#include "libs/segment.h"

/*!
\class Box2 box.h
\brief %Axis aligned box in the plane.

\sa Box
\ingroup PlanarGroup
*/

const double Box2::epsilon = 1.0e-5; //!< Epsilon value used to check intersections and some round off errors.

const Box2 Box2::Infinity(Vector2(-Math::Infinity), Vector2(Math::Infinity)); //!< Huge bounding box, which should enclose any other 

const Box2 Box2::Null(0.0); //!< Null box, equivalent to \code Box2(Vector2(0.0)); \endcode

const Box2 Box2::Unit(Vector2(0.0), Vector2(1.0)); //!< Unit box, defined as \code Box2(Vector2(0.0),Vector2(1.0)); \endcode

/*!
\brief Compute the bounding box of a set of points.
\param p Set of points (should not be empty).
*/
Box2::Box2(const QVector<Vector2>& p)
{
  a = p.at(0);
  b = p.at(0);
  for (int i = 1; i < p.size(); i++)
  {
    a = Vector2::Min(a, p.at(i));
    b = Vector2::Max(b, p.at(i));
  }
}

/*!
\brief Create a box embedding two boxes.
\param x,y Argument boxes.
*/
Box2::Box2(const Box2& x, const Box2& y)
{
  a = Vector2::Min(x.a, y.a);
  b = Vector2::Max(x.b, y.b);
}

/*!
\brief Creates an axis aligned bounding box from a box and a transformation matrix.
\param box The box.
\param t Transformation matrix.
*/
Box2::Box2(const Box2& box, const Matrix2& t)
{
  // Set starting vertex by transformation
  a = b = t * box.Vertex(0);

  // Frame other vertices and update box
  for (int i = 1; i < 4; i++)
  {
    Vector2 p = t * box.Vertex(i);
    a = Vector2::Min(a, p);
    b = Vector2::Max(b, p);
  }
}

/*!
\brief Creates an axis aligned bounding box from a box and a frame.
\param box The box.
\param frame Transformation.
*/
Box2::Box2(const Box2& box, const Frame2& frame)
{
  // Set starting vertex by transformation
  a = b = frame.Transform(box.Vertex(0));

  // Frame other vertices and update box
  for (int i = 1; i < 4; i++)
  {
    Vector2 p = frame.Transform(box.Vertex(i));
    a = Vector2::Min(a, p);
    b = Vector2::Max(b, p);
  }
}

/*!
\brief Create a box from a Qt size.
\param size Size.
*/
Box2::Box2(const QSize& size)
{
  a = Vector2::Null;
  b = Vector2(size.width(), size.height());
}

/*!
\brief Creates the tightest embedding cube from an arbitrarilly shaped box.

This function creates a cube located at the same center point, and its side length
equal to the maximum side of the argument box.
*/
void Box2::SetCubic()
{
  Vector2 c = 0.5 * (a + b);
  Vector2 r = 0.5 * (b - a);
  r = Vector2(Math::Max(r[0], r[1]));
  a = c - r;
  b = c + r;
}

/*!
\brief Return the tightest embedding cube from an arbitrarilly shaped box.
\sa SetSubic()
*/
Box2 Box2::Cube() const
{
  Vector2 c = 0.5 * (a + b);
  Vector2 r = 0.5 * (b - a);
  r = Vector2(Math::Max(r[0], r[1]));
  return Box2(c - r, c + r);
}

/*!
\brief Creates the biggest cube inscribed in the box.

This function creates a cube located at the same center point,
and its side length equal to the minimum side of the argument box.

\sa SetCubic()
*/
void Box2::SetInscribedCubic()
{
  Vector2 c = 0.5 * (a + b);
  Vector2 r = 0.5 * (b - a);
  r = Vector(Math::Min(r[0], r[1]));
  a = c - r;
  b = c + r;
}

/*!
\brief Extend the limits of the box by a given distance.

Note that this is the same as performing the Minkowski sum with a cubic box of size r.
\param r Range.
*/
void Box2::Extend(const double& r)
{
  a -= Vector2(r);
  b += Vector2(r);
}

/*!
\brief Extend the limits of the box by a given distance.

Note that this is the same as performing the Minkowski sum with a cubic box of size r.
\param r Range.
*/
Box2 Box2::Extended(const double& r) const
{
  return Box2(a - Vector2(r), b + Vector2(r));
}

/*!
\brief Extend the limits of the box given a point.

If the point lies inside the box, the vertices of the box are unchanged.
\param p Point.
*/
void Box2::Extend(const Vector2& p)
{
  a = Vector2::Min(a, p);
  b = Vector2::Max(b, p);
}

/*!
\brief Compute the horizontal or vertical segment dividing the box into u and 1-u parts.
\param u Real.
\param horizontal Boolean, set to true for horizontal split, false for vertical.
*/
Segment2 Box2::GetSegment(const double& u, bool horizontal) const
{
  if (horizontal == true)
  {
    double y = Math::Lerp(a[1], b[1], u);
    return Segment2(Vector2(a[0], y), Vector2(b[0], y));
  }
  else
  {
    double x = Math::Lerp(a[0], b[0], u);
    return Segment2(Vector2(x, a[1]), Vector2(x, b[1]));
  }
}

/*!
\brief Test if a point is inside the box.
\param p Point.
*/
bool Box2::Inside(const Vector2& p) const
{
  if ((p[0] < a[0]) || (p[0] > b[0]) || (p[1] < a[1]) || (p[1] > b[1]))
    return false;
  else
    return true;
}

/*!
\brief Computes the intersection between two boxes.

Note that if the intersection is empty, the resulting box is invalid.

\param x Argument box.
*/
Box2 Box2::Intersection(const Box2& x) const
{
  return Box2(Vector2::Max(a, x.a), Vector2::Min(b, x.b));
}

/*!
\brief Test if a point lies withing a given range of the box.

The following two pieces of code are almost equivalent:
\code
bool a=Box2(1.0).Inside(Vector(1.5),1.0); // Do not modify the box and test range directly
\endcode
In this case we modify the geometry:
\code
Box2 box(1.0);
box.Extend(1.0); // We change the box geometry
// Note that we could have used a Minkowski sum: box=box+Box(1.0);
bool a=box.Inside(Vector(1.5));
\endcode

\param p Point.
\param r Range.
*/
bool Box2::Inside(const Vector2& p, const double& r) const
{
  if (p[0]<a[0] - r || p[0]>b[0] + r || p[1]<a[1] - r || p[1]>b[1] + r)
    return false;
  else
    return true;
}

/*!
\brief Computes the intersection between a box and a ray.

\param ray The ray.
*/
bool Box2::Intersect(const Ray2& ray) const
{
  double ta, tb;

  return Intersect(ray, ta, tb);
}

/*!
\brief Computes the intersection between a box and a ray.

\param ray The ray.
\param tmin, tmax Intersection depths
*/
bool Box2::Intersect(const Ray2& ray, double& tmin, double& tmax) const
{
  tmin = -1e16;
  tmax = 1e16;

  Vector2 p = ray.Origin();
  Vector2 d = ray.Direction();

  double t;
  // Ox
  if (d[0] < -epsilon)
  {
    t = (a[0] - p[0]) / d[0];
    if (t < tmin)
      return false;
    if (t <= tmax)
      tmax = t;
    t = (b[0] - p[0]) / d[0];
    if (t >= tmin)
    {
      if (t > tmax)
        return false;
      tmin = t;
    }
  }
  else if (d[0] > epsilon)
  {
    t = (b[0] - p[0]) / d[0];
    if (t < tmin)
      return false;
    if (t <= tmax)
      tmax = t;
    t = (a[0] - p[0]) / d[0];
    if (t >= tmin)
    {
      if (t > tmax)
        return false;
      tmin = t;
    }
  }
  else if (p[0]<a[0] || p[0]>b[0])
    return false;

  // Oy
  if (d[1] < -epsilon)
  {
    t = (a[1] - p[1]) / d[1];
    if (t < tmin)
      return false;
    if (t <= tmax)
      tmax = t;
    t = (b[1] - p[1]) / d[1];
    if (t >= tmin)
    {
      if (t > tmax)
        return false;
      tmin = t;
    }
  }
  else if (d[1] > epsilon)
  {
    t = (b[1] - p[1]) / d[1];
    if (t < tmin)
      return false;
    if (t <= tmax)
      tmax = t;
    t = (a[1] - p[1]) / d[1];
    if (t >= tmin)
    {
      if (t > tmax)
        return false;
      tmin = t;
    }
  }
  else if (p[1]<a[1] || p[1]>b[1])
    return false;

  return true;
}

/*!
\brief Computes the intersection between a box and a line.

Sorted intersection depths are returned if intersection occurs.
\param s The line.
\param tmin, tmax Intersection depths
*/
bool Box2::Intersect(const Line2& s, double& tmin, double& tmax) const
{
  tmin = -1e16;
  tmax = 1e16;

  Vector2 p = s.Vertex(0);
  Vector2 d = s.Vertex(1) - s.Vertex(0);

  double t;
  // Ox
  if (d[0] < -epsilon)
  {
    t = (a[0] - p[0]) / d[0];
    if (t < tmin)
      return false;
    if (t <= tmax)
      tmax = t;
    t = (b[0] - p[0]) / d[0];
    if (t >= tmin)
    {
      if (t > tmax)
        return false;
      tmin = t;
    }
  }
  else if (d[0] > epsilon)
  {
    t = (b[0] - p[0]) / d[0];
    if (t < tmin)
      return false;
    if (t <= tmax)
      tmax = t;
    t = (a[0] - p[0]) / d[0];
    if (t >= tmin)
    {
      if (t > tmax)
        return false;
      tmin = t;
    }
  }
  else if (p[0]<a[0] || p[0]>b[0])
    return false;

  // Oy
  if (d[1] < -epsilon)
  {
    t = (a[1] - p[1]) / d[1];
    if (t < tmin)
      return false;
    if (t <= tmax)
      tmax = t;
    t = (b[1] - p[1]) / d[1];
    if (t >= tmin)
    {
      if (t > tmax)
        return false;
      tmin = t;
    }
  }
  else if (d[1] > epsilon)
  {
    t = (b[1] - p[1]) / d[1];
    if (t < tmin)
      return false;
    if (t <= tmax)
      tmax = t;
    t = (a[1] - p[1]) / d[1];
    if (t >= tmin)
    {
      if (t > tmax)
        return false;
      tmin = t;
    }
  }
  else if (p[1]<a[1] || p[1]>b[1])
    return false;

  return true;
}

/*!
\brief Convert a planar Box2 to a Box.
\param a,b Lower and upper coordinates (note that a should be smaller than b).
*/
Box Box2::ToBox(const double& a, const double& b) const
{
  return Box(Box2::a.ToVector(a), Box2::b.ToVector(b));
}

/*!
\brief Compute the intersection between a box and a segment.

\param s The line.
\param tmin, tmax Intersection depths
\sa Box2::Intersect(const Line2&,double&,double&) const
*/
bool Box2::Intersect(const Segment2& s, double& tmin, double& tmax) const
{
  double a, b;

  // Intersect with line
  if (Intersect(Line2(s), a, b))
  {
    // Check intersection values
    tmin = Math::Max(a, 0.0);
    tmax = Math::Min(b, 1.0);
    if (tmin < tmax)
    {
      return true;
    }
    else
    {
      return false;
    }
  }
  return false;
}

/*!
\brief Check the intersection between a box and a segment.

\param s The segment.
*/
bool Box2::Intersect(const Segment2& s) const
{
  double a, b;

  // Intersect with segment
  return Intersect(s, a, b);
}

/*!
\brief Check the intersection between a box and a line.

This function is a convenience fuction for:
\code
double a,b;
// Intersect with segment
return Intersect(l,a,b);
\endcode
\param l The line.
*/
bool Box2::Intersect(const Line2& l) const
{
  double a, b;

  // Intersect with segment
  return Intersect(l, a, b);
}

/*!
\brief Compute the squared distance between the box and a point.
\param p Point.
*/
double Box2::R(const Vector2& p) const
{
  double r = 0.0;

  for (int i = 0; i < 2; i++)
  {
    if (p[i] < a[i])
    {
      double s = p[i] - a[i];
      r += s * s;
    }
    else if (p[i] > b[i])
    {
      double s = p[i] - b[i];
      r += s * s;
    }
  }
  return r;
}

/*!
\brief Computes the signed distance between the box and a point.
\param p Point.
*/
double Box2::Signed(const Vector2& p) const
{
  Vector2 c = 0.5 * (a + b);
  Vector2 d = 0.5 * (b - a);

  // To center
  Vector2 pc = p - c;

  // Symmetry
  Vector2 q = Abs(pc) - d;

  // Exterior distance
  double r = Norm(Vector2::Max(q, Vector2::Null));

  // Interior distance
  double i = Math::Min(q.Max(), 0.0);

  return r + i;
}

/*!
\brief Compute the squared Euclidean distance between two boxes.

This function computes the squared distance to avoid the computation of a square root.

\param y The box.
*/
double Box2::R(const Box2& y) const
{
  double r = 0.0;
  for (int i = 0; i < 2; i++)
  {
    if (a[i] > y.b[i])
    {
      r += (a[i] - y.b[i]) * (a[i] - y.b[i]);
    }
    else if (b[i] < y.a[i])
    {
      r += (y.a[i] - b[i]) * (y.a[i] - b[i]);
    }
    else
    {
    }
  }
  return r;
}

/*!
\brief Overloaded.
\param s Stream.
\param box The box.
*/
std::ostream& operator<<(std::ostream& s, const Box2& box)
{
  s << "Box2(" << box.a << ',' << box.b << ")";
  return s;
}

/*!
\brief Inflates a box so that its dimensions should be a fraction of its maximum side length.

\param n Fraction.
\param x,y Two integers.
*/
void Box2::SetParallelepipedic(int n, int& x, int& y)
{
  // Diagonal
  Vector2 d = (b - a);

  // Maximum side length
  double e = Math::Max(d[0], d[1]);

  double size = e / n;

  SetParallelepipedic(size, x, y);
}

/*!
\brief Creates a parallelepipedic box whose dimensions are integer
multiples of a given input reference size.

\param size Reference size, the dimension of the box will be a multiple of this size.
\param x,y Two integers.
*/
void Box2::SetParallelepipedic(const double& size, int& x, int& y)
{
  // Diagonal
  Vector2 d = (b - a);

  // Integer sizes
  // Bug tracking: adding 0.99 avoids keeping track of which indexes are the maxima 
  x = int(d[0] / size + 0.99);
  y = int(d[1] / size + 0.99);

  // Expand if necessary
  if (x == 0) { x++; }
  if (y == 0) { y++; }

  // Center
  Vector2 c = 0.5 * (a + b);

  // Diagonal
  Vector2 e = Vector2(x, y) * size / 2.0;
  a = c - e;
  b = c + e;
}

/*!
\brief Compute the coordinates of a grid aligned point.

This function computes the coordinates of a point inside the box as if the box was decomposed into a regular grid.

\param i,j Integer coordinates.
\param x,y Virtual grid size.
*/
Vector2 Box2::Vertex(int i, int j, int x, int y) const
{
  return Vector2(a[0] + i * (b[0] - a[0]) / (x - 1), a[1] + j * (b[1] - a[1]) / (y - 1));
}

/*!
\brief Generate a random vector inside the box.
\param random %Random number generator.
*/
Vector2 Box2::RandomInside(Random& random) const
{
  double r = Math::Max(b[0] - a[0], b[1] - a[1]);
  Vector2 p;
  while (true)
  {
    p = a + r * Vector2(random.Uniform(), random.Uniform());
    if (Inside(p)) break;
  }
  return p;
}

/*!
\brief Generate a random vector on the perimeter of the box.
\param random %Random number generator.
*/
Vector2 Box2::RandomOn(Random& random) const
{
  double x = b[0] - a[0];
  double y = b[1] - a[1];

  double u = 2.0 * (x + y) * random.Uniform();

  // Bottom horizontal segment
  if (u < x)
  {
    return Vector2::Lerp(a, Vertex(1), (u - x) / x);
  }

  u -= x;
  // Right vertical segment
  if (u < y)
  {
    return Vector2::Lerp(Vertex(1), b, (u - y) / y);
  }

  u -= y;
  // Top horizontal segment
  if (u < x)
  {
    return Vector2::Lerp(b, Vertex(2), (u - x) / x);
  }

  u -= x;
  // Left vertical segment
  return Vector2::Lerp(Vertex(2), a, (u - y) / y);
}

/*!
\brief Translates a box.

\param t Translation vector.
*/
void Box2::Translate(const Vector2& t)
{
  a += t;
  b += t;
}

/*!
\brief Translated box.

\param t Translation vector.
*/
Box2 Box2::Translated(const Vector2& t) const
{
  return Box2(a + t, b + t);
}

/*!
\brief Compute the box translated to origin.

This is the same as:
\code
Box centered=box;
centered.Translate(-box.Center());
\endcode
*/
Box2 Box2::Centered() const
{
  Vector2 t = 0.5 * (a + b);
  return Box2(a - t, b - t);
}

/*!
\brief Scales a box.

Note that this function handles negative coefficients in
the scaling vector.
\param s Scaling vector.
*/
void Box2::Scale(const Vector2& s)
{
  a *= s;
  b *= s;
  // Swap coordinates for negative coefficients 
  for (int i = 0; i < 2; i++)
  {
    if (s[i] < 0.0)
    {
      Math::Swap(a[i], b[i]);
    }
  }
}

/*!
\brief Scales a box.

Note that this function handles negative coefficients.
\param s Scaling factor.
*/
void Box2::Scale(const double& s)
{
  Box2::Scale(Vector2(s));
}

/*!
\brief Scales a box and return the scaled box.

\param s Scaling factor.
\sa Box2::Scaled(const Vector2&)
*/
Box2 Box2::Scaled(const double& s) const
{
  if (s > 0.0)
  {
    return Box2(a * s, b * s);
  }
  else
  {
    return Box2(b * s, a * s);
  }
}

/*!
\brief Scales a box and return the scaled box.

\param s Scaling vector.
\sa Box2::Scale
*/
Box2 Box2::Scaled(const Vector2& s) const
{
  Box2 box(a.Scaled(s), b.Scaled(s));

  // Swap coordinates for negative coefficients 
  for (int i = 0; i < 2; i++)
  {
    if (s[i] < 0.0)
    {
      Math::Swap(box[0][i], box[1][i]);
    }
  }
  return box;
}

/*!
\brief Scales a box according to a Qt size.

This function preserves the width of the box, and scales its height according to the ratio.
\param s Size.
\sa Box2::Scaled
*/
Box2 Box2::Scaled(const QSize& s) const
{
  double r = s.height() / double(s.width());
  return Box2(Vector2(a[0], r * a[1]), Vector2(b[0], r * b[1]));
}

/*!
\brief Scale the box so that the largest side should equal the argument value.

\param s Target largest size.
*/
Box2 Box2::ScaledTo(const double& s) const
{
  Vector2 c = b - a;
  if (c[0] >= c[1])
  {
    return Box2(a, b).Scaled(Vector2(s / c[0]));
  }
  else
  {
    return Box2(a, b).Scaled(Vector2(s / c[1]));
  }
}

/*!
\brief Compute the box embedding the rotated box.
\param a Rotation angle.
*/
Box2 Box2::Rotated(const double& a) const
{
  return Rotated(Matrix2::Rotation(a));
}

/*!
\brief Compute the box embedding the rotated box.
\param r Rotation matrix.
*/
Box2 Box2::Rotated(const Matrix2& r) const
{
  return Box2(*this, r);
}

/*!
\brief Computes quadrant index of a vertex with respect to the box center.

\sa Vector2::Quadrant()
\param p Point.
*/
int Box2::Quadrant(const Vector2& p) const
{
  Vector2 c = Center();
  return c.Quadrant(p);
}

/*!
\brief Draw a rectangle.
\param scene Graphics scene.
\param pen The pen.
\param brush The brush, should the box be filled.
*/
void Box2::Draw(QGraphicsScene& scene, const QPen& pen, const QBrush& brush) const
{
  scene.addRect(a[0], a[1], b[0] - a[0], b[1] - a[1], pen, brush);
}

/*!
\brief Create the Qt rectangle.
*/
QRectF Box2::GetQtRect() const
{
  return QRectF(a[0], a[1], b[0] - a[0], b[1] - a[1]);
}

/*!
\brief Compute the range of index for tiling the argument box so that it covers the box.
\param t Tiling box.
\sa Tile(const Box2&), Tile(const QRect&)
*/
QRect Box2::TileRange(const Box2& t) const
{
  Vector2 d = t.Diagonal();

  // Integer starting point

  int x = int((a[0] - t.a[0]) / d[0]) - 1;
  int y = int((a[1] - t.a[1]) / d[1]) - 1;

  int sx = int((b[0] - (t.a[0] + x * d[0])) / d[0]) + 1;
  int sy = int((b[1] - (t.a[1] + y * d[1])) / d[1]) + 1;

  return QRect(x, y, sx, sy);
}

/*!
\brief Return the tiled the box using integer coordinates.
\param x,y Integer coordinates.
\sa TileRange(const Box2&), Tile(const QRect&)
*/
Box2 Box2::Tile(int x, int y) const
{
  // Translation vector
  Vector2 t = (b - a).Scaled(Vector2(x, y));

  return Translated(t);
}

/*!
\brief Return the tiled the box using integer coordinates.
\param r Rectangle defining tiling coordinates.
\sa TileRange(const Box2&), Tile(int,int)
*/
Box2 Box2::Tile(const QRect& r) const
{
  Vector2 u = (b - a).Scaled(Vector2(r.x(), r.y()));
  Vector2 v = (b - a).Scaled(Vector2(r.x() + r.width() - 1, r.y() + r.height() - 1));
  return Box2(a + u, b + v);
}

/*!
\brief Compute a Poisson disc distribution inside a box.

This function uses a simple dart throwing algorithm.

\sa SphereTile

\param r Radius of the discs.
\param n Number of candidate points.
\param s Set of points.
\param a Add set of points flag: set to true it the set points should be added to the sampling, set to false to define constraints (typically for borders).
\param random %Random number generator.

\sa Box2::Poisson(const double& , int , Random& )
*/
QVector<Vector2> Box2::Poisson(const double& r, int n, const QVector<Vector2>& s, bool a, Random& random) const
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
\brief Compute a Poisson disc distribution inside a box.

This function uses a simple dart throwing algorithm.

\sa SphereTile

\param r Radius of the discs.
\param n Number of candidate points.
\param random %Random number generator.
*/
QVector<Vector2> Box2::Poisson(const double& r, int n, Random& random) const
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