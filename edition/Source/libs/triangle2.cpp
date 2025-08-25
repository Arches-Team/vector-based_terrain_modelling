// Triangle

#include <QtWidgets/QGraphicsScene>

#include "libs/triangle.h"
#include "libs/segment.h"

/*!
\class Triangle2 triangle.h
\brief Base minimum storage triangle class in the plane.

\ingroup PlanarGroup
*/

double Triangle2::epsilon = 1.0e-7;

/*!
\brief Create a triangle.
\param t Triangle.
*/
Triangle2::Triangle2(const Triangle& t)
{
  p[0] = t[0];
  p[1] = t[1];
  p[2] = t[2];
}

/*!
\brief Computes the axis aligned box enclosing the triangle.
*/
Box2 Triangle2::GetBox() const
{
  return Box2(Vector2::Min(Vector2::Min(p[0], p[1]), p[2]), Vector2::Max(Vector2::Max(p[0], p[1]), p[2]));
}

/*!
\brief Compute the barycenter of the triangle.
*/
Vector2 Triangle2::Center() const
{
  return (p[0] + p[1] + p[2]) / 3.0;
}

/*!
\brief Compute the area of the triangle.
\sa Triangle2::SignedArea
*/
double Triangle2::Area() const
{
  return 0.5 * Math::Abs((p[0] - p[1]) / (p[2] - p[0]));
}

/*!
\brief Compute the area of the triangle.
\sa Triangle2::Area
*/
double Triangle2::SignedArea() const
{
  return 0.5 * ((p[0] - p[1]) / (p[2] - p[0]));
}

/*!
\brief Compute the perimeter of the triangle.
 */
double Triangle2::Perimeter() const
{
  return Norm(p[1] - p[0]) + Norm(p[2] - p[1]) + Norm(p[0] - p[2]);
}

/*!
\brief Compute the radius of the circle inscribed in the
triangle.

Some algebra will show that it is half the ratio
of the half perimeter and the surface of the triangle.
*/
double Triangle2::InscribedRadius() const
{
  Vector2 u = p[0] - p[1];
  Vector2 v = p[2] - p[0];
  Vector2 w = p[1] - p[2];
  double a = Norm(u) + Norm(v) + Norm(w);
  double s = Math::Abs(u / v);
  return s / a;
}

/*!
\brief Compute the barycentric coordinates of an argument vector with respect to the triangle.

If the point lies within the triangle, the coordinates are positive and form a partition of the unity.
This no longer holds if the point is outside of the triangle.

Note that this function is not the most efficient, as it relies on the evaluation of the area
formed by the triangles pab, pbc, and pca.

\sa Triangle::BarycentricCoordinates()

\param p Point.
*/
Vector Triangle2::BarycentricCoordinates(const Vector2& p) const
{
  double a = Area();
  double aa = Triangle2(p, Triangle2::p[1], Triangle2::p[2]).Area();
  double ab = Triangle2(p, Triangle2::p[0], Triangle2::p[2]).Area();
  double ac = Triangle2(p, Triangle2::p[0], Triangle2::p[1]).Area();

  return Vector(aa, ab, ac) / a;
}

/*!
\brief Compute the barycenter given barycentric coordinates.

\sa Triangle2::BarycentricCoordinates()

\param weights Weights.
*/
Vector2 Triangle2::BaryCenter(const Vector& weights) const
{
  return weights[0] * p[0] + weights[0] * p[1] + weights[0] * p[2];
}

/*!
\brief Compute the orthocenter.
*/
Vector2 Triangle2::OrthoCenter() const
{
  // Lengths
  double c = SquaredNorm(p[0] - p[1]);
  double a = SquaredNorm(p[1] - p[2]);
  double b = SquaredNorm(p[2] - p[0]);

  return BaryCenter(Vector(1.0 / (b + c - a), 1.0 / (c + a - b), 1.0 / (a + b - c)));
}

/*!
\brief Compute the circle inscribed in the triangle.

The center of the circle
is the center of (A,a), (B,b) and (C,c) where A, B and C are the vertices of
the triangle and coefficients a, b and c represent the length of their facing
edge.
*/
Circle2 Triangle2::Inscribed() const
{
  Vector2 u = p[0] - p[1];
  Vector2 v = p[2] - p[0];
  Vector2 w = p[1] - p[2];

  double a = Norm(w);
  double b = Norm(v);
  double c = Norm(u);

  double l = a + b + c;
  return Circle2((p[0] * a + p[1] * b + p[2] * c) / l, fabs(u / v) / l);
}

/*!
\brief Compute the radius of the circumscribed circle of the triangle.
*/
double Triangle2::CircumscribedRadius() const
{
  double u = Norm(p[0] - p[1]);
  double v = Norm(p[1] - p[2]);
  double w = Norm(p[2] - p[0]);
  return u * v * w / sqrt((u + v + w) * (-u + v + w) * (u - v + w) * (u + v - w));
}

/*!
\brief Overloaded.
\param t %Triangle.
\param s Stream.
*/
std::ostream& operator<<(std::ostream& s, const Triangle2& t)
{
  s << "Triangle2(" << t.p[0] << ',' << t.p[1] << ',' << t.p[2] << ')';
  return s;
}

/*!
\brief Computes the aspect ratio of the triangle.

The aspect ratio is defined as twice
the radius of the inscribed circle divided by the radius of the
circumscribing circle.

\sa Triangle::Aspect()
*/
double Triangle2::Aspect() const
{
  double ab = Norm(p[1] - p[0]);
  double bc = Norm(p[2] - p[1]);
  double ca = Norm(p[0] - p[2]);

  double s = 0.5 * (ab + bc + ca);

  double u = (s - ab) * (s - bc) * (s - ca);

  return 8.0 * u / (ab * bc * ca);
}

/*!
\brief Draw a triangle.
\param scene Graphics scene.
\param pen The pen.
\param brush The brush, should the triangle be filled.
*/
void Triangle2::Draw(QGraphicsScene& scene, const QPen& pen, const QBrush& brush) const
{
  QPolygonF polygon(3);

  // Fill with three vertices
  polygon[0] = QPointF(p[0][0], p[0][1]);
  polygon[1] = QPointF(p[1][0], p[1][1]);
  polygon[2] = QPointF(p[2][0], p[2][1]);

  scene.addPolygon(polygon, pen, brush);
}

/*!
\brief Test if a box lies inside the triangle.

\param box The box.
*/
bool Triangle2::Inside(const Box2& box) const
{
  for (int i = 0; i < 4; i++)
  {
    if (!Inside(box.Vertex(i)))
    {
      return false;
    }
  }
  return true;
}

/*!
\brief Test if a circle lies inside the triangle.

\param circle The circle.
*/
bool Triangle2::Inside(const Circle2& circle) const
{
  const Vector2& c = circle.Center();
  const double& r = circle.Radius();
  for (int i = 0; i < 3; i++)
  {
    Vector2 n = Normalized(p[(i + 1) % 3] - p[i]).Orthogonal();
    double d = (c - p[i]) * n;
    if (d >= -r)
    {
      return false;
    }
  }
  return true;
}

/*!
\brief Test if a point lies inside the triangle.

This is an optimized version with some early tests to avoid some computations.

\sa Triangle2::Inside(const Vector2&,double&,double&)
\param p Point.
*/
bool Triangle2::Inside(const Vector2& p) const
{
  Vector2 ab = Triangle2::p[1] - Triangle2::p[0];
  Vector2 ac = Triangle2::p[2] - Triangle2::p[0];
  Vector2 ap = p - Triangle2::p[0];

  // Compute dot products
  double abac = ab * ac;
  double abap = ab * ap;
  double acac = ac * ac;
  double acap = ac * ap;

  // Compute barycentric coordinates
  double u = (acac * abap - abac * acap);
  if (u < 0.0) return false;

  // Compute last dot product
  double abab = ab * ab;
  double v = (abab * acap - abac * abap);
  if (v < 0.0) return false;

  double d = (abab * acac - abac * abac);

  // Check if point is in triangle
  return (u + v <= d);
}

/*!
\brief Test if a point lies inside the triangle and computes its uv-coordinates.
\param p Point.
\param u, v UV-coordinates.
*/
bool Triangle2::Inside(const Vector2& p, double& u, double& v) const
{
  Vector2 ab = Triangle2::p[1] - Triangle2::p[0];
  Vector2 ac = Triangle2::p[2] - Triangle2::p[0];
  Vector2 ap = p - Triangle2::p[0];

  // Compute dot products
  double abab = ab * ab;
  double abac = ab * ac;
  double abap = ab * ap;
  double acac = ac * ac;
  double acap = ac * ap;
  double d = 1.0 / (abab * acac - abac * abac);

  // Compute barycentric coordinates
  u = (acac * abap - abac * acap) * d;
  v = (abab * acap - abac * abap) * d;

  // Check if point is in triangle
  return (u >= 0) && (v >= 0) && (u + v <= 1);
}

/*!
\brief Generate a random vector inside the triangle.

\sa Triangle::RandomInside

\param random %Random number generator.
*/
Vector2 Triangle2::RandomInside(Random& random) const
{
  double u = sqrt(random.Uniform());
  double v = random.Uniform();

  return (1.0 - u) * p[0] + u * ((1.0 - v) * p[1] + v * p[2]);
}

/*!
\brief Computes the normal vector to a triangle.
\param q Point.
*/
Vector2 Triangle2::Normal(const Vector2& q) const
{
  // Edges
  const Vector2 e[3] = { Normalized(p[1] - p[0]), Normalized(p[2] - p[1]), Normalized(p[0] - p[2]) };

  for (int i = 0; i < 3; i++)
  {
    Vector2 pa = q - p[i];

    // Scalar product to tell side
    double s = pa * e[i];
    if (s <= 0.0)
    {
      if (pa * e[(i + 2) % 3] >= 0.0)
      {
        // Vertex A
        return pa;
      }
    }
    else
    {
      pa = q - p[(i + 1) % 3];
      s = pa * e[i];

      if (s < 0.0)
      {
        if (pa * e[i].Orthogonal() >= 0.0)
        {
          // Edge [AB]
          return pa - s * e[i];
        }
      }
    }
  }

  // Inside
  return Vector2::Null;
}

/*!
\brief Computes the squared distance between a point and a triangle.
\param q Point.
*/
double Triangle2::R(const Vector2& q) const
{
  Vector2 n = Normal(q);
  return n * n;
}

/*!
\brief Check the intersection with a circle.
\param circle The circle.
*/
bool Triangle2::Intersect(const Circle2& circle) const
{
  if (R(circle.Center()) < Math::Sqr(circle.Radius()))
  {
    return true;
  }
  else
  {
    return false;
  }
}

/*!
\brief Check the intersection with a segment.
\param segment The segment.
*/
bool Triangle2::Intersect(const Segment2& segment) const
{
  if (Segment2(p[0], p[1]).Intersect(segment)) return true;
  if (Segment2(p[1], p[2]).Intersect(segment)) return true;
  if (Segment2(p[2], p[0]).Intersect(segment)) return true;

  return false;
}

/*!
\brief Check the intersection with a box.
\param box The box.
*/
bool Triangle2::Intersect(const Box2& box) const
{
  // Inclusion of only one point means intersection
  if (Inside(box.Vertex(0))) return true;

  // For points inside box, we still test three as a speed up
  if (box.Inside(p[0])) return true;
  if (box.Inside(p[1])) return true;
  if (box.Inside(p[2])) return true;

  // Box intersecting segments of the triangle
  if (box.Intersect(Segment2(p[0], p[1]))) return true;
  if (box.Intersect(Segment2(p[1], p[2]))) return true;
  if (box.Intersect(Segment2(p[2], p[0]))) return true;

  // Box intersecting segments of the triangle
  for (int i = 0; i < 4; i++)
  {
    if (Intersect(Segment2(box.Vertex(i), box.Vertex((i + 1) % 4)))) return true;
  }
  return false;
}

/*!
\brief Compute the intersection with a ray.

\param ray The ray.
\param ta, tb Intersection depths
*/
bool Triangle2::Intersect(const Ray2& ray, double& ta, double& tb) const
{
  ta = -Math::Infinity;
  tb = Math::Infinity;

  // Half space intersection
  for (int i = 0; i < 3; i++)
  {
    Vector2 pi = p[i];
    Vector2 ni = (p[(i + 1) % 3] - p[i]).Orthogonal();

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
\brief Rotate the triangle.
\param r Rotation matrix.
*/
void Triangle2::Rotate(const Matrix2& r)
{
  for (int i = 0; i < 3; i++)
  {
    p[i] = r * p[i];
  }
}

/*!
\brief Translate the triangle.
\param t Translation vector.
*/
void Triangle2::Translate(const Vector2& t)
{
  for (int i = 0; i < 3; i++)
  {
    p[i] += t;
  }
}

/*!
\brief Scale the triangle.
\param s Scaling factor.
*/
void Triangle2::Scale(const double& s)
{
  for (int i = 0; i < 3; i++)
  {
    p[i] *= s;
  }
}

/*!
\brief Scale the triangle.
\param s Scaling vector.
*/
void Triangle2::Scale(const Vector2& s)
{
  for (int i = 0; i < 3; i++)
  {
    p[i] = p[i].Scaled(s);
  }
}

/*!
\brief Shrinks a triangle.

\param e Erosion radius.
*/
void Triangle2::Shrink(const double& e)
{
  // Planes orthogonal to the edges
  Line2 ab(p[0], p[1]);
  Line2 bc(p[0], p[1]);
  Line2 ca(p[0], p[1]);

  // Move
  ab.Translate(Normalized(ab.GetAxis()).Orthogonal() * e);
  bc.Translate(Normalized(bc.GetAxis()).Orthogonal() * e);
  ca.Translate(Normalized(ca.GetAxis()).Orthogonal() * e);

  // Points
  ca.Intersection(ca, p[0]);
  ab.Intersection(ab, p[1]);
  bc.Intersection(bc, p[2]);
}

/*!
\brief Test if the edges of two triangles intersect.
\param t %Triangle.
*/
bool Triangle2::IntersectEdge(const Triangle2& t) const
{
  if (WhichSide(p[1], t[2], t[0]) >= 0.0)
  {
    if (WhichSide(p[1], p[0], t[0]) >= 0.0)
    {
      if (WhichSide(t[2], p[0], p[1]) >= 0.0) return true;
      else return false;
    }
    else
    {
      if (WhichSide(t[0], p[1], p[2]) >= 0.0)
      {
        if (WhichSide(t[0], p[2], p[0]) >= 0.0) return true;
        else return false;
      }
      else return false;
    }
  }
  else
  {
    if (WhichSide(p[2], t[2], t[0]) >= 0.0)
    {
      if (WhichSide(p[2], p[0], t[0]) >= 0.0)
      {
        if (WhichSide(t[2], p[0], p[2]) >= 0.0) return true;
        else
        {
          if (WhichSide(t[2], p[1], p[2]) >= 0.0) return true;
          else return false;
        }
      }
      else  return false;
    }
    else return false;
  }
}

/*!
\brief Test if two triangles are intersecting by performing a vertex intersection test.
\param t %Triangle.

Triangles should have their vertices in trigonometric order.
*/
bool Triangle2::IntersectVertex(const Triangle2& t) const
{
  if (WhichSide(p[1], t[2], t[0]) >= 0.0)
    if (WhichSide(p[1], t[2], t[1]) <= 0.0)
      if (WhichSide(p[1], p[0], t[0]) > 0.0) {
        if (WhichSide(p[1], p[0], t[1]) <= 0.0) return true;
        else return false;
      }
      else {
        if (WhichSide(p[2], p[0], t[0]) >= 0.0)
          if (WhichSide(t[0], p[1], p[2]) >= 0.0) return true;
          else return false;
        else return false;
      }
    else
      if (WhichSide(p[1], p[0], t[1]) <= 0.0)
        if (WhichSide(p[2], t[2], t[1]) <= 0.0)
          if (WhichSide(t[1], p[1], p[2]) >= 0.0) return true;
          else return false;
        else return false;
      else return false;
  else
    if (WhichSide(p[2], t[2], t[0]) >= 0.0)
      if (WhichSide(t[2], p[1], p[2]) >= 0.0)
        if (WhichSide(p[2], p[0], t[0]) >= 0.0) return true;
        else return false;
      else
        if (WhichSide(t[1], p[1], p[2]) >= 0.0) {
          if (WhichSide(t[1], t[2], p[2]) >= 0.0) return true;
          else return false;
        }
        else return false;
    else  return false;
}

/*!
\brief Test if two triangles intersect.
\param t %Triangle.

Triangles should have their vertices in trigonometric order.
*/
bool Triangle2::Intersect(const Triangle2& t) const
{
  if (WhichSide(p[0], t[0], t[1]) >= 0.0)
  {
    if (WhichSide(p[0], t[1], t[2]) >= 0.0)
    {
      if (WhichSide(p[0], t[2], t[0]) >= 0.0) return true;
      else return Triangle2(p[0], p[1], p[2]).IntersectEdge(Triangle2(t[0], t[1], t[2]));
    }
    else {
      if (WhichSide(p[0], t[2], t[0]) >= 0.0)
        return Triangle2(p[0], p[1], p[2]).IntersectEdge(Triangle2(t[2], t[0], t[1]));
      else return Intersect(Triangle2(t[0], t[1], t[2]));
    }
  }
  else {
    if (WhichSide(p[0], t[1], t[2]) >= 0.0)
    {
      if (WhichSide(p[0], t[2], t[0]) >= 0.0)
        return Triangle2(p[0], p[1], p[2]).IntersectEdge(Triangle2(t[1], t[2], t[0]));
      else return Intersect(Triangle2(t[1], t[2], t[0]));
    }
    else return Intersect(Triangle2(t[2], t[0], t[1]));
  }
}

/*!
\brief Test if two triangles overlap.
\param t %Triangle.
*/
bool Triangle2::Overlap(const Triangle2& t) const
{
  if (WhichSide(p[2], p[0], p[1]) < 0.0)
    if (WhichSide(t[2], t[0], t[1]) < 0.0)
      return Triangle2(p[0], p[2], p[1]).Intersect(Triangle2(t[0], t[2], t[1]));
    else
      return Triangle2(p[0], p[2], p[1]).Intersect(Triangle2(t[0], t[1], t[2]));
  else
    if (WhichSide(t[2], t[0], t[1]) < 0.0)
      return Triangle2(p[0], p[1], p[2]).Intersect(Triangle2(t[0], t[2], t[1]));
    else
      return Triangle2(p[0], p[1], p[2]).Intersect(Triangle2(t[0], t[1], t[2]));
}

/*!
\brief Compute a Poisson disc distribution inside the triangle.

This function uses a simple dart throwing algorithm.

\param r Radius of the discs.
\param n Number of candidate points.
\param random %Random number generator.
*/
QVector<Vector2> Triangle2::Poisson(const double& r, int n, Random& random) const
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

/*!
\brief Create an equilateral triangle.
\param c Center.
\param r Radius.
\param a Rotation angle.
*/
Triangle2 Triangle2::Equilateral(const Vector2& c, const double& r, const double& a)
{
  Triangle2 t(Vector2(r, 0.0), Vector2(-0.5 * r, 0.5 * Math::Sqrt3 * r), Vector2(-0.5 * r, -0.5 * Math::Sqrt3 * r));

  if (a != 0.0)
  {
    t.Rotate(Matrix2::Rotation(a));
  }
  t.Translate(c);
  return t;

}
