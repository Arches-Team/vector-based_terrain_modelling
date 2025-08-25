// Segment

#include "libs/segment.h"

#include "libs/ray.h"

/*!
\class Segment segment.h
\brief A segment defined by its end vertices.

Axis inherits Segment and extends the structure by storing the
normalized axis vector and its length.
%Lines are implemented by the class Line.

<P><I>How do I compute the intersection of two segments?</I>
<BR>In general, segments do not intersect, use the distance function.
\code
Segment a(Vector(0,1,-1),Vector(-1,0,3));
Segment b(Vector(0,1,-1),Vector(-1,0,3));
double r=a.R(b);
\endcode
Alternatively, it is possible to compute the intersection between two
segments in the plane, using Segment2::Intersect(const Segment2&).

\ingroup KernelGroup
*/

/*
\brief Create a segment from a segment in the plane.
\param s Plane segment.
\param u,v Elevations.
*/
Segment::Segment(const Segment2& s, const double& u, const double& v) :Segment(s.Vertex(0).ToVector(u), s.Vertex(1).ToVector(v))
{
}

/*!
\brief Test if two segments are almost equal.

This function computes the distance between the end vertices of the
two argument segments. If the norm infinity of this difference is
greater than the epsilon threshold value, then the segments are different.

\param s Other segment.
\param epsilon Precision.
This is a convenience function which is the same as:
\code
Segment a,b; // Two segments
double e; // Epsilon value
bool t=Vector::Equal(a.Vertex(0), b.Vertex(0), e) && Vector::Equal(a.Vertex(1), b.Vertex(1), e);
\endcode
*/
bool Segment::Equal(const Segment& s, const double& epsilon) const
{
  if (Vector::Equal(a, s.a, epsilon) && Vector::Equal(b, s.b, epsilon))
    return true;
  return false;
}

/*!
\brief Compute the bounding box of the segment.
*/
Box Segment::GetBox() const
{
  return Box(Vector::Min(a, b), Vector::Max(a, b));
}

/*!
\brief Rotates a segment.

Only the vertices and the axis are modified, the length is preserved.

\param r Rotation matrix.
*/
void Segment::Rotate(const Matrix& r)
{
  a = r * a;
  b = r * b;
}

/*!
\brief Rotates a segment.

\param r Rotation matrix.
*/
Segment Segment::Rotated(const Matrix& r) const
{
  return Segment(r * a, r * b);
}

/*!
\brief Inverse transformation.

\param f The frame.
*/
Segment Segment::InverseTransformed(const Frame& f) const
{
  return Segment(f.InverseTransform(a), f.InverseTransform(b));
}

/*!
\brief Translates a segment.

\param t Translation vector.
*/
void Segment::Translate(const Vector& t)
{
  a += t;
  b += t;
}

/*!
\brief Translates a segment.

\param t Translation vector.
*/
Segment Segment::Translated(const Vector& t) const
{
  return Segment(a + t, b + t);
}

/*!
\brief Uniformly scales a segment.

\param s Scaling factor.
*/
void Segment::Scale(const double& s)
{
  a *= s;
  b *= s;
}

/*!
\brief Scales a segment.

\param s Scaling vector.
*/
Segment Segment::Scaled(const Vector& s) const
{
  return Segment(a.Scaled(s), b.Scaled(s));
}

/*!
\brief Uniformly scales a segment.

\param s Scaling factor.
*/
Segment Segment::Scaled(const double& s) const
{
  return Segment(a * s, b * s);
}

/*!
\brief Overloaded.
\param s Stream.
\param segment The segment.
*/
std::ostream& operator<<(std::ostream& s, const Segment& segment)
{
  s << "Segment(" << segment.a << ',' << segment.b << ")";
  return s;
}

/*!
\brief Compute the squared distance to the segment.
\param p Point.
*/
double Segment::R(const Vector& p) const
{
  Vector n = p - a;

  // Do not normalize there
  Vector axis = b - a;
  double s = n * axis;

  // Squared distance 
  double r = n * n;

  // Not vertex a
  if (s > 0.0)
  {
    // Normalize only there : compute length and scale scalar product
    double length = Norm(axis);

    // Dividing s by length is more efficient than normalizing the axis
    s /= length;

    r -= s * s;

    // Vertex b
    if (s > length)
    {
      s -= length;
      r += s * s;
    }
  }
  return r;
}

/*!
\brief Compute the normal vector between a point and its projection onto the segment.
\param p Point.
*/
Vector Segment::Normal(const Vector& p) const
{
  Vector n = p - b;

  // Do not normalize there
  Vector axis = b - a;

  // Not that vertex
  if (axis * n <= 0.0)
  {
    n = p - a;
    double s = axis * n;
    // Not that one either
    if (s >= 0.0)
    {
      // Normalize only there : compute length and scale scalar product
      double length = axis * axis;

      // Dividing s by length is more efficient than normalizing the axis
      s /= length;

      n -= s * axis;
    }
  }
  return n;
}

/*!
\brief Compute the squared distance to the segment.
\param p Point.
\param s Linear coefficient corresponding to the projection of argument point p onto the segment.
*/
double Segment::R(const Vector& p, double& s) const
{
  Vector n = p - a;
  // Do not normalize there
  Vector axis = b - a;

  s = n * axis;

  // Squared distance 
  double r = n * n;

  // Not vertex a
  if (s > 0.0)
  {
    // Normalize only there : compute length and scale scalar product
    double length = Norm(axis);
    // Dividing s by length is more efficient than normalizing the axis
    s /= length;

    r -= s * s;

    // Vertex b
    if (s > length)
    {
      s -= length;
      r += s * s;
      s = 1.0;
    }
    else
    {
      s /= length;
    }
  }
  // Vertex a
  else
  {
    s = 0.0;
  }
  return r;
}

/*!
\brief Compute the squared distance between two segments.

\param segment %The segment.
*/
double Segment::R(const Segment& segment) const
{
  Vector u = b - a;
  Vector v = segment.b - segment.a;
  Vector w = a - segment.a;
  double uu = u * u;
  double uv = u * v;
  double vv = v * v;
  double uw = u * w;
  double vw = v * w;
  double e = uu * vv - uv * uv;

  double sN, sD = e;       // sc = sN / sD, default sD = D >= 0
  double tN, tD = e;       // tc = tN / tD, default tD = D >= 0

  // Compute the line parameters of the two closest points
  // Lines are almost parallel
  if (e < 1e-6)
  {
    // force using point P0 on segment S1 to prevent possible division by 0.0 later
    sN = 0.0;
    sD = 1.0;
    tN = vw;
    tD = vv;
  }
  // Get the closest points on the infinite lines
  else
  {
    sN = uv * vw - vv * uw;
    tN = uu * vw - uv * uw;
    if (sN < 0.0)
    {        // sc < 0 => the s=0 edge is visible
      sN = 0.0;
      tN = vw;
      tD = vv;
    }
    else if (sN > sD)
    {  // sc > 1  => the s=1 edge is visible
      sN = sD;
      tN = vw + uv;
      tD = vv;
    }
  }

  if (tN < 0.0)
  {            // tc < 0 => the t=0 edge is visible
    tN = 0.0;
    // recompute sc for this edge
    if (-uw < 0.0)
    {
      sN = 0.0;
    }
    else if (-uw > uu)
    {
      sN = sD;
    }
    else
    {
      sN = -uw;
      sD = uu;
    }
  }
  else if (tN > tD)
  {      // tc > 1  => the t=1 edge is visible
    tN = tD;
    // recompute sc for this edge
    if ((-uw + uv) < 0.0)
    {
      sN = 0;
    }
    else if ((-uw + uv) > uu)
    {
      sN = sD;
    }
    else
    {
      sN = (-uw + uv);
      sD = uu;
    }
  }

  double s = (fabs(sN) < 1e-6 ? 0.0 : sN / sD);
  double t = (fabs(tN) < 1e-6 ? 0.0 : tN / tD);

  // Difference of the two closest points
  Vector st = w + (s * u) - (t * v);

  return st * st;
}

/*!
\brief Computes the polynomial equation of the Euclidean distance between a
ray and the line corresponding to the segment.
\param ray The ray.
*/
Quadric Segment::Equation(const Ray& ray) const
{
  Vector axis = (b - a) / Norm(b - a);
  Vector pa = a - ray.Origin();

  double dx = ray.Direction() * axis;
  double pax = pa * axis;
  double dpa = ray.Direction() * pa;

  return Quadric(ray.Direction() * ray.Direction() - dx * dx, 2.0 * (dx * pax - dpa), pa * pa - pax * pax);
}

/*!
\brief Check if a segment intersects a box.
\param box The box.
*/
bool Segment::Intersect(const Box& box) const
{
  Vector ba = b - a;

  Vector d = 0.5 * ba;
  Vector c = 0.5 * (a + b);

  double fd[3];
  Vector cc = c - box.Center();

  fd[0] = fabs(d[0]);
  if (fabs(cc[0]) > ba[0] + fd[0])
    return false;

  fd[1] = fabs(d[1]);
  if (fabs(cc[1]) > ba[1] + fd[1])
    return false;

  fd[2] = fabs(d[2]);
  if (fabs(cc[2]) > ba[2] + fd[2])
    return false;

  if (fabs(d[1] * cc[2] - d[2] * cc[1]) > ba[1] * fd[2] + ba[2] * fd[1])
    return false;

  if (fabs(d[2] * cc[0] - d[0] * cc[2]) > ba[0] * fd[2] + ba[2] * fd[0])
    return false;

  if (fabs(d[0] * cc[1] - d[1] * cc[0]) > ba[0] * fd[1] + ba[1] * fd[0])
    return false;

  return true;
}

/*!
\brief Compute the polynomial equation of the distance function along the ray.
\param ray The ray.
\param a First vertex of the segment.
\param axis %Axis.
*/
Quadric Segment::EdgeEquation(const Ray& ray, const Vector& a, const Vector&, const Vector& axis)
{
  Vector pa = a - ray.Origin();

  double dx = ray.Direction() * axis;
  double pax = pa * axis;
  double dpa = ray.Direction() * pa;

  return Quadric(ray.Direction() * ray.Direction() - dx * dx, 2.0 * (dx * pax - dpa), pa * pa - pax * pax);
}

/*!
\brief Compute the intersection between the line f(x)=y and a line such that f(a)=va and f(b)=vb, on the segment <B>ab</B>.
\param a,b Points.
\param va, vb Values.
\param y Threshold.
*/
Vector Segment::Intersect(const Vector& a, const Vector& b, const double va, const double vb, const double y)
{
  double t = (y - va) / (vb - va);
  return a + t * (b - a);
}