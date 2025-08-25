// Cylinder

#include "libs/cylinder.h"
#include "libs/quadric.h"
#include "libs/circle.h"

const double Cylinder::epsilon = 1.0e-4;

const Cylinder Cylinder::Unit(Vector::Null, Vector::Z, 1.0);

/*!
\class Cylinder cylinder.h
\brief A cylinder characterized by its end vertices and radius.

Note that only the radius is stored, the squared radius value is
computed on the fly whenever needed.

\ingroup KernelGroup
*/

/*!
\brief Creates a cylinder.
\param a, b End vertices of the axis.
\param r Radius of the cylinder.
*/
Cylinder::Cylinder(const Vector& a, const Vector& b, const double& r) :Axis(a, b), r(r)
{
}

/*!
\brief Creates a vertical cylinder.
\param a, b Height of apexes.
\param r Radius.
*/
Cylinder::Cylinder(const double& a, const double& b, const double& r) :Cylinder(Vector(0.0, 0.0, a), Vector(0.0, 0.0, b), r)
{
}

/*!
\brief Creates a vertical cylinder.
\param c Base center vertex.
\param h Height.
\param r Radius.

Same as:
\code
Cylinder cylinder(c, c + Vector(0.0, 0.0, h), r);
\endcode
*/
Cylinder::Cylinder(const Vector& c, const double& h, const double& r) :Cylinder(c, c + Vector(0.0, 0.0, h), r)
{
}

/*!
\brief Create a vertical cylinder as an extruded circle.
\param c %Circle.
\param h Height.
\param symmetric Extrude in both directions, set to false by default.

Same as:
\code
Circle c;
double h;
  Cylinder cylinder(c.Center(), c.Center() + c.Axis() * h, c.Radius());
\endcode
*/
Cylinder::Cylinder(const Circle& c, const double& h, bool symmetric) :Cylinder((symmetric == false) ? c.Center() : c.Center() - c.Axis() * h, c.Center() + c.Axis() * h, c.Radius())
{
}

/*!
\brief Creates the boudning cylinder embedding a box rotated around an axis.
\param box The box.
\param ax The axis.
*/
Cylinder::Cylinder(const Box& box, const Axis& ax) : Axis(ax)
{
  double lengths[2];

  // Radial coordinates of the first vertex of box
  Vector2 ra = Radial(box.Vertex(0));

  // Radius
  r = ra[0];

  // Lengths
  lengths[0] = ra[1];
  lengths[1] = ra[1];

  // Incrementally create the bounding box
  for (int i = 1; i < 8; i++)
  {
    // Vertex of box
    Vector2 ri = Axis::Radial(box.Vertex(i));

    // Radius
    r = Math::Max(r, ri[0]);

    lengths[0] = Math::Min(lengths[0], ri[1]);
    lengths[1] = Math::Max(lengths[0], ri[1]);
  }

  // Keep normalized axis from input axis, but change vertexes and length
  a = ax.Vertex(lengths[0]);
  b = ax.Vertex(lengths[1]);
  length = lengths[1] - lengths[0];
}

/*!
\brief Squared distance between a point and the cylinder.

From A. Barbier and E. Galin. Fast distance
Computation between a Point and Cylinders, Cones, Line-swept-Spheres and
Cone-Spheres. <I>Journal of Graphics Tools</I>, <B>9</B>(2), 11-19, 2004.

Cost in worst case: 11 +, 9 *, 1 sqrt, 3 if.

\param p Point.

\code
double Cylinder::R(const Vector& p) const
{
Vector n = p - a;  // Cost: 3+
double s = axis*n; // Cost: 3* 2+
double ra = n*n;   // Cost: 3* 2+

// Vertex A
if (s < 0.0)
{
double z = s*s;    // Cost: 1*
double u = ra - z; // Cost: 1+
// Inside disc at vertex A
if (u < r*r)       // Cost: 1?
{
return z;
}
// Circle A
else
{
u = sqrt(u);    // Cost: 1 sqrt
u -= r;         // Cost: 1+
return z + u*u; // Cost: 1* 1+
}
}
else
{
// Vertex B
if (s > length)
{
double z = s - length; // 1+
z *= z;                // 1*
double u = ra - s*s;   // 1+ 1*
// Inside disc at vertex B
if (u < r*r)
{
return z;
}
// Circle A
else
{
u = sqrt(u);     // 1 sqrt
u -= r;          // 1+
return z + u*u;  // 1* 1+
}
}
// Cylinder
else
{
double z = s*s;    // 1*
double u = ra - z; // 1+
if (u < r*r)
{
return 0;
}
else
{
u = sqrt(u);   // 1 sqrt
u -= r;        // 1+
return u*u;    // 1*
}
}
}
}
\endcode
*/
double Cylinder::R(const Vector& p) const
{
  Vector n = p - a;  // Cost: 3+
  double s = axis * n; // Cost: 3* 2+
  double ra = n * n;   // Cost: 3* 2+

  // Vertex A
  if (s < 0.0)
  {
    double z = s * s;    // Cost: 1*
    double u = ra - z; // Cost: 1+
    // Inside disc at vertex A
    if (u < r * r)       // Cost: 1?
    {
      return z;
    }
    // Circle A
    else
    {
      u = sqrt(u);    // Cost: 1 sqrt
      u -= r;         // Cost: 1+
      return z + u * u; // Cost: 1* 1+
    }
  }
  else
  {
    // Vertex B
    if (s > length)
    {
      double z = s - length; // 1+
      z *= z;                // 1*
      double u = ra - s * s;   // 1+ 1*
      // Inside disc at vertex B
      if (u < r * r)
      {
        return z;
      }
      // Circle A
      else
      {
        u = sqrt(u);     // 1 sqrt
        u -= r;          // 1+
        return z + u * u;  // 1* 1+
      }
    }
    // Cylinder
    else
    {
      double z = s * s;    // 1*
      double u = ra - z; // 1+
      if (u < r * r)
      {
        return 0.0;
      }
      else
      {
        u = sqrt(u);   // 1 sqrt
        u -= r;        // 1+
        return u * u;    // 1*
      }
    }
  }
}

/*!
\brief Computes the squared distance between a point and the cylinder.

\param p Point.
\param s Returned position of the projected point onto the axis of the cylinder.

\sa Cylinder::R(const Vector&);
*/
double Cylinder::R(const Vector& p, double& s) const
{
  Vector n = p - a;   // Cost: 3 +
  s = axis * n;       // Cost: 3*2 +
  double ra = n * n;  // Cost: 3*2 +

  // Vertex A
  if (s < 0.0)
  {
    double z = s * s;  // Cost: 1 *
    double u = ra - z; // Cost: 1 +

    // Set s to 0.0 before return
    s = 0.0;

    // Inside disc at vertex A
    if (u < r * r)    // Cost: 1 ?
    {
      return z;
    }
    // Circle A
    else
    {
      u = sqrt(u);   // Cost: 1 sqrt
      u -= r;        // Cost: 1 +
      return z + u * u;// Cost: 1*1 +
    }
  }
  else
  {
    // Vertex B
    if (s > length)
    {
      double z = s - length;// 1+
      z *= z;// 1*
      double u = ra - s * s;// 1+ 1*

      // Set s to 1.0 before return
      s = 1.0;

      // Inside disc at vertex B
      if (u < r * r)
      {
        return z;
      }
      else
        // Circle A
      {
        u = sqrt(u);// 1sqrt
        u -= r;// 1+
        return z + u * u;// 1* 1+
      }
    }
    // Cylinder
    else
    {
      double z = s * s; // 1*
      double u = ra - z;// 1+

      // Normalize s to unit interval before return
      s /= length;

      if (u < r * r)
      {
        return 0;
      }
      else
      {
        u = sqrt(u);// 1sqrt
        u -= r; // 1+
        return u * u;// 1*
      }
    }
  }
}

/*!
\brief Compute the normal vector to a cylinder.

Project the point in space onto the axis of the cylinder. If projection
lies inside the segment, then return the normal to the segment minus
the right radius direction, otherwise perform computations for the
end caps.

Return Vector::Null if the point is inside the cylinder.
\param p Point.
*/
Vector Cylinder::Normal(const Vector& p) const
{
  Vector n = p - b;
  double s = axis * n;

  // Not that vertex
  if (s <= 0.0)
  {
    n = p - a;
    s = axis * n;
    // Not that one either
    if (s >= 0.0)
    {
      n -= s * axis;
      double u = n * n;
      // Inside cylinder
      if (u < r * r)
      {
        return Vector::Null;
      }
      // Return scaled vector
      else
      {
        u = sqrt(u);
        return n * ((1.0 - r / u));
      }
    }
  }

  // Vertex cases
  Vector t = s * axis;
  Vector q = n - t;
  double u = q * q;
  // Inside cylinder
  if (u < r * r)
  {
    return t;
  }
  // Circle
  else
  {
    u = sqrt(u);
    return t + q * (1.0 - r / u);
  }
}


/*!
\brief Compute the normal vector to the surface of the cylinder.

\param p Point.
*/
Vector Cylinder::SignedNormal(const Vector& p) const
{
  Vector n = p - b;
  double sb = axis * n;
  double s = sb;

  // Not that vertex
  if (sb <= 0.0)
  {
    n = p - a;
    double sa = axis * n;
    s = sa;
    // Not that one either
    if (sa >= 0.0)
    {
      n -= sa * axis;
      double u = n * n;
      // Inside cylinder
      if (u < r * r)
      {
        // Check caps, with sb < 0.0 and sa > 0.0
        Vector nc = axis;
        if (sa > -sb)
        {
          s = -sb;
          nc *= sb;
        }
        else
        {
          s = sa;
          nc *= sa;
        }
        // Radial distance
        u = sqrt(u);

        // Distance to cylinder
        if (-u < s)
        {
          n *= r / u - 1.0;
        }
        // Caps
        else
        {
          n = nc;
        }
        return n;
      }
      // Return scaled vector
      else
      {
        u = sqrt(u);
        return n * (1.0 - r / u);
      }
    }
  }

  // Vertex cases
  Vector t = s * axis;
  Vector q = n - t;
  double u = q * q;
  // Inside infinite cylinder
  if (u < r * r)
  {
    return t;
  }
  // Circle
  else
  {
    u = sqrt(u);
    return t + q * (1.0 - r / u);
  }
}

/*!
\brief Signed Euclidean distance between a point and a cylinder.

\param p Point.
*/
double Cylinder::Signed(const Vector& p) const
{
  Vector n = p - a;   // Cost: 3 +
  double s = axis * n;// Cost: 3*2 +
  double ra = n * n;  // Cost: 3*2 +

  // Vertex A
  if (s < 0.0)
  {
    double z = s * s;  // Cost: 1 *
    double u = ra - z; // Cost: 1 +

    // Inside disc at vertex A
    if (u < r * r)    // Cost: 1 ?
    {
      return sqrt(z); // 1 Sqrt
    }
    // Circle A
    else
    {
      u = sqrt(u);   // Cost: 1 sqrt
      u -= r;     // Cost: 1 +
      return sqrt(z + u * u);// Cost: 1*1 +  1 Sqrt
    }
  }
  else
  {
    // Vertex B
    if (s > length)
    {
      double z = s - length;// 1+
      z *= z;// 1*
      double u = ra - s * s;// 1+ 1*

      // Inside disc at vertex B
      if (u < r * r)
      {
        return sqrt(z); //  1 Sqrt
      }
      else
        // Circle A
      {
        u = sqrt(u);// 1sqrt
        u -= r;// 1+
        return sqrt(z + u * u);// 1* 1+  1 Sqrt
      }
    }
    // Cylinder
    else
    {
      double z = s * s; // 1*
      double u = ra - z;// 1+

      u = sqrt(u) - r;// 1 Sqrt 1+

      // Inside
      if (u < 0.0)
      {
        return Math::Max(u, -s, s - length);
      }
      else
      {
        return u;
      }
    }
  }
}

/*!
\brief Check if a point is inside or outside the cylinder.
\param p Point.
*/
bool Cylinder::Inside(const Vector& p) const
{
  Vector n = p - b;
  // Not that vertex
  if (axis * n <= 0.0)
  {
    n = p - a;
    double s = axis * n;
    // Not that one either
    if (s >= 0.0)
    {
      n -= s * axis;
      if (n * n < r * r)
      {
        return true;
      }
    }
  }
  return false;
}

/*!
\brief Compute the intersection between a cylinder and a ray.

Intersections are sorted according to intersection depth.

\param ray The ray.
\param ta, tb Intersection depths.
*/
int Cylinder::Intersect(const Ray& ray, double& ta, double& tb) const
{
  // Origin to base vertex
  Vector pa = a - ray.Origin();
  Vector pb = b - ray.Origin();

  double dx = ray.Direction() * axis;

  double pax = pa * axis;
  double pbx = pb * axis;

  double u[2];
  u[0] = pax / dx;
  u[1] = pbx / dx;

  Math::Sort(u[0], u[1]);

  double dpa = ray.Direction() * pa;

  Quadric q(ray.Direction() * ray.Direction() - dx * dx, 2.0 * (dx * pax - dpa), pa * pa - pax * pax - r * r);

  double v[2];
  // Escape if does not intersect infinite cylinder
  if (!q.Solve(v[0], v[1]))
    return 0;

  // Intersection
  if (u[0] > v[0])
  {
    ta = u[0];
  }
  else
  {
    ta = v[0];
  }

  if (u[1] < v[1])
  {
    tb = u[1];
  }
  else
  {
    tb = v[1];
  }

  // No intersections at all
  if (ta > tb)
    return 0;

  return 2;
}

/*!
\brief Compute the intersections between a cylinder and a ray.

Intersections are sorted according to intersection depth.

\param ray The ray.
\param ta, tb Intersection depths.
\param na, nb Normals at intersection points.
*/
int Cylinder::Intersect(const Ray& ray, double& ta, double& tb, Vector& na, Vector& nb) const
{
  // Origin to base vertex
  Vector pa = a - ray.Origin();
  Vector pb = b - ray.Origin();

  double dx = ray.Direction() * axis;

  double pax = pa * axis;
  double pbx = pb * axis;

  double u[2];
  u[0] = pax / dx;
  u[1] = pbx / dx;
  if (u[1] < u[0])
  {
    double temp = u[0];
    u[0] = u[1];
    u[1] = temp;
    na = axis;
    nb = -axis;
  }
  else
  {
    na = -axis;
    nb = axis;
  }

  double dpa = ray.Direction() * pa;

  Quadric q(ray.Direction() * ray.Direction() - dx * dx, 2.0 * (dx * pax - dpa), pa * pa - pax * pax - r * r);

  double v[2];
  // Escape if does not intersect infinite cylinder
  if (!q.Solve(v[0], v[1]))
    return 0;

  // Intersection
  if (u[0] > v[0])
  {
    ta = u[0];
  }
  else
  {
    ta = v[0];
    na = ray(ta) - a;
    na -= (axis * na) * axis;
  }

  if (u[1] < v[1])
  {
    tb = u[1];
  }
  else
  {
    tb = v[1];
    nb = ray(tb) - a;
    nb -= (axis * nb) * axis;
  }

  // No intersections at all
  if (ta > tb)
    return 0;

  return 2;
}

/*!
\brief Check the intersection between a cylinder and a ray.
\param ray The ray.
*/
bool Cylinder::Intersect(const Ray& ray) const
{
  double ta, tb;

  if (Cylinder::Intersect(ray, ta, tb))
  {
    if (ta < 0.0)
    {
      return false;
    }
    return true;
  }
  return false;
}

/*!
\brief Compute the first positive intersection between a cylinder and a ray.
\param ray The ray.
\param t Intersection depth.
\param n Normal at intersection point.
*/
int Cylinder::Intersect(const Ray& ray, double& t, Vector& n) const
{
  double ta, tb;
  Vector na, nb;

  if (Cylinder::Intersect(ray, ta, tb, na, nb))
  {
    if (ta > 0.0)
    {
      t = ta;
      n = na;
      return 1;
    }
    else if (tb > 0.0)
    {
      t = tb;
      n = nb;
      return 1;
    }
  }
  return 0;
}

/*!
\brief Computes the axis-aligned bounding box of a cylinder.

\sa Circle::GetBox()
*/
Box Cylinder::GetBox() const
{
  Vector s = Axis::BoxVector(axis);

  return Box(Vector::Min(a, b) - r * s, Vector::Max(a, b) + r * s);
}

/*!
\brief Computes the inverse mapping of a point.

Projects the point onto the z-axis unit cylinder, uv-coordinates are set to unit interval [0,1].
\param p Point.
\return Inverse mapping coordinates.
*/
Vector2 Cylinder::Cast(const Vector& p)
{
  double u = 0.0;
  // Radius
  double d = sqrt(p[0] * p[0] + p[1] * p[1]);
  if (d != 0.0)
  {
    if (p[1] == 0.0)
    {
      if (p[0] > 0.0)
        u = 0.0;
      else
        u = Math::Pi;
    }
    else
    {
      u = acos(p[0] / d);
      if (p[1] < 0.0)
      {
        u = Math::TwoPi - u;
      }
    }
    u /= Math::TwoPi;
  }
  return Vector2(u, Math::Mod(p[2], 1.0));
}

/*!
\brief Uniformly scales a cylinder.

\param s Scaling factor.
*/
void Cylinder::Scale(const double& s)
{
  Axis::Scale(s);
  r *= fabs(s);
}

/*!
\brief Transformed cylinder.

\param frame %Frame.
*/
Cylinder Cylinder::Transformed(const FrameScaled& frame) const
{
  return Cylinder(frame.Transform(a), frame.Transform(b), Abs(frame.S()).Max());
}

/*!
\brief Overloaded.
\param s Stream.
\param cylinder The cylinder.
*/
std::ostream& operator<<(std::ostream& s, const Cylinder& cylinder)
{
  s << "Cylinder(" << cylinder.a << ',' << cylinder.b << ',' << cylinder.r << ')';
  return s;
}

/*!
\brief Generate a random vector inside the cylinder.
\param random %Random number generator.
*/
Vector Cylinder::RandomInside(Random& random) const
{
  // Size of the embedding cube
  double c = Math::Max(2.0 * r, length);

  // Perform selection of parameters
  Vector p;
  while (true)
  {
    p = c * Vector(random.Uniform(), random.Uniform(), random.Uniform()) - Vector(c / 2.0, c / 2.0, 0.0);
    if ((SquaredNorm(Vector2(p)) < r * r) && (p[2] < length)) break;
  }

  Vector x, y;
  axis.Orthonormal(x, y);

  return a + p[0] * x + p[1] * y + p[2] * axis;
}

/*!
\brief Return a cylinder rotated by a given matrix.
\param r Rotation matrix.
*/
Cylinder Cylinder::Rotated(const Matrix& r) const
{
  return Cylinder(r * a, r * b, Cylinder::r);
}

/*!
\brief Return a cylinder translated by a given vector.
\param t Translation vector.
*/
Cylinder Cylinder::Translated(const Vector& t) const
{
  return Cylinder(a + t, b + t, r);
}

/*!
\brief Return a scaled cylinder.
\param s Scaling factor.
*/
Cylinder Cylinder::Scaled(const double& s) const
{
  return Cylinder(s * a, s * b, s * r);
}