// Cone

#include "libs/cone.h"
#include "libs/quadric.h"

/*!
\class Cone cone.h
\brief Cones defined as truncated cones.

\ingroup KernelGroup

*/

const double Cone::epsilon = 1.0e-4;

const Cone Cone::Unit(1.0, 1.0, 0.0);

/*!
\brief Creates a cone characterized by its end-vertices and radii.

\param a, b End vertices of the axis.
\param ra, rb Radii.
*/
Cone::Cone(const Vector& a, const Vector& b, const double& ra, const double& rb) :Axis(a, b), ra(ra), rb(rb)
{
  // Normalize
  Cone::rlength = (rb - ra) / length;

  // Compute the length of the side of cone, i.e. its slant height
  Cone::conelength = sqrt((rb - ra) * (rb - ra) + length * length);

  // Line segment
  Cone::side = Vector2(rb - ra, length);
  Cone::side /= conelength;
}

/*!
\brief Creates a vertical cone characterized height and radii.

\param ab Height.
\param ra, rb Radii.
*/
Cone::Cone(const double& ab, const double& ra, const double& rb) :Cone(Vector::Null, Vector(0.0, 0.0, ab), ra, rb)
{
}

/*!
\brief Check if a point is inside or outside the cone.
\param p Point.
*/
bool Cone::Inside(const Vector& p) const
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
      double r = ra + rlength * s;
      if (n * n < r * r)
      {
        return true;
      }
    }
  }
  return false;
}

/*!
\brief Compute the intersections between the cone and a ray.

\param ray The ray.
\param ta, tb Intersection depths.
*/
int Cone::Intersect(const Ray& ray, double& ta, double& tb) const
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

  Quadric d(ray.Direction() * ray.Direction() - dx * dx, 2.0 * (dx * pax - dpa), pa * pa - pax * pax);
  Quadric r(rlength * rlength * dx * dx, 2.0 * (ra * rlength * dx + rlength * rlength * pax * dx), ra * ra + rlength * rlength * pax * pax + 2.0 * ra * rlength * pax);

  Quadric q = d - r;
  double v[2];

  // Escape if does not intersect infinite cone
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
\brief Compute the intersections between the cone and a ray.

\param ray The ray.
\param ta, tb Intersection depths.
\param na, nb Normals at intersection points.
*/
int Cone::Intersect(const Ray& ray, double& ta, double& tb, Vector& na, Vector& nb) const
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

  Quadric d(ray.Direction() * ray.Direction() - dx * dx, 2.0 * (dx * pax - dpa), pa * pa - pax * pax);
  Quadric r(rlength * rlength * dx * dx, 2.0 * (ra * rlength * dx + rlength * rlength * pax * dx), ra * ra + rlength * rlength * pax * pax + 2.0 * ra * rlength * pax);

  Quadric q = d - r;
  double v[2];

  // Escape if does not intersect infinite cone
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
    na = Normal(ray(ta));
  }

  if (u[1] < v[1])
  {
    tb = u[1];
  }
  else
  {
    tb = v[1];
    nb = Normal(ray(tb));
  }

  // No intersections at all
  if (ta > tb)
    return 0;

  return 2;
}

/*!
\brief Compute the first positive intersection between the cone and a ray.
\param ray The ray.
\param t Intersection depth.
\param n Normal at intersection point.
*/
int Cone::Intersect(const Ray& ray, double& t, Vector& n) const
{
  double ta, tb;
  Vector na, nb;

  if (Cone::Intersect(ray, ta, tb, na, nb))
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
\brief Compute the axis-aligned bounding box of the cone.
*/
Box Cone::GetBox() const
{
  Vector s = Axis::BoxVector(axis);

  return Box(Box(a - ra * s, a + ra * s), Box(b - rb * s, b + rb * s));
}

/*!
\brief Compute the total surface area of the cone.

The total surface area is the sum of the area of the cap
discs, and the lateral surface area.
*/
double Cone::Area() const
{
  return Math::Pi * (ra * (ra + conelength) + rb * (rb + conelength));
}

/*!
\brief Rotates the cone.

\param r Rotation matrix.
*/
void Cone::Rotate(const Matrix& r)
{
  Axis::Rotate(r);
}

/*!
\brief Translates the cone.

\param t Translation vector.
*/
void Cone::Translate(const Vector& t)
{
  Axis::Translate(t);
}

/*!
\brief Uniformly scales the cone.

\param s Scaling factor.
*/
void Cone::Scale(const double& s)
{
  Axis::Scale(s);

  ra *= fabs(s);
  rb *= fabs(s);

  // Normalize
  Cone::rlength = (rb - ra) / length;

  // Compute the length of side of cone, i.e. its slant height
  Cone::conelength = sqrt((rb - ra) * (rb - ra) + length * length);
}

/*!
\brief Computes the vector distance between the cone and a point.

\param p Point.
*/
Vector Cone::Normal(const Vector& p) const
{
  // Compute coordinates 
  // Where y derives from scalar product
  Vector n = (p - a);
  double y = n * axis;

  // And x is derived from vector algebra
  Vector j = y * axis;
  Vector i = n - j;
  double x = Norm(i);

  if (ra > rb)
  {
    // Small cylinder test 
    if (x < rb)
    {
      // Cap
      if (y < 0.0)
      {
        return j;
      }
      // Other
      else if (y > length)
      {
        return j - (b - a);
      }
      // Inside cone
      else
      {
        return Vector::Null;
      }
    }

    // Change frame, so that point is now on large cap
    i *= (1.0 - ra / x);
    x -= ra;

    // At this point,
    // Distance to large cap
    if (y < 0.0)
    {
      // Distance to plane cap
      if (x < 0.0)
      {
        return j;
      }
      // Distance to plane circle
      else
      {
        return i + j;
      }
    }

    // Coordinates of p
    Vector u = i + j;

    // At this step of the algorithm, reference point is on the large cap
    // Need to use another reference frame

    // Normalize vector i
    Vector ip = i / x;

    // Rotate vector j
    Vector jp = (b - a) - (ra - rb) * ip;

    // Rotate vector i
    ip = length * ip + (ra - rb) * axis;

    // Note that although ip and jp are orthogonal,
    // they are not unit, and their squared norm is (ra-rb)^2+length^2

    y = u * jp;
    // Large cap
    if (y < 0.0)
    {
      return u;
    }
    else
    {
      // Top cap
      if (y > conelength * conelength)
      {
        return u - jp;
      }
      // Cone 
      else
      {
        x = u * ip;
        if (x > 0.0)
        {
          return (x / (conelength * conelength)) * ip;
        }
        else
        {
          return Vector::Null;
        }
      }
    }
  }
  else
  {
    // Small cylinder test 
    if (x < ra)
    {
      // Cap
      if (y < 0.0)
      {
        return j;
      }
      // Other
      else if (y > length)
      {
        return j - (b - a);
      }
      // Inside cone
      else
      {
        return Vector::Null;
      }
    }

    // Change frame, so that point is now on large cap
    i *= (1.0 - ra / x);
    x -= ra;

    // Distance to large cap
    if (y > length)
    {
      // Distance to plane cap
      if (x < rb - ra)
      {
        return j - (b - a);
      }
      // Distance to plane circle
      else
      {
        return i * (1.0 - (rb - ra) / x) + j - (b - a);
      }
    }

    // Coordinates of p
    Vector u = i + j;

    // At this step of the algorithm, reference point is on the large cap
    // Need to use another reference frame

    // Normalize vector i
    Vector ip = i / x;

    // Rotate vector j
    Vector jp = (b - a) - (ra - rb) * ip;

    // Rotate vector i
    ip = length * ip + (ra - rb) * axis;

    // Note that although ip and jp are orthogonal,
    // they are not unit, and their squared norm is (ra-rb)^2+length^2

    y = u * jp;
    // Large cap
    if (y < 0.0)
    {
      return u;
    }
    else
    {
      // Top cap
      if (y > conelength * conelength)
      {
        return u - jp;
      }
      // Cone 
      else
      {
        x = u * ip;
        if (x > 0.0)
        {
          return ((u * ip) / (conelength * conelength)) * ip;
        }
        else
        {
          return Vector::Null;
        }
      }
    }
  }
}

/*!
\brief Compute the distance between a point and the cone.

Details appear in A. Barbier and E. Galin. Fast distance Computation between a Point and Cylinders, Cones, Line-swept-Spheres and Cone-Spheres.
Journal of Graphics Tools. <b>9</b>(2), 11-19 (2004).

\param p Point.
*/
double Cone::R(const Vector& p) const
{
  // Compute coordinates 
  // Where y derives from scalar product
  Vector n = p - a;      // Cost: 3 +
  double y = n * axis;   // Cost: 3*2+

  // And x is derived from vector algebra
  double x = sqrt(n * n - y * y);// Norm(i);  // Cost: 4* 3+ 1 sqrt
  double e = 0.0;
  if (ra > rb)
  {
    // Small cylinder test 
    if (x < rb)        // Cost: 1 ?
    {
      // Cap
      if (y < 0.0)     // Cost: 1 ?
      {
        e = y * y;       // Cost: 1 *
      }
      // Other
      else if (y > length) // Cost: 1 ? 
      {
        e = y - length;      // Cost: 1 -
        e *= e;            // Cost: 1 *
      }
      // Inside cone
      //else
      //{
      //  e=0.0;
      //}
      return e; // Total for those cases: 
    }
    else
    {
      // Change frame, so that point is now on large cap
      x -= ra;   // Cost: 1 -

      // At this point,
      // Distance to large cap
      if (y < 0.0)  // Cost: 1 ?
      {
        // Distance to cap
        if (x < 0.0) // Cost: 1 ?
        {
          e = y * y; // Cost: 1 *
        }
        // Distance to circle
        else
        {
          e = x * x + y * y; // Cost: 1+2 *
        }
      }
      else
      {
        double ry = x * side[0] + y * side[1]; // Cost : 1+2 *
        double rx = x * side[1] - y * side[0]; // Cost : 1+2 *

        if (ry < 0.0) // Cost: 1 ?
        {
          e = x * x + y * y; // Cost: 1+ 2 *
        }
        else if (ry > conelength) // Cost: 1 ?
        {
          ry -= conelength; // Cost: 1 +
          e = rx * rx + ry * ry;  // Cost: 1+2 *
        }
        else
        {
          if (rx > 0.0)  // Cost: 1 ?
          {
            e = rx * rx;  // Cost: 1 *
          }
          //else
          //{
          //  e=0.0;
          //}
        }
      }
    }
  }
  else
  {
    // Small cylinder test 
    if (x < ra)
    {
      // Cap
      if (y < 0.0)
      {
        e = y * y;
      }
      // Other
      else if (y > length)
      {
        e = y - length;
        e *= e;
      }
      // Inside cone
      else
      {
        e = 0.0;
      }
      return e;
    }
    else
    {
      // Change frame, so that point is now on large cap
      x -= rb;
      y -= length;

      // Distance to large cap
      if (y > 0.0)
      {
        // Distance to cap
        if (x < 0.0)
        {
          e = y * y;
        }
        // Distance to circle
        else
        {
          e = x * x + y * y;
        }
      }
      else
      {
        double ry = x * side[0] + y * side[1];
        double rx = x * side[1] - y * side[0];

        if (ry > 0.0)
        {
          e = x * x + y * y;
        }
        else if (ry < -conelength)
        {
          ry += conelength;
          e = rx * rx + ry * ry;
        }
        else
        {
          if (rx < 0.0)
          {
            e = 0.0;
          }
          else
          {
            e = rx * rx;
          }
        }
      }
    }
  }
  return e;
}

/*!
\brief Compute the signed distance between a point and the cone.
\param p Point.
*/
double Cone::Signed(const Vector& p) const
{
  // Compute coordinates 
  // Where y derives from scalar product
  Vector n = p - a;
  double y = n * axis;

  // And x is derived from vector algebra
  double x = sqrt(n * n - y * y);

  if (ra > rb)
  {
    // Small cylinder test 
    if (x < rb)
    {
      // Cap
      if (y < 0.0)
      {
        return -y;
      }
      // Other
      else if (y > length)
      {
        return y - length;
      }
      // Inside cone
      else
      {
        x -= ra;
        double rx = x * side[1] - y * side[0];
        return Math::Max(-y, y - length, rx);
      }
    }
    else
    {
      // Change frame, so that point is now on large cap
      x -= ra;

      // Distance to large cap
      if (y < 0.0)
      {
        // Distance to cap
        if (x < 0.0)
        {
          return -y;
        }
        // Distance to circle
        else
        {
          return sqrt(x * x + y * y);
        }
      }
      else
      {
        double ry = x * side[0] + y * side[1];
        double rx = x * side[1] - y * side[0];

        if (ry < 0.0)
        {
          return sqrt(x * x + y * y);
        }
        else if (ry > conelength)
        {
          ry -= conelength;
          return sqrt(rx * rx + ry * ry);
        }
        else
        {
          if (rx > 0.0)
          {
            return rx;
          }
          else
          {
            // Inside
            return Math::Max(-y, y - length, rx);
          }
        }
      }
    }
  }
  else
  {
    // Small cylinder test 
    if (x < ra)
    {
      // Cap
      if (y < 0.0)
      {
        return -y;
      }
      // Other
      else if (y > length)
      {
        return y - length;
      }
      // Inside cone
      else
      {
        x -= rb;
        y -= length;
        double rx = x * side[1] - y * side[0];
        return Math::Max(-y, y - length, rx);
      }
    }
    else
    {
      // Change frame, so that point is now on large cap
      x -= rb;
      y -= length;

      // At this point,
      // Distance to large cap
      if (y > 0.0)
      {
        // Distance to cap
        if (x < 0.0)
        {
          return y;
        }
        // Distance to circle
        else
        {
          return sqrt(x * x + y * y);
        }
      }
      else
      {
        double ry = x * side[0] + y * side[1];
        double rx = x * side[1] - y * side[0];

        if (ry > 0.0)
        {
          return sqrt(x * x + y * y);
        }
        else if (ry < -conelength)
        {
          ry += conelength;
          return sqrt(rx * rx + ry * ry);
        }
        else
        {
          if (rx < 0.0)
          {
            // Inside
            return Math::Max(-y, y - length, rx);
          }
          else
          {
            return rx;
          }
        }
      }
    }
  }
}

/*!
\brief Overloaded.
\param s Stream.
\param cone The cone.
*/
std::ostream& operator<<(std::ostream& s, const Cone& cone)
{
  s << "Cone(" << cone.a << ',' << cone.b << ',' << cone.ra << ',' << cone.rb << ')';
  return s;
}

/*!
\brief Compute the volume of the cone.
*/
double Cone::Volume() const
{
  return Math::Pi * (ra * ra + ra * rb + rb * rb) * length / 3.0;
}

/*!
\brief Generate a random vector inside the cone.
\param random %Random number generator.
*/
Vector Cone::RandomInside(Random& random) const
{
  Box box = GetBox();

  while (true)
  {
    Vector p = box.RandomInside(random);
    if (Inside(p))
      return p;
  }

  // Should never happen
  return Vector::Null;
}

/*!
\brief Generate a random direction.
\param axis Axis.
\param alpha Angle.
\param random %Random number generator.
*/
Vector Cone::RandomDirection(const Vector& axis, const double& alpha, Random& random)
{
  // Make an orthonormal basis
  Vector x, y;
  axis.Orthonormal(x, y);

  const Vector& z = axis;

  // Pick randomly around (0, 0, 1)
  double rz = random.Uniform(cos(alpha), 1);
  double r = sqrt(1.0 - rz * rz);

  // Angle
  double t = random.Uniform(-Math::Pi, Math::Pi);

  double rx = r * cos(t);
  double ry = r * sin(t);

  return rx * x + ry * y + rz * z;
}
