// Capsule

#include "libs/capsule.h"
#include "libs/sphere.h"

/*!
\class Capsule capsule.h
\brief A line-swept-sphere.

\image html capsule.png

Capsules, which are line-swept-spheres, make good bounding volumes. The function
Intersect(const Capsule&, const Capsule&) implements an efficient algorithm to
detect the intersection between two capsules.

Note that capsule-ray intersection member function assumes that the ray direction
is normalized.

\ingroup KernelGroup

*/

const Capsule Unit(-Vector::Z, Vector::Z, 1.0);

/*!
\brief Creates a capsule.
\param a, b Vertexes
\param r Radius.
*/
Capsule::Capsule(const Vector& a, const Vector& b, const double& r) :Cylinder(a, b, r)
{
}

/*!
\brief Computes the tight bounding box of the capsule.

\sa Circle::GetBox()
*/
Box Capsule::GetBox() const
{
  return Box(Vector::Min(a, b), Vector::Max(a, b)).Extended(r);
}

/*!
\brief Compute the volume.
*/
double Capsule::Volume() const
{
  return Math::Pi * r * r * (length + 4.0 / 3.0 * r);
}

/*!
\brief Compute the surface area.
*/
double Capsule::Area() const
{
  return Math::TwoPi * r * (length + 2.0 * r);
}

/*!
\brief Computes the squared distance between a point and the capsule.
\param p Point.
*/
double Capsule::R(const Vector& p) const
{
  Vector n = p - a;
  double s = axis * n;
  double ra = n * n;

  // Vertex A
  if (s < 0.0)
  {
    // Inside sphere
    if (ra < r * r)
    {
      return 0.0;
    }
    else
    {
      ra = sqrt(ra);
      ra -= r;
      return ra * ra;
    }
  }
  else
  {
    // Vertex B
    if (s > length)
    {
      double y = s - length;
      ra += y * y - s * s;
      // Inside sphere
      if (ra < r * r)
      {
        return 0.0;
      }
      else
      {
        ra = sqrt(ra);
        ra -= r;
        return ra * ra;
      }
    }
    // Cylinder
    else
    {
      ra -= s * s;
      if (ra < r * r)
      {
        return 0.0;
      }
      else
      {
        ra = sqrt(ra);
        ra -= r;
        return ra * ra;
      }
    }
  }
}

/*!
\brief Computes the signed distance between a point and the capsule.
\param p Point.
*/
double Capsule::Signed(const Vector& p) const
{
  Vector n = p - a;
  double s = axis * n;
  double ra = n * n;

  // Vertex A
  if (s < 0.0)
  {
  }
  // Vertex B
  else if (s > length)
  {
    double y = s - length;
    ra += y * y - s * s;
  }
  // Cylinder
  else
  {
    // Bug tracking: rounding leaks might occur which means that ra could become ~ <0.0
    ra -= s * s;
  }

  // Computing fabs avoids rounding leaks (see above)
  ra = sqrt(fabs(ra));
  ra -= r;
  return ra;
}

/*!
\brief Check if a point is inside the capsule.
\param p Point.
*/
bool Capsule::Inside(const Vector& p) const
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
    else
    {
      // Half sphere
      if (n * n < r * r)
      {
        return true;
      }
    }
  }
  else
  {
    // Half sphere
    if (n * n < r * r)
    {
      return true;
    }
  }
  return false;
}

/*!
\brief Check if two capsules intersect.

The algorithm computes the squared distance and compares it with the squared sum of radii.
\param cs The other capsule.
\param u,v The location of the closests points on the two capsules where minimum distance occurs.
\sa Capsule::R()
*/
bool Capsule::Intersect(const Capsule& cs, double& u, double& v) const
{
  double d = R(cs, u, v);
  double e = r + cs.Radius();
  e *= e;
  if (d > e)
  {
    return false;
  }
  else
  {
    return true;
  }
}

/*!
\brief Computes the squared distance between two capsules.
\param cs The line-swept-sphere.
\param u,v The location of the closests points on the segments of the capsules where minimum distance occurs.
*/
double Capsule::R(const Capsule& cs, double& u, double& v) const
{
  Vector aa = a - cs.a;
  const double a00 = 1.0;
  double a01 = -axis * cs.axis;
  const double a11 = 1.0;
  double aaxa = aa * axis;
  double aaxb = -aa * cs.axis;
  double ra = aa * aa;
  double fDet = fabs(a00 * a11 - a01 * a01);
  double e, t;

  if (fDet >= Cylinder::epsilon)
  {
    // segments are not parallel
    u = a01 * aaxb - a11 * aaxa;
    v = a01 * aaxa - a00 * aaxb;

    if (u >= 0.0)
    {
      if (u <= fDet)
      {
        if (v >= 0.0)
        {
          if (v <= fDet)  // region 0 (interior)
          {
            // minimum at two interior points of 3D lines
            double fInvDet = (1.0) / fDet;
            u *= fInvDet;
            v *= fInvDet;
            e = u * (a00 * u + a01 * v + (2.0) * aaxa) + v * (a01 * u + a11 * v + (2.0) * aaxb) + ra;
          }
          else  // region 3 (side)
          {
            v = 1.0;
            t = a01 + aaxa;
            if (t >= 0.0)
            {
              u = 0.0;
              e = a11 + (2.0) * aaxb + ra;
            }
            else if (-t >= a00)
            {
              u = 1.0;
              e = a00 + a11 + ra + (2.0) * (aaxb + t);
            }
            else
            {
              u = -t / a00;
              e = t * u + a11 + (2.0) * aaxb + ra;
            }
          }
        }
        else  // region 7 (side)
        {
          v = 0.0;
          if (aaxa >= 0.0)
          {
            u = 0.0;
            e = ra;
          }
          else if (-aaxa >= a00)
          {
            u = 1.0;
            e = a00 + (2.0) * aaxa + ra;
          }
          else
          {
            u = -aaxa / a00;
            e = aaxa * u + ra;
          }
        }
      }
      else
      {
        if (v >= 0.0)
        {
          if (v <= fDet)  // region 1 (side)
          {
            u = 1.0;
            t = a01 + aaxb;
            if (t >= 0.0)
            {
              v = 0.0;
              e = a00 + (2.0) * aaxa + ra;
            }
            else if (-t >= a11)
            {
              v = 1.0;
              e = a00 + a11 + ra + (2.0) * (aaxa + t);
            }
            else
            {
              v = -t / a11;
              e = t * v + a00 + (2.0) * aaxa + ra;
            }
          }
          else  // region 2 (corner)
          {
            t = a01 + aaxa;
            if (-t <= a00)
            {
              v = 1.0;
              if (t >= 0.0)
              {
                u = 0.0;
                e = a11 + (2.0) * aaxb + ra;
              }
              else
              {
                u = -t / a00;
                e = t * u + a11 + (2.0) * aaxb + ra;
              }
            }
            else
            {
              u = 1.0;
              t = a01 + aaxb;
              if (t >= 0.0)
              {
                v = 0.0;
                e = a00 + (2.0) * aaxa + ra;
              }
              else if (-t >= a11)
              {
                v = 1.0;
                e = a00 + a11 + ra + (2.0) * (aaxa + t);
              }
              else
              {
                v = -t / a11;
                e = t * v + a00 + (2.0) * aaxa + ra;
              }
            }
          }
        }
        else  // region 8 (corner)
        {
          if (-aaxa < a00)
          {
            v = 0.0;
            if (aaxa >= 0.0)
            {
              u = 0.0;
              e = ra;
            }
            else
            {
              u = -aaxa / a00;
              e = aaxa * u + ra;
            }
          }
          else
          {
            u = 1.0;
            t = a01 + aaxb;
            if (t >= 0.0)
            {
              v = 0.0;
              e = a00 + (2.0) * aaxa + ra;
            }
            else if (-t >= a11)
            {
              v = 1.0;
              e = a00 + a11 + ra + (2.0) * (aaxa + t);
            }
            else
            {
              v = -t / a11;
              e = t * v + a00 + (2.0) * aaxa + ra;
            }
          }
        }
      }
    }
    else
    {
      if (v >= 0.0)
      {
        if (v <= fDet)  // region 5 (side)
        {
          u = 0.0;
          if (aaxb >= 0.0)
          {
            v = 0.0;
            e = ra;
          }
          else if (-aaxb >= a11)
          {
            v = 1.0;
            e = a11 + (2.0) * aaxb + ra;
          }
          else
          {
            v = -aaxb / a11;
            e = aaxb * v + ra;
          }
        }
        else  // region 4 (corner)
        {
          t = a01 + aaxa;
          if (t < 0.0)
          {
            v = 1.0;
            if (-t >= a00)
            {
              u = 1.0;
              e = a00 + a11 + ra + (2.0) * (aaxb + t);
            }
            else
            {
              u = -t / a00;
              e = t * u + a11 + (2.0) * aaxb + ra;
            }
          }
          else
          {
            u = 0.0;
            if (aaxb >= 0.0)
            {
              v = 0.0;
              e = ra;
            }
            else if (-aaxb >= a11)
            {
              v = 1.0;
              e = a11 + (2.0) * aaxb + ra;
            }
            else
            {
              v = -aaxb / a11;
              e = aaxb * v + ra;
            }
          }
        }
      }
      else   // region 6 (corner)
      {
        if (aaxa < 0.0)
        {
          v = 0.0;
          if (-aaxa >= a00)
          {
            u = 1.0;
            e = a00 + (2.0) * aaxa + ra;
          }
          else
          {
            u = -aaxa / a00;
            e = aaxa * u + ra;
          }
        }
        else
        {
          u = 0.0;
          if (aaxb >= 0.0)
          {
            v = 0.0;
            e = ra;
          }
          else if (-aaxb >= a11)
          {
            v = 1.0;
            e = a11 + (2.0) * aaxb + ra;
          }
          else
          {
            v = -aaxb / a11;
            e = aaxb * v + ra;
          }
        }
      }
    }
  }
  else
  {
    // segments are parallel
    if (a01 > 0.0)
    {
      // direction vectors form an obtuse angle
      if (aaxa >= 0.0)
      {
        u = 0.0;
        v = 0.0;
        e = ra;
      }
      else if (-aaxa <= a00)
      {
        u = -aaxa / a00;
        v = 0.0;
        e = aaxa * u + ra;
      }
      else
      {
        u = 1.0;
        t = a00 + aaxa;
        if (-t >= a01)
        {
          v = 1.0;
          e = a00 + a11 + ra + (2.0) * (a01 + aaxa + aaxb);
        }
        else
        {
          v = -t / a01;
          e = a00 + (2.0) * aaxa + ra + v * (a11 * v + (2.0) * (a01 + aaxb));
        }
      }
    }
    else
    {
      // direction vectors form an acute angle
      if (-aaxa >= a00)
      {
        u = 1.0;
        v = 0.0;
        e = a00 + (2.0) * aaxa + ra;
      }
      else if (aaxa <= 0.0)
      {
        u = -aaxa / a00;
        v = 0.0;
        e = aaxa * u + ra;
      }
      else
      {
        u = 0.0;
        if (aaxa >= -a01)
        {
          v = 1.0;
          e = a11 + (2.0) * aaxb + ra;
        }
        else
        {
          v = -aaxa / a01;
          e = ra + v * ((2.0) * aaxb + a11 * v);
        }
      }
    }
  }

  return fabs(e);
}

/*!
\brief Compute the intersections between a capsule and a ray.

Note that this member function assumes that the ray is normalized, i.e. has a unit direction vector.

Intersections are sorted according to intersection depth.

\param ray The normalized ray.
\param ta, tb Intersection depths.
*/
int Capsule::Intersect(const Ray& ray, double& ta, double& tb) const
{
  int n = 0;

  // Cylinder
  if (Cylinder::Intersect(ray, ta, tb))
  {
    n = 1;
  }

  double u[2];

  // Sphere
  if (Sphere(a, r).Intersect(ray, u[0], u[1]))
  {
    if (n == 1)
    {
      if (u[0] < ta)
      {
        ta = u[0];
      }
      if (u[1] > tb)
      {
        tb = u[1];
      }
    }
    else
    {
      ta = u[0];
      tb = u[1];
      n = 1;
    }
  }

  // Sphere
  if (Sphere(b, r).Intersect(ray, u[0], u[1]))
  {
    if (n == 1)
    {
      if (u[0] < ta)
      {
        ta = u[0];
      }
      if (u[1] > tb)
      {
        tb = u[1];
      }
    }
    else
    {
      ta = u[0];
      tb = u[1];
      n = 1;
    }
  }

  return n;
}

/*!
\brief Check the intersection between a capsule and a ray.

Note the function assumes that the ray is normalized, i.e. has a unit direction vector.

\param ray The normalized ray.
*/
bool Capsule::Intersect(const Ray& ray) const
{
  double ta, tb;
  // Cylinder
  if (Cylinder::Intersect(ray, ta, tb) || Sphere(a, r).Intersect(ray) || Sphere(b, r).Intersect(ray))
  {
    return true;
  }

  return false;
}

/*!
\brief Compute the intersections between a capsule and a ray.

Intersections are sorted according to intersection depth.

Note that this function is not fully optimized, as it relies on the
member functions Sphere::Intersect() and Cylinder::Intersect() to compute
the intersection with the union of two spheres and a cylinder. Several
enhancements could be made by inlining the functions with the corresponding
code and removing redundant computations.

\param ray The ray.
\param ta, tb Intersection depths.
\param na, nb Normals at intersection points.
*/
int Capsule::Intersect(const Ray& ray, double& ta, double& tb, Vector& na, Vector& nb) const
{
  int n = 0;

  // Cylinder
  if (Cylinder::Intersect(ray, ta, tb, na, nb))
  {
    n = 1;
  }

  double u[2];
  Vector nu[2];

  // Sphere
  if (Sphere(a, r).Intersect(ray, u[0], u[1], nu[0], nu[1]))
  {
    if (n == 1)
    {
      if (u[0] < ta)
      {
        ta = u[0];
        na = nu[0];
      }
      if (u[1] > tb)
      {
        tb = u[1];
        nb = nu[1];
      }
    }
    else
    {
      ta = u[0];
      na = nu[0];
      tb = u[1];
      nb = nu[1];
      n = 1;
    }
  }

  // Sphere
  if (Sphere(b, r).Intersect(ray, u[0], u[1], nu[0], nu[1]))
  {
    if (n == 1)
    {
      if (u[0] < ta)
      {
        ta = u[0];
        na = nu[0];
      }
      if (u[1] > tb)
      {
        tb = u[1];
        nb = nu[1];
      }
    }
    else
    {
      ta = u[0];
      na = nu[0];
      tb = u[1];
      nb = nu[1];
      n = 1;
    }
  }

  return n;
}

/*!
\brief Compute the first positive intersection between a capsule and a ray.
\param ray The ray.
\param t Intersection depth.
\param n Normal at intersection point.
*/
int Capsule::Intersect(const Ray& ray, double& t, Vector& n) const
{
  double ta, tb;
  Vector na, nb;

  if (Capsule::Intersect(ray, ta, tb, na, nb))
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
\brief Return a capsule rotated by a given matrix.
\param r Rotation matrix.
*/
Capsule Capsule::Rotated(const Matrix& r) const
{
  return Capsule(r * a, r * b, Capsule::r);
}

/*!
\brief Return a capsule translated by a given vector.
\param t Translation vector.
*/
Capsule Capsule::Translated(const Vector& t) const
{
  return Capsule(a + t, b + t, r);
}

/*!
\brief Return a scaled capsule.
\param s Scaling factor.
*/
Capsule Capsule::Scaled(const double& s) const
{
  return Capsule(s * a, s * b, s * r);
}

/*!
\brief Return a capsule transformed by a frame.
\param frame The frame.
*/
Capsule Capsule::Transformed(const Frame& frame) const
{
  return Capsule(frame.Transform(a), frame.Transform(b), r);
}

/*!
\brief Overloaded.
\param s Stream.
\param c %Capsule.
*/
std::ostream& operator<<(std::ostream& s, const Capsule& c)
{
  s << "Capsule(" << c.a << ',' << c.b << ',' << c.r << ')';
  return s;
}
