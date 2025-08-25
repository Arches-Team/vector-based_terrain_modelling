// Plane

#include "libs/plane.h"
#include "libs/matrix.h"
#include "libs/segment.h"
#include "libs/quadric.h"
#include "libs/sphere.h"
#include "libs/axis.h"

/*!
\class Plane plane.h
\brief A plane defined by its analytic equation.

This class serves as a base core class for binary space patitionning
or constructive solid geometry in ray-tracing.

\ingroup KernelGroup

*/

const double Plane::epsilon = 1e-06;
const Plane Plane::XY(Vector::Z, 0.0);

/*!
\brief Creates a plane given normal and vertex.

\param n %Plane normal, should be unit.
\param p Any point on the plane.
*/
Plane::Plane(const Vector& n, const Vector& p) :n(n), c(p* n)
{
}

/*!
\brief Creates a plane given normal passing through the origin.

\param n %Plane normal, should be unit.
*/
Plane::Plane(const Vector& n) :Plane(n, Vector::Null)
{
}

/*!
\brief Check if a point is inside or outside the plane.

This function compares the sign strictly to 0.0, use Plane::Side() for a more fuzzy implementation.
\param p Point.
\sa Plane::Side().
*/
bool Plane::Inside(const Vector& p) const
{
  if (n * p - c < 0.0)
  {
    return true;
  }
  return false;
}

/*!
\brief Check if a point lies inside, outside or even on the plane.

The epsilon tolerance used for testing if the point is on the plane
is 10<SUP>-6</SUP>.
\param p Point.
*/
int Plane::Side(const Vector& p) const
{
  double r = n * p - c;

  // Epsilon test
  if (r > epsilon)
    return 1;
  else if (r < -epsilon)
    return -1;
  else
    return 0;
}

/*!
\brief Finds a vertex on the plane.

Basically computes the point on the plane which is nearest to the origin.
*/
Vector Plane::Vertex() const
{
  return -c * n;
}

/*!
\brief Compute the intersection between a plane and a ray.

The intersection depth is returned if intersection occurs.
\param ray The ray.
\param t Intersection depth.
*/
bool Plane::Intersect(const Ray& ray, double& t) const
{
  // Check if parallel
  double x = n * ray.Direction();

  if ((x < epsilon) && (x > -epsilon))
    return false;

  double y = c - n * ray.Origin();

  // Depth
  t = y / x;

  return true;
}

/*!
\brief Compute the intersection between a plane and a line.

The line should not be parallel to the plane.
\param l The line.
*/
Vector Plane::Intersect(const Line& l) const
{
  // We do not check if parallel
  double x = n * (l.Vertex(1) - l.Vertex(0));

  double y = c - n * l.Vertex(0);

  // Depth
  double t = y / x;

  return l.Vertex(0) + (l.Vertex(1) - l.Vertex(0)) * t;
}

/*!
\brief Compute the intersection between a plane and a segment.

The intersection point can be computed using Segment::Vertex(const double&) if intersection occurs.
\sa Plane::Intersect(const Line&) const
\sa Plane::Intersect(const Ray&) const
\param s The segment.
\param t Position of the point if intersection occurs.
*/
bool Plane::Intersect(const Segment& s, double& t) const
{
  const Vector& a = s.Vertex(0);
  const Vector& b = s.Vertex(1);

  // Check if parallel
  double x = n * (b - a);
  if ((x < epsilon) && (x > -epsilon))
    return false;

  double y = c - n * a;

  // Position
  t = y / x;

  // Outside of segment
  if ((t < 0.0) || (t > 1.0))
    return false;
  return true;
}

/*!
\brief Transforms a plane given a homogeneous transformation matrix.
\param M %Matrix.
*/
void Plane::Transform(const Matrix4& M)
{
  // Find a point on the plane
  Vector x = Vertex();

  // Find a point in normal direction 
  Vector y = x + n;

  // Transform
  Vector tx = M * x;
  Vector ty = M * y;
  Vector tn = Normalized(ty - tx);

  n = tn;
  c = -tx * tn;
}

/*!
\brief Translates a plane in the direction of its normal.
\param t Translation distance.
*/
void Plane::Translate(const double& t)
{
  // p'=p+(n t) so c'=(p'.n)=(p.n)+(n.n) t=c+t since n is normalized
  c += t;
}

/*!
\brief Rotate a plane.
\param r Rotation matrix.
*/
void Plane::Rotate(const Matrix& r)
{
  Vector p = r * Vertex();
  n = r * n;
  c = p * n;
}

/*!
\brief Compute the intersection between three planes.

Returns Vector::Null if two planes are complanar.
\param a, b, c Three planes.
*/
Vector Plane::Intersection(const Plane& a, const Plane& b, const Plane& c)
{
  double e = Matrix4(Matrix(a.Normal(), b.Normal(), c.Normal())).Determinant();
  if (e < Plane::epsilon)
  {
    return Vector::Null;
  }

  Vector p = (a.Vertex() * a.Normal()) * (b.Normal() / c.Normal()) + (b.Vertex() * b.Normal()) * (c.Normal() / a.Normal()) + (c.Vertex() * c.Normal()) * (a.Normal() / b.Normal());
  return p / (-e);
}


/*!
\brief Compute the intesection between two planes.

Planes should not be coplanar.

\param plane Argument plane.
*/
Axis Plane::Intersect(const Plane& plane) const
{
  Vector u = plane.n / n;
  Normalize(u);

  // Vertex on plane A
  Vector q = Vertex();

  // Move it onto plane B, by the distance to the plane B and in the direction of the normal of plane B
  q += plane.Signed(q) * plane.n;

  return Axis(q, q + u);
}

/*!
\brief Compute the intersection of a set of half spaces.

This is a O(n^4) algorithm

A better version in O(n^3) could be implemented: see Finding the Intersection
of Half-Spaces in Time O(n ln n). Preparata and Muller, Theoretical Computer Science 8, 45-55, 1979.

Returns a set of point representing the minimal convex polygon embedding all half spaces.
The function doesn't check if the intersection is bounded or not.

\param planes Set of planes.
*/
QVector<Vector> Plane::ConvexPoints(const QVector<Plane>& planes)
{
  QVector<Vector> pts;
  for (int i = 0; i < planes.size(); i++)
  {
    for (int j = i + 1; j < planes.size(); j++)
    {
      for (int k = j + 1; k < planes.size(); k++)
      {
        Vector p = Intersection(planes[i], planes[j], planes[k]);
        if (p != Vector::Null)
        {
          bool isInside = true;
          for (int l = 0; l < planes.size(); l++)
          {
            // Do not check point ijk with one of its generating plane
            if (l == i || l == j || l == k)
              continue;
            int s = planes[l].Side(p);
            if (s > 0)
            {
              isInside = false;
              break;
            }
          }
          if (isInside)
            pts.append(p);
        }
      }
    }
  }
  return pts;
}

/*!
\brief Compute the inverse mapping coordinates of a point.
\param p Point that will be projected on the plane.
\param u, v Returned coordinates.
*/
void Plane::Cast(const Vector& p, double& u, double& v)
{
  u = Math::Mod(p[0], 1.0);
  v = Math::Mod(p[1], 1.0);
}

/*!
\brief Compute the point symmetric to the argument point.
\param p Point.
*/
Vector Plane::Symmetry(const Vector& p) const
{
  return p - 2.0 * Math::Abs(Eval(p)) * n;
}

/*!
\brief Compute the squared distance to the half space delimited by the plane.
\param p Point.

\sa Plane::Eval()
*/
double Plane::R(const Vector& p) const
{
  double r = n * p - c;
  if (r < 0.0)
  {
    return 0.0;
  }
  else
  {
    return r * r;
  }
}

/*!
\brief Compute the signed distance between a point and the plane.
\param p The point.
\sa Plane::Eval(const Vector&) const
*/
double Plane::Signed(const Vector& p) const
{
  return n * p - c;
}

/*!
\brief Return a Frame attached to the plane.

\sa Axis::GetFrame()
*/
Frame Plane::GetFrame() const
{
  Vector x, y;
  n.Orthonormal(x, y);
  return Frame(Vertex(), x, y, n);
}

/*!
\brief Compute the quadric equation of the Euclidean distance to the plane.
\param ray The ray.
*/
Quadric Plane::Equation(const Ray& ray) const
{
  double s = ray.Origin() * n - c;
  double t = ray.Direction() * n;
  return Quadric(t * t, 2.0 * s * t, s * s);
}

/*!
\brief Overloaded.

\param s Stream.
\param plane The plane.
*/
std::ostream& operator<<(std::ostream& s, const Plane& plane)
{
  s << "Plane(" << plane.n << ',' << plane.c << ')';
  return s;
}

/*!
\brief Compute the boundix box of the symmetric argument box.
\param box %The box.
*/
Box Plane::Symmetric(const Box& box) const
{
  Vector v[8];
  for (int i = 0; i < 8; i++)
  {
    v[i] = Symmetry(box.Vertex(i));
  }

  return Box(v, 8);
}

/*!
\brief Compute the symmetric segment.
\param segment %The segment.
*/
Segment Plane::Symmetric(const Segment& segment) const
{
  return Segment(Symmetry(segment.Vertex(0)), Symmetry(segment.Vertex(1)));
}

/*!
\brief Compute the symmetric sphere.
\param sphere %The sphere.
*/
Sphere Plane::Symmetric(const Sphere& sphere) const
{
  return Sphere(Symmetry(sphere.Center()), sphere.Radius());
}

/*!
\brief Compute the refracted direction.

\param d Direction.
\param eta Index of incoming ray medium divided by the index of the other material.
\param r Refracted direction.
*/
bool Plane::Refract(const Vector& d, const double& eta, Vector& r) const
{
  double ni = n * d;
  double k = 1.0 - eta * eta * (1.0 - ni * ni);
  if (k < 0.0)
  {
    return false;
  }
  else
  {
    r = eta * d - (eta * ni + sqrt(k)) * n;
    return true;
  }
}

/*!
\brief Compute the refracted direction.

\param d Direction.
\param n Normal.
\param eta Index of incoming ray medium divided by the index of the other material.
\param r Refracted direction.
*/
bool Plane::Refract(const Vector& d, const Vector& n, const double& eta, Vector& r)
{
  double ni = n * d;
  double k = 1.0 - eta * eta * (1.0 - ni * ni);
  if (k < 0.0)
  {
    return false;
  }
  else
  {
    r = eta * d - (eta * ni + sqrt(k)) * n;
    return true;
  }
}
