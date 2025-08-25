
// Sphere

#include "libs/sphere.h"

const double Sphere::epsilon = 1.0e-4;

/*!
\class Sphere sphere.h
\brief Spheres.

Spheres are defined by their center and radius. The
squared radius is not stored in the data-structure, and is
computed on the fly whenever needed.

This class implements several high level functions to compute
the minimal bounding sphere of a set of points.
It also implements the Minkowski sum of two spheres.

The sphere-ray intersection functions assume that the ray
is normalized, i.e. has a unit direction vector.

\ingroup KernelGroup

*/

const Sphere Sphere::Null(0.0);

const Sphere Sphere::Infinity(Math::Infinity);

const Sphere Sphere::Unit(1.0);

/*!
\brief Creates a sphere centered at origin with specified radius.
\param r Radius.
*/
Sphere::Sphere(const double& r):c(Vector::Null),r(r)
{
}

/*!
\brief Creates a sphere given center and radius.
\param c Center.
\param r Radius.
*/
Sphere::Sphere(const Vector& c, const double& r):c(c),r(r)
{
}

/*!
\brief Creates a sphere given two points.

Basically set the center
to the middle of the segment and set the radius to the half
length of the segment.
\param a, b The two vertices.
*/
Sphere::Sphere(const Vector& a, const Vector& b)
{
  c = (a + b) * 0.5;
  r = 0.5 * Norm(b - a);
}

/*!
\brief Creates a sphere given three vertices.

Set the center to the center of the triangle, and compute
the corresponding radius. Note that the argument vertices should not
be aligned.
\param x,y,z The three vertices.
*/
Sphere::Sphere(const Vector& x, const Vector& y, const Vector& z)
{
  Vector a = x - z;
  Vector b = y - z;

  double na = a * a;
  double nb = b * b;
  double ab = a * b;

  double d = na * nb - ab * ab;

  if (fabs(d) > 1e-06)
  {
    double e = 0.5 / d;
    double u = e * nb * (na - ab);
    double v = e * na * (nb - ab);

    c = u * x + v * y + (1.0 - u - v) * z;
    r = Norm(c - x);
  }
  else
  {
    c = Vector(0.0);
    r = 0.0;
  }
}

/*!
\brief Create a sphere circumsizing four vertices.

This is the same as a sphere bounding a tetrahedron.

\sa Tetrahedra::GetSphere()

\param p0, p1, p2, p3 Argument vertices
*/
Sphere::Sphere(const Vector& p0, const Vector& p1, const Vector& p2, const Vector& p3)
{
  Vector e10 = p1 - p0;
  Vector e20 = p2 - p0;
  Vector e30 = p3 - p0;

  // Column matrix
  Matrix A(e10, e20, e30);

  // Transpose
  A = A.T();

  Vector b = 0.5 * Vector(e10 * e10, e20 * e20, e30 * e30);

  // Inverse
  Matrix AI = Inverse(A);

  Vector x = AI * b;

  c = p0 + x;
  r = Norm(x);
}

/*!
\brief Compute the squared distance between a point and the sphere.
\param p The point.
*/
double Sphere::R(const Vector& p) const
{
  double a = (p - c) * (p - c);
  if (a < r * r)
  {
    return 0.0;
  }
  a = sqrt(a) - r;
  a *= a;
  return a;
}

/*!
\brief Compute the signed distance between a point and the sphere.
\param p The point.
*/
double Sphere::Signed(const Vector& p) const
{
  double a = Norm(p - c) - r;

  return a;
}

/*!
\brief Compute the signed distance between two spheres.

If the spheres intersect, the result is negative.
\param s The other Sphere.
*/
double Sphere::R(const Sphere& s) const
{
  return Norm(s.c - c) - s.r - r;
}

/*!
\brief Check if a point is inside or outside the sphere.
\param p The point.
*/
bool Sphere::Inside(const Vector& p) const
{
  Vector q = c - p;
  return (q * q < r * r) ? true : false;
}

/*!
\brief Check the intersection between a sphere and a ray.

Note that intersections are sorted.

This function assumes that the ray is normalized, i.e. has a unit direction vector.

\param ray The (normalized) ray.

\sa Sphere::Intersect(const Ray&,double&,double&) const
*/
bool Sphere::Intersect(const Ray& ray) const
{
  // Center
  Vector n = c - ray.Origin();

  // Distance to center
  double k = n * n;
  double t = n * ray.Direction();
  if (!(k < r * r) && (t < epsilon))
    return false;

  double h = r * r - k + t * t;

  if (h < epsilon)
    return false;

  return true;
}

/*!
\brief Compute the intersection between a sphere and a ray.

Note that intersections are sorted.

This function assumes that the ray is normalized, i.e. has a unit direction vector.

\param ray The (normalized) ray.
\param ta, tb Intersection depths.
*/
int Sphere::Intersect(const Ray& ray, double& ta, double& tb) const
{
  // Center
  Vector n = c - ray.Origin();

  // Distance to center
  double k = n * n;
  double t = n * ray.Direction();
  if (!(k < r * r) && (t < epsilon))
    return 0;

  double h = r * r - k + t * t;

  if (h < epsilon)
    return 0;

  // Intersection depths
  h = sqrt(h);
  ta = t - h;
  tb = t + h;

  return 1;
}

/*!
\brief This function computes the intersections between a sphere and a ray.

The intersections are sorted.

This function assumes that the ray is normalized, i.e. has a unit direction vector.

\param ray The (normalized) ray.
\param ta, tb Intersection depths.
\param na, nb Normals at intersection points.
*/
int Sphere::Intersect(const Ray& ray, double& ta, double& tb, Vector& na, Vector& nb) const
{
  // Center
  Vector n = c - ray.Origin();

  // Distance to center
  double k = n * n;
  double t = n * ray.Direction();

  if (!(k < r * r) && (t < epsilon))
    return 0;

  double h = r * r - k + t * t;

  if (h < epsilon)
    return 0;

  // Intersection depths
  h = sqrt(h);
  ta = t - h;
  tb = t + h;

  na = (ray(ta) - c) / r;
  nb = (ray(tb) - c) / r;

  return 1;
}

/*!
\brief This function computes the first intersection between a sphere and a ray.

This function assumes that the ray is normalized, i.e., has a unit direction vector.

\param ray The (normalized) ray.
\param t Returned intersection depth.
*/
bool Sphere::Intersect(const Ray& ray, double& t) const
{
  double tb;
  if (Sphere::Intersect(ray, t, tb))
  {
    if (t > 0.0)
    {
      return true;
    }
    if (tb > 0.0)
    {
      t = tb;
      return true;
    }
  }
  return false;
}

/*!
\brief Computes the normal vector between a point and
the sphere.

Simply project point onto the sphere, and return the corresponding
Euclidean distance vector.
\param p Point.
*/
Vector Sphere::Normal(const Vector& p) const
{
  Vector n = p - c;

  double t = n * n;

  if (t < r * r)
    return Vector(0.0);

  t = sqrt(t);
  n *= (1.0 - r / t);

  return n;
}

/*!
\brief Rotates a sphere.

\param r Rotation matrix.
*/
void Sphere::Rotate(const Matrix& r)
{
  c = r * c;
}

/*!
\brief Rotates a sphere.

\param r Rotation matrix.
*/
Sphere Sphere::Rotated(const Matrix& r) const
{
  return Sphere(r * c, Sphere::r);
}

/*!
\brief Transforms a sphere.

\param t Transformation.
*/
Sphere Sphere::Transformed(const Frame& t) const
{
  return Sphere(t.Transform(c), r);
}

/*!
\brief Inverse transforms a sphere.

\param t Transformation.
*/
Sphere Sphere::InverseTransformed(const Frame& t) const
{
  return Sphere(t.InverseTransform(c), r);
}

/*!
\brief Translates a sphere.

\param t Translation vector.
*/
void Sphere::Translate(const Vector& t)
{
  c += t;
}

/*!
\brief Translate a sphere.

\param t Translation vector.
*/
Sphere Sphere::Translated(const Vector& t) const
{
  return Sphere(c + t, r);
}

/*!
\brief Scales a sphere by a given vector.

Since the sphere may become an ellipsoid, compute the enclosing sphere.

\param s Scaling vector.
*/
Sphere Sphere::Scaled(const Vector& s) const
{
  return Sphere(c.Scaled(s), r * NormInfinity(s));
}

/*!
\brief Uniformly scales a sphere.

\param s Scaling factor.
*/
void Sphere::Scale(const double& s)
{
  // Center
  c *= s;

  // Radius
  r *= fabs(s);
}

/*!
\brief Extend the sphere, i.e. increase the radius of the sphere.

Contrary to Scale, this function does not change the center.

\param e Radius change.
*/
void Sphere::Extend(const double& e)
{
  r += e;
}

/*!
\brief Extend the sphere, i.e. increase the radius of the sphere.

\sa Extend

\param e Radius change.
*/
Sphere Sphere::Extended(const double& e) const
{
  return Sphere(c, r + e);
}

/*!
\brief Overloaded.

\param s Stream.
\param sphere The sphere.
*/
std::ostream& operator<<(std::ostream& s, const Sphere& sphere)
{
  s << "Sphere(" << sphere.c << ',' << sphere.r << ')';
  return s;
}

/*!
\brief Generate a random point on the sphere.
\param random %Random number generator.
*/
Vector Sphere::RandomSurface(Random& random) const
{
  return c + r * RandomNormal(random);
}

/*!
\brief Generate a random unit vector orthonormal to the sphere.
\param random %Random number generator.
*/
Vector Sphere::RandomNormal(Random& random)
{
  double z = 2.0 * random.Uniform() - 1.0;
  double w = sqrt(1.0 - z * z);

  Vector2 xy = Circle2(1.0).RandomOn(random);
  xy *= w;
  return xy.ToVector(z);
}

/*!
Compute the i-th point of the Fibonacci distribution of points on the sphere.

In other words, phyllotaxic samples for nearly optimal distribution over a sphere.

\param i Point index.
\param n Number of samples on the sphere.
*/
Vector Sphere::Fibonacci(int i, int n) const
{
  int u = 1;
  double offset = 2.0 / n;
  double increment = Math::Pi * (3.0 - Math::Sqrt5);

  //  for i in range(samples):
  double y = ((i * offset) - 1.0) + (offset / 2.0);
  double ry = sqrt(1.0 - y * y);

  double phi = ((i + u) % n) * increment;

  double x = cos(phi) * ry;
  double z = sin(phi) * ry;

  return c + r * Vector(x, y, z);
}

/*!
\brief Create a %Poisson %Disc sampling of the sphere.

Simple dart throwing algorithm.

\param r Radius of the discs.
\param n Maximum number of darts thrown on the sphere.
\param random %Random number generator.
*/
QVector<Vector> Sphere::Poisson(const double& r, int n, Random& random) const
{
  QVector<Vector> p;

  // Collision radius
  double c = 4.0 * r * r;

  // Create instances
  for (int i = 0; i < n; i++)
  {
    Vector t = RandomSurface(random);
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
    {
      p.append(t);
    }
  }
  return p;
}

/*!
\brief Generate a random vector inside the sphere.
\param random %Random number generator.
*/
Vector Sphere::RandomInside(Random& random) const
{
  Vector p;
  // Compute random vector inside unit sphere
  while (true)
  {
    p = Vector(random.Uniform(-1.0, 1.0), random.Uniform(-1.0, 1.0), random.Uniform(-1.0, 1.0));
    if (p * p <= 1.0) break;
  }
  return c + r * p;
}

/*!
\brief %Box-sphere intersection test.

The algorithm calculates the square of the distance from the box to the
sphere by analyzing the orientation of the sphere relative to the box in
a single loop. If the box is not axis aligned, transform the center
of the sphere to the box's local coordinate frame.

After J. Arvo, <I>A simple method for box-sphere intersection testing</I>,
In <B>A. Glassner, Graphics Gems</B>: 335-339, Academic Press 1990.

\sa Box::R()
*/
bool Sphere::Intersect(const Box& box) const
{
  double d = 0.0;

  // find the square of the distance from the sphere to the box
  for (int i = 0; i < 3; i++)
  {
    if (c[i] < box[0][i])
    {
      double s = c[i] - box[0][i];
      d += s * s;
    }
    else if (c[i] > box[1][i])
    {
      double s = c[i] - box[1][i];
      d += s * s;
    }
  }
  return d <= r * r;
}

/*!
\brief Check if two spheres intersect.

If the centers are identical, then the spheres intersect if and only
if their radii are equal. Otherwise, there exists a circle of intersection,
the normal of the plane of intersection is derived from the centers
of the spheres.

\param sphere %Sphere.
\param circle Returned %Circle.
*/
bool Sphere::Intersect(const Sphere& sphere, Circle& circle) const
{
  circle = Circle(Vector::Null, Vector::Z,0.0);

  Vector n = c - sphere.c;
  double length = n * n;
  double s = r + sphere.r;
  if (length > s * s)
    return false;

  double rb = sphere.r * sphere.r;

  double t = 0.5 * (1.0 + (r * r - rb) / length);
  if (t < 0.0 || t>1.0)
    return false;

  double rc = r * r - t * t * length;
  if (rc < 0.0)
    return false;

  rc = sqrt(rc);

  circle = Circle(c + t * n, n/length, rc);
  return true;
}

/*!
\brief Check if two spheres intersect.

\param sphere Sphere.
*/
bool Sphere::Intersect(const Sphere& sphere) const
{
  Vector n = c - sphere.c;

  if (n * n > Math::Sqr(r + sphere.r))
    return false;

  return true;
}

/*!
\brief Compute the minimal bounding sphere of a set of points.

Error around 1%~2%. Based on Bo Tian, <I>Bouncing bubble: a
fast algorithm for minimal enclosing pall problem</I>, 2012.
\param p Set of points.
*/
Sphere::Sphere(const QVector<Vector>& p)
{
  // Starting center
  c = p[0];
  r = 0.0001;

  QVector<Vector>::const_iterator it;

  for (int i = 0; i < 2; i++)
  {
    for (it = p.begin(); it != p.end(); it++)
    {
      Vector q = *it;
      Vector cq = q - c;
      double len = Norm(cq);
      if (len > r)
      {
        double alpha = len / r;
        double a = alpha * alpha;
        r = 0.5 * (alpha + 1 / alpha) * r;
        c = 0.5 * ((1 + 1 / a) * c + (1 - 1 / a) * q);
      }
    }
  }

  for (it = p.begin(); it != p.end(); it++)
  {
    Vector q = *it;
    Vector cq = q - c;
    double len = Norm(cq);
    if (len > r)
    {
      r = (r + len) / 2.0;
      c = c + ((len - r) / len * cq);
    }
  }
}

/*!
\brief Compute the great-circle or orthodromic distance.

It is the geodesic distance between two points on the surface of a sphere.

This function is computationnally intensive, requiring several square root and inverse
trigonometric function evaluations.

\param a,b Points on the sphere.
*/
double Sphere::R(const Vector& a, const Vector& b) const
{
  Vector na = Normalized(a - c);
  Vector nb = Normalized(b - c);
  return r * atan2(Norm(na / nb), na * nb);
}

/*!
\brief Compute the volume of the intersection of two spheres.

If the spheres do not intersect, the function returns 0; if the smallest
is included in the largest, the function returns the volume of the smallest,
otherwise complex computations are performed.

\sa Circle::Area(const Circle2&)

\param sphere The other sphere.
*/
double Sphere::Volume(const Sphere& sphere) const
{
  // Check if spheres do not intersect
  double cc = SquaredNorm(c - sphere.c);
  double rr = Math::Sqr(r + sphere.r);
  if (cc >= rr)
  {
    return 0.0;
  }

  // Sort radii so that ra is the large radius and rb is the small one
  double ra = r;
  double rb = sphere.r;
  Math::Sort(rb, ra);

  // Distance between centers
  double scc = sqrt(cc);

  // Small sphere is fully inside the large one
  if (scc + rb <= ra)
  {
    return Math::Pi * rb * rb;
  }

  return Math::Pi * Math::Sqr(ra + rb - cc) * (scc + 2.0 * cc * rb - 3.0 * rb * rb + 2.0 * cc * ra + 6.0 * ra * rb - 3.0 * ra * ra) / (12.0 * cc);
}

/*!
\brief Compute the Euler coordinates of a point.
\param p Point.
\return A vector containing the theta and phi angles.
*/
Vector2 Sphere::Euler(const Vector& p) const
{
  // Unit
  Vector u = Normalized(p - c);

  // Theta
  double theta = acos(p[2]);

  // Phi
  double phi = atan2(p[0], p[1]);

  return Vector2(theta, phi);
}

/*!
\brief Compute the Euler coordinates of a point defined in a rectangle map.
\param i, j Point coordinates.
\param x, y Size of the map.
\return A vector containing the theta and phi angles.
*/
Vector2 Sphere::EquiRectangular(int i, int j, int x, int y)
{
  // Equi-rectangular projection
  double theta = (double(i) / double(x)) * 2.0 * Math::Pi;
  double phi = Math::Pi * (double(y - 1 - j) / double(y - 1) - 0.5);
  return Vector2(theta, phi);
}

/*!
\brief Compute the intersection of three spheres.
\param sa, sb, sc Spheres.
\param u, v Returned points.
*/
bool Sphere::Intersection(const Sphere& sa, const Sphere& sb, const Sphere& sc, Vector& u, Vector& v)
{
  const Vector a = sa.Center();
  const Vector b = sb.Center();
  const Vector c = sc.Center();

  const double ra = sa.Radius();
  const double rb = sb.Radius();
  const double rc = sc.Radius();

  Vector ab = b - a;
  Vector e_x = ab / Norm(ab);
  Vector ac = c - a;
  double i = e_x * ac;
  Vector temp3 = ac - i * e_x;
  Vector e_y = temp3 / Norm(temp3);
  Vector e_z = e_x / e_y;

  double d = Norm(b - a);
  double j = e_y * ac;
  double x = (ra * ra - rb * rb + d * d) / (2.0 * d);
  double y = (ra * ra - rc * rc - 2 * i * x + i * i + j * j) / (2.0 * j);
  double temp4 = ra * ra - x * x - y * y;
  if (temp4 < 0.0)
    return false;
  double z = sqrt(temp4);
  u = a + x * e_x + y * e_y + z * e_z;
  v = a + x * e_x + y * e_y - z * e_z;
  return true;
}