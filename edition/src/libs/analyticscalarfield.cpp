// Analytic fields

#include "libs/scalarfield.h"
#include "libs/sphere.h"
#include "libs/segment.h"
#include "libs/color.h"
#include "libs/quadrangle.h"
#include "libs/rectangle.h"
#include "libs/voxel.h"

/*!
\class AnalyticScalarField scalarfield.h
\brief A core analytic three-dimensional scalar field.

\ingroup StructureGroup
*/

const double AnalyticScalarField::Epsilon = 1e-6;

int AnalyticScalarField::ncubes = 0;

/*!
\brief Get the number of cubes processed.
*/
int AnalyticScalarField::Cubes()
{
  return ncubes;
}

/*!
\brief Compute the color.
\param p Point.
\param n Normal at point.
\return %Color.
*/
#pragma warning(push)
#pragma warning(disable: 4100)  
Color AnalyticScalarField::GetMaterial(const Vector& p, const Vector& n) const
{
  return Color::Black;
}
#pragma warning(pop)

/*!
\brief Compute the value of the field.
\param p Point.
*/
#pragma warning(push)
#pragma warning(disable: 4100)  
double AnalyticScalarField::Value(const Vector& p) const
{
  return 0.0;
}
#pragma warning(pop)

/*!
\brief Compute the bounding box of the scalar field.

Simply return Box::Infinity.
*/
Box AnalyticScalarField::GetBox() const
{
  return Box::Infinity;
}

/*!
\brief Compute the  Lipschitz constant.

By default, return 1.0.
*/
double AnalyticScalarField::K() const
{
  return 1.0;
}

/*!
\brief Compute the local Lipschitz constant inside a box.
\param box %Box.
*/
#pragma warning(push)
#pragma warning(disable: 4100) 
double AnalyticScalarField::K(const Box& box) const
{
  return K();
}
#pragma warning(pop)

/*!
\brief Compute the local Lipschitz constant inside a box.
\param sphere %Sphere.
*/
double AnalyticScalarField::K(const Sphere& sphere) const
{
  return K(sphere.GetBox());

}

/*!
\brief Compute the local Lipschitz constant along a segment.
\param segment The segment.
*/
double AnalyticScalarField::K(const Segment& segment) const
{
  return K(segment.GetBox());
}

/*!
\brief Compute the gradient of the field.
\param p Point.
*/
Vector AnalyticScalarField::Gradient(const Vector& p) const
{
  double x = Value(Vector(p[0] + Epsilon, p[1], p[2])) - Value(Vector(p[0] - Epsilon, p[1], p[2]));
  double y = Value(Vector(p[0], p[1] + Epsilon, p[2])) - Value(Vector(p[0], p[1] - Epsilon, p[2]));
  double z = Value(Vector(p[0], p[1], p[2] + Epsilon)) - Value(Vector(p[0], p[1], p[2] - Epsilon));

  return Vector(x, y, z) * (0.5 / Epsilon);
}

/*!
\brief Compute the normal to the surface.

\sa AnalyticScalarField::Gradient(const Vector&) const

\param p Point (should be on the surface).
*/
Vector AnalyticScalarField::Normal(const Vector& p) const
{
  Vector normal = Normalized(Gradient(p));
  if (sign)
  {
    return normal;
  }
  else
  {
    return -normal;
  }
}

/*!
\brief Compute the Hessian symmetric matrix of the field function.

This method uses a numerical approximation of the successive derivatives
in each direction.

\param p Point.
\return Hessian matrix.
*/
Matrix AnalyticScalarField::Hessian(const Vector& p) const
{
  double a[3][3][3];
  Matrix H;

  for (int i = -1; i < 2; i++)
  {
    for (int j = -1; j < 2; j++)
    {
      for (int k = -1; k < 2; k++)
      {
        a[1 + i][1 + j][1 + k] = Value(p + Vector(i, j, k) * Epsilon);
      }
    }
  }
  H[0] = (a[2][1][1] - 2.0 * a[1][1][1] + a[0][1][1]) / (Epsilon * Epsilon);
  H[1] = (a[2][2][1] - a[2][0][1] - a[0][2][1] + a[0][0][1]) / (2.0 * Epsilon * Epsilon);
  H[2] = (a[2][1][2] - a[2][1][0] - a[0][1][2] + a[0][1][0]) / (2.0 * Epsilon * Epsilon);

  H[3] = H[1];
  H[4] = (a[1][2][1] - 2.0 * a[1][1][1] + a[1][0][1]) / (Epsilon * Epsilon);
  H[5] = (a[1][2][2] - a[1][2][0] - a[1][0][2] + a[1][0][0]) / (2.0 * Epsilon * Epsilon);

  H[6] = H[2];
  H[7] = H[5];
  H[8] = (a[1][1][2] - 2.0 * a[1][1][1] + a[1][1][0]) / (Epsilon * Epsilon);

  return H;
}

/*!
\brief Compute the intersection between a segment and an implicit surface.

This method iteratively converges to the solution by performing a bisection on
the segment. It avoids the internal computation if the length of the segment
if it happens to be known already.

The segment should straddle the implicit surface only once, otherwise
should several intersections exist, there is no guarantee that the bissection will
converge to an intersection.

\param a,b End vertices of the segment straddling the surface.
\param va,vb Field function value at those end vertices.
\param length Distance between vertices.
\param epsilon Precision.
\return Point on the implicit surface.
*/
Vector AnalyticScalarField::Dichotomy(Vector a, Vector b, double va, double vb, double length, const double& epsilon) const
{
  int ia = va > 0.0 ? 1 : -1;

  // Get an accurate first guess
  Vector c = Vector::Solve(a, b, va, vb);

  while (length > epsilon)
  {
    double vc = Value(c);
    int ic = vc > 0.0 ? 1 : -1;
    if (ia + ic == 0)
    {
      b = c;
    }
    else
    {
      ia = ic;
      a = c;
    }
    length *= 0.5;
    c = 0.5 * (a + b);
  }
  return c;
}

/*!
\brief Find a random sample point inside or outside the surface.
\param p Returned point.
\param s Prescribed position with respect to the surface, true if inside, false if outside.
\param box Box domain where the point is to be found.
\param n Maximum number of random points being evaluated.
\param random %Random number generator.
\return Return true if a sample has been found, false otherwise.
*/
bool AnalyticScalarField::Find(Vector& p, bool s, const Box& box, int n, Random& random) const
{
  // Search randomly
  for (int i = 0; i < n; i++)
  {
    // Set vertex
    p = box.RandomInside(random);
    double v = Value(p);
    if (s == (v > 0.0))
    {
      return true;
    }
  }
  return false;
}

/*!
\brief Find a point sample on the implicit surface.

Simply find two random points inside and outside of
the surface, and perform bisection to converge to
the surface.

\param box Sampling box.
\param p Point sample on the surface, if succeded in finding one.
\param random %Random number generator.
*/
bool AnalyticScalarField::GetSample(Vector& p, const Box& box, Random& random) const
{
  Vector a, b;

  if (!Find(a, true, box, 10000, random))
    return false;
  if (!Find(b, false, box, 10000, random))
    return false;

  p = Dichotomy(a, b, Value(a), Value(b), Norm(b - a), 1e-4);
  return true;
}

/*!
\brief Find a set of random sample points on the implicit surface.

Simply call AnalyticScalarField::GetSample() several times.

\param box Sampling box.
\param s Set of sample points.
\param n Number of (desired) sampes.
\param random %Random number generator.
\return The number of samples found on the surface.
*/
int AnalyticScalarField::GetSamples(QVector<Vector>& s, const Box& box, int n, Random& random) const
{
  for (int i = 0; i < n; i++)
  {
    Vector p;
    if (GetSample(p, box, random))
    {
      s.append(p);
    }
  }
  return s.size();
}

/*!
\brief Compute a Poisson sphere distribution on an implicit surface.

This function uses a simple O(n<SUP>2</SUP>) dart throwing algorithm.

\param r Radius of the sphere.
\param n Number of candidate points.
\param rng %Random number generator.
*/
QVector<Vector> AnalyticScalarField::Poisson(double r, int n, Random& rng) const
{
  Box box = GetBox();
  QVector<Vector> p;

  // Collision radius
  double c = 4.0 * r * r;

  // Create instances
  for (int i = 0; i < n; i++)
  {
    Vector t;
    if (!GetSample(t, box, rng))
      continue;

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
\brief Compute the field function at the vertices a grid.

\param r Rectangle.
\param x,y Discretization.
*/
ScalarField2 AnalyticScalarField::Sample(const Rectangles& r, int x, int y) const
{
  Box2 box(1.0);
  box.Scale(Vector2(r.Width(), r.Height()));
  ScalarField2 s(box, x, y);

  for (int i = 0; i < x; i++)
  {
    double u = Math::Unit(i, x);
    for (int j = 0; j < y; j++)
    {
      double v = Math::Unit(j, y);
      Vector p = r.Vertex(u, v);
      s(i, j) = Value(p);
    }
  }
  return s;
}

/*!
\brief Sample the scalar field on a quadrangle.

Returns the 2D scalarfield.

\param quad %Quadrangle
\param n discretization
*/
ScalarField2 AnalyticScalarField::Sample(const Quadrangle& quad, int n) const
{
  ScalarField2 f = ScalarField2(Box2(quad.GetBox()), n, n, 0.0);
  for (int i = 0; i < n; i++)
  {
    double v = Math::Unit(i, n);
    for (int j = 0; j < n; j++)
    {
      double u = Math::Unit(j, n);
      Vector p = quad.Vertex(u, v);
      f(i, j) = Value(p);
    }
  }
  return f;
}

/*!
\brief Given a starting vertex x, project it onto the implicit surface
following the gradient.

This function returns a vertex on the implicit surface.
It updates a pair of points that are moved so that
they straddle the surface, before invoking a bisection step.

\param p Starting vertex.
\param step Stepping distance.
\param epsilon Precision.
\param n Maximum number of iterations when marching in the direction of the gradient.
\return Point on the implicit surface.
*/
Vector AnalyticScalarField::Cast(const Vector& p, const double& step, const double& epsilon, int n) const
{
  Vector a = p;
  double va = Value(a);
  // QVector<double> list;
  // list << va;
   // On the surface
  if (fabs(va) < epsilon)
  {
    return a;
  }

  Vector b;
  double vb;

  // Loop
  for (int j = 0; j < n; j++)
  {
    // Include the sign of the function to get a vector in the direction of the surface
    Vector g = -Normalized(va * Gradient(a));

    b = a + step * g;
    vb = Value(b);

    // Escape if straddling the surface
    if ((va * vb) <= 0.0)
    {
      return Dichotomy(a, b, va, vb, Norm(b - a), epsilon);
    }
    a = b;
    va = vb;
  }
  // Should not get there
  std::cout << "Should not get there" << std::endl;

  return 0.5 * (a + b);
}

/*!
\brief Compute the gaussian and mean curvatures.
\param p Point.
\param gaussian Returned gaussian curvature.
\param mean Returned mean curvature.
*/
void AnalyticScalarField::Curvature(const Vector& p, double& gaussian, double& mean) const
{
  // Degenerate surface point default
  mean = gaussian = 0.0;
  Vector g = Gradient(p);

  if (g == Vector::Null)
    return;

  // Hessian
  Matrix H = Hessian(p);

  // Temporary
  double fxfx = g[0] * g[0];
  double fxfy = g[0] * g[1];
  double fxfz = g[0] * g[2];
  double fyfy = g[1] * g[1];
  double fyfz = g[1] * g[2];
  double fzfz = g[2] * g[2];

  double n = fxfx + fyfy + fzfz;

  // Gaussian curvature
  gaussian = fxfx * (H[4] * H[8] - H[5] * H[5]) + fyfy * (H[0] * H[8] - H[2] * H[2]) + fzfz * (H[0] * H[4] - H[1] * H[1])
    + 2.0 * (fxfy * (H[2] * H[5] - H[1] * H[8]) + fxfz * (H[1] * H[5] - H[2] * H[4]) + fyfz * (H[1] * H[2] - H[0] * H[5]));
  gaussian /= n * n;

  // Mean curvature
  mean = (H[0] * (fyfy + fzfz) + H[4] * (fxfx + fzfz) + H[8] * (fxfx + fyfy) - 2.0 * (H[1] * fxfy + H[2] * fxfz + H[5] * fyfz)) / 2.0;
  mean /= n * sqrt(n);
}

/*!
\brief Find the intersections between a ray and an implicit surface on a given interval.
\param ray The ray.
\param a,b Interval of analysis.
\param k Lipschitz constant.
\param s Array for storing roots.
\param n Search for n roots at most.
\param epsilon Precision parameter.
*/
int AnalyticScalarField::Roots(const Ray& ray, const double& a, const double& b, const double& k, double* s, int n, const double& epsilon) const
{
  double ya = Value(ray(a));
  double yb = Value(ray(b));

  return Roots(ray, a, b, k, ya, yb, s, n, epsilon);
}

static double stack_t1[256];
static double stack_t2[256];
static int stack_s1[256];
static int stack_s2[256];

static int stack_n;

/*!
\brief This functions searches the intersections between a ray and an implicit surface
on a given interval.

The original algorithm recursively isolate the roots by applying the Lipschitz
criteria over the interval. It has been rewritten into an iterative form
which is more efficient.

\param ray The ray.
\param a,b Interval of analysis.
\param k Lipschitz constant.
\param ya,yb Field function values at the end-points of the interval of analysis.
\param x Array for storing roots.
\param n Search for at most n roots.
\param epsilon Precision parameter.
*/
int AnalyticScalarField::Roots(const Ray& ray, const double& a, const double& b, const double& k, const double& ya, const double& yb, double* x, int n, const double& epsilon) const
{
  // Number of roots
  int np = 0;
  double epsilonI = 1e-4;

  // Initialize stack 
  stack_t1[0] = a;
  stack_t2[0] = b;
  stack_s1[0] = ya > 0 ? 1 : -1;
  stack_s2[0] = yb > 0 ? 1 : -1;
  stack_n = 1;

  // Escape if reached maximum expected roots
  while (np < n)
  {
    stack_n--;
    if (stack_n == -1) return np;

    double T1 = stack_t1[stack_n];
    double T2 = stack_t2[stack_n];
    int S1 = stack_s1[stack_n];
    int S2 = stack_s2[stack_n];

    int S3 = 1;

    double T3 = 0.5 * (T1 + T2);
    double dt = T2 - T3;

    double f3 = Value(ray(T3));
    double ff3 = f3;

    // Compute the absolute value and store the sign
    if (f3 < 0.0)
    {
      ff3 = -f3;
      S3 = -1;
    }

    // Lipschitz test
    if (ff3 > dt * k)
      continue;

    // Small interval check
    if ((dt < epsilon) && ((S1 == S2) || (ff3 < epsilonI)))
    {
      // Root found
      if (S1 != S2)
      {
        x[np] = T3;
        np++;
        //continue;
      }
      // No root
      else
      {
        //continue;
      }
    }
    // Subdivide interval
    else
    {
      double r0 = ff3 / k;

      // Stack last interval first
      stack_t1[stack_n] = T3 + r0;
      stack_t2[stack_n] = T2;
      stack_s1[stack_n] = S3;
      stack_s2[stack_n] = S2;
      stack_n++;

      // Stack first interval last
      stack_t1[stack_n] = T1;
      stack_t2[stack_n] = T3 - r0;
      stack_s1[stack_n] = S1;
      stack_s2[stack_n] = S3;
      stack_n++;
    }
  }
  return np;
}

/*!
\brief Sample the implicit surface to get a scalar field.

\param n Subdivision.
\param box The box.
*/
ScalarField AnalyticScalarField::Sample(const Box& box, int n) const
{
  Box clipped = box;
  int nx, ny, nz;
  clipped.SetParallelepipedic(n, nx, ny, nz);
  Array ar(box, nx, ny, nz);

  ScalarField field = ScalarField(ar, 0.0, sign);

  {
#pragma omp parallel for shared(field) num_threads(8)
    for (int i = 0; i < nx; i++)
    {
      for (int j = 0; j < ny; j++)
      {
        for (int k = 0; k < nz; k++)
        {
          field(i, j, k) = Value(ar.Vertex(i, j, k));
        }
      }
    }
  }
  return field;
}

/*!
\brief Dual polygonization.

Compute a scalar field with function values sampled at the vertices of the grid,
and invoque the dual scalar field polygonization algorithm.
\sa ScalarField::Dual
\param n Subdivision.
\param g %Mesh.
\param box The box.
*/
void AnalyticScalarField::Dual(int n, Mesh& g, const Box& box) const
{
  ScalarField field = Sample(box, n);

  // Number of processed cubes
  ncubes = field.CellSize();

  field.Dual(g);
}

/*!
\brief Find the first intersection between a ray and an implicit surface
on a given interval, using sphere tracing.

\param ray The ray.
\param k Lipschitz constant.
\param a,b Interval where intersection is tracked.
\param t Nearest root, if it exists.
\param s Return number of steps.
\param epsilon Precision on field function values for detecting ray-surface intersection.

Simply march along the ray
with variable steps depending on the magnitude of the field function at
point samples. Intersection occurs if the absolute value of the field
function get lower than epsilon.
*/
bool AnalyticScalarField::SphereTrace(const Ray& ray, const double& k, const double& a, const double& b, double& t, int& s, const double& epsilon) const
{
  t = a;
  s = 0;
  while (t < b)
  {
    Vector p = ray(t);
    s++;
    double i = Value(p);
    // Got inside
    if (i > 0.0)
    {
      return true;
    }
    t += Math::Min(fabs(i) / k, epsilon);
  }
  return false;
}

/*!
\brief Polygonize a cube.

The algorithm uses either a look-up table disambiguation scheme,
or a decomposition into six tetrahedra.

The fundamental problem is to form a facet approximation
to an isosurface through a scalar field sampled on a rectangular
grid. Given one grid cell defined by its vertices and scalar
values at each vertex, it is necessary to create planar facets
that best represent the isosurface through that grid cell.

The isosurface may not be pass through the grid cell, it may cut off
any one of the vertices, or it may pass through in any one
of a number of more complicated ways. Each possibility will be
characterised by the number of vertices that have values above
or below the isosurface. If one vertex is above the
isosurface say and an adjacent vertex is below the isosurface
then we know the isosurface cuts the edge between these
two vertices.

The position that it cuts the edge will be
linearly interpolated, the ratio of the length between the two vertices
will be the same as the ratio of the isosurface value to the
values at the vertices of the grid cell.

\param box The box.
\param triangles Set of triangles.
\param epsilon Precision.
\param tetrahedra Boolean specifying whether tetrahedral decomposition should be used.
*/
void AnalyticScalarField::Polygonize(const Box& box, QVector<Triangle>& triangles, const double& epsilon, bool tetrahedra) const
{
  Vector vertex[12];

  // Use cube lookup table
  if (tetrahedra == false)
  {
    double x[8];
    for (int i = 0; i < 8; i++)
    {
      x[i] = Value(box.Vertex(i));
    }
    //  Determine the index into the edge table which
    //  tells us which vertices are inside of the surface
    int cubeindex = 0;
    if (x[0] < 0.0) cubeindex |= 1;
    if (x[1] < 0.0) cubeindex |= 2;
    if (x[2] < 0.0) cubeindex |= 4;
    if (x[3] < 0.0) cubeindex |= 8;
    if (x[4] < 0.0) cubeindex |= 16;
    if (x[5] < 0.0) cubeindex |= 32;
    if (x[6] < 0.0) cubeindex |= 64;
    if (x[7] < 0.0) cubeindex |= 128;

    // Cube is entirely in/out of the surface
    if (ScalarField::edgeTable[cubeindex] == 0)
      return;

    // Find the vertices where the surface intersects the cube 
    if (ScalarField::edgeTable[cubeindex] & 1)
      vertex[0] = Dichotomy(box.Vertex(0), box.Vertex(1), x[0], x[1], epsilon);
    if (ScalarField::edgeTable[cubeindex] & 2)
      vertex[1] = Dichotomy(box.Vertex(2), box.Vertex(3), x[2], x[3], epsilon);
    if (ScalarField::edgeTable[cubeindex] & 4)
      vertex[2] = Dichotomy(box.Vertex(4), box.Vertex(5), x[4], x[5], epsilon);
    if (ScalarField::edgeTable[cubeindex] & 8)
      vertex[3] = Dichotomy(box.Vertex(6), box.Vertex(7), x[6], x[7], epsilon);
    if (ScalarField::edgeTable[cubeindex] & 16)
      vertex[4] = Dichotomy(box.Vertex(0), box.Vertex(2), x[0], x[2], epsilon);
    if (ScalarField::edgeTable[cubeindex] & 32)
      vertex[5] = Dichotomy(box.Vertex(1), box.Vertex(3), x[1], x[3], epsilon);
    if (ScalarField::edgeTable[cubeindex] & 64)
      vertex[6] = Dichotomy(box.Vertex(4), box.Vertex(6), x[4], x[6], epsilon);
    if (ScalarField::edgeTable[cubeindex] & 128)
      vertex[7] = Dichotomy(box.Vertex(5), box.Vertex(7), x[5], x[7], epsilon);
    if (ScalarField::edgeTable[cubeindex] & 256)
      vertex[8] = Dichotomy(box.Vertex(0), box.Vertex(4), x[0], x[4], epsilon);
    if (ScalarField::edgeTable[cubeindex] & 512)
      vertex[9] = Dichotomy(box.Vertex(1), box.Vertex(5), x[1], x[5], epsilon);
    if (ScalarField::edgeTable[cubeindex] & 1024)
      vertex[10] = Dichotomy(box.Vertex(2), box.Vertex(6), x[2], x[6], epsilon);
    if (ScalarField::edgeTable[cubeindex] & 2048)
      vertex[11] = Dichotomy(box.Vertex(3), box.Vertex(7), x[3], x[7], epsilon);

    // Create the triangle 
    for (int i = 0; ScalarField::TriangleTable[cubeindex][i] != -1; i += 3)
    {
      triangles.append(Triangle(vertex[ScalarField::TriangleTable[cubeindex][i]],
        vertex[ScalarField::TriangleTable[cubeindex][i + 1]],
        vertex[ScalarField::TriangleTable[cubeindex][i + 2]]));
    }
  }
  // Decompose cube into tetrahedra
  else
  {
    Polygonize(Tetrahedra(box, 0), triangles, epsilon);
    Polygonize(Tetrahedra(box, 1), triangles, epsilon);
    Polygonize(Tetrahedra(box, 2), triangles, epsilon);
    Polygonize(Tetrahedra(box, 3), triangles, epsilon);
    Polygonize(Tetrahedra(box, 4), triangles, epsilon);
    Polygonize(Tetrahedra(box, 5), triangles, epsilon);
  }
}

/*!
\brief Polygonize a tetrahedral cell.
\param tetra Tetrahedral cell.
\param triangles Set of triangles.
\param epsilon Precision.
*/
void AnalyticScalarField::Polygonize(const Tetrahedra& tetra, QVector<Triangle>& triangles, const double& epsilon) const
{
  double x[4];
  for (int i = 0; i < 4; i++)
  {
    x[i] = Value(tetra.Vertex(i));
  }
  // Sign at vertices
  int s[4];

  // Edge locations
  Vector e[6];

  // Determine which of the 16 cases we have given 
  // are above or below the isosurface
  int index = 0;
  s[0] = (x[0] > 0.0);
  s[1] = (x[1] > 0.0);
  s[2] = (x[2] > 0.0);
  s[3] = (x[3] > 0.0);

  if (s[0]) index |= 1;
  if (s[1]) index |= 2;
  if (s[2]) index |= 4;
  if (s[3]) index |= 8;

  // Compute edge vertices and normals
  if (s[0] != s[1])
  {
    e[5] = Dichotomy(tetra.Vertex(0), tetra.Vertex(1), x[0], x[1], epsilon);
  }
  if (s[0] != s[2])
  {
    e[4] = Dichotomy(tetra.Vertex(0), tetra.Vertex(2), x[0], x[2], epsilon);
  }
  if (s[0] != s[3])
  {
    e[2] = Dichotomy(tetra.Vertex(0), tetra.Vertex(3), x[0], x[3], epsilon);
  }
  if (s[1] != s[2])
  {
    e[3] = Dichotomy(tetra.Vertex(1), tetra.Vertex(2), x[1], x[2], epsilon);
  }
  if (s[1] != s[3])
  {
    e[1] = Dichotomy(tetra.Vertex(1), tetra.Vertex(3), x[1], x[3], epsilon);
  }
  if (s[2] != s[3])
  {
    e[0] = Dichotomy(tetra.Vertex(2), tetra.Vertex(3), x[2], x[3], epsilon);
  }

  // Form the vertices of the triangles for each case
  // 14 productive cases (0000 and 1111 do not yield triangles)

  switch (index) {
  case 1:
    triangles.append(Triangle(e[2], e[4], e[5]));
    break;
  case 2:
    triangles.append(Triangle(e[1], e[5], e[3]));
    break;
  case 3:
    triangles.append(Triangle(e[1], e[2], e[4]));
    triangles.append(Triangle(e[1], e[4], e[3]));
    break;
  case 4:
    triangles.append(Triangle(e[0], e[3], e[4]));
    break;
  case 5:
    triangles.append(Triangle(e[0], e[3], e[2]));
    triangles.append(Triangle(e[2], e[3], e[5]));
    break;
  case 6:
    triangles.append(Triangle(e[0], e[1], e[4]));
    triangles.append(Triangle(e[4], e[1], e[5]));
    break;
  case 7:
    triangles.append(Triangle(e[0], e[1], e[2]));
    break;
  case 8:
    triangles.append(Triangle(e[2], e[1], e[0]));
    break;
  case 9:
    triangles.append(Triangle(e[1], e[4], e[5]));
    triangles.append(Triangle(e[4], e[1], e[0]));
    break;
  case 10:
    triangles.append(Triangle(e[0], e[2], e[3]));
    triangles.append(Triangle(e[3], e[2], e[5]));
    break;
  case 11:
    triangles.append(Triangle(e[0], e[4], e[3]));
    break;
  case 12:
    triangles.append(Triangle(e[1], e[3], e[4]));
    triangles.append(Triangle(e[2], e[1], e[4]));
    break;
  case 13:
    triangles.append(Triangle(e[1], e[3], e[5]));
    break;
  case 14:
    triangles.append(Triangle(e[2], e[5], e[4]));
    break;
  }
}

/*!
\brief Voxelize the scalar field.

\param box %The box.
\param n Subdivision.
\param cell Returned reference cube.
\param cubes Set of vectors defining the position of the instantiated reference cube.
*/
void AnalyticScalarField::Voxelize(const Box& box, int n, Box& cell, QVector<Vector>& cubes) const
{
  int nx, ny, nz;

  Box clipped = box;
  clipped.SetParallelepipedic(n, nx, ny, nz);

  // Diagonal of a cell
  Vector d = Array(clipped, nx, ny, nz).CellDiagonal();

  cell = Box(-d, d);

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      for (int k = 0; k < nz; k++)
      {
        Vector p = d / 2.0 + clipped.Vertex(i, j, k, nx, ny, nz);
        if (Value(p) > 0.0)
        {
          cubes.append(p);
        }
      }
    }
  }
}


/*!
\brief Voxelization.

\param n Subdivision.
\param voxel %Voxel.
\param box The box.
*/
void AnalyticScalarField::Voxelize(int n, Voxel& voxel, const Box& box) const
{
  Box clipped = box;
  int nx, ny, nz;
  clipped.SetParallelepipedic(n, nx, ny, nz);

  voxel = Voxel(box, nx, ny, nz);

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      for (int k = 0; k < nz; k++)
      {
        voxel(i, j, k) = Value(voxel.Vertex(i, j, k)) > 0.0 ? 1 : 0;
      }
    }
  }
}


/*!
\brief Computes the volume of the %BlobTree.

Create an octree subdivision of space and check the differents cells
of the octree using a Lipschitz criterion to determine whether cells
are inside, outside or straddling the object.

\param depth Recursion depth in the creation of the octree.
*/
double AnalyticScalarField::Volume(int depth) const
{
  Box box = GetBox();
  box.SetCubic();

  return Volume(box, depth);
}

/*!
\brief Compute the volume of an implicit surface.
\param box Cube where volume computations will be performed.
\param depth Recursion depth.
\param e Error, <I>i.e.</I>, amount of volume that could not be classified as inside or outside.

\sa BlobTree::Volume(const Box&,int) const
*/
double AnalyticScalarField::Volume(const Box& box, int depth, double& e) const
{
  double k = K(box);
  double r = box.Radius();
  double s = Value(box.Center());

  // Lipschitz criterion
  if (fabs(s) > k * r)
  {
    // Box completely inside
    if (Inside(s))
    {
      return box.Volume();
    }
    // Box outside
    else
    {
      return 0.0;
    }
  }
  // Box may be straddling the surface or at least the Lipschitz criterion was not sufficient to classify
  else
  {
    // Terminal case
    if (depth == 0)
    {
      // We are inside the volume, we can at least use the Lipschitz radius to add a little volume
      if (Inside(s))
      {
        double a = Math::Min(r / Math::Sqrt3, fabs(s) / k);
        double v = Sphere::Volume(a);
        e += box.Volume() - v;
        return v;
      }
      else
      {
        e += box.Volume();
        return 0.0;
      }
    }
    // Recurse
    else
    {
      double v = 0.0;
      for (int i = 0; i < 8; i++)
      {
        v += Volume(box.Sub(i), depth - 1, e);
      }
      return v;
    }
  }
}

/*!
\brief Compute the volume of a BlobTree object.
\param box Cube where volume computations will be performed.
\param depth Recursion depth.
*/
double AnalyticScalarField::Volume(const Box& box, int depth) const
{
  double e;
  return Volume(box, depth, e);
}

/*!
\brief Computes the gravity center of the %BlobTree.

This function returns a sphere, the center is the center of gravity, and the radius is the volume.

\param depth Recursion depth.
*/
Sphere AnalyticScalarField::Center(int depth) const
{
  Box box = GetBox();
  box.SetCubic();

  // Use sphere to keep track of center and volume  
  Sphere s = Center(box, depth);

  return s;
}

/*!
\brief Compute the center of gravity of the %BlobTree.
\param box Initial box.
\param depth Recursion depth.
\sa BlobTree::Center(int)
*/
Sphere AnalyticScalarField::Center(const Box& box, int depth) const
{
  Vector c = box.Center();
  double k = K(box);
  double s = Value(c);
  double r = box.Radius();

  // Lipschitz criterion
  if (fabs(s) > k * r)
  {
    // Box completely inside
    if (Inside(s))
    {
      double v = box.Volume();
      return Sphere(c * v, v);
    }
    else
    {
      return Sphere(0.0);
    }
  }
  else
  {
    if (depth == 0)
    {
      if (Inside(s))
      {
        double a = Math::Min(r / Math::Sqrt3, fabs(s) / k);
        double v = Sphere::Volume(a);
        return Sphere(c * v, v);
      }
      else
      {
        return Sphere(0.0);
      }
    }
    else
    {
      Sphere e(0.0);
      for (int i = 0; i < 8; i++)
      {
        e = e + Center(box.Sub(i), depth - 1);
      }
      return e;
    }
  }
}
