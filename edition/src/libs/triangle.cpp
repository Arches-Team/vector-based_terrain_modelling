// Triangle

#include "libs/triangle.h"
#include "libs/plane.h"
#include "libs/segment.h"

// Used for drawing
#include <QtGui/QPolygonF>


/*!
\class Triangle triangle.h
\brief Base minimum storage triangle class.

The data-structure of the triangle does not include the normalized normal vector.

\ingroup KernelGroup
*/

/*!
\brief Compute a point in the triangle, given uv-coordinates.
\param u,v Coordinates.
*/
Vector Triangle::Vertex(const double& u, const double& v) const
{
  return p[0] + u * (p[1] - p[0]) + v * (p[2] - p[0]);
}

/*!
\brief Shrinks a triangle.

This function erodes a triangle by moving its edges by a given distance.
\param e Erosion radius.
*/
void Triangle::Shrink(const double& e)
{
  Vector normal = Normal();

  // Plane of the triangle
  Plane plane(normal, p[0]);

  // Planes orthogonal to the edges
  Plane ab(normal / Normalized((p[1] - p[0])), p[0]);
  Plane bc(normal / Normalized((p[2] - p[1])), p[1]);
  Plane ca(normal / Normalized((p[0] - p[2])), p[2]);

  // Move
  ab.Translate(e);
  bc.Translate(e);
  ca.Translate(e);

  // Points
  p[0] = Plane::Intersection(plane, ca, ab);
  p[1] = Plane::Intersection(plane, ab, bc);
  p[2] = Plane::Intersection(plane, bc, ca);
}

/*!
\brief Compute the normal vector of the triangle.

This function is expensive as it requires normalizing the cross product of the edge vectors.
\sa TriangleEdge
*/
Vector Triangle::Normal() const
{
  return Normalized((p[1] - p[0]) / (p[2] - p[0]));
}

/*!
\brief Compute the normal vector of the triangle, and scale the normal using its area.

This function is less expensive than Triangle::Normal() as it does not require normalizing the cross product of the edge vectors.

IIt is useful for computing the vertex normals of a triangle mesh by averaging the face normals of the triangles and weighting
the influence of the different triangles with their area.
\sa Triangle::Normal()
*/
Vector Triangle::AreaNormal() const
{
  return 0.5 * ((p[1] - p[0]) / (p[2] - p[0]));
}

/*!
\brief Compute the intersection between a ray and a triangle.

Return the barycentric coordinates of the intersection
if within the triangle.

After Tomas Moller and Ben Trumbore,
<I>Fast, minimum storage ray-triangle intersection</I>,
<B>Journal of graphics tools</B>, 2(1):21-28, 1997.

\param ray The ray (direction should be of unit length).
\param t Intersection depth.
\param u,v Parametric coordinates of the intersection depth in the triangle.
*/
bool Triangle::Intersect(const Ray& ray, double& t, double& u, double& v) const
{
  Vector e[2];

  // Find edge vectors
  e[0] = p[1] - p[0];
  e[1] = p[2] - p[0];

  Vector pvec = ray.Direction() / e[1];

  // Determinant, if determinant is near zero, ray lies in plane of triangle
  double det = e[0] * pvec;

  if ((det > -epsilon) && (det < epsilon))
    return false;
  det = 1.0 / det;

  // calculate distance from first vertex to ray origin
  Vector tvec = ray.Origin() - p[0];

  // calculate U parameter and test bounds
  u = (tvec * pvec) * det;
  if ((u < 0.0) || (u > 1.0))
    return false;

  // calculate V parameter and test bounds
  Vector qvec = tvec / e[0];

  v = (ray.Direction() * qvec) * det;
  if ((v < 0.0) || (u + v) > 1.0)
    return false;

  // Ray intersects triangle
  t = (e[1] * qvec) * det;
  return true;
}

/*!
\brief Compute the intersection between a ray
and a triangle.

After Tomas Moller and Ben Trumbore,
<I>Fast, minimum storage ray-triangle intersection</I>,
<B>Journal of graphics tools</B>, 2(1):21-28, 1997.
\param ray The ray (direction should be of unit length).
\param t Intersection depth.
*/
bool Triangle::Intersect(const Ray& ray, double& t) const
{
  Vector e[2];

  // Find edge vectors
  e[0] = p[1] - p[0];
  e[1] = p[2] - p[0];

  // Determinant
  Vector pvec = ray.Direction() / e[1];

  // If determinant is near zero, ray lies in plane of triangle
  double det = e[0] * pvec;

  if ((det > -epsilon) && (det < epsilon))
    return false;
  det = 1.0 / det;

  // Distance from first vertex to ray origin
  Vector tvec = ray.Origin() - p[0];

  // Calculate U parameter and test bounds
  double u = (tvec * pvec) * det;
  if ((u < 0.0) || (u > 1.0))
    return false;

  // Calculate V parameter and test bounds
  Vector qvec = tvec / e[0];

  double v = (ray.Direction() * qvec) * det;
  if ((v < 0.0) || (u + v) > 1.0)
    return false;

  // Ray intersects triangle
  t = (e[1] * qvec) * det;

  return true;
}

/*!
\brief Translates a triangle by a given vector.

\param u Translation vector.
*/
void Triangle::Translate(const Vector& u)
{
  for (int i = 0; i < 3; i++)
  {
    p[i] += u;
  }
}

/*!
\brief Rotates a triangle.

\param r Rotation matrix.
*/
void Triangle::Rotate(const Matrix& r)
{
  for (int i = 0; i < 3; i++)
  {
    p[i] = r * p[i];
  }
}

/*!
\brief Scale a triangle.

\param u Scaling vector.
*/
void Triangle::Scale(const Vector& u)
{
  for (int i = 0; i < 3; i++)
  {
    p[i] *= u;
  }
}

/*!
\brief Transforms a triangle given a frame transformation.
\param t Transformation.
*/
void Triangle::Transform(const Frame& t)
{
  for (int i = 0; i < 3; i++)
  {
    p[i] = t.Transform(p[i]);
  }
}

/*!
\brief Transforms a triangle given a frame transformation.
\param t Transformation.
*/
Triangle Triangle::Transformed(const FrameScaled& t) const
{
  return Triangle(t.Transform(p[0]), t.Transform(p[1]), t.Transform(p[2]));
}

/*!
\brief Computes the axis aligned box enclosing the triangle.
*/
Box Triangle::GetBox() const
{
  return Box(Vector::Min(Vector::Min(p[0], p[1]), p[2]), Vector::Max(Vector::Max(p[0], p[1]), p[2]));
}

/*!
\brief Compute the radius of the circle inscribed in the
triangle.

Some algebra will show that it is half the ratio
of the half perimeter and the surface of the triangle.
*/
double Triangle::InscribedRadius() const
{
  Vector u = p[0] - p[1];
  Vector v = p[2] - p[0];
  Vector w = p[1] - p[2];
  double a = Norm(u) + Norm(v) + Norm(w);
  double s = Norm(u / v);
  return s / a;
}

/*!
\brief Compute the radius of the circumscribed circle of the triangle.
*/
double Triangle::CircumscribedRadius() const
{
  double u = Norm(p[0] - p[1]);
  double v = Norm(p[1] - p[2]);
  double w = Norm(p[2] - p[0]);
  return u * v * w / sqrt((u + v + w) * (-u + v + w) * (u - v + w) * (u + v - w));
}

/*!
\brief Compute the sphere inscribed in the triangle.

The center of the sphere
is the center of (A,a), (B,b) and (C,c) where A, B and C are the vertices of
the triangle and coefficients a, b and c represent the length of their facing
edge.
*/
Sphere Triangle::Inscribed() const
{
  Vector u = p[0] - p[1];
  Vector v = p[2] - p[0];
  Vector w = p[1] - p[2];
  double a = Norm(w);
  double b = Norm(v);
  double c = Norm(u);
  double l = a + b + c;
  return Sphere((p[0] * a + p[1] * b + p[2] * c) / l, Norm(u / v) / l);
}

/*!
\brief Compute the sphere circumscribing the triangle.

This function directly relies on the constructor of the class
Sphere.

\sa Sphere::Sphere(const Vector&, const Vector&, const Vector&);
*/
Sphere Triangle::Circumscribed() const
{
  return Sphere(p[0], p[1], p[2]);
}

/*!
\brief Computes the aspect ratio of the triangle.

It is defined as the radius of the inscribed circle divided by the radius of the
circumscribing circle. This function is faster than:

\code
Triangle t;
double a=2.0*t.InscribedRadius()/t.CircumscribedRadius();
\endcode
*/
double Triangle::Aspect() const
{
  double ab = Norm(p[1] - p[0]);
  double bc = Norm(p[2] - p[1]);
  double ca = Norm(p[0] - p[2]);

  double s = 0.5 * (ab + bc + ca);

  double u = (s - ab) * (s - bc) * (s - ca);

  return 8.0 * u / (ab * bc * ca);
}

/*!
\brief Overloaded.
\param t %Triangle.
\param s Stream.
*/
std::ostream& operator<<(std::ostream& s, const Triangle& t)
{
  s << "Triangle(" << t.p[0] << ',' << t.p[1] << ',' << t.p[2] << ')';
  return s;
}

/*!
\brief Computes the squared distance between a point and a triangle.
\param q Point.
*/
double Triangle::R(const Vector& q) const
{
  // Those parameters are constant for a given triangle, and could be preprocessed.
  Vector ab = p[1] - p[0];
  Vector ac = p[2] - p[0];
  double abab = ab * ab;
  double abac = ab * ac;
  double acac = ac * ac;
  double delta = fabs(abab * acac - abac * abac);

  Vector pa = p[0] - q;
  double paab = pa * ab;
  double paac = pa * ac;
  double papa = pa * pa;
  double fS = abac * paac - acac * paab;
  double fT = abac * paab - abab * paac;
  double r;

  if (fS + fT <= delta)
  {
    if (fS < 0.0)
    {
      if (fT < 0.0)  // region 4
      {
        if (paab < 0.0)
        {
          fT = 0.0;
          if (-paab >= abab)
          {
            fS = 1.0;
            r = abab + 2.0 * paab + papa;
          }
          else
          {
            fS = -paab / abab;
            r = paab * fS + papa;
          }
        }
        else
        {
          fS = 0.0;
          if (paac >= 0.0)
          {
            fT = 0.0;
            r = papa;
          }
          else if (-paac >= acac)
          {
            fT = 1.0;
            r = acac + 2.0 * paac + papa;
          }
          else
          {
            fT = -paac / acac;
            r = paac * fT + papa;
          }
        }
      }
      else  // region 3
      {
        fS = 0.0;
        if (paac >= 0.0)
        {
          fT = 0.0;
          r = papa;
        }
        else if (-paac >= acac)
        {
          fT = 1.0;
          r = acac + 2.0 * paac + papa;
        }
        else
        {
          fT = -paac / acac;
          r = paac * fT + papa;
        }
      }
    }
    else if (fT < 0.0)  // region 5
    {
      fT = 0.0;
      if (paab >= 0.0)
      {
        fS = 0.0;
        r = papa;
      }
      else if (-paab >= abab)
      {
        fS = 1.0;
        r = abab + 2.0 * paab + papa;
      }
      else
      {
        fS = -paab / abab;
        r = paab * fS + papa;
      }
    }
    else  // region 0
    {
      // minimum at interior point
      double fInvDet = 1.0 / delta;
      fS *= fInvDet;
      fT *= fInvDet;
      r = fS * (abab * fS + abac * fT + 2.0 * paab) + fT * (abac * fS + acac * fT + 2.0 * paac) + papa;
    }
  }
  else
  {
    double fTmp0, fTmp1, nu, de;

    if (fS < 0.0)  // region 2
    {
      fTmp0 = abac + paab;
      fTmp1 = acac + paac;
      if (fTmp1 > fTmp0)
      {
        nu = fTmp1 - fTmp0;
        de = abab - 2.0 * abac + acac;
        if (nu >= de)
        {
          fS = 1.0;
          fT = 0.0;
          r = abab + 2.0 * paab + papa;
        }
        else
        {
          fS = nu / de;
          fT = 1.0 - fS;
          r = fS * (abab * fS + abac * fT + 2.0 * paab) + fT * (abac * fS + acac * fT + 2.0 * paac) + papa;
        }
      }
      else
      {
        fS = 0.0;
        if (fTmp1 <= 0.0)
        {
          fT = 1.0;
          r = acac + 2.0 * paac + papa;
        }
        else if (paac >= 0.0)
        {
          fT = 0.0;
          r = papa;
        }
        else
        {
          fT = -paac / acac;
          r = paac * fT + papa;
        }
      }
    }
    else if (fT < 0.0)  // region 6
    {
      fTmp0 = abac + paac;
      fTmp1 = abab + paab;
      if (fTmp1 > fTmp0)
      {
        nu = fTmp1 - fTmp0;
        de = abab - 2.0 * abac + acac;
        if (nu >= de)
        {
          fT = 1.0;
          fS = 0.0;
          r = acac + 2.0 * paac + papa;
        }
        else
        {
          fT = nu / de;
          fS = 1.0 - fT;
          r = fS * (abab * fS + abac * fT + 2.0 * paab) + fT * (abac * fS + acac * fT + 2.0 * paac) + papa;
        }
      }
      else
      {
        fT = 0.0;
        if (fTmp1 <= 0.0)
        {
          fS = 1.0;
          r = abab + 2.0 * paab + papa;
        }
        else if (paab >= 0.0)
        {
          fS = 0.0;
          r = papa;
        }
        else
        {
          fS = -paab / abab;
          r = paab * fS + papa;
        }
      }
    }
    else  // region 1
    {
      nu = acac + paac - abac - paab;
      if (nu <= 0.0)
      {
        fS = 0.0;
        fT = 1.0;
        r = acac + 2.0 * paac + papa;
      }
      else
      {
        de = abab - 2.0 * abac + acac;
        if (nu >= de)
        {
          fS = 1.0;
          fT = 0.0;
          r = abab + 2.0 * paab + papa;
        }
        else
        {
          fS = nu / de;
          fT = 1.0 - fS;
          r = fS * (abab * fS + abac * fT + 2.0 * paab) + fT * (abac * fS + acac * fT + 2.0 * paac) + papa;
        }
      }
    }
  }

  return fabs(r);
}

/*!
\brief Computes the normal vector to a triangle.
\param q Point.
*/
Vector Triangle::Normal(const Vector& q) const
{
  // Those parameters are constant for a given triangle, and could be preprocessed.
  Vector ab = p[1] - p[0];
  Vector ac = p[2] - p[0];
  double abab = ab * ab;
  double abac = ab * ac;
  double acac = ac * ac;
  double delta = fabs(abab * acac - abac * abac);

  Vector pa = p[0] - q;
  double paab = pa * ab;
  double paac = pa * ac;
  //double papa = pa * pa;
  double fS = abac * paac - acac * paab;
  double fT = abac * paab - abab * paac;

  if (fS + fT <= delta)
  {
    if (fS < 0.0)
    {
      if (fT < 0.0)  // region 4
      {
        if (paab < 0.0)
        {
          fT = 0.0;
          if (-paab >= abab)
          {
            fS = 1.0;
          }
          else
          {
            fS = -paab / abab;
          }
        }
        else
        {
          fS = 0.0;
          if (paac >= 0.0)
          {
            fT = 0.0;
          }
          else if (-paac >= acac)
          {
            fT = 1.0;
          }
          else
          {
            fT = -paac / acac;
          }
        }
      }
      else  // region 3
      {
        fS = 0.0;
        if (paac >= 0.0)
        {
          fT = 0.0;
        }
        else if (-paac >= acac)
        {
          fT = 1.0;
        }
        else
        {
          fT = -paac / acac;
        }
      }
    }
    else if (fT < 0.0)  // region 5
    {
      fT = 0.0;
      if (paab >= 0.0)
      {
        fS = 0.0;
      }
      else if (-paab >= abab)
      {
        fS = 1.0;
      }
      else
      {
        fS = -paab / abab;
      }
    }
    else  // region 0
    {
      // minimum at interior point
      double fInvDet = 1.0 / delta;
      fS *= fInvDet;
      fT *= fInvDet;
    }
  }
  else
  {
    double fTmp0, fTmp1, nu, de;

    if (fS < 0.0)  // region 2
    {
      fTmp0 = abac + paab;
      fTmp1 = acac + paac;
      if (fTmp1 > fTmp0)
      {
        nu = fTmp1 - fTmp0;
        de = abab - 2.0 * abac + acac;
        if (nu >= de)
        {
          fS = 1.0;
          fT = 0.0;
        }
        else
        {
          fS = nu / de;
          fT = 1.0 - fS;
        }
      }
      else
      {
        fS = 0.0;
        if (fTmp1 <= 0.0)
        {
          fT = 1.0;
        }
        else if (paac >= 0.0)
        {
          fT = 0.0;
        }
        else
        {
          fT = -paac / acac;
        }
      }
    }
    else if (fT < 0.0)  // region 6
    {
      fTmp0 = abac + paac;
      fTmp1 = abab + paab;
      if (fTmp1 > fTmp0)
      {
        nu = fTmp1 - fTmp0;
        de = abab - 2.0 * abac + acac;
        if (nu >= de)
        {
          fT = 1.0;
          fS = 0.0;
        }
        else
        {
          fT = nu / de;
          fS = 1.0 - fT;
        }
      }
      else
      {
        fT = 0.0;
        if (fTmp1 <= 0.0)
        {
          fS = 1.0;
        }
        else if (paab >= 0.0)
        {
          fS = 0.0;
        }
        else
        {
          fS = -paab / abab;
        }
      }
    }
    else  // region 1
    {
      nu = acac + paac - abac - paab;
      if (nu <= 0.0)
      {
        fS = 0.0;
        fT = 1.0;
      }
      else
      {
        de = abab - 2.0 * abac + acac;
        if (nu >= de)
        {
          fS = 1.0;
          fT = 0.0;
        }
        else
        {
          fS = nu / de;
          fT = 1.0 - fS;
        }
      }
    }
  }
  return q - (p[0] + fS * ab + fT * ac);
}

/*!
\brief Create an n-adic subdivision of a triangle.

\param n Subdivision level.
\param vertex Array of points.
\param index Array of indexes defining the triangles.
*/
void Triangle::Subdivide(int n, QVector<Vector>& vertex, QVector<int>& index) const
{
  // The triangle will be split into n strips
  // Let e a given strip, every strip is composed of 2e-1 triangles
  Vector e01 = (p[1] - p[0]) / n;
  Vector e12 = (p[2] - p[1]) / n;

  // Starting point.
  vertex.append(p[0]);
  int k = 0;

  //pour chaque Ã©tage sauf le dernier
  for (int j = 1; j <= n; j++)
  {
    for (int i = 0; i <= j; i++)
    {
      vertex.append(p[0] + e01 * j + e12 * i);
      k++;

      if (i != j)
      {
        // Counter clockwize triangle
        index.append(k);
        index.append(k + 1);
        index.append(k - j);

        if (j != n)
        {
          // Clockwize triangle
          index.append(k);
          index.append(k + j + 2);
          index.append(k + 1);
        }
      }
    }
  }
}

/*!
\brief Compute the barycentric coordinates of the projection of the point onto the plane of the triangle.

\param p Point.

\sa Triangle2::BarycentricCoordinates
*/
Vector Triangle::BarycentricCoordinates(const Vector& p) const
{
  // Create a frame
  Frame frame = Frame::Orthonormal(Center(), Normal());

  // Transform the triangle, compute its projection 
  Triangle2 t=Triangle2(Transformed(frame.Inverse()));

  // Apply the same inverse transformation to the point
  Vector2 q = Vector2(frame.InverseTransform(p));

  return t.BarycentricCoordinates(q);
}

/*!
\brief Compute the intersection between a triangle and a plane.

\param plane The plane.
\param segment Returned segment if intersection occurs.
\return The number of intersections, denoted as n. If n=1 then only one vertex of the triangle is one the plane, whereas if n=3 then the triangle is on the plane.
\sa int Triangle::Intersect(const double&,Segment&)
*/
int Triangle::Intersect(const Plane& plane, Segment& segment) const
{
  // Side
  int sa = plane.Side(p[0]);
  int sb = plane.Side(p[1]);
  int sc = plane.Side(p[2]);

  int s = abs(sa + sb + sc);

  // All points are on the same side of the plane
  if (s == 3)
  {
    return 0;
  }
  // One point lies on the plane
  else if (s == 2)
  {
    if (sa == 0)
    {
      segment = Segment(p[0], p[0]);
    }
    else if (sb == 0)
    {
      segment = Segment(p[1], p[1]);
    }
    else
    {
      segment = Segment(p[2], p[2]);
    }
    return 1;
  }
  // Either two points are on and one on one side (0 0 1) which means one edge of the triangle is on the plane, either two on one side and one on the other (1 1 -1) which means true intersection
  else if (s == 1)
  {
    if ((sa == 0) && (sb == 0))
    {
      segment = Segment(p[0], p[1]);
    }
    else if ((sa == 0) && (sc == 0))
    {
      segment = Segment(p[0], p[2]);
    }
    else if ((sb == 0) && (sc == 0))
    {
      segment = Segment(p[1], p[2]);
    }
    // True intersection
    else
    {
      if (sa == sb)
      {
        segment = Segment(plane.Intersect(Line(p[0], p[2])), plane.Intersect(Line(p[1], p[2])));
      }
      else if (sa == sc)
      {
        segment = Segment(plane.Intersect(Line(p[0], p[1])), plane.Intersect(Line(p[2], p[1])));
      }
      else //if (sb==sc)
      {
        segment = Segment(plane.Intersect(Line(p[1], p[0])), plane.Intersect(Line(p[2], p[0])));
      }
    }
    return 2;
  }
  // Either all points are on the plane, or one is on the plane and two on opposite sides
  else
  {
    if ((sa == 0) && (sb == 0) && (sc == 0))
    {
      return 3;
    }
    else
    {
      if (sa == 0)
      {
        segment = Segment(plane.Intersect(Line(p[1], p[2])), p[0]);
      }
      else if (sb == 0)
      {
        segment = Segment(plane.Intersect(Line(p[0], p[2])), p[1]);
      }
      else
      {
        segment = Segment(plane.Intersect(Line(p[0], p[1])), p[2]);
      }
      return 2;
    }
  }
}

/*!
\brief Compute the intersection between a triangle and a horizontal plane.

\param z Height of horizontal plane.
\param segment Returned segment if intersection occurs.
\return The number of intersections, denoted as n. If n=1 then only one vertex of the triangle is one the plane, whereas if n=3 then the triangle is on the plane.
\sa int Triangle::Intersect(const double&,Segment&)
*/
int Triangle::Intersect(const double& z, Segment& segment) const
{
  return Intersect(Plane(Vector::Z, z), segment);
}

/*!
\brief Check if a sphere intersects or contains a triangle.

It does not attempt to compute the set of intersection.

If a triangle vertex lies exactly on the sphere, then no intersection occurs.
Apart from this limit case, the algorithm proceeds as follows.
<BR><B>1.</B> All three triangle vertices are contained
in the sphere. The sphere completely contains the triangle.
<BR><B>2.</B> At least one vertex is inside the sphere and at least one vertex is
outside the sphere. The sphere and triangle intersect.
<BR><B>3.</B> All three vertices are outside the sphere. The sphere and triangle
intersect when the distance from sphere center to triangle is less
than the sphere radius.

\param sphere The sphere.
*/
bool Triangle::Intersect(const Sphere& sphere) const
{
  double r = sphere.Radius() * sphere.Radius();
  int inside = 0;
  Vector TS;

  // Test if vertices are inside the sphere
  for (int i = 0; i < 3; i++)
  {
    TS = p[i] - sphere.Center();
    if (TS * TS <= r)
    {
      inside++;
    }
  }

  // Triangle does not intersect sphere
  if (inside == 3)
  {
    return false;
  }

  // Triangle transversely intersects sphere
  if (inside > 0)
  {
    return true;
  }

  // All vertices are outside the sphere. Determine the minimum squared
  // distance between the sphere center and the triangle. 

  Vector edge0 = p[1] - p[0];
  Vector edge1 = p[2] - p[0];
  double A = edge0 * edge0;
  double B = edge0 * edge1;
  double C = edge1 * edge1;
  double D = edge0 * TS;
  double E = edge1 * TS;
  double F = TS * TS;
  double det = A * C - B * B;  // Determinant should be not null for triangles
  double invdet = 1.0 / det;
  double s = (B * E - C * D) * invdet;
  double t = (B * D - A * E) * invdet;

  if (s + t <= 1.0)
  {
    if (s < 0.0)
    {
      if (t < 0.0)  // region 4
      {
        if (D < 0)
        {
          t = 0.0;
          s = -D / A;
          if (s > 1.0) { s = 1.0; }
        }
        else if (E < 0)
        {
          s = 0.0;
          t = -E / C;
          if (t > 1.0) { t = 1.0; }
        }
        else
        {
          s = 0.0;
          t = 0.0;
        }
      }
      else  // region 3
      {
        s = 0.0;
        t = -E / C;
        if (t < 0.0) { t = 0.0; }
        else if (t > 1.0) { t = 1.0; }
      }
    }
    else if (t < 0.0)  // region 5
    {
      t = 0.0;
      s = -D / A;
      if (s < 0.0) { s = 0.0; }
      else if (s > 1.0) { s = 1.0; }
    }
    else  // region 0
    {
      // minimum at interior point
    }
  }
  else
  {
    if (s < 0.0)  // region 2
    {
      if (B - C + D - E < 0.0)
      {
        s = -(B - C + D - E) / (A - 2 * B + C);
        if (s < 0.0) { s = 0.0; }
        else if (s > 1.0) { s = 1.0; }
        t = 1.0 - s;
      }
      else if (C + E > 0.0)
      {
        s = 0.0;
        t = -E / C;
        if (t < 0.0) { t = 0.0; }
        else if (t > 1.0) { t = 1.0; }
      }
      else
      {
        s = 0.0;
        t = 1.0;
      }
    }
    else if (t < 0.0)  // region 6
    {
      if (A - B + D - E > 0.0)
      {
        t = (A - B + D - E) / (A - 2 * B + C);
        if (t < 0.0) { t = 0.0; }
        else if (t > 1.0) { t = 1.0; }
        s = 1.0 - t;
      }
      else if (A + D > 0.0)
      {
        t = 0.0;
        s = -D / A;
        if (s < 0.0) { s = 0.0; }
        else if (s > 1.0) { s = 1.0; }
      }
      else
      {
        s = 1.0;
        t = 0.0;
      }
    }
    else  // region 1
    {
      s = -(B - C + D - E) / (A - 2 * B + C);
      if (s < 0.0) { s = 0.0; }
      else if (s > 1.0) { s = 1.0; }
      t = 1.0 - s;
    }
  }
  double qz = s * (A * s + B * t + 2 * D) + t * (B * s + C * t + 2 * E) + F;

  if (qz < r)
  {
    return true;
  }
  else
  {
    return false;
  }
}

/*!
\brief Generate a random vector inside the triangle.

After R. Osada, T. Funkhouser, B. Chazelle, D. Dobkin. Shape Distributions,
<I>ACM Transactions on Graphics</I>, <B>21</B>(4), 807-832, 2002.

\param random %Random number generator.
*/
Vector Triangle::RandomInside(Random& random) const
{
  double u = sqrt(random.Uniform());
  double v = random.Uniform();

  return (1.0 - u) * p[0] + u * ((1.0 - v) * p[1] + v * p[2]);
}