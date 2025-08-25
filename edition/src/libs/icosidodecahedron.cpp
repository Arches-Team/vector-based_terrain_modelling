// Icosidodecahedron

// Self include
#include "libs/icosidodecahedron.h"
#include "libs/plane.h"

/*!
\class Icosidodecahedron dodecahedron.h
\brief An icosidodecahedron.

\sa Dodecahedron, Icosahedron
\ingroup ExtendedKernelGroup
*/
const double Icosidodecahedron::Phi = (1 + sqrt(5.0)) / 2.0;

const Vector Icosidodecahedron::vertex[30] = {
  Vector(Phi, 0, 0), Vector(-Phi, 0, 0), Vector(0, Phi, 0), Vector(0, -Phi, 0), Vector(0, 0, Phi), Vector(0, 0, -Phi),
  Vector(Phi / 2, Phi * Phi / 2, 0.5), Vector(-Phi / 2, Phi * Phi / 2, 0.5), Vector(Phi / 2, -Phi * Phi / 2, 0.5), Vector(-Phi / 2, -Phi * Phi / 2, 0.5),
  Vector(Phi / 2, Phi * Phi / 2, -0.5), Vector(-Phi / 2, Phi * Phi / 2, -0.5), Vector(Phi / 2, -Phi * Phi / 2, -0.5), Vector(-Phi / 2, -Phi * Phi / 2, -0.5),
  Vector(Phi * Phi / 2, 0.5, Phi / 2), Vector(-Phi * Phi / 2, 0.5, Phi / 2), Vector(Phi * Phi / 2, -0.5, Phi / 2), Vector(-Phi * Phi / 2, -0.5, Phi / 2),
  Vector(Phi * Phi / 2, 0.5, -Phi / 2), Vector(-Phi * Phi / 2, 0.5, -Phi / 2), Vector(Phi * Phi / 2, -0.5, -Phi / 2), Vector(-Phi * Phi / 2, -0.5, -Phi / 2),
  Vector(0.5, Phi / 2, Phi * Phi / 2), Vector(-0.5, Phi / 2, Phi * Phi / 2), Vector(0.5, -Phi / 2, Phi * Phi / 2), Vector(-0.5, -Phi / 2, Phi * Phi / 2),
  Vector(0.5, Phi / 2, -Phi * Phi / 2), Vector(-0.5, Phi / 2, -Phi * Phi / 2), Vector(0.5, -Phi / 2, -Phi * Phi / 2), Vector(-0.5, -Phi / 2, -Phi * Phi / 2)
};

const Vector Icosidodecahedron::normal[32] = {
  Vector(1.0 / Phi, 0, 1) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(0, 1, 1.0 / Phi) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(-1.0 / Phi, 0, 1) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(0, -1, 1.0 / Phi) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(1.0 / Phi, 0, -1) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(0, -1, -1.0 / Phi) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(-1.0 / Phi, 0, -1) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(0, 1, -1.0 / Phi) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(1.0 / Phi, 0, 1) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(0, 1, 1.0 / Phi) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(-1.0 / Phi, 0, 1) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(0, -1, 1.0 / Phi) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(1.0 / Phi, 0, -1) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(0, -1, -1.0 / Phi) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(-1.0 / Phi, 0, -1) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(0, 1, -1.0 / Phi) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(1.0 / Phi, 0, 1) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(0, 1, 1.0 / Phi) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(-1.0 / Phi, 0, 1) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(0, -1, 1.0 / Phi) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(1.0 / Phi, 0, -1) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(0, -1, -1.0 / Phi) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(-1.0 / Phi, 0, -1) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(0, 1, -1.0 / Phi) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(1.0 / Phi, 0, 1) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(0, 1, 1.0 / Phi) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(-1.0 / Phi, 0, 1) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(0, -1, 1.0 / Phi) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(1.0 / Phi, 0, -1) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(0, -1, -1.0 / Phi) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(-1.0 / Phi, 0, -1) / sqrt(1 + 1.0 / (Phi * Phi)),
  Vector(0, 1, -1.0 / Phi) / sqrt(1 + 1.0 / (Phi * Phi))
};

const int Icosidodecahedron::pentagon[12][5] = {
  { 0, 1, 15, 4, 12 }, { 0, 12, 8, 9, 13 }, { 0, 13, 5, 14, 1 }, { 1, 14, 10, 11, 15 }, { 2, 3, 17, 7, 16 },
  { 2, 16, 11, 10, 19 }, { 2, 19, 6, 18, 3 }, { 18, 9, 8, 17, 3 }, { 15, 11, 16, 7, 4 }, { 4, 7, 17, 8, 12 },
  { 13, 9, 18, 6, 5 }, { 5, 6, 19, 10, 14 }
};
/*
0 14 16
3 24 25
1 19 21
0 18 14
19 27 5
1 15 19
18 10 14
2 10 27
0 16 20
0 20 18
28 18 20
5 26 28
28 26 18
27 26 5
18 26 10
10 26 27
16 12 20
3 28 12
28 20 12
3 29 28
19 29 21
19 5 29
5 28 29
3 8 24
16 24 8
16 8 12
3 12 8
3 25 9
10 6 14
2 6 10
2 23 6
23 22 6
14 6 22
14 4 16
14 22 4
16 4 24
24 4 25
23 25 4
23 4 22
15 17 23
23 17 25
1 21 17
9 17 21
1 17 15
25 17 9
19 11 27
2 27 11
15 11 19
3 13 29
29 13 21
3 9 13
9 21 13
2 7 23
15 23 7
2 11 7
15 7 11
*/
const int Icosidodecahedron::triangle[20][3] = {
  { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
  { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
  { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
  { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }
};

const double Icosidodecahedron::R2 = sqrt(3.0);

/*!
\brief Create an icosidodecahedron given its center and radius.
\param c Center.
\param r Radius.
*/
Icosidodecahedron::Icosidodecahedron(const Vector& c, const double& r) :c(c), r(r)
{
}

/*!
\brief Create an icosidodecahedron given its radius.
\param r Radius.
*/
Icosidodecahedron::Icosidodecahedron(const double& r) :Icosidodecahedron(Vector::Null, r)
{
}

/*!
\brief Get all vertexes.
*/
QVector<Vector> Icosidodecahedron::Vertexes() const
{
  QVector<Vector> v(30);
  for (int i = 0; i < 30; i++)
  {
    v[i] = Vertex(i);
  }
  return v;
}

/*!
\brief Return the plane of the k-th face.
\param k Face, pentagon if k in [0,11], and triangle if k in [12,31]
*/
Plane Icosidodecahedron::GetPlane(int k) const
{
  // Pentagons
  if (k < 12)
  {
    return Plane(normal[k], c + r * vertex[pentagon[k][0]]);
  }
  // Triangles
  else
  {
    return Plane(normal[k], c + r * vertex[triangle[k - 12][0]]);
  }
}

/*!
\brief Compute the volume.
*/
double Icosidodecahedron::Volume() const
{
  // Edge length
  double u = r / R2;
  return (1.0 / 6.0) * (45.0 + 17.0 * sqrt(5.0)) * u * u * u;
}

/*!
\brief Compute the area.
*/
double Icosidodecahedron::Area() const
{
  double u = r / R2;
  return (5.0 * sqrt(3.0) + 3.0 * sqrt(25.0 + 10.0 * sqrt(5.0))) * u * u;
}


/*!
\brief Compute the intersection with a ray.

\param ray The ray.
\param ta, tb Intersection depths.
\param fa, fb Index of the intersected faces.
*/
bool Icosidodecahedron::Intersect(const Ray& ray, double& ta, double& tb, int& fa, int& fb) const
{
  // Intersection depths with planes
  ta = -Math::Infinity;
  tb = Math::Infinity;

  // Indexes
  fa = -1;
  fb = -1;

  for (int i = 0; i < 32; i++)
  {
    double t;
    if (GetPlane(i).Intersect(ray, t))
    {
      if ((Normal(i) * ray.Direction()) < 0.0)
      {
        if (ta > t)
        {
          t = ta;
          fa = i;
        }
      }
      else
      {
        if (tb < t)
        {
          t = tb;
          fb = i;
        }
      }
    }
    if (ta > tb) return false;
  }
  return true;
}

/*!
\brief Compute the approximate signed distance to the icosahedron.
\param p Point.
*/
double Icosidodecahedron::Signed(const Vector& p) const
{
  double s = GetPlane(0).Signed(p);
  for (int i = 1; i < 32; i++)
  {
    s = Math::Max(s, GetPlane(i).Signed(p));
  }
  return s;
}

/*!
\brief Check if a point is inside.

\param p The point.
*/
bool Icosidodecahedron::Inside(const Vector& p) const
{
  for (int i = 0; i < 32; i++)
  {
    if (!Plane(Normal(i), Vertex(i, 0)).Inside(p))
    {
      return false;
    }
  }
  return true;
}

/*!
\brief Overloaded.
\param s Stream.
\param icosi The icosidodecahedron.
*/
std::ostream& operator<<(std::ostream& s, const Icosidodecahedron& icosi)
{
  s << "Box(" << icosi.c << ',' << icosi.r << ")";
  return s;
}