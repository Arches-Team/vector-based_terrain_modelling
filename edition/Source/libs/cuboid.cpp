// Cuboid

// Self include
#include "libs/cuboid.h"

#include "libs/quadrangle.h"
#include "libs/mesh.h"

/*!
\class Cuboid cuboid.h
\brief Cuboids.

Cuboids implement convex hexahedra, i.e. polyhedra with six faces.

\ingroup ExtendedKernelGroup
*/

const int Cuboid::edge[24] =
{
  0,1,2,3,4,5,6,7,
  0,2,1,3,4,6,5,7,
  0,4,1,5,2,6,3,7
};

const int Cuboid::quadrangle[24] =
{
  0, 2, 3, 1,
  4, 5, 6, 7,
  0, 4, 6, 2,
  1, 3, 5, 7,
  0, 1, 4, 5,
  3, 2, 7, 6
};

/*!
\brief Create a cuboid.
\param box The box.
*/
Cuboid::Cuboid(const Box& box)
{
  for (int i = 0; i < 8; i++)
  {
    a[i] = box.Vertex(i);
  }
}

/*!
\brief Create a cuboid given vertexes.
\param a000, a100, a010, a110, a001, a101, a011, a111 Vertexes.
*/
Cuboid::Cuboid(const Vector& a000, const Vector& a100, const Vector& a010, const Vector& a110, const Vector& a001, const Vector& a101, const Vector& a011, const Vector& a111)
{
  a[0] = a000;
  a[1] = a100;
  a[2] = a010;
  a[3] = a110;

  a[4] = a001;
  a[5] = a101;
  a[6] = a011;
  a[7] = a111;
}

/*!
\brief Translates a cuboid.

\param t Translation vector.
*/
void Cuboid::Translate(const Vector& t)
{
  for (int i = 0; i < 8; i++)
  {
    a[i] += t;
  }
}

/*!
\brief Rotates a cuboid.

\param r Rotation matrix.
*/
void Cuboid::Rotate(const Matrix& r)
{
  for (int i = 0; i < 8; i++)
  {
    a[i] = r * a[i];
  }
}

/*!
\brief Scales a cuboid.

Note that this function handles negative coefficients in
the scaling vector (by swapping coordinates if need be).
\param s Scaling vector.
*/
void Cuboid::Scale(const Vector& s)
{
  for (int i = 0; i < 8; i++)
  {
    a[i] *= s;
  }
}

/*!
\brief Compute the area.
*/
double Cuboid::Area() const
{
  double area = 0.0;

  // Sum quadrangles
  for (int i = 0; i < 6; i++)
  {
    area += Quadrangle(Vertex(i, 0), Vertex(i, 1), Vertex(i, 2), Vertex(i, 3)).Area();
  }
  return area;
}
