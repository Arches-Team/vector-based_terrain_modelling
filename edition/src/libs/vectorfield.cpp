// Vector fields

#include "libs/vectorfield.h"
#include "libs/curvepoint.h"

/*!
\class VectorField vectorfield.h
\brief A base three-dimensional field of Vector values.

\ingroup StructureGroup
*/

/*!
\brief Create the field structure.
\param box The box.
\param x,y,z Size of the array.
\param v Default value of field.
*/
VectorField::VectorField(const Box& box, int x, int y, int z, const Vector& v) :Array(box, x, y, z)
{
  field.fill(v, nx * ny * nz);
}

/*!
\brief Create the field structure from an analytic vector field.
\param box The box.
\param x,y,z Size of the array.
\param a The analytic vector field.
*/
VectorField::VectorField(const Box& box, int x, int y, int z, const AnalyticVectorField& a) :Array(box, x, y, z)
{
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      for (int k = 0; k < nz; k++)
      {
        field[VertexIndex(i, j, k)] = a.Value(Vertex(i, j, k));
      }
    }
  }
}

/*!
\brief Return the field vector at a given array vertex.
\param i,j,k Integer coordinates of the vertex.
\sa at(int,int,int)
*/
Vector VectorField::Value(int i, int j, int k) const
{
  return field.at(VertexIndex(i, j, k));
}

/*!
\brief Get the field value at a given point.

This function relies on a bi-linear interpolation of the vectors.

\sa Math::Bilinear
\param p Point.
*/
Vector VectorField::Value(const Vector& p) const
{
  double u, v, w;
  int i, j, k;
  CellInteger(p, i, j, k, u, v, w);

  // Test domain
  if (!InsideCellIndex(i, j, k))
    return Vector::Null;

  return Vector::Trilinear(at(i, j, k), at(i + 1, j, k), at(i + 1, j + 1, k), at(i, j + 1, k), at(i, j, k + 1), at(i + 1, j, k + 1), at(i + 1, j + 1, k + 1), at(i, j + 1, k + 1), u, v, w);
}

/*!
\brief Set the whole vector field to a constant vector.
\param v %Vector.
*/
void VectorField::Set(const Vector& v)
{
  field.fill(v, nx * ny * nz);
}

/*!
\brief Compute the new position of a point in the vector field using forward Euler integration.

Direct forward Euler integration simply returns <B>p</B> + t v(<B>p</B>), error is bounded by O(t<SUP>2</SUP>).

Improved Eulerâ€™s method use a prediction-correction scheme, error is bounded by O(t<SUP>3</SUP>).
\param p Point.
\param t Time step.
\param s Scheme: false is simple Euler method, true for prediction-correction.
*/
Vector VectorField::Euler(const Vector& p, const double& t, bool s) const
{
  if (s == false)
  {
    return p + t * Value(p);
  }
  else
  {
    Vector vp = Value(p);
    Vector q = p + 0.5 * t * vp;
    return p + 0.5 * t * (vp + Value(q));
  }
}

/*!
\brief Compute the new position of a point in the vector field using Runge Kutta 4 integration.

Error is bounded by O(t<SUP>4</SUP>).

\param p Point.
\param t Time step.
*/
Vector VectorField::RungeKutta(const Vector& p, const double& t) const
{
  Vector k1 = t * Value(p);
  Vector v1 = p + 0.5 * k1;

  Vector k2 = t * Value(v1);
  Vector v2 = p + 0.5 * k2;

  Vector k3 = t * Value(v2);
  Vector v3 = p + k3;

  Vector k4 = t * Value(v3);

  return p + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
}


/*!
\brief Compute trajectory of a point using Euler integration.

\param p Point.
\param t Time step.
\param n Number of steps.
\param s Scheme: false is simple Euler method, true for prediction-correction.
*/
PointCurve VectorField::EulerSteps(const Vector& p, const double& t, int n, bool s) const
{
  PointCurve c;

  Vector q = p;

  c.Append(q);
  for (int i = 1; i < n; i++)
  {
    q = Euler(q, t, s);
    if (!Inside(q))
    {
      break;
    }
    c.Append(q);
  }
  return c;
}

/*!
\class AnalyticVectorField vectorfield.h
\brief A core analytic three-dimensional vector field.

\ingroup StructureGroup
*/

const double AnalyticVectorField::epsilon = 1e-6;

/*!
\brief Compute the value of the field.
\param p Point.
*/
#pragma warning(push)
#pragma warning(disable: 4100)  
Vector AnalyticVectorField::Value(const Vector& p) const
{
  return Vector(0.0);
}
#pragma warning(pop)

/*!
\brief Compute the divergence of the field.
\param p Point.
*/
double AnalyticVectorField::Divergence(const Vector& p) const
{
  const double x = Value(Vector(p[0] + epsilon, p[1], p[2]))[0] - Value(Vector(p[0] - epsilon, p[1], p[2]))[0];
  const double y = Value(Vector(p[0], p[1] + epsilon, p[2]))[1] - Value(Vector(p[0], p[1] - epsilon, p[2]))[1];
  const double z = Value(Vector(p[0], p[1], p[2] + epsilon))[2] - Value(Vector(p[0], p[1], p[2] - epsilon))[2];

  return (x + y + z) * (0.5 / epsilon);
}

/*!
\brief Compute the curl of the field.
\param p Point.
*/
Vector AnalyticVectorField::Curl(const Vector& p) const
{
  const Vector dx = Value(Vector(p[0] + epsilon, p[1], p[2])) - Value(Vector(p[0] - epsilon, p[1], p[2]));
  const Vector dy = Value(Vector(p[0], p[1] + epsilon, p[2])) - Value(Vector(p[0], p[1] - epsilon, p[2]));
  const Vector dz = Value(Vector(p[0], p[1], p[2] + epsilon)) - Value(Vector(p[0], p[1], p[2] - epsilon));

  const double x = dy[2] - dz[1];
  const double y = dz[0] - dx[2];
  const double z = dx[1] - dy[0];

  return Vector(x, y, z) * (0.5 / epsilon);
}

/*!
\class AnalyticVectorField2 vectorfield.h
\brief A core analytic two-dimensional vector field.

\ingroup StructureGroup
*/

const double AnalyticVectorField2::epsilon = 1e-6;

/*!
\brief Compute the value of the field.
\param p Point.
*/
#pragma warning(push)
#pragma warning(disable: 4100)  
Vector2 AnalyticVectorField2::Value(const Vector2& p) const
{
  return Vector2(0.0);
}
#pragma warning(pop)

/*!
\brief Compute the divergence of the field.
\param p Point.
*/
double AnalyticVectorField2::Divergence(const Vector2& p) const
{
  const double x = Value(Vector2(p[0] + epsilon, p[1]))[0] - Value(Vector2(p[0] - epsilon, p[1]))[0];
  const double y = Value(Vector2(p[0], p[1] + epsilon))[1] - Value(Vector2(p[0], p[1] - epsilon))[1];

  return (x + y) * (0.5 / epsilon);
}

/*!
\brief Compute the curl of the field.
\param p Point.
*/
double AnalyticVectorField2::Curl(const Vector2& p) const
{
  const Vector2 dx = Value(Vector2(p[0] + epsilon, p[1])) - Value(Vector2(p[0] - epsilon, p[1]));
  const Vector2 dy = Value(Vector2(p[0], p[1] + epsilon)) - Value(Vector2(p[0], p[1] - epsilon));

  return dx[1] - dy[0];
}

