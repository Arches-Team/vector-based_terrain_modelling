// Quad

#include "libs/quadrangle.h"

/*!
\class Quadrangle quadrangle.h
\brief Quadrangles.

Quadrangles should be planar and convex.

\ingroup ExtendedKernelGroup
*/

double Quadrangle::epsilon = 1.0e-10;

/*!
\brief Create a quadrangle.
\param a, b, c, d Vertices.
*/
Quadrangle::Quadrangle(const Vector& a, const Vector& b, const Vector& c, const Vector& d)
{
  p[0] = a;
  p[1] = b;
  p[2] = c;
  p[3] = d;
}

/*!
\brief Create a horizontal quadrangle.
\param x,y Half size, the quadrangle vertexes will have (&plusmn;x &plusmn;y) coordinates.
*/
Quadrangle::Quadrangle(double x, double y) :Quadrangle(Vector(-x, -y, 0.0), Vector(-x, y, 0.0), Vector(x, y, 0.0), Vector(x, -y, 0.0))
{
}

/*!
\brief Create a horizontal quadrangle.
\param r Radius, i.e., half size.
\sa Quadrangle::Quadrangle(double, double)
*/
Quadrangle::Quadrangle(double r) :Quadrangle(Vector(-r, -r, 0.0), Vector(-r, r, 0.0), Vector(r, r, 0.0), Vector(r, -r, 0.0))
{
}

/*!
\brief Compute the normal of the quadrangle.

The quadrangle should be flat.
*/
Vector Quadrangle::Normal() const
{
  return Normalized((p[1] - p[0]) / (p[2] - p[0]));
}

/*!
\brief Compute the position of a vertex inside the quadrangle.
\param u,v Parametric coordinates of the vertex.
*/
Vector Quadrangle::Vertex(const double& u, const double& v) const
{
  return Vector::Bilinear(p[0], p[1], p[2], p[3], u, v);
}

/*!
\brief Compute the normal of a vertex inside the quadrangle.
\param u,v Parametric coordinates of the vertex.
*/
Vector Quadrangle::Normal(const double& u, const double& v) const
{
  Vector tu = -(1 - v) * p[0] - (v)*p[3] + (1 - v) * p[1] + (v)*p[2];
  Vector tv = -(1 - u) * p[0] + (1 - u) * p[3] - (u)*p[1] + (u)*p[2];

  return Normalized(tu / tv);
}

/*!
\brief Rotate the quadrangle.
\param r Rotation matrix.
*/
void Quadrangle::Rotate(const Matrix& r)
{
  for (int i = 0; i < 4; i++)
  {
    p[i] = r * p[i];
  }
}

/*!
\brief Translate the quadrangle.
\param t Translation vector.
*/
void Quadrangle::Translate(const Vector& t)
{
  for (int i = 0; i < 4; i++)
  {
    p[i] += t;
  }
}

/*!
\brief Return a translated quadrangle.
\param t Translation vector.
*/
Quadrangle Quadrangle::Translated(const Vector& t) const
{
  return Quadrangle(p[0] + t, p[1] + t, p[2] + t, p[3] + t);
}

/*!
\brief Scale the quadrangle.
\param s Scaling factor.
*/
void Quadrangle::Scale(const double& s)
{
  for (int i = 0; i < 4; i++)
  {
    p[i] *= s;
  }
}

/*!
\brief Transform the quadrangle.
\param frame Transformation.
*/
void Quadrangle::Transform(const FrameScaled& frame)
{
  for (int i = 0; i < 4; i++)
  {
    p[i] = frame.Transform(p[i]);
  }
}

/*!
\brief Transform the quadrangle.
\param frame Transformation.
*/
Quadrangle Quadrangle::Transformed(const FrameScaled& frame) const
{
  return Quadrangle(frame.Transform(p[0]), frame.Transform(p[1]), frame.Transform(p[2]), frame.Transform(p[3]));
}

/*
\brief Compute the box embedding the shape.
*/
Box Quadrangle::GetBox() const
{
  return Box(Vector::Min(Vector::Min(p[0], p[1]), Vector::Min(p[2], p[3])), Vector::Max(Vector::Max(p[0], p[1]), Vector::Max(p[2], p[3])));
}

/*
\brief Overloaded.
\param s Stream.
\param q The Quadrangle.
*/
std::ostream& operator<<(std::ostream& s, const Quadrangle& q)
{
  s << "Quadrangle(" << q.p[0] << ',' << q.p[1] << ',' << q.p[2] << ',' << q.p[3] << ')';
  return s;
}

/*
\brief If quadrangle is planar, compute its area.
*/
double Quadrangle::Area() const
{
  return 0.5 * (Norm((p[1] - p[0]) / (p[2] - p[0])) + (Norm((p[2] - p[0]) / (p[3] - p[0]))));
}

