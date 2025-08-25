// Rectangle

#include "libs/rectangle.h"

/*!
\class Rectangles rectangle.h
\brief %Rectangles.
*/

double Rectangles::epsilon = 1.0e-10;

/*!
\brief Create a horizontal rectangle in the xy-plane, centered at the origin.
\param a, b Size.
*/
Rectangles::Rectangles(const double& a, const double& b) :a(a), b(b), c(Vector::Null), x(Vector::X), y(Vector::Y), z(Vector::Z)
{
}

/*!
\brief Create a rectangle given its origin and two orthogonal axes.
\param c Origin
\param x, y Axes.
*/
Rectangles::Rectangles(const Vector& c, const Vector& x, const Vector& y) :Rectangles(c, Normalized(x), Normalized(y), Norm(x), Norm(y))
{
}

/*!
\brief Create a rectangle given its origin and two orthogonal axes.
\param c Origin
\param x, y Unit axes (normalized).
\param a,b %Size.
*/
Rectangles::Rectangles(const Vector& c, const Vector& x, const Vector& y, const double& a, const double& b) :c(c), x(x), y(y), a(a), b(b)
{
  z = Rectangles::x / Rectangles::y;
}

/*!
\brief Compute the squared distance between a point and the rectangle.
\param p Point.
*/
double Rectangles::R(const Vector& p) const
{
  //return SquaredNorm(Normal(p));

  Vector q = p - c;

  Vector n(q * x, q * y, q * z);

  // Distance
  double d = Math::Sqr(n[2]);

  // Symmetry
  Vector2 nxy = Abs(Vector2(n)) - Vector2(a, b);

  if (nxy[0] > 0.0)
  {
    d += Math::Sqr(nxy[0]);
  }

  if (nxy[1] > 0.0)
  {
    d += Math::Sqr(nxy[1]);
  }

  return d;
}

/*!
\brief Computes the vector distance between the box and a point.
\param p Point.
*/
Vector Rectangles::Normal(const Vector& p) const
{
  Vector q = p - c;

  Vector n(q * x, q * y, q * z);

  if (n[0] < -a)
  {
    n[0] += a;
  }
  else if (n[0] > a)
  {
    n[0] -= a;
  }
  else
  {
    n[0] = 0.0;
  }

  if (n[1] < -b)
  {
    n[1] += b;
  }
  else if (n[1] > b)
  {
    n[1] -= b;
  }
  else
  {
    n[1] = 0.0;
  }

  return n[0] * x + n[1] * y + n[2] * z;
}

/*!
\brief Rotate the rectangle.
\param r Rotation matrix.
*/
void Rectangles::Rotate(const Matrix& r)
{
  c = r * c;
  x = r * x;
  y = r * y;
  z = r * z;
}

/*!
\brief Translate the rectangle.
\param t Translation vector.
*/
void Rectangles::Translate(const Vector& t)
{
  c += t;
}

/*!
\brief Scale the rectangle.
\param s Scaling factor.
*/
void Rectangles::Scale(const double& s)
{
  a *= s;
  b *= s;
}

/*!
\brief Compute the box embedding the shape.
*/
Box Rectangles::GetBox() const
{
  Vector p[4];
  p[0] = c + x * a + y * b;
  p[1] = c - x * a + y * b;
  p[2] = c + x * a - y * b;
  p[3] = c - x * a - y * b;

  return Box(Vector::Min(Vector::Min(p[0], p[1]), Vector::Min(p[2], p[3])), Vector::Max(Vector::Max(p[0], p[1]), Vector::Max(p[2], p[3])));
}

/*!
\brief Overloaded.
\param s Stream.
\param rectangle The rectangle.
*/
std::ostream& operator<<(std::ostream& s, const Rectangles& rectangle)
{
  s << "Rectangle(" << rectangle.c << ',' << rectangle.a * rectangle.x << ',' << rectangle.b * rectangle.y << ')';
  return s;

}

/*!
\brief Compute the coordinates of a vertex in the rectangle.
\param u,v Coordinates of the point.
*/
Vector Rectangles::Vertex(const double& u, const double& v) const
{
  return c + (2.0 * u - 1.0) * a * x + (2.0 * v - 1.0) * b * y;
}
