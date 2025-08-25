// Quadrangle

#include <QtWidgets/QGraphicsScene>

#include "libs/quadrangle.h"

/*!
\class Quadrangle2 quadrangle.h
\brief Convex quadrangles in the plane.

Note that the quadrangle should be convex, otherwise some functions will not work properly.

Some functions such as the squared distance R() and normal Normal() are computationally intensive.

\image html quadrangle.png

\ingroup PlanarGroup
*/

double Quadrangle2::epsilon = 1.0e-8;

/*!
\brief Create a quadrangle.
\param a, b, c, d Vertices (given in trigonometric order).
*/
Quadrangle2::Quadrangle2(const Vector2& a, const Vector2& b, const Vector2& c, const Vector2& d) :p{ a,b,c,d }
{
}

/*!
\brief Create a planar quadrangle from a quadrangle.
\param q %Quadrangle.
*/
Quadrangle2::Quadrangle2(const Quadrangle& q)
{
  p[0] = q.Vertex(0);
  p[1] = q.Vertex(1);
  p[2] = q.Vertex(2);
  p[3] = q.Vertex(3);
}

/*!
\brief Create a horizontal quadrangle.
\param r Radius, i.e., half size.
\sa Quadrangle::Quadrangle(double, double)
*/
Quadrangle2::Quadrangle2(const double& r) :Quadrangle2(Vector2(-r, -r), Vector2(-r, r), Vector2(r, r), Vector2(r, -r))
{
}

/*!
\brief Compute the position of a vertex inside the quadrangle.
\param u,v Parametric coordinates of the vertex.
*/
Vector2 Quadrangle2::Vertex(const double& u, const double& v) const
{
  return Vector2::Bilinear(p[0], p[1], p[2], p[3], u, v);
}

/*!
\brief Translate the quadrangle.
\param t Translation vector.
*/
void Quadrangle2::Translate(const Vector2& t)
{
  for (int i = 0; i < 4; i++)
  {
    p[i] += t;
  }
}

/*!
\brief Scale the quadrangle.
\param s Scaling factor.
*/
void Quadrangle2::Scale(const double& s)
{
  for (int i = 0; i < 4; i++)
  {
    p[i] *= s;
  }
}

/*!
\brief Rotate the quadrangle.
\param r Rotation matrix.
*/
void Quadrangle2::Rotate(const Matrix2& r)
{
  for (int i = 0; i < 4; i++)
  {
    p[i] = r * p[i];
  }
}

/*!
\brief Transform a quadrangle.
\param f Transform.
*/
void Quadrangle2::Transform(const Frame2& f)
{
  for (int i = 0; i < 4; i++)
  {
    p[i] = f.Transform(p[i]);
  }
}

/*!
\brief Compute the box embedding the shape.
*/
Box2 Quadrangle2::GetBox() const
{
  return Box2(Vector2::Min(Vector2::Min(p[0], p[1]), Vector2::Min(p[2], p[3])), Vector2::Max(Vector2::Max(p[0], p[1]), Vector2::Max(p[2], p[3])));
}

/*!
\brief Check if a point is inside the quadrangle.
\param q Point.
*/
bool Quadrangle2::Inside(const Vector2& q) const
{
  for (int i = 0; i < 4; i++)
  {
    if (IsRight(q, p[i], p[(i + 1) % 4]))
    {
      return false;
    }
  }
  return true;
}

/*!
\brief Compute the squared distance between a point and the quadrangle.
\param q Point.
\sa Normal()
*/
double Quadrangle2::R(const Vector2& q) const
{
  Vector2 n = Normal(q);
  return n * n;
}

/*!
\brief Computes the normal vector between a point and a quadrangle.

Let <b>q</b> the projection of <b>p</b> onto the box, the normal vector
is defined as <b>n</b>=<b>p</b>-<b>q</b>. The returned vector is null
if the point is inside the quadrangle.

\param q Point.
*/
Vector2 Quadrangle2::Normal(const Vector2& q) const
{
  // Edges and edge normals
  Vector2 e[4];
  Vector2 ne[4];

  // Point to vertices
  Vector2 pq[4];

  // Pre-compute some edge and vertex data
  for (int i = 0; i < 4; i++)
  {
    e[i] = Normalized(p[(i + 1) % 4] - p[i]);
    ne[i] = e[i].Orthogonal();
    pq[i] = q - p[i];
  }

  for (int i = 0; i < 4; i++)
  {
    double neq = ne[i] * pq[i];
    // On one side of the quadrangle : may be we can compute the normal
    if (neq > 0.0)
    {
      // Check vertex i 
      if (pq[i] * e[i] < 0.0)
      {
        // Inside sector : return normal as distance to vertex i
        if (pq[i] * e[(i + 3) % 4] > 0.0)
        {
          return pq[i];
        }
      }
      else
      {
        // Check vertex i+1 
        if (pq[(i + 1) % 4] * e[i] > 0.0)
        {
          // Inside sector : return normal as distance to vertex i+1
          if (pq[(i + 1) % 4] * e[(i + 1) % 4] < 0.0)
          {
            return pq[(i + 1) % 4];
          }
        }
        // Edge : return normal
        else
        {
          return neq * ne[i];
        }
      }
    }
  }
  return Vector::Null;
}

/*!
\brief Overloaded.
\param s Stream.
\param q The Quadrangle.
*/
std::ostream& operator<<(std::ostream& s, const Quadrangle2& q)
{
  s << "Quadrangle2(" << q.p[0] << ',' << q.p[1] << ',' << q.p[2] << ',' << q.p[3] << ')';
  return s;
}

/*!
\brief Compute the inverse bi-linear coordinates of a point in the quadrangle.
\param q Point.
\param u,v Bi-linear coordinates.
\return Boolean that checks whether a solution has been found.
*/
bool Quadrangle2::InverseBilinear(const Vector2& q, double& u, double& v) const
{
  Vector2 e = p[1] - p[0];
  Vector2 f = p[3] - p[0];
  Vector2 g = p[0] - p[1] + p[2] - p[3];
  Vector2 h = q - p[0];

  double k2 = g / f;
  double k1 = e / f + h / g;
  double k0 = h / e;

  if (Math::Abs(k2) < 0.00000001)
  {
    // Linear Case 
    v = -k0 / k1;
  }
  else
  {
    double delta = k1 * k1 - 4.0 * k0 * k2;
    // Non Linear Case 
    if (delta >= 0.0)
    {
      v = (-k1 + sqrt(delta)) / (2.0 * k2);
      // Test that v is between 0 and 1 (this applies only if we query p inside the quad which is true when delta>=0)
      /*
      if (v < 0.0 || v>1.0)
        v = (-k1 - sqrt(k1 * k1 - 4 * k0 * k2)) / (2 * k2);
        */
    }
    else
    {
      // We could handle the else case (p is not inside the quad)
      //U = (H[0] - F[0]*V)/(E[0] + G[0]*V);
      return false;
    }
  }
  double div = (e[0] + g[0] * v);
  if (Math::Abs(div) < Quadrangle2::epsilon)
    u = (h[1] - f[1] * v) / (e[1] + g[1] * v);
  else
    u = (h[0] - f[0] * v) / div;
  return true;
}

/*!
\brief Draw the quadrangle.
\param scene Graphics scene.
\param pen The pen.
\param brush The brush, should the triangle be filled.
*/
void Quadrangle2::Draw(QGraphicsScene& scene, const QPen& pen, const QBrush& brush) const
{
  QPolygonF polygon;

  for (int i = 0; i < 4; i++)
  {
    polygon << QPointF(p[i][0], p[i][1]);
  }

  scene.addPolygon(polygon, pen, brush);
}

/*!
\brief Compute the perimeter of the quadrangle.
*/
double Quadrangle2::Perimeter() const
{
  return Norm(p[1] - p[0]) + Norm(p[2] - p[1]) + Norm(p[3] - p[2]) + Norm(p[0] - p[3]);
}

/*!
\brief Compute the area of the quadrangle.

The area is defined by Bretschneider's formula.
The trigonometric adjustment in Bretschneider's formula for non-cyclicality
of the quadrilateral can be rewritten non-trigonometrically in terms of the
sides and the diagonals:

A<SUP>2</SUP>=(<I>s-a</I>)(<I>s-b</I>)(<I>s-c</I>)(<I>s-d</I>)-1/4(<I>ac</I>+<I>bd</I>+<I>pq</I>)(<I>ac</I>+<I>bd</I>-<I>pq</I>).
*/
double Quadrangle2::Area() const
{
  double a = Norm(p[1] - p[0]);
  double b = Norm(p[2] - p[1]);
  double c = Norm(p[3] - p[2]);
  double d = Norm(p[0] - p[3]);

  // Half perimeter
  double s = 0.5 * (a + b + c + d);

  // Diagonals
  double pq = Norm(p[2] - p[0]) * Norm(p[3] - p[1]);

  // Squared area
  double area = (s - a) * (s - b) * (s - c) * (s - d) - 0.25 * (a * c + b * d + pq) * (a * c + b * d - pq);

  return sqrt(area);
}
