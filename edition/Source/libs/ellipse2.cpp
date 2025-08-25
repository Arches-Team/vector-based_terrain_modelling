// Ellipse

#include <QtWidgets/QGraphicsScene>

#include "libs/ellipse.h"
#include "libs/polygon.h"

/*!
\class Ellipse2 ellipse.h
\brief %Axis aligned ellipses.

\sa Circle2
\ingroup PlanarGroup
*/

/*!
\brief Create an ellipse.
\param c Center.
\param a,b Axes lengths.
\param u Major axis direction.
*/
Ellipse2::Ellipse2(const Vector2& c, const double& a, const double& b, const Vector2& u) :c(c), a(a), b(b), u(u)
{
}

/*!
\brief Create an ellipse.
\param a,b Axes lengths.
*/
Ellipse2::Ellipse2(const double& a, const double& b) :Ellipse2(Vector2::Null, Math::Max(a, b), Math::Min(a, b), a >= b ? Vector2::X : Vector2::Y)
{
}


/*!
\brief Return the half distance between the focus points of the ellipse.

The focus points can be derived as:
\code
Ellipse e;
Vector2 f=e.Center()+e.C()*Vector2::X; // And e.Center()-e.C()*Vector2::X;
\endcode
*/
double Ellipse2::C() const
{
  return sqrt(a * a - b * b);
}

/*!
\brief Return the focus points of the ellipse.

\param i Index of the focus point, true for focus in the direction of the axis, false for the point in the opposite direction.

\sa Ellipse2::C()
*/
Vector2 Ellipse2::Focus(bool i) const
{
  double f = sqrt(a * a - b * b);

  if (i == true)
  {
    return c + f * u;
  }
  else
  {
    return c - f * u;
  }
}

/*!
\brief Return the parameter of the ellipse.
*/
double Ellipse2::P() const
{
  return b * b / a;
}

/*
\brief Length.
*/
double Ellipse2::Length() const
{
  // Equal to 2 sqrt(b*b+c*c) where c is half the distance between the focus points of the ellipse, thus equal to 2 a.
  return 2.0 * a;
}

/*!
\brief Compute the field function value of an ellipse.
\param p Point.
*/
double Ellipse2::Value(const Vector2& p) const
{
  Vector2 q = p - c;
  q = Vector2(u * q, u.Orthogonal() * q);

  double aa = a * a;
  double bb = b * b;

  return sqrt(q[0] * q[0] / aa + q[1] * q[1] / bb) - 1.0;
}

/*!
\brief Compute the gradient of the field function.
\param p Point.
*/
Vector2 Ellipse2::Gradient(const Vector2& p) const
{
  Vector2 q = p - c;
  if (q == Vector2::Null)
  {
    return Vector2::X;
  }
  q = Vector2(u * q, u.Orthogonal() * q);

  double aa = a * a;
  double bb = b * b;

  return (1.0 / sqrt(q[0] * q[0] / bb + q[1] * q[1] / bb)) * Vector2(q[0] / aa, q[1] / bb);
}

/*!
\brief Check if a point is inside the circle.
\param p Point.
*/
bool Ellipse2::Inside(const Vector2& p) const
{
  if (Value(p) < 0.0)
  {
    return true;
  }
  else
  {
    return false;
  }
}


/*!
\brief Draw a circle.
\param scene Graphics scene.
\param pen The pen.
\param brush The brush, should the ellipse be filled.
*/
void Ellipse2::Draw(QGraphicsScene& scene, const QPen& pen, const QBrush& brush) const
{
  QPen thin;
  thin.setWidthF(0.01);
  Circle2(c, 0.01 * (a + b)).Draw(scene, thin);
  Segment2(c, c + u).DrawArrow(scene, 0.25, thin);
  Segment2(c, c + u.Orthogonal()).DrawArrow(scene, 0.25, thin);
  Polygon2(*this).Draw(scene, pen, brush);
  // scene.addEllipse(c[0] - a, c[1] - b, a + a, b + b, pen, brush);

  /*
  QPainterPath ellipse;
  //ellipse.moveTo(c.ToQt());
  ellipse.arcTo(QRectF((c-Vector2(a,b))[0]-a, (c - Vector2(a, b))[1] - b,2.0*a,2.0*b), 0, 360);


   double theta = atan2(u[0],u[1]) * (180 / M_PI);

   QPainter painter(&scene);

  ellipse->save();
  painter->translate(p1.x(), p1.y());
  painter->rotate(theta);
  painter->translate(-p1.x(), -p1.y());
  painter->drawPath(ellipse);
  painter->restore();
  */
}

/*!
\brief Compute the bounding box.
*/
Box2 Ellipse2::GetBox() const
{
  // Axis aligned ellipse
  // return Box2(c - Vector2(a, b), c + Vector2(a, b));

  // General case:  x = +/- sqrt (a^2 cos^2 t + b^2 sin^2 t) and y = +/- sqrt (a^2 sin^2 t + b^2 cos^2 t)
  // If u is normalized, then ux = cos t and uy = sin t and since we do not care about the sign then

  double cu = u[0] * u[0];
  double su = u[1] * u[1];
  double x = sqrt(a * a * cu + b * b * su);
  double y = sqrt(a * a * su + b * b * cu);

  return Box2(c - Vector2(x, y), c + Vector2(x, y));
}

/*!
\brief Generate a random vector inside the circle.
\param random %Random number generator.
*/
Vector2 Ellipse2::RandomInside(Random& random) const
{
  // Radius
  double ru;

  // Squared relative lenths
  double iab = 1.0 / Math::Sqr(Math::Min(a, b) / Math::Max(a, b));

  // Random samples
  double rx, ry;
  do
  {
    rx = 2.0 * random.Uniform() - 1.0;
    ry = 2.0 * random.Uniform() - 1.0;
    ru = rx * rx + iab * ry * ry;
  } while (ru >= 1.0);

  return c + Math::Max(a, b) * Vector2(rx, ry);
}

/*!
\brief Compute a vertex on the ellipse.

\param t Angle.
*/
Vector2 Ellipse2::Vertex(const double& t) const
{
  return c + a * u * cos(t) + b * u.Orthogonal() * sin(t);
}

/*!
\brief Compute the curvature at a vertex on the ellipse.

\param t Angle.
*/
double Ellipse2::Curvature(const double& t) const
{
  double a4 = a * a * a * a;
  double b4 = b * b * b * b;

  // Vertex
  Vector2 p(a * cos(t), b * sin(t));

  // Square components
  p *= p;

  return Math::Sqrt32(a4 * p[1] + b4 * p[0]) / (a4 * b4);
}

/*!
\brief Compute the squared distance between a point and the ellipse.

\param p Point.
*/
double Ellipse2::R(const Vector2& p) const
{
  // Symmetry
  Vector2 qq = p - c;
  qq = Vector2(u * qq, u.Orthogonal() * qq);
  qq = Abs(qq);

  double aa = a;
  double bb = b;

  // Quadrant symmetry
  if (qq[0] > qq[1])
  {
    Math::Swap(qq[0], qq[1]);
    Math::Swap(aa, bb);
  }

  double l = bb * bb - aa * aa;

  double pax = aa * qq[0] / l;
  double pby = bb * qq[1] / l;
  double m2 = pax * pax;
  double n2 = pby * pby;

  double c = (m2 + n2 - 1.0) / 3.0;
  double c3 = c * c * c;

  double q = c3 + m2 * n2 * 2.0;
  double d = c3 + m2 * n2;
  double g = pax + pax * n2;

  double co;

  if (d < 0.0)
  {
    double h = acos(q / c3) / 3.0;
    double s = cos(h);
    double t = sin(h) * Math::Sqrt3;
    double rx = sqrt(-c * (s + t + 2.0) + m2);
    double ry = sqrt(-c * (s - t + 2.0) + m2);
    co = (ry + Math::CopySign(rx, l) + abs(g) / (rx * ry) - pax) / 2.0;
  }
  else
  {
    double h = 2.0 * pax * pby * sqrt(d);
    double s = Math::CopySign(cbrt(fabs(q + h)), q + h);
    double u = Math::CopySign(cbrt(fabs(q - h)), q - h);
    double rx = -s - u - c * 4.0 + 2.0 * m2;
    double ry = (s - u) * Math::Sqrt3;
    double rm = sqrt(rx * rx + ry * ry);
    co = (ry / sqrt(rm - rx) + 2.0 * g / rm - pax) / 2.0;
  }

  double si = sqrt(1.0 - co * co);

  Vector2 pe(aa * co, bb * si);

  return SquaredNorm(qq - pe);
}

/*!
\brief Compute the signed distance between a point and the ellipse.

\param p Point.
*/
double Ellipse2::Signed(const Vector2& p) const
{
  // Symmetry
  Vector2 pp = Abs(p - c);
  double aa = a;
  double bb = b;

  // Quadrant symmetry
  if (pp[0] > pp[1])
  {
    Math::Swap(pp[0], pp[1]);
    aa = b;
    bb = a;
  }

  double l = bb * bb - aa * aa;

  double pax = aa * pp[0] / l;
  double pby = bb * pp[1] / l;
  double m2 = pax * pax;
  double n2 = pby * pby;

  double c = (m2 + n2 - 1.0) / 3.0;
  double c3 = c * c * c;

  double q = c3 + m2 * n2 * 2.0;
  double d = c3 + m2 * n2;
  double g = pax + pax * n2;

  double co;

  if (d < 0.0)
  {
    double h = acos(q / c3) / 3.0;
    double s = cos(h);
    double t = sin(h) * Math::Sqrt3;
    double rx = sqrt(-c * (s + t + 2.0) + m2);
    double ry = sqrt(-c * (s - t + 2.0) + m2);
    co = (ry + Math::CopySign(rx, l) + abs(g) / (rx * ry) - pax) / 2.0;
  }
  else
  {
    double h = 2.0 * pax * pby * sqrt(d);
    double s = Math::CopySign(cbrt(fabs(q + h)), q + h);
    double u = Math::CopySign(cbrt(fabs(q - h)), q - h);
    double rx = -s - u - c * 4.0 + 2.0 * m2;
    double ry = (s - u) * Math::Sqrt3;
    double rm = sqrt(rx * rx + ry * ry);
    co = (ry / sqrt(rm - rx) + 2.0 * g / rm - pax) / 2.0;
  }

  double si = sqrt(1.0 - co * co);

  Vector2 pe(aa * co, bb * si);

  return Math::CopySign(Norm(pp - pe), pp[1] - pe[1]);
}

/*!
\brief Translate an ellipse.
\param t Translation vector.
*/
Ellipse2 Ellipse2::Translated(const Vector2& t) const
{
  return Ellipse2(c + t, a, b, u);
}

/*!
\brief Scales an ellipse.

\param s Scaling factor.
*/
Ellipse2 Ellipse2::Scaled(const double& s) const
{
  return Ellipse2(c * s, a * s, b * s, u);
}

/*!
\brief Compute the matrix form of the ellipse.
*/
Matrix2 Ellipse2::MatrixForm()const
{
  // Cosine and sine are the main axis normalized coordinates
  double c = u[0];
  double s = u[1];

  double A = c * c / (a * a) + s * s / (b * b);
  double B = c * s / (a * a) - s * c / (b * b);
  double C = s * s / (a * a) + c * c / (b * b);;
  return Matrix2(A, B, B, C);
}

/*!
\brief Scales an ellipse.

\sa Ellipse2::Scaled(const Vector2&) const

\param s Scaling vector.
\param local Boolean, set to true if scaling occurs in the local frame of the ellipse (default), otherwise computations are performed in the world coordinate system.
*/
Ellipse2 Ellipse2::Scaled(const Vector2& s,bool local) const
{
  // Local frame: easy case
  if (local == true)
  {
    return Ellipse2(c.Scaled(s), a * s[0], b * s[1], u);
  }
  
  // Apply scaling
  Matrix2 M = MatrixForm().Scaled(Matrix2(s[0] * s[0], s[1] * s[0], s[1] * s[0], s[1] * s[1]));

  // Get parameters from eigen values and vectors.
    double e[2];
    Vector2 v[2];
    M.EigenSolveSymmetric(e,v);
    std::cout << e[0] << std::endl;
    std::cout << e[1] << std::endl;
    std::cout << v[0] << std::endl;
    std::cout << v[1] << std::endl;


    return Ellipse2(c.Scaled(s), 1.0/sqrt(e[0]),1.0/sqrt(e[1]),v[0]);
}

/*!
\brief Rotates an ellipse.
\param a Rotation angle.
*/
Ellipse2 Ellipse2::Rotated(const double& a) const
{
  return Rotated(Matrix2::Rotation(a));
}

/*!
\brief Rotates an ellipse.
\param r Rotation matrix.
*/
Ellipse2 Ellipse2::Rotated(const Matrix2& r) const
{
  /*
  std::cout << u << std::endl;
  std::cout << r  << std::endl;
  std::cout << r * u << std::endl;
  */
  return Ellipse2(r * c, a, b, r * u);
}


/*!
\brief Return a ellipse transformed by a frame.
\param f Transformation.
*/
Ellipse2 Ellipse2::Transformed(const Frame2& f) const
{
  return Ellipse2(f.Transform(c), a, b, f.TransformDirection(u));
}

/*!
\brief Overloaded.
\param s Stream.
\param ellipse The ellipse.
*/
std::ostream& operator<<(std::ostream& s, const Ellipse2& ellipse)
{
  s << "Ellipse2(" << ellipse.c << ',' << ellipse.a << ',' << ellipse.b << ',' << ellipse.u<<')';
  return s;
}

/*!
\brief Perimeter of an ellipse.

The exact perimeter of an ellipse is 4 a E(e) where e is the eccentricity, and the function E is the complete elliptic integral of the second kind.

Here we use an approximation proposed by Srinivasa Ramanujan, <I>Modular Equations and Approximations to Pi</I>.
*/
double Ellipse2::Perimeter() const
{
  return Math::Pi * ((a + b) - sqrt(10.0 * a * b + 3.0 * (a * a + b * b)));
}

/*!
\brief Linear interpolation between two ellipses.

Linear interpolation of centers, radii, and main axis.
*/
Ellipse2 Ellipse2::Lerp(const double& t, const Ellipse2& ea, const Ellipse2& eb)
{
  return Ellipse2((1.0 - t) * ea.c + t * eb.c, (1.0 - t) * ea.a + t * eb.a, (1.0 - t) * ea.b + t * eb.b, (1.0 - t) * ea.u + t * eb.u);
}
