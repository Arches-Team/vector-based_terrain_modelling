// Frames 

#include "libs/frame.h"

/*!
\class Frame frame.h
\brief Solid transformations.

The class groups the rotation matrix and the translation vector.

The orientation of the frame should be defined by the appropriate rotation matrix.

The inverse transformation is implicitly characterized as the transposed matrix and
the opposite translation vector.

Conversion constructors allow the following:
\code
  Frame f = Frame2::Rotation(Math::Pi/4.0); // Rotation around z-axis
\endcode
\ingroup KernelGroup
*/

const Frame Frame::Id(Matrix::Identity, Vector::Null);

/*!
\brief Creates a frame given a rotation matrix and a translation vector.
\param r Frame matrix.
\param t Translation vector.
*/
Frame::Frame(const Matrix& r, const Vector& t) :r(r), t(t)
{
}

/*!
\brief Creates a frame given the origin and its orthogonal unit vectors.
\param c Origin.
\param x, y, z Orthonormal vectors.
*/
Frame::Frame(const Vector& c, const Vector& x, const Vector& y, const Vector& z) :r(Matrix(x, y, z)), t(c)
{
}

/*!
\brief Creates a frame from another planar frame.
\param frame %Frame
*/
Frame::Frame(const Frame2& frame):r(Matrix(frame.R())),t(frame.T().ToVector(0.0))
{
}

/*!
\brief Creates a translation transformation.

%Matrix is identity whereas the translation vector is
computed from the argument vector.
\param t Translation vector.
*/
Frame Frame::Translation(const Vector& t)
{
  return Frame(Matrix::Identity, t);
}

/*!
\brief Creates a rotation frame.
\param a A vector of angles (should be in radians) that defines the rotation x-y-z axes.
*/
Frame Frame::Rotation(const Vector& a)
{
  return Frame(Matrix::Rotation(a), Vector::Null);
}

/*!
\brief Create a rotation frame about an arbitrary axis.

\param axis %Axis.
\param angle Angle (should be in radian).
*/
Frame Frame::Rotation(const Vector& axis, const double& angle)
{
  return Frame(Matrix::Rotation(axis, angle), Vector::Null);
}

/*!
\brief Create a frame that rotates a normalized vector into another one.
\param a, b Initial and final vectors (which should be normalized).
*/
Frame Frame::Rotation(const Vector& a, const Vector& b)
{
  return Frame(Matrix::Rotation(a, b), Vector(0.0));
}

/*!
\brief Given a point and a direction, compute the frame that brings these into a canonical coordinate system.

This is useful for computing the canonical coordinate system
of revolution objects such as cylinders or cones.

\param p Origin of the frame.
\param axis Vertical z-axis of the frame.
*/
Frame Frame::Canonical(const Vector& p, const Vector& axis)
{
  Vector a;
  if (fabs(axis[2]) == 1.0)
  {
    a = Vector(1.0, 0.0, 0.0);
  }
  else
  {
    a = Vector(-axis[1], axis[0], 0.0);
  }

  return Frame(Matrix::Rotation(a, acos(axis[2])), p);
}

/*!
\brief Compute a frame given an origin and direction vector.

\sa Vector::Orthonormal(Vector&,Vector&)
\param c Origin.
\param n Direction.
*/
Frame Frame::Orthonormal(const Vector& c, const Vector& n)
{
  Vector x, y;
  Vector z = Normalized(n);
  z.Orthonormal(x, y);
  return Frame(c, x, y, z);
}

/*!
\brief Overloaded.
\param s Stream.
\param frame The frame.
*/
std::ostream& operator<<(std::ostream& s, const Frame& frame)
{
  s << "Frame(" << frame.r << ',' << frame.t << ')';
  return s;
}

/*!
\brief Compute the coordinates of a point on a circle inside the frame;
\param theta Euler angle.
\param i,j Integers representing which column vector of the matrix will be used as basis.
*/
Vector Frame::CircleVertex(const double& theta, int i, int j) const
{
  return t + cos(theta) * r.C(i) + sin(theta) * r.C(j);
}

/*!
\brief Compute the coordinates of the normal a point on a circle inside the frame;
\param theta Euler angle.
\param i,j Integers representing which column vector of the matrix will be used as basis.
*/
Vector Frame::CircleNormal(const double& theta, int i, int j) const
{
  return cos(theta) * r.C(i) + sin(theta) * r.C(j);
}

/*!
\brief Compute the coordinates of a point on a sphere inside the frame.
\param r Radius.
\param theta, phi Polar coordinates.
\param i,j,k Integers representing which column vector of the matrix will be used as basis.
*/
Vector Frame::SphereVertex(const double& r, const double& theta, const double& phi, int i, int j, int k) const
{
  return t + r * (cos(theta) * cos(phi) * Frame::r.C(i) + sin(theta) * cos(phi) * Frame::r.C(j) + sin(phi) * Frame::r.C(k));
}

/*!
\brief Compute the coordinates of the normal of a point on a sphere inside the frame.
\param theta, phi Polar coordinates.
\param i,j,k Integers representing which column vector of the matrix will be used as basis.
*/
Vector Frame::SphereNormal(const double& theta, const double& phi, int i, int j, int k) const
{
  return cos(theta) * cos(phi) * r.C(i) + sin(theta) * cos(phi) * r.C(j) + sin(phi) * r.C(k);
}

/*!
\brief Compute the inverse transformation.
*/
Frame Frame::Inverse() const
{
  return Frame(r.T(), -(r.T() * t));
}

/*!
\brief Compose the frame with another one.

\param frame The frame.
*/
void Frame::Compose(const Frame& frame)
{
  t = frame.r * t + frame.t;
  r = frame.r * r;
}

/*!
\brief Compose the frame with another one.

\param frame The frame.
*/
Frame Frame::Composed(const Frame& frame) const
{
  return Frame(frame.r * r, frame.r * t + frame.t);
}

/*!
\brief Transform a ray out of the frame coordinate system.
\param ray The ray.
*/
Ray Frame::Transform(const Ray& ray) const
{
  return Ray(r * ray.Origin() + t, r * ray.Direction());
}

/*!
\brief Transform a ray into the frame coordinate system.
\param ray The ray.
*/
Ray Frame::InverseTransform(const Ray& ray) const
{
  return Ray(r.T() * (ray.Origin() - t), r.T() * ray.Direction());
}

