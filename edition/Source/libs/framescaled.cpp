// Frames 

#include "libs/frame.h"

/*!
\class FrameScaled frame.h
\brief A transformation defined as the composition of a rotation, a translation and a scale.

The data-structure groups the rotation matrix, the translation and scaling vectors.

The orientation of the frame should be defined by the appropriate rotation matrix.

In general, a FrameScaled is defined by its rotation matrix, and a translation and scaling vector:
\code
FrameScaled frame(Matrix::Rotation(Vector( Math::Pi/4.0,0.0, 0.0)),Vector(0.0),Vector(1.0,1.0,2.0));
\endcode

When there is no rotation, the rotation matrix can be defined using the identity matrix:
\code
FrameScaled frame(Matrix::Identity,Vector(-2.0,4.0,0.0),Vector(2.0,1.0,1.0));
\endcode
It is also possible and simpler to use:
\code
FrameScaled frame(Vector(-2.0,4.0,0.0),Vector(2.0,1.0,1.0));
\endcode

\ingroup ExtendedKernelGroup
*/

const FrameScaled FrameScaled::Id(Matrix(Vector(1.0)));

/*!
\brief Creates a frame.
\param r Rotation matrix.
\param t Translation vector.
\param s Scaling vector.
*/
FrameScaled::FrameScaled(const Matrix& r, const Vector& t, const Vector& s) :Frame(r, t), s(s)
{
}

/*!
\brief Creates a frame given the origin of the frame and orthogonal unit vectors.
\param c Origin.
\param x, y, z Orthogonal unit vectors.
*/
FrameScaled::FrameScaled(const Vector& c, const Vector& x, const Vector& y, const Vector& z) :Frame(c, x, y, z), s(Vector(1.0))
{
}

/*!
\brief Creates a translation frame.
\param t Translation.
*/
FrameScaled::FrameScaled(const Vector& t) :FrameScaled(Matrix::Identity, t, Vector(1.0))
{
}

/*!
\brief Creates a frame given the translation and scaling vectors.
\param t Translation.
\param s Scaling vector.
*/
FrameScaled::FrameScaled(const Vector& t, const Vector& s) :FrameScaled(Matrix::Identity, t, s)
{
}

/*!
\brief Creates a frame given the translation and scale.
\param t Translation.
\param s Uniform scaling factor.
*/
FrameScaled::FrameScaled(const Vector& t, const double& s) :FrameScaled(t, Vector(s))
{
}

/*!
\brief Creates a frame.
\param frame The frame.
\param s Scaling.
*/
FrameScaled::FrameScaled(const Frame& frame, const Vector& s) :Frame(frame), s(s)
{
}

/*!
\brief Creates a frame.
\param frame The frame.
\param s Scaling.
*/
FrameScaled::FrameScaled(const Frame2& frame, const Vector& s) :Frame(frame), s(s)
{
}

/*!
\brief Compute and return the inverse transformation.
*/
FrameScaled FrameScaled::Inverse() const
{
  return FrameScaled(Frame::Inverse(), s.Inverse());
}

/*!
\brief Rotate the shape.
\param r Rotation vector representing the rotation angles around x, y and z axes.
*/
void FrameScaled::Rotate(const Vector& r)
{
  Matrix a = Matrix::Rotation(r);

  Frame::t = a * Frame::t;
  Frame::r = a * Frame::r;
}

/*!
\brief Rotate the shape in object coordinates (does not modify the translation vector).
\param r Rotation vector representing the rotation angles around x, y and z axes.
*/
void FrameScaled::ObjectRotate(const Vector& r)
{
  Matrix a = Matrix::Rotation(r);

  Frame::r = Frame::r * a;
}

/*!
\brief Translate the shape.
\param t Translation vector.
*/
void FrameScaled::Translate(const Vector& t)
{
  Frame::t += t;
}

/*!
\brief Scale the shape.
\param u Scaling vector.
*/
void FrameScaled::Scale(const Vector& u)
{
  FrameScaled::s *= u;
}

/*!
\brief Scale the shape.
\param u Scaling factor.
*/
void FrameScaled::Scale(const double& u)
{
  FrameScaled::s *= u;
}

/*!
\brief Get scaling vector from transformation.
*/
Vector FrameScaled::S() const
{
  return s;
}

/*!
\brief Transforms a given input point.
\param p Point.
*/
Vector FrameScaled::Transform(const Vector& p) const
{
  return Frame::Transform(p.Scaled(s));
}

/*!
\brief Transforms a given input normal.
\param n Normal vector.
*/
Vector FrameScaled::TransformDirection(const Vector& n) const
{
  // Apply rotation
  return Frame::TransformDirection(n.Scaled(s.Inverse()));
}

/*!
\brief Inverse transformation of a given input point.
\param p Point.
*/
Vector FrameScaled::InverseTransform(const Vector& p) const
{
  return Frame::InverseTransform(p).Scaled(s.Inverse());
}

/*!
\brief Transforms a given input normal.
\param n Normal vector.
*/
Vector FrameScaled::InverseTransformDirection(const Vector& n) const
{
  // Apply rotation
  return Frame::InverseTransformDirection(n.Scaled(s.Inverse()));
}

/*!
\brief Linear interpolation of two frames.
\param a, b The two frames.
\param alpha Interpolation parameter.
*/
void FrameScaled::Lerp(const FrameScaled& a, const FrameScaled& b, const double& alpha)
{
  Vector ar = a.r.GetRotationAngles();
  Vector br = b.r.GetRotationAngles();
  Vector cr = (1.0 - alpha) * ar + alpha * br;

  r = Matrix::Rotation(cr);
  //r = Matrix::Lerp(alpha, a.r, b.r);
  t = (1.0 - alpha) * a.t + alpha * b.t;
  s = (1.0 - alpha) * a.s + alpha * b.s;
}

/*!
\brief Compose two frames.
\param frame The frame.
*/
FrameScaled FrameScaled::Composed(const FrameScaled& frame) const
{
  return FrameScaled(frame.r * r, frame.r * Matrix(frame.s) * t + frame.t, s.Scaled(frame.s));
}

/*!
\brief Compose two frames.
\param frame The frame.
*/
void FrameScaled::Compose(const FrameScaled& frame)
{
  t = frame.r * Matrix(frame.s) * t + frame.t;
  r = frame.r * r;
  s *= frame.s;
}

/*!
\brief Get the homogeneous matrix out of the frame.
*/
Matrix4 FrameScaled::GetMatrix4() const
{
  return Matrix4(r * Matrix(s), t);
}

/*!
\brief Overloaded.
\param s Stream.
\param frame The frame.
*/
std::ostream& operator<<(std::ostream& s, const FrameScaled& frame)
{
  s << "FrameScaled(" << frame.r << ',' << frame.t << ',' << frame.s << ')';
  return s;
}

/*!
\brief Creates a translation transformation.

\param t Translation vector.
*/
FrameScaled FrameScaled::Translation(const Vector& t)
{
  return FrameScaled(Matrix::Identity, t);
}

/*!
\brief Creates a rotation frame.
\param a A vector of angles (expressed in radians) that defines the rotation around x-y-z axes.
*/
FrameScaled FrameScaled::Rotation(const Vector& a)
{
  return FrameScaled(Matrix::Rotation(a), Vector(0.0));
}

/*!
\brief Create a rotation frame about an arbitrary axis.

\param axis %Axis.
\param angle Angle (should be in radian).
*/
FrameScaled FrameScaled::Rotation(const Vector& axis, const double& angle)
{
  return FrameScaled(Matrix::Rotation(axis, angle), Vector(0.0));
}

/*!
\brief Create a frame that rotates a normalized vector into another one.
\param a, b Initial and final vectors (should be normalized).
*/
FrameScaled FrameScaled::Rotation(const Vector& a, const Vector& b)
{
  return FrameScaled(Matrix::Rotation(a, b), Vector(0.0));
}

/*!
\brief Create a scaling frame.
\param s Scaling vector.
*/
FrameScaled FrameScaled::Scaling(const Vector& s)
{
  return FrameScaled(Matrix::Identity, Vector::Null, s);
}

/*!
\brief Create a uniform scaling frame.
\param s Scaling.
*/
FrameScaled FrameScaled::Scaling(const double& s)
{
  return FrameScaled(Matrix::Identity, Vector::Null, Vector(s));
}
