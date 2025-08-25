// Frames 

#include "libs/frame.h"
#include "libs/circle.h"
#include "libs/segment.h"

/*!
\class Frame2 frame.h
\brief Solid transformations in the plane.

\sa Frame
\ingroup PlanarGroup
*/

const Frame2 Frame2::Id(Matrix2::Identity, Vector2::Null);

/*!
\brief Creates a frame given a rotation angle and a translation vector.
\param a Angle.
\param t Translation vector.
*/
Frame2::Frame2(const double& a, const Vector2& t) :r(Matrix2::Rotation(a)), t(t)
{
}

/*!
\brief Creates a frame given a rotation matrix and a translation vector.
\param r Frame matrix.
\param t Translation vector.
*/
Frame2::Frame2(const Matrix2& r, const Vector2& t) :r(r), t(t)
{
}

/*!
\brief Creates a frame given the origin and its orthogonal unit vectors.
\param c Origin.
\param x Axis, the orthonormal vectors will be set as (x,x&#10178;).
*/
Frame2::Frame2(const Vector2& c, const Vector2& x) :Frame2(Matrix2(x, x.Orthogonal()), c)
{
}

/*!
\brief Creates a translation transformation.

%Matrix is identity whereas the translation vector is
computed from the argument vector.
\param t Translation vector.
*/
Frame2 Frame2::Translation(const Vector2& t)
{
  return Frame2(Matrix2::Identity, t);
}

/*!
\brief Create a rotation frame.

\param angle Angle (should be in radian).
*/
Frame2 Frame2::Rotation(const double& angle)
{
  return Frame2(Matrix2::Rotation(angle), Vector(0.0));
}

/*!
\brief Overloaded.
\param s Stream.
\param frame The frame.
*/
std::ostream& operator<<(std::ostream& s, const Frame2& frame)
{
  s << "Frame2(" << frame.r << ',' << frame.t << ')';
  return s;
}

/*!
\brief Compute the inverse transformation.
*/
Frame2 Frame2::Inverse() const
{
  return Frame2(r.T(), -(r.T() * t));
}

/*!
\brief Compose the frame with another one.

\param frame The frame.
*/
void Frame2::Compose(const Frame2& frame)
{
  t = frame.r * t + frame.t;
  r = frame.r * r;
}

/*!
\brief Compose the frame with another one.

\param frame The frame.
*/
Frame2 Frame2::Composed(const Frame2& frame) const
{
  return Frame2(frame.r * r, frame.r * t + frame.t);
}


/*!
\brief Draw a circle.
\param scene Graphics scene.
\param pen The pen.
*/
void Frame2::Draw(QGraphicsScene& scene, const QPen& pen) const
{
  Segment2(t, t + r * Vector::X).DrawArrow(scene, 0.1, pen, QBrush(pen.color()));
  Segment2(t, t + r * Vector::Y).DrawArrow(scene, 0.1, pen, QBrush(pen.color()));

  Circle2(t, 0.1).Draw(scene, pen, QBrush(QColor(255, 255, 255)));
}
