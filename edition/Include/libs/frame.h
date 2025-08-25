// Frames 

#pragma once

#include "libs/matrix.h"
#include "libs/ray.h"
#include <QtWidgets/QGraphicsScene>

// Frame class 
class Frame2
{
protected:
  Matrix2 r;  //!< Rotation matrix.
  Vector2 t;  //!< Translation vector.
public:
  Frame2(const Matrix2 & = Matrix2::Identity, const Vector2 & = Vector2::Null);

  explicit Frame2(const double&, const Vector2 & = Vector2::Null);
  explicit Frame2(const Vector2&, const Vector2&);

  //! Empty
  ~Frame2() {}

  // Translation transformation
  static Frame2 Translation(const Vector2&);

  // Rotation transformations
  static Frame2 Rotation(const double&);

  // Access to rotation matrix and translation vector
  Matrix2 R() const;
  Vector2 T() const;

  // Composition
  void Compose(const Frame2&);
  Frame2 Composed(const Frame2&) const;

  // Inverse
  Frame2 Inverse() const;

  Vector2 Transform(const Vector2&) const;
  Vector2 InverseTransform(const Vector2&) const;
  Vector2 TransformDirection(const Vector2&) const;
  Vector2 InverseTransformDirection(const Vector2&) const;

  // Rays
  friend std::ostream& operator<<(std::ostream&, const Frame2&);

  void Draw(QGraphicsScene&, const QPen&) const;

public:
  static const Frame2 Id; //!< Identity.
};

/*!
\brief Returns the rotation matrix of the frame.
*/
inline Matrix2 Frame2::R() const
{
  return r;
}

/*!
\brief Returns the translation vector of the frame.
*/
inline Vector2 Frame2::T() const
{
  return t;
}

/*!
\brief Transform a point out of the frame coordinate system.
\param p Point.
*/
inline Vector2 Frame2::Transform(const Vector2& p) const
{
  return r * p + t;
}

/*!
\brief Transform a point into the frame coordinate system.
\param p Point.
*/
inline Vector2 Frame2::InverseTransform(const Vector2& p) const
{
  return r.T() * (p - t);
}

/*!
\brief Transform a direction vector out of the frame coordinate system.
\param n %Vector.
*/
inline Vector2 Frame2::TransformDirection(const Vector2& n) const
{
  return r * n;
}

//! Transform a direction vector into the frame coordinate system.
inline Vector2 Frame2::InverseTransformDirection(const Vector2& n) const
{
  return r.T() * n;
}




// Frame class for instances
class FrameScaled2 :protected Frame2
{
protected:
  Vector2 s; //!< Scaling vector
public:
  explicit FrameScaled2(const Matrix2 & = Matrix2::Identity, const Vector2 & = Vector2::Null, const Vector2 & = Vector2(1.0));
  explicit FrameScaled2(const Vector2&, const Vector2&);
  explicit FrameScaled2(const Vector2&, const double&);

  // Conversion constructor
  FrameScaled2(const Frame2&, const Vector2 & = Vector2(1.0));

  //! Empty
  ~FrameScaled2() {}

  // Access to rotation matrix, translation vector, scaling vector
  Matrix2 R() const;
  Vector2 T() const;
  Vector2 S() const;

  FrameScaled2 Inverse() const;

  void Rotate(const double&);
  void Scale(const Vector2&);
  void Scale(const double&);
  void Translate(const Vector2&);

  // Scaling
  static FrameScaled2 Scaling(const Vector2&);
  static FrameScaled2 Scaling(const double&);

  Vector2 Transform(const Vector2&) const;
  Vector2 InverseTransform(const Vector2&) const;
  Vector2 TransformDirection(const Vector2&) const;
  Vector2 InverseTransformDirection(const Vector2&) const;

  void Compose(const FrameScaled2&);
  FrameScaled2 Composed(const FrameScaled2&) const;

  friend std::ostream& operator<<(std::ostream&, const FrameScaled2&);
public:
  static const FrameScaled2 Id; //!< Identity.
};

/*!
\brief Returns the rotation matrix of the frame.
*/
inline Matrix2 FrameScaled2::R() const
{
  return r;
}

/*!
\brief Returns the translation vector of the frame.
*/
inline Vector2 FrameScaled2::T() const
{
  return t;
}

// Frame class
class Frame
{
protected:
  Matrix r;  //!< Rotation matrix.
  Vector t;  //!< Translation vector.
public:
  Frame(const Matrix & = Matrix::Identity, const Vector & = Vector::Null);
  explicit Frame(const Vector&, const Vector&, const Vector&, const Vector&);

  // Conversion constructor
  Frame(const Frame2&);

  //! Empty
  ~Frame() {}

  // Translation transformation
  static Frame Translation(const Vector&);

  // Rotation transformations
  static Frame Rotation(const Vector&);
  static Frame Rotation(const Vector&, const double&);
  static Frame Rotation(const Vector&, const Vector&);
  static Frame Canonical(const Vector&, const Vector&);
  static Frame Orthonormal(const Vector&, const Vector&);

  // Access to rotation matrix and translation vector
  Matrix R() const;
  Vector T() const;

  // Composition
  void Compose(const Frame&);
  Frame Composed(const Frame&) const;

  // Inverse
  Frame Inverse() const;

  Vector Transform(const Vector&) const;
  Vector InverseTransform(const Vector&) const;
  Vector TransformDirection(const Vector&) const;
  Vector InverseTransformDirection(const Vector&) const;

  // Rays
  Ray Transform(const Ray&) const;
  Ray InverseTransform(const Ray&) const;

  friend std::ostream& operator<<(std::ostream&, const Frame&);

  // Coordinates of vertices in a given frame
  Vector CircleVertex(const double&, int = 1, int = 2) const;
  Vector CircleNormal(const double&, int = 1, int = 2) const;
  Vector SphereVertex(const double&, const double&, const double&, int = 1, int = 2, int = 0) const;
  Vector SphereNormal(const double&, const double&, int = 1, int = 2, int = 0) const;

public:
  static const Frame Id; //!< Identity.
};

/*!
\brief Returns the rotation matrix of the frame.

To get the the i-th basis vector of the frame, simply code:
\code
Frame frame;
Vector y=frame.R().C(1); // Get rotation matrix, and then column vector.
\endcode

\sa Frame2::R() const
*/
inline Matrix Frame::R() const
{
  return r;
}

/*!
\brief Returns the translation vector of the frame.
*/
inline Vector Frame::T() const
{
  return t;
}

/*!
\brief Transform a point out of the frame coordinate system.
\param p Point.
*/
inline Vector Frame::Transform(const Vector& p) const
{
  return r * p + t;
}

/*!
\brief Transform a point into the frame coordinate system.
\param p Point.
*/
inline Vector Frame::InverseTransform(const Vector& p) const
{
  return r.T() * (p - t);
}

/*!
\brief Transform a direction vector out of the frame coordinate system.
\param n %Vector.
*/
inline Vector Frame::TransformDirection(const Vector& n) const
{
  return r * n;
}

//! Transform a direction vector into the frame coordinate system.
inline Vector Frame::InverseTransformDirection(const Vector& n) const
{
  return r.T() * n;
}

// Frame class for instances
class FrameScaled :protected Frame
{
protected:
  Vector s; //!< Scaling vector
public:
  explicit FrameScaled(const Matrix & = Matrix::Identity, const Vector & = Vector::Null, const Vector & = Vector(1.0));
  explicit FrameScaled(const Vector&);
  explicit FrameScaled(const Vector&, const Vector&);
  explicit FrameScaled(const Vector&, const double&);
  explicit FrameScaled(const Vector&, const Vector&, const Vector&, const Vector&);

  // Conversion constructor
  FrameScaled(const Frame&, const Vector & = Vector(1.0));
  FrameScaled(const Frame2&, const Vector & = Vector(1.0));

  //! Empty
  ~FrameScaled() {}

  // Access to rotation matrix, translation vector, scaling vector
  Matrix R() const;
  Vector T() const;
  Vector S() const;

  void Rotate(const Vector&);
  void ObjectRotate(const Vector&);
  void Scale(const Vector&);
  void Scale(const double&);
  void Translate(const Vector&);

  FrameScaled Inverse() const;

  // Translations 
  static FrameScaled Translation(const Vector&);

  // Rotations 
  static FrameScaled Rotation(const Vector&);
  static FrameScaled Rotation(const Vector&, const double&);
  static FrameScaled Rotation(const Vector&, const Vector&);

  // Scaling
  static FrameScaled Scaling(const Vector&);
  static FrameScaled Scaling(const double&);

  // Compute matrix
  Matrix4 GetMatrix4() const;

  Vector Transform(const Vector&) const;
  Vector InverseTransform(const Vector&) const;
  Vector TransformDirection(const Vector&) const;
  Vector InverseTransformDirection(const Vector&) const;

  void Compose(const FrameScaled&);
  FrameScaled Composed(const FrameScaled&) const;

  void Lerp(const FrameScaled&, const FrameScaled&, const double&);

  friend std::ostream& operator<<(std::ostream&, const FrameScaled&);
public:
  static const FrameScaled Id; //!< Identity.
};

/*!
\brief Returns the rotation matrix of the frame.
*/
inline Matrix FrameScaled::R() const
{
  return r;
}

/*!
\brief Returns the translation vector of the frame.
*/
inline Vector FrameScaled::T() const
{
  return t;
}
