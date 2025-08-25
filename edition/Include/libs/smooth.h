// Smooth 

#pragma once

#include "libs/sphere.h"
#include "libs/circle.h"
#include "libs/segment.h"
#include "libs/ellipse.h"

// Smooth vector 
class SmoothVertex
{
protected:
  Vector c; //!< Center.
  double r; //!< Falloff radius.
  double s; //!< Intensity.
  bool cubic; //!< Cubic or quadric falloff.
public:
  explicit SmoothVertex(const Vector&, const double&, const double& = 1.0, bool = true);
  double Value(const Vector&) const;

  Vector Center() const;
  double Radius() const;
  Box GetBox() const;

  double K() const;
};

/*!
\brief Compute the bounding box.
*/
inline Box SmoothVertex::GetBox() const
{
  return Box(c, r);
}

/*!
\brief Return the radius.
*/
inline double SmoothVertex::Radius() const
{
  return r;
}

/*!
\brief Return the center.
*/
inline Vector SmoothVertex::Center() const
{
  return c;
}

/*!
\brief Compute the Lipschitz constant.
*/
inline double SmoothVertex::K() const
{
  return 1.72 * fabs(s) / r;
}


// Smooth segment 
class SmoothSegment :protected Segment
{
protected:
  double r = 0.0; //!< Falloff radius.
  double s = 1.0; //!< Intensity.
  bool cubic = true; //!< Cubic falloff, or quadric falloff
public:
  explicit SmoothSegment(const Vector&, const Vector&, const double&, const double& = 1.0, bool = true);
  double Value(const Vector&) const;

  double Radius() const;
  Box GetBox() const;
};

/*!
\brief Compute the bounding box.
*/
inline Box SmoothSegment::GetBox() const
{
  return Segment::GetBox().Extended(r);
}

/*!
\brief Return the radius.
*/
inline double SmoothSegment::Radius() const
{
  return r;
}

// Smooth sphere class
class SmoothSphere :protected Sphere
{
protected:
  double re = 0.0; //!< Falloff radius.
public:
  explicit SmoothSphere(const Vector&, const double&, const double&);
  double Value(const Vector&) const;
  Box GetBox() const;
};

/*!
\brief Compute the bounding box.
*/
inline Box SmoothSphere::GetBox() const
{
  return Sphere::GetBox().Extended(re);
}

// Smooth vector class
class SmoothVertex2
{
protected:
  Vector2 c; //!< Center.
  double r = 0.0; //!< Falloff radius.
public:
  explicit SmoothVertex2(const Vector2&, const double&);
  double Value(const Vector2&) const;
  Box2 GetBox() const;
};

/*!
\brief Compute the bounding box.
*/
inline Box2 SmoothVertex2::GetBox() const
{
  return Box2(c, r);
}

// Smooth disc class
class SmoothDisc2 :protected Disc2
{
protected:
  double re = 0.0; //!< Falloff radius.
public:
  explicit SmoothDisc2(const Vector2&, const double&, const double&);
  double Value(const Vector2&) const;
  Vector2 Gradient(const Vector2&) const;
  Box2 GetBox() const;
  Circle2 GetCircle() const;

  using Disc2::Center;
};

/*!
\brief Compute the bounding box.
*/
inline Box2 SmoothDisc2::GetBox() const
{
  return Disc2::GetBox().Extended(re);
}

/*!
\brief Compute the bounding circle.
*/
inline Circle2 SmoothDisc2::GetCircle() const
{
  return Circle2(c, r + re);
}

// Smooth box class
class SmoothBox2 :protected Box2
{
protected:
  double re; //!< Falloff radius.
public:
  explicit SmoothBox2(const Box2&, const double&);
  double Value(const Vector2&) const;
  Box2 GetBox() const;
};

/*!
\brief Compute the bounding box.
*/
inline Box2 SmoothBox2::GetBox() const
{
  return Box2(a, b).Extended(re);
}

class SmoothEllipse2 :protected Ellipse2
{
protected:
public:
  explicit SmoothEllipse2(const Ellipse2&);
  double Value(const Vector2&) const;

  using Ellipse2::Center;
  using Ellipse2::GetBox;
};
