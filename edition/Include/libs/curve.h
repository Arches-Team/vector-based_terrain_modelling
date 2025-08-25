// Curve
#pragma once


#include "libs/matrix.h"
#include "libs/quartic.h"
#include "libs/box.h"
#include "libs/ray.h"

// Quadric curve class
class QuadricCurve
{
protected:
  Quadric x, y, z; //!< %Quadric functions for every coordinate.
public:
  //! Empty.
  QuadricCurve() {}
  explicit QuadricCurve(const Quadric&, const Quadric&, const Quadric&);
  //! Empty.
  ~QuadricCurve() {}

  // Access class components
  Quadric& operator[] (int);
  Quadric operator[] (int) const;

  // Unary operators
  QuadricCurve operator+ () const { return *this; }
  QuadricCurve operator- () const;

  // Evaluates curve
  Vector operator()(const double&) const;
  Vector Eval(const double&) const;

  // Compute curvature 
  double Curvature(const double&) const;

  // Compute primary and secondary derivatives
  Vector Tangent(const double&) const;
  Vector Normal(const double&) const;

  // Get bounding box
  Box GetBox() const;

  // Compute distance between point and curve
  double R(const Vector&, double&) const;
  double Signed(const Vector&) const;

  // Compute the curvilign abscissa of a point on the curve
  double S(const double&, int = 256) const;
  double U(const double&, int = 256) const;

  // Computes the length of the quadric curve
  double LengthIntegral(const double&, const double&) const;
  double Length(int = 256) const;
  double Length(const double&, const double&, int = 256) const;

  // Frenet Frame
  Matrix GetMatrix(const double&, bool = false) const;
  Frame GetFrame(const double&, bool = false) const;
  QVector<Frame> GetFrames(const double&, bool = false) const;

  static QuadricCurve Bezier(const Vector&, const Vector&, const Vector&);
  static QVector<QuadricCurve> Bezier(const QVector<Vector>&);

  friend std::ostream& operator<<(std::ostream&, const QuadricCurve&);
};

//! Access curve components
inline Quadric& QuadricCurve::operator[] (int i)
{
  if (i == 0) { return x; }
  else if (i == 1) { return y; }
  else { return z; }
}

//! Access curve components
inline Quadric QuadricCurve::operator[] (int i) const
{
  if (i == 0) { return x; }
  else if (i == 1) { return y; }
  else { return z; }
}

// Cubic curve class
class CubicCurve
{
protected:
  Cubic x, y, z; //!< %Cubic polynomial functions for every coordinate.
public:
  //! Empty.
  CubicCurve() {}
  CubicCurve(const Cubic&, const Cubic&, const Cubic&);
  //! Empty.
  ~CubicCurve() {}

  // Access class components
  Cubic& operator[] (int);
  Cubic operator[] (int) const;

  // Unary operators

  //! Overloaded out of consistency.
  CubicCurve operator+ () const { return *this; }
  CubicCurve operator- () const;

  CubicCurve Translated(const Vector&) const;

  // Computes the tangent curve equation
  QuadricCurve Tangent() const;
  Vector Tangent(const double&) const;
  Vector Normal(const double&) const;

  double Curvature(const double&) const;

  // Compute the curvilign abscissa of a point on the curve
  double S(const double&, int = 256) const;
  double U(const double&, int = 256) const;

  // Computes the length of the cubic curve
  double Length(int = 256) const;
  double Length(const double&, const double&, int = 256) const;
  double InverseArcLength(const double&, int = 256) const;

  // Evaluates curve
  Vector operator()(const double&) const;
  Vector Eval(const double&) const;

  // Get bounding box
  Box GetBox() const;

  // Compute distance between point and curve
  double R(const Vector&, double&) const;

  // Frenet Frame
  Matrix GetMatrix(const double&, bool = false) const;
  Frame GetFrame(const double&, bool = false) const;
  QVector<Frame> GetFrames(const double&, bool = false) const;

  // Approximation
  Vector BezierControl(int) const;
  void Approximate(QuadricCurve&, QuadricCurve&, const double& = 0.5) const;

  friend std::ostream& operator<<(std::ostream&, const CubicCurve&);
public:
  static CubicCurve Hermite(const Vector&, const Vector&, const Vector&, const Vector&);
  static CubicCurve Bezier(const Vector&, const Vector&, const Vector&, const Vector&);
};

//! Access curve components
inline Cubic& CubicCurve::operator[] (int i)
{
  if (i == 0) { return x; }
  else if (i == 1) { return y; }
  else { return z; }
}

//! Access curve components
inline Cubic CubicCurve::operator[] (int i) const
{
  if (i == 0) { return x; }
  else if (i == 1) { return y; }
  else { return z; }
}

// Quadric curve class
class QuadricCurve2
{
protected:
  Quadric x, y; //!< Quadric functions for every coordinate.
public:
  //! Empty.
  QuadricCurve2() {}
  explicit QuadricCurve2(const Quadric&, const Quadric&);
  explicit QuadricCurve2(const QuadricCurve&);

  //! Empty.
  ~QuadricCurve2() {}

  // Access class components
  Quadric& operator[] (int);
  Quadric operator[] (int) const;

  // Unary operators
  QuadricCurve2 operator+ () const { return *this; }
  QuadricCurve2 operator- () const;

  // Evaluates curve
  Vector2 operator()(const double&) const;
  Vector2 Eval(const double&) const;

  // Curvature 
  double Curvature(const double&) const;

  // Computes the tangent curve equation
  Vector2 Tangent(const double&) const;
  Vector2 Normal(const double&) const;

  // Get bounding box
  Box2 GetBox() const;
  void Translate(const Vector2&);

  // Compute distance between point and curve
  double R(const Vector2&, double&) const;

  // Compute the curvilign abscissa of a point on the curve
  double S(const double&, int = 256) const;
  double U(const double&, int = 256) const;

  // Computes the length of the quadric curve
  double LengthIntegral(const double& = 0.0, const double& = 1.0) const;
  double Length(int = 256) const;
  double Length(const double&, const double&, int = 256) const;

  double UV(const Vector2&, double&, double&) const;
  Vector2 Normal(const Vector2&, double&) const;

  // Frenet Frame
  Matrix2 GetMatrix(const double&) const;
  Frame2 GetFrame(const double&) const;
  QVector<Frame2> GetFrames(const double&) const;

  int Intersect(const Ray2&, double[2]) const;
  static QuadricCurve2 Bezier(const Vector2&, const Vector2&, const Vector2&);
  static QuadricCurve2 ThroughPoints(const Vector2&, const Vector2&, const Vector2&, const double& = 0.5);
  friend std::ostream& operator<<(std::ostream&, const QuadricCurve2&);

  void Draw(QGraphicsScene&, const QPen & = QPen(), const QBrush & = QBrush()) const;
public:
  static QVector<QuadricCurve2> Bezier(const QVector<Vector2>&);
protected:
  // Compute tdistance between point and curve
  double R(const Vector2&, double*, int) const;
};

//! Access curve components
inline Quadric& QuadricCurve2::operator[] (int i)
{
  if (i == 0) { return x; }
  else { return y; }
}

//! Access curve components
inline Quadric QuadricCurve2::operator[] (int i) const
{
  if (i == 0) { return x; }
  else { return y; }
}

// Quadric curve class
class CubicCurve2
{
protected:
  Cubic x, y; //!< %Cubic polynomial functions for every coordinate.
public:
  //! Empty.
  CubicCurve2() {}
  explicit CubicCurve2(const Cubic&, const Cubic&);
  explicit CubicCurve2(const CubicCurve&);

  //! Empty.
  ~CubicCurve2() {}

  // Access class components
  Cubic& operator[] (int);
  Cubic operator[] (int) const;

  // Conversion
  CubicCurve ToCubicCurve(const Cubic & = Cubic(0.0, 0.0, 0.0, 0.0)) const;

  // Unary operators

  // Overloaded out of consistency
  CubicCurve2 operator+ () const { return *this; }
  CubicCurve2 operator- () const;

  // Computes the tangent curve equation
  QuadricCurve2 Tangent() const;
  Vector2 Tangent(const double&) const;
  Vector2 Normal(const double&) const;

  // Compute the curvilign abscissa of a point on the curve
  double S(const double&, int = 256) const;
  double U(const double&, int = 256) const;

  int Intersect(const Ray2&, double[3]) const;

  // Computes the length of the cubic curve
  double Length(int = 256) const;
  double Length(const double&, const double&, int = 256) const;
  double InverseArcLength(const double&, int = 256) const;

  // Evaluates curve
  Vector2 operator()(const double&) const;
  Vector2 Eval(const double&) const;

  // Curvature 
  double Curvature(const double&) const;
  double Curvature(const double&, const double&) const;

  // Get bounding box
  Box2 GetBox() const;

  // Compute distance between point and curve
  double R(const Vector2&, double&) const;
  Vector2 Normal(const Vector2&, double&) const;
  double UV(const Vector2&, double&, double&) const;
  Vector2 UV(const Vector2&) const;

  // Frenet Frame
  Matrix2 GetMatrix(const double&) const;
  Frame2 GetFrame(const double&) const;
  QVector<Frame2> GetFrames(const double&) const;

  Vector2 BezierControl(int) const;
  void Approximate(QuadricCurve2&, QuadricCurve2&, const double& =0.5) const;
  
  // Drawing
  void Draw(QGraphicsScene&, const QPen & = QPen(), const QBrush & = QBrush()) const;

  friend std::ostream& operator<<(std::ostream&, const CubicCurve2&);

public:
  static CubicCurve2 Hermite(const Vector2&, const Vector2&, const Vector2&, const Vector2&);
  static CubicCurve2 Bezier(const Vector2&, const Vector2&, const Vector2&, const Vector2&);
};

//! Access curve components
inline Cubic& CubicCurve2::operator[] (int i)
{
  if (i == 0) { return x; }
  else { return y; }
}

//! Access curve components
inline Cubic CubicCurve2::operator[] (int i) const
{
  if (i == 0) { return x; }
  else { return y; }
}

/*!
\brief Convert a planar cubic curve to a cubic curve.

This is a convenience function, which is basically the same as:
\code
CubicCurve2 c; // Curve in the plane
Cubic z;       // Cubic for the z coordinate
CubicCurve3 c3=CubicCurve(c[0],c[1],z);
\endcode
\param z Missing cubic.
*/
inline CubicCurve CubicCurve2::ToCubicCurve(const Cubic& z) const
{
  return CubicCurve(x, y, z);
}

/*!
\brief Compute the i-th BÃ©zier control point of the curve.

Computations have been optimized.
\param i Index.
*/
inline Vector2 CubicCurve2::BezierControl(int i) const
{
  if (i == 0)
  {
    // Equivalent to c(0.0)
    return Vector2(x[0], y[0]);
  }
  else if (i == 1)
  {
    return Vector2(x[0] + x[1] / 3.0, y[0] + y[1] / 3.0);
  }
  else if (i == 2)
  {
    return Vector2(x[0] + (2.0 * x[1] + x[2]) / 3.0, y[0] + (2.0 * y[1] + y[2]) / 3.0);
  }
  else
  {
    // Equivalent to c(1.0)
    return Vector2(x[0] + x[1] + x[2] + x[3], y[0] + y[1] + y[2] + y[3]);
  }
}

class FrameCurve {
protected:
  QVector<Frame> frame; //!< Array of frames
  QVector<bool> smooth; //!< Array of flags defining whether a frame in the curve should be considered a smooth or not
  bool closed = false; //!< Closed curve flag.
public:
  //! Empty.
  FrameCurve() {}
  FrameCurve(const QVector<Frame>&, bool = false);
  FrameCurve(const QVector<Frame>&, const QVector<bool>&, bool = false);
  //! Empty.
  ~FrameCurve() {}

  // Access 
  Vector Vertex(int) const;

  Frame& GetFrame(int);
  Frame GetFrame(int) const;

  bool& Smooth(int);
  bool Smooth(int) const;

  // Get bounding box
  Box GetBox() const;

  // Computes the length of the quadric curve
  double Length() const;
  int Size() const;

  void Translate(const Vector&);
  void Rotate(const Matrix&);

  bool IsClosed() const;
};

//! Access to curve vertices
inline Vector FrameCurve::Vertex(int i) const
{
  return frame.at(i).T();
}

//! Access to curve frames
inline Frame& FrameCurve::GetFrame(int i)
{
  return frame[i];
}

//! Access to curve frames
inline Frame FrameCurve::GetFrame(int i) const
{
  return frame.at(i);
}

//! Access to curve vertices
inline bool& FrameCurve::Smooth(int i)
{
  return smooth[i];
}

//! Access to curve vertices
inline bool FrameCurve::Smooth(int i) const
{
  return smooth.at(i);
}

//! Check if the curve is closed
inline bool FrameCurve::IsClosed() const
{
  return closed;
}

/*!
\brief Return the size of the sampled curve.
*/
inline int FrameCurve::Size() const
{
  return frame.size();
}

