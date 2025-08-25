// Curve
#pragma once


#include "libs/vectorset.h"
#include "libs/curve.h"

class QuadricCurve2Set;
class CubicCurve2Set;

class PointCurve :public VectorSet {
protected:
  bool closed = false; //!< Closed curve flag.
public:
  //! Empty.
  PointCurve() { }
  explicit PointCurve(const CubicCurve&, int, const double& = 0.0, const double& = 1.0);
  explicit PointCurve(const QuadricCurve&, int, const double& = 0.0, const double& = 1.0);
  explicit PointCurve(const QVector<Vector>&, bool = false);

  //! Empty.
  ~PointCurve() {}

  int Size() const;
  Vector Tangent(int) const;

  bool IsClosed() const;

  // Projection
  double R(const Vector&, double&, int&) const;


  void Rotate(const Matrix&);
  void Translate(const Vector&);
  void Scale(const double&);
  void Transform(const Frame&);
  PointCurve Transformed(const Frame&) const;

  // Compute the length of the curve
  double Length() const;
  static double Length(const QVector<Vector>&, bool = false);
};

/*!
\brief Check if the curve is closed
*/
inline bool PointCurve::IsClosed() const
{
  return closed;
}

/*!
\brief Return the size of the sampled curve.
*/
inline int PointCurve::Size() const
{
  return v.size();
}

class PointCurve2 :public VectorSet2 {
protected:
  bool closed = false; //!< Closed curve flag.
public:
  //! Empty.
  PointCurve2() {}
  explicit PointCurve2(const VectorSet2&, bool = false);
  explicit PointCurve2(const CubicCurve2&, int, const double& = 0.0, const double& = 1.0);
  explicit PointCurve2(const QuadricCurve2&, int, const double& = 0.0, const double& = 1.0);
  explicit PointCurve2(const QVector<Vector2>&, bool = false);
  explicit PointCurve2(const QVector<Vector2>&, const QVector<int>&, bool = false);

  //! Empty.
  ~PointCurve2() {}

  int Size() const;
  Vector2 Tangent(int) const;
  double Curvature(int) const;

  PointCurve Transform(const Frame&) const;

  bool IsClosed() const;
  void Close();

  // Computes the length of the quadric curve
  double Length() const;
  static double Length(const QVector<Vector2>&, bool = false);

  double Sinuosity() const;
  Polygon2 GetPolygon() const;
  CubicCurve2Set ToCubicCurve() const;
  QuadricCurve2Set ToQuadricCurve() const;

  PointCurve2 Translated(const Vector2&) const;
  PointCurve2 Transformed(const Frame2&) const;
  PointCurve2 Rotated(const double&) const;
  PointCurve2 Scaled(const Vector2&) const;

  // Projection
  double R(const Vector2&, double&, int&) const;
  double R(const Vector2&, double&, int&, int&) const;

  // Draw
  void Draw(QGraphicsScene&, const QPen & = QPen()) const;
};

/*!
\brief Check if the curve is closed
*/
inline bool PointCurve2::IsClosed() const
{
  return closed;
}

/*!
\brief Define the curve as closed.
*/
inline void PointCurve2::Close()
{
  closed = true;
}

/*!
\brief Return the size of the sampled curve.
*/
inline int PointCurve2::Size() const
{
  return v.size();
}
