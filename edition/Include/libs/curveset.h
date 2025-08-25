// Curve

#pragma once

#include "libs/curve.h"

class QuadricCurveSet
{
protected:
  QVector<QuadricCurve> curve; //!< Set of quadric curves.
  QVector<double> lengths;     //!< Length of every curve (internal optimization)
  double length = 0.0;       //!< Total length (internal optimization)
public:
  QuadricCurveSet();
  explicit QuadricCurveSet(const QVector<Vector>&);
  explicit QuadricCurveSet(const QVector<QuadricCurve>&);
  //! Empty
  ~QuadricCurveSet() {}

  // Frenet Frame
  Matrix GetMatrix(const double&) const;
  Frame GetFrame(const double&) const;

  // Access to members
  int Size() const;
  QuadricCurve operator()(int) const;

  QVector<Vector> GetDiscretisation(const double&) const;
  int U(const double&, double&) const;
  double R(const Vector&, double&, int&) const;
  double R(const Vector&) const;
  double Signed(const Vector&) const;

  // Get bounding box
  Box GetBox() const;

  // Get length
  double GetLength() const;
  double GetLength(int) const;

  QuadricCurveSet& operator+=(const QuadricCurveSet&);
};

/*!
\brief Returns the number of elements of the piecewise quadric curve.
*/
inline int QuadricCurveSet::Size() const
{
  return curve.size();
}

/*!
\brief Read only access to the elements of the piecewise quadric curve.
*/
inline QuadricCurve QuadricCurveSet::operator()(int k) const
{
  return curve[k];
}

/*!
\brief Return the length of the curve.
*/
inline double QuadricCurveSet::GetLength() const
{
  return length;
}

/*!
\brief Return the length of the sub-curve.
\param k Curve index.
*/
inline double QuadricCurveSet::GetLength(int k) const
{
  return lengths.at(k);
}

class CubicCurveSet
{
protected:
  QVector<CubicCurve> curve; //!< Set of cubic (spline) curves.
  QVector<double> lengths; //!< Length of every curve (internal optimization).
  double length = 0.0; //!< Total length (internal optimization).
public:
  CubicCurveSet();
  explicit CubicCurveSet(const QVector<CubicCurve>&);
  explicit CubicCurveSet(const QVector<Vector>&, const Vector & = Vector::Null, const Vector & = Vector::Null);
  //! Empty
  ~CubicCurveSet() {}

  // Frenet Frame
  Matrix GetMatrix(const double&) const;
  Frame GetFrame(const double&) const;

  // Access to members
  int Size() const;
  CubicCurve operator()(int) const;

  // Get bounding box
  Box GetBox() const;

  // Quadric approximation
  QuadricCurveSet Approximate(double) const;

  // Get length
  double GetLength() const;
  double GetLength(int) const;

  QVector<Vector> GetDiscretisation(const double&) const;
  QVector<Vector> GetDiscretisation(const double&, QVector<Vector>&) const;
  int U(const double&, double&) const;
  double R(const Vector&, double&, int&) const;
};

inline CubicCurve CubicCurveSet::operator()(int i) const
{
  return curve[i];
}

/*!
\brief Return the number of elements of the piecewise cubic curve.
*/
inline int CubicCurveSet::Size() const
{
  return curve.size();
}

/*!
\brief Return the length of the curve.
*/
inline double CubicCurveSet::GetLength() const
{
  return length;
}

/*!
\brief Return the length of the sub-curve.
\param k Curve index.
*/
inline double CubicCurveSet::GetLength(int k) const
{
  return lengths.at(k);
}

class CubicCurve2Set
{
protected:
  QVector<CubicCurve2> curve; //!< Set of 2D cubic curves.
  QVector<double> lengths; //!< Length of every curve.
  double length = 0.0; //!< Total length.
public:
  CubicCurve2Set();
  explicit CubicCurve2Set(const QVector<CubicCurve2>&);
  explicit CubicCurve2Set(const QVector<Vector2>&, const Vector2 & = Vector2::Null, const Vector2 & = Vector2::Null);
  explicit CubicCurve2Set(const CubicCurveSet&);

  //! Empty
  ~CubicCurve2Set() {}

  // Access to members
  int Size() const;
  CubicCurve2 operator()(int) const;

  QVector<Vector2> GetDiscretisation(const double&) const;
  QVector<Vector2> GetDiscretisation(const double&, QVector<Vector2>&) const;

  // Get bounding box
  Box2 GetBox() const;

  // Get length
  double GetLength() const;
  double GetLength(int) const;

  double Sinuosity() const;

  bool Inside(const Vector2&) const;

  int U(const double&, double&) const;
  double R(const Vector2&, double&, int&) const;
  double Signed(const Vector2&) const;
  int Intersect(const Ray2& ray, int i, double t[3]) const;

  void UV(const Vector2&, double&, double&, int&) const;
  void Draw(QGraphicsScene&, const QPen & = QPen(), const QBrush & = QBrush()) const;
};

/*!
\brief Return the i-th sub-curve.
\param i Index.
*/
inline CubicCurve2 CubicCurve2Set::operator()(int i) const
{
  return curve.at(i);
}

/*!
\brief Return the number of elements of the piecewise cubic curve.
*/
inline int CubicCurve2Set::Size() const
{
  return curve.size();
}

/*!
\brief Return the length of the curve.
*/
inline double CubicCurve2Set::GetLength() const
{
  return length;
}

/*!
\brief Return the length of the sub-curve.
\param k Curve index.
*/
inline double CubicCurve2Set::GetLength(int k) const
{
  return lengths.at(k);
}

class QuadricCurve2Set
{
protected:
  QVector<QuadricCurve2> curve; //!< Set of quadric curves.
  QVector<double> lengths;     //!< Length of every curve (internal optimization)
  double length = 0.0;       //!< Total length (internal optimization)
public:
  QuadricCurve2Set();
  explicit QuadricCurve2Set(const QVector<Vector2>&);
  explicit QuadricCurve2Set(const QVector<QuadricCurve2>&);
  explicit QuadricCurve2Set(const QuadricCurveSet&);

  ~QuadricCurve2Set();

  // Access to members
  int Size() const;
  QuadricCurve2 operator()(int) const;

  void Translate(const Vector2&);

  // Get bounding box
  Box2 GetBox() const;

  // Get length
  double GetLength() const;
  double GetLength(int) const;

  double Sinuosity() const;

  QVector<Vector2> GetDiscretisation(const double&) const;
  int U(const double&, double&) const;
  double R(const Vector2&, double&, int&) const;
  double UV(const Vector2&, double&, double&, int&) const;

  bool Inside(const Vector2&) const;
  bool Intersect(const Ray2&, double&) const;
  int Intersections(const Ray2&) const;

  QuadricCurve2Set& operator+=(const QuadricCurve2Set&);
  void Draw(QGraphicsScene&, const QPen & = QPen(), const QBrush & = QBrush()) const;
};

/*!
\brief Returns the number of elements of the piecewise quadric curve.
*/
inline int QuadricCurve2Set::Size() const
{
  return curve.size();
}

/*!
\brief Read only access to the elements of the piecewise quadric curve.
*/
inline QuadricCurve2 QuadricCurve2Set::operator()(int k) const
{
  return curve[k];
}

/*!
\brief Return the length of the curve.
*/
inline double QuadricCurve2Set::GetLength() const
{
  return length;
}

/*!
\brief Return the length of the sub-curve.
\param k Curve index.
*/
inline double QuadricCurve2Set::GetLength(int k) const
{
  return lengths.at(k);
}
