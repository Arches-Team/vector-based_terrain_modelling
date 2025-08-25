// Fundamentals

#pragma once


#include <QtGui/QImage>

#include "libs/array.h"
#include "libs/scalarfield.h"

class ScalarField2;
class PointCurve2;

class AnalyticVectorField2
{
public:
  virtual Vector2 Value(const Vector2&) const;
  virtual double Divergence(const Vector2&) const;
  virtual double Curl(const Vector2&) const;

protected:
  static const double epsilon; //!< Epsilon value for partial derivatives
};

class VectorField2 :public Array2
{
protected:
  QVector<Vector2> field; //!< Field samples.
public:
  //! Empty.
  VectorField2() {}
  explicit VectorField2(const Box2&, int, int, const Vector2 & = Vector2::Null);
  explicit VectorField2(const Box2&, int, int, const AnalyticVectorField2&);

  //! Empty
  ~VectorField2() {}

  void Set(const Vector2&);

  // Access to elements
  Vector2 at(int, int) const;
  Vector2 at(const QPoint&) const;
  Vector2& operator()(int, int);
  Vector2& operator()(const QPoint&);

  Vector2 at(int) const;
  Vector2& operator[](int);

  VectorField2& operator+= (const VectorField2&);
  VectorField2& operator-= (const VectorField2&);

  ScalarField2 GetNorm() const;
  ScalarField2 GetAngle() const;
  ScalarField2 GetComponent(int) const;

  QImage CreateAngleImage() const;

  void Normalize();
  void Unit();
  VectorField2 Orthogonal() const;
  void GetRange(Vector2&, Vector2&) const;
  void GetNormRange(double&, double&) const;
  double AverageNorm() const;

  virtual Vector2 Value(const Vector2&) const;
  virtual Vector2 Value(int, int) const;

  // Divergence
  virtual ScalarField2 Divergence() const;
  virtual double Divergence(int, int) const;

  void Add(const VectorField2&);
  void Sub(const VectorField2&);
  void Mul(const double&);

  // Curl
  virtual ScalarField2 Curl() const;
  virtual double Curl(int, int) const;

  void Rotate(const Matrix2&);

  void Lerp(const VectorField2&, const VectorField2&, const double&);

  // modify the VectorField
  void Multiply(const Vector2&, const double&, const double&);
  void Radial(const Vector2&, const double&, const double&);
  void Clone(const Vector2&, const Vector2&, const Vector2&, const double&);
  void CloneZone(const Vector2&, const QPolygon&, int, const VectorField2*  = nullptr);
  void CloneBrokenLine(const QVector<Vector2>&, const QVector<Vector2>&, int , const VectorField2*  = nullptr);

  QImage CreateImage() const;
  QImage CreateBlueRedImage() const;
  QImage CreateGradientImage(const double &) const;
  QImage CreateStreamLines() const;
  QImage CreateImageParticles(const double&) const;
  void Draw(QGraphicsScene&, const QPen & = QPen(), const QBrush & = QBrush()) const;

  Vector2 Euler(const Vector2&, const double&, bool = false) const;
  Vector2 RungeKutta(const Vector2&, const double&) const;
  PointCurve2 EulerSteps(const Vector2&, const double&, int, bool = false) const;
  ScalarField2 LineIntegral(int) const;

  VectorField2 GaussianBlur(const double&) const;
  VectorField2 GaussianBlur(int) const;
};

/*!
\brief Return the field value at a given array vertex.
\param i,j Integer coordinates of the vertex.
*/
inline Vector2 VectorField2::at(int i, int j) const
{
  return field.at(VertexIndex(i, j));
}

/*!
\brief Return the field value at a given array vertex.
\param p Point with its coordinates.
*/
inline Vector2 VectorField2::at(const QPoint& p) const
{
  return field.at(VertexIndex(p.x(), p.y()));
}

/*!
\brief Return the field value at a given array vertex.
\param i,j Integer coordinates of the vertex.
*/
inline Vector2& VectorField2::operator()(int i, int j)
{
  return field[VertexIndex(i, j)];
}

/*!
\brief Return the field value at a given array vertex.
\param p Point.
*/
inline Vector2& VectorField2::operator()(const QPoint& p)
{
  return field[VertexIndex(p.x(), p.y())];
}

/*!
\brief Return the data in the field.
\param c Index.
*/
inline Vector2 VectorField2::at(int c) const
{
  return field.at(c);
}

/*!
\brief Return the data in the field.
\param c Index.
*/
inline Vector2& VectorField2::operator[](int c)
{
  return field[c];
}

class AnalyticVectorField
{
public:
  virtual Vector Value(const Vector&) const;
  virtual double Divergence(const Vector&) const;
  virtual Vector Curl(const Vector&) const;

protected:
  static const double epsilon; //!< Epsilon value for partial derivatives
};

class PointCurve;

class VectorField :public Array
{
protected:
  QVector<Vector> field; //!< Field samples.
public:
  //! Empty.
  VectorField() {}
  explicit VectorField(const Box&, int, int, int, const Vector & = Vector::Null);
  explicit VectorField(const Box&, int, int, int, const AnalyticVectorField&);

  //! Empty
  ~VectorField() {}

  void Set(const Vector&);

  // Access to elements
  Vector at(int, int, int) const;
  Vector& operator()(int, int, int);

  Vector at(int) const;
  Vector& operator[](int);

  virtual Vector Value(const Vector&) const;
  virtual Vector Value(int, int, int) const;

  Vector Euler(const Vector&, const double&, bool = false) const;
  Vector RungeKutta(const Vector&, const double&) const;
  PointCurve EulerSteps(const Vector&, const double&, int, bool = false) const;
};

/*!
\brief Return the field value at a given array vertex.
\param i,j,k Integer coordinates of the vertex.
*/
inline Vector VectorField::at(int i, int j, int k) const
{
  return field.at(VertexIndex(i, j, k));
}

/*!
\brief Return the field value at a given array vertex.
\param i,j,k Integer coordinates of the vertex.
*/
inline Vector& VectorField::operator()(int i, int j, int k)
{
  return field[VertexIndex(i, j, k)];
}

/*!
\brief Return the data in the field.
\param c Index.
*/
inline Vector VectorField::at(int c) const
{
  return field.at(c);
}

/*!
\brief Return the data in the field.
\param c Index.
*/
inline Vector& VectorField::operator[](int c)
{
  return field[c];
}
