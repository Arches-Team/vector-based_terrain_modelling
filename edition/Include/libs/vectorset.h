// Vector sets

#pragma once

#include <QtCore/QVector> 

#include "libs/evector.h"
#include "libs/box.h"
#include "libs/convex.h"
#include "libs/octogon.h"

class VectorSet2;

// Sphere class
class VectorSet
{
protected:
  QVector<Vector> v; //!< Set of vectors. 
public:
  //! Empty.
  VectorSet() {}
  explicit VectorSet(const VectorSet2&);
  explicit VectorSet(const QVector<Vector>&);
  explicit VectorSet(Vector*, unsigned int);

  VectorSet VectorToPoint() const;
  void Reverse();

  //! Empty.
  ~VectorSet() {}

  Vector At(int) const;
  Vector Vertex(int) const;
  Vector& operator[](int);
  int Size() const;

  // Parameters
  Vector Barycenter() const;
  void Append(const Vector&);

  Box GetBox() const;

  void Rotate(const Matrix&);
  void Translate(const Vector&);
  void Scale(const double&);
  void Scale(const Vector&);
  void Transform(const Frame&);
  void Transform(const FrameScaled&);

  QVector<Vector>& Get();
  const QVector<Vector>& Get() const;
  void Clean(const double&);

  void Load(const QString&);

  int Nearest(const Vector&) const;
  QVector<Vector> Nearest(const Vector&, int) const;
  QVector<int> NearestIndexes(const Vector&, int) const;
  QVector<int> NearestIndexes(int, int, const double& = Math::Infinity) const;
  friend std::ostream& operator<<(std::ostream&, const VectorSet&);

public:
  // Provide range-based for loops
  auto begin() { return v.begin(); }
  auto end() { return v.end(); }
  auto cbegin() const { return v.begin(); }
  auto cend() const { return v.end(); }
  auto begin() const { return v.begin(); }
  auto end() const { return v.end(); }
};

/*!
\brief Access to the internal structure.
*/
inline QVector<Vector>& VectorSet::Get()
{
  return v;
}

/*!
\brief Access to the internal structure.
*/
inline const QVector<Vector>& VectorSet::Get() const
{
  return v;
}

/*!
\brief Access to element.
*/
inline Vector VectorSet::At(int i) const
{
  return v.at(i);
}

/*!
\brief Access to element.
*/
inline Vector VectorSet::Vertex(int i) const
{
  return v.at(i);
}

/*!
\brief Access to element.
*/
inline Vector& VectorSet::operator[](int i)
{
  return v[i];
}

/*!
\brief Return the size of the set.
*/
inline int VectorSet::Size() const
{
  return v.size();
}

// Sphere class
class VectorSet2
{
protected:
  QVector<Vector2> v; //!< Vectors.
public:
  //! Empty.
  VectorSet2() {}
  VectorSet2(const QVector<Vector2>&);

  //! Empty.
  ~VectorSet2() {}

  VectorSet2 VectorToPoint() const;
  void Reverse();
  
  int Size() const;

  Vector2 At(int) const;
  Vector2& operator[](int);

  // Parameters
  Vector2 Barycenter() const;
  void Append(const Vector2&);
  void Append(const QVector<Vector2>&);
  void Append(const VectorSet2&);

  void Vibrate(const double&);

  Box2 GetBox() const;
  Convex2 GetHull() const;

  void Rotate(const Matrix2&);
  void Rotate(const double&);
  void Translate(const Vector2&);
  void Scale(const double&);
  void Scale(const Vector2&);

  VectorSet2 Transformed(const Frame2&) const;
  VectorSet2 Translated(const Vector2&) const;
  VectorSet2 Rotated(const double&) const;
  VectorSet2 Scaled(const Vector2&) const;

  VectorSet2 Cut(const Box2&) const;
  VectorSet2 Cut(const Octogon2&) const;

  QVector<Vector2>& Get();
  const QVector<Vector2>& Get() const;

  int Nearest(const Vector2&) const;
  QVector<Vector2> Nearest(const Vector2&, int) const;
  QVector<int> NearestIndexes(const Vector2&, int) const;
  QVector<int> NearestIndexes(int, int, const double& = Math::Infinity) const;

  friend std::ostream& operator<<(std::ostream&, const VectorSet2&);

  VectorSet2 operator+(const VectorSet2&) const;
  void operator+=(const VectorSet2&);
public:
  // Provide range-based for loops
  auto begin() { return v.begin(); }
  auto end() { return v.end(); }
  auto cbegin() const { return v.begin(); }
  auto cend() const { return v.end(); }
  auto begin() const { return v.begin(); }
  auto end() const { return v.end(); }
};

/*!
\brief Access to element.
*/
inline Vector2 VectorSet2::At(int i) const
{
  return v.at(i);
}

/*!
\brief Access to element.
*/
inline Vector2& VectorSet2::operator[](int i)
{
  return v[i];
}

/*!
\brief Access to the internal structure.
*/
inline QVector<Vector2>& VectorSet2::Get()
{
  return v;
}
/*!
\brief Access to the internal structure.
*/
inline const QVector<Vector2>& VectorSet2::Get() const
{
  return v;
}

/*!
\brief Return the size of the set.
*/
inline int VectorSet2::Size() const
{
  return v.size();
}
