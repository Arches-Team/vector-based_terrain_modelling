// Fundamentals

#pragma once

#include "libs/box.h"
#include "libs/circle.h"
#include "libs/vectorset.h"

class ParticleSet2 { // GUNDO TODO BUG :protected VectorSet2 {
protected:
  QVector<Vector2> points; //!< Set of particles
  double r = 0.0; //!< Radius of particles.
public:
  ParticleSet2();
  explicit ParticleSet2(const double&);
  explicit ParticleSet2(const Vector2&, const double&);
  explicit ParticleSet2(const QVector<Vector2>&, const double&);

  int Size() const;

  Circle2 At(int) const;
  double Radius() const;
  QVector<Vector2> GetCenters() const;

  Circle2 GetCircle() const;
  Box2 GetBox() const;

  void Append(const Vector2&);

  bool Intersect(const Circle2&) const;
  bool Intersect(const Box2&) const;

  void Draw(QGraphicsScene&, const QPen & = QPen(QColor(150, 150, 200), 0.025), const QBrush & = QBrush(QColor(200, 200, 250))) const;
};

/*!
\brief Empty constructor.
*/
inline ParticleSet2::ParticleSet2() :r(1.0) {}

/*!
\brief Return the radius of the particles.
*/
inline double ParticleSet2::Radius() const
{
  return r;
}

/*!
\brief Return the number of particles.
*/
inline int ParticleSet2::Size() const
{
  return points.size();
}

/*!
\brief Return the centers of the particles.
*/
inline QVector<Vector2> ParticleSet2::GetCenters() const
{
  return points;
}