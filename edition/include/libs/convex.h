// Convex

#pragma once

#include "libs/polygon.h"

class Convex2 :public Polygon2
{
public:
  //! Empty
  Convex2() {}
  explicit Convex2(Vector2*, int);
  explicit Convex2(const Vector2&, const Vector2&, const Vector2&);
  explicit Convex2(const QVector<Vector>&);
  explicit Convex2(const QVector<Vector2>&);

  explicit Convex2(const Box2&);
  explicit Convex2(const Hexagon2&);

  //! Empty
  ~Convex2() {}

  bool Cut(const Line2&);
  bool AddIntersection(const Segment2&);

  double R(const Vector2&) const;
  Vector2 Normal(const Vector2&) const;

  bool Inside(const Convex2&) const;

  static Convex2 Hull(QVector<Vector2>);

  static QVector<Convex2> VoronoiCells(const QVector<Vector2>&, const Box2&);
  static Convex2 VoronoiCell(const QVector<Vector2>&, const Box2&, const Vector2&);

  Convex2 Minkowski(const Convex2&) const;
};
