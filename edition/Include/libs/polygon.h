#pragma once

#include <QtCore/QVector>

#include "libs/segment.h"
#include "libs/triangle.h"
#include "libs/ellipse.h"

class Pentagon2;
class Hexagon2;

class Polygonal
{
protected:
  QVector<Vector> vertex; //!< Array of vertices.
public:
  Polygonal();
  explicit Polygonal(Vector*, int);
  explicit Polygonal(const QVector<Vector>&);
  explicit Polygonal(const QVector<Vector2>&);
  explicit Polygonal(const QVector<Segment>&);
  explicit Polygonal(const Box2&);

  Box GetBox() const;

  Vector Normal() const;

  void Translate(const Vector&);
  void Scale(const double&);
  void Rotate(const Matrix&);

  // Add vertex
  void Append(const Vector&);

  // Access to data
  Vector& Vertex(int);
  Vector Vertex(int) const;

  int Size() const;
  double Length() const;
  double PlanarLength() const;

  Vector PointAtLength(const double&) const;
  Vector NormalAtLength(const double&) const;

  Vector Normal(const Vector&) const;
  double R(const Vector&) const;

  // Area
  double Area() const;

  // Centers
  Vector Centroid() const;
  Vector Center() const;
  void Subdivide(int, QVector<Vector>&, QVector<int>&) const;

  friend std::ostream& operator<<(std::ostream&, const Polygonal&);
};

/*!
\brief Add a vertex to the polygon.
\param p Point.
*/
inline void Polygonal::Append(const Vector& p)
{
  vertex.append(p);
}
/*!
\brief Create an empty polygon.
*/
inline Polygonal::Polygonal()
{
}

//! Read write access to the i-th point.
inline Vector& Polygonal::Vertex(int i)
{
  return vertex[i];
}

//! Read only access to the i-th point.
inline Vector Polygonal::Vertex(int i) const
{
  return vertex.at(i);
}

//! Return the number of vertices of the polygon.
inline int Polygonal::Size() const
{
  return vertex.size();
}

// Forward declare classes: Quadrangle2
class Quadrangle2;

class Polygon2
{
protected:
  QVector<Vector2> q; //!< Array of vertices.
public:
  Polygon2();
  explicit Polygon2(Vector2*, int);
  explicit Polygon2(const Vector2&, const Vector2&, const Vector2&);
  explicit Polygon2(const Vector2&, const Vector2&, const Vector2&, const Vector2&);
  explicit Polygon2(const QVector<Vector>&);
  explicit Polygon2(const QVector<Vector2>&);
  explicit Polygon2(const QVector<Vector2>&, const QVector<int>&);

  explicit Polygon2(const Triangle2&);
  explicit Polygon2(const Box2&);
  explicit Polygon2(const Hexagon2&);
  explicit Polygon2(const Pentagon2&);
  explicit Polygon2(const Quadrangle2&);
  explicit Polygon2(const Polygonal&);
  explicit Polygon2(const Ellipse2&, int = 72);

  double Hausdorff(const Polygon2&, bool = false) const;

  Box2 GetBox() const;

  void Translate(const Vector2&);
  Polygon2 Translated(const Vector2&) const;
  void Scale(const double&);
  void Rotate(const Matrix2&);

  // Add vertex
  void Append(const Vector2&);

  // Access to vertices
  Vector2& Vertex(int);
  Vector2 Vertex(int) const;

  QVector<Vector2> Vertices();
  const QVector<Vector2>& Vertices() const;

  Vector2 Edge(int) const;

  int Size() const;
  double Length() const;
  bool IsConvex() const;

  Vector2 PointAtLength(const double&) const;
  Vector2 NormalAtLength(const double&) const;

  // Inside-outside
  bool Inside(const Vector2&) const;

  double R(const Vector2&) const;
  double RC(const Vector2&) const;
  double R(const Line2&) const;
  double Signed(const Vector2&) const;

  bool Intersect(const Circle2&) const;
  int Where(const Circle2&) const;

  // Static functions that avoid constructors
  static bool Inside(const QVector<Polygon2>&, const Vector2&);
  static double R(const QVector<Polygon2>&, const Vector2&);

  bool IntersectSegment(const Segment2&) const;

  // Area
  double Area() const;

  // Centers
  Vector2 Centroid() const;
  Vector2 Center() const;
  void Subdivide(int, QVector<Vector2>&, QVector<int>&) const;
  void Expand(const double&);

  QPolygonF GetQt() const;

  friend std::ostream& operator<<(std::ostream&, const Polygon2&);
  void Draw(QGraphicsScene&, const QPen & = QPen(), const QBrush & = QBrush()) const;

  Polygon2 Resampled(const double&) const;
  QVector<Vector2> Poisson(const double&, int, bool = false, Random & = Random::R239) const;
  Vector2 RandomInside(Random & = Random::R239) const;

  // Triangulate
  QVector<int> EarClip() const;
public:
  // Provide range-based for loops
  auto begin() { return q.begin(); }
  auto end() { return q.end(); }
  auto cbegin() const { return q.begin(); }
  auto cend() const { return q.end(); }
  auto begin() const { return q.begin(); }
  auto end() const { return q.end(); }
};

/*!
\brief Get the array of vertices.
*/
inline QVector<Vector2> Polygon2::Vertices()
{
  return q;
}

/*!
\brief Get the array of vertices.
*/
inline const QVector<Vector2>& Polygon2::Vertices() const
{
  return q;
}

/*!
\brief Add a vertex to the polygon.
\param p Point.
*/
inline void Polygon2::Append(const Vector2& p)
{
  q.append(p.ToVector());
}

/*!
\brief Create an empty polygon.
*/
inline Polygon2::Polygon2()
{
}

/*!
\brief Read write access to the i-th point.
\param i Index.
*/
inline Vector2& Polygon2::Vertex(int i)
{
  return q[i];
}

/*!
\brief Read only access to the i-th point.
\param i Index.
*/
inline Vector2 Polygon2::Vertex(int i) const
{
  return q.at(i);
}

/*!
\brief Return the i-th edge, starting from the i-th vertex.
\param i Index.
*/
inline Vector2 Polygon2::Edge(int i) const
{
  return q.at((i + 1) % q.size()) - q.at(i);
}

/*!
\brief Return the number of vertices of the polygon.
*/
inline int Polygon2::Size() const
{
  return q.size();
}

class Polygons2 {
protected:
  QVector<Polygon2> poly; //!< %Set of polygons.
public:
  Polygons2();
  Polygons2(const QVector<Polygon2>&);
  Polygons2(const QVector<Polygonal>&);

  int Size() const;
  //! Empty
  ~Polygons2() {}

  bool Inside(const Vector2&) const;
  double RC(const Vector2&) const;
  const Polygon2& At(int) const;
  void Add(const Polygon2&);
  void Draw(QGraphicsScene&, const QPen & = QPen(), const QBrush & = QBrush()) const;
};

/*!
\brief Return the size of the set.
*/
inline int Polygons2::Size() const
{
  return poly.size();
}

/*!
\brief Access.
\param i Index.
*/
inline const Polygon2& Polygons2::At(int i) const
{
  return poly.at(i);
}