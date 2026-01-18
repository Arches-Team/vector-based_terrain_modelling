// Segment
#pragma once

// Forward declare seldom used pre-processing intensive Qt class
class QGraphicsScene;

#include "libs/matrix.h"
#include "libs/box.h"

// Forward declare classes: Quadric
class Quadric;

// Segment class
class Segment
{
protected:
  Vector a = Vector::Null, b = Vector::Z;   //!< End vertices of the segment.
public:
  //! Empty.
  Segment() {}
  explicit Segment(const Vector&, const Vector&);
  explicit Segment(const Segment2&, const double&, const double&);

  //! Empty.
  ~Segment() {}

  void Rotate(const Matrix&);
  void Translate(const Vector&);
  void Scale(const double&);

  Segment Translated(const Vector&) const;
  Segment Scaled(const Vector&) const;
  Segment Scaled(const double&) const;
  Segment Rotated(const Matrix&) const;
  Segment InverseTransformed(const Frame&) const;

  Quadric Equation(const Ray&) const;

  // End vertices
  Vector Vertex(int) const;
  Vector& Vertex(int);

  // Vertex along segment
  Vector VertexAt(const double&) const;
  Vector Center() const;

  bool Intersect(const Box&) const;

  Vector GetAxis() const;
  Box GetBox() const;
  double Length() const;
  friend std::ostream& operator<<(std::ostream&, const Segment&);

  // Distance
  double R(const Vector&) const;
  double R(const Vector&, double&) const;
  double R(const Segment&) const;
  Vector Normal(const Vector&) const;

  bool Equal(const Segment&, const double&) const;
  static Quadric EdgeEquation(const Ray&, const Vector&, const Vector&, const Vector&);
  static Vector Intersect(const Vector&, const Vector&, const double, const double, double = 0.0);
};

/*!
\brief Creates a segment given its end vertices.
\param a, b End vertices of the segment.
*/
inline Segment::Segment(const Vector& a, const Vector& b) :a(a), b(b)
{
}

//! Return one of the end vertex of the axis.
inline Vector Segment::Vertex(int i) const
{
  if (i == 0) return a;
  else return b;
}

//! Return one of the end vertex of the axis.
inline Vector& Segment::Vertex(int i)
{
  if (i == 0) return a;
  else return b;
}

//! Return axis length.
inline double Segment::Length() const
{
  return Norm(b - a);
}

//! Returns the normalized axis vector. 
inline Vector Segment::GetAxis() const
{
  return (b - a) / Norm(b - a);
}

/*!
\brief Compute a point on the segment.
\param t Interpolation parameter.
*/
inline Vector Segment::VertexAt(const double& t) const
{
  return (1.0 - t) * a + t * b;
}

/*!
\brief Compute the center of the segment.
*/
inline Vector Segment::Center() const
{
  return 0.5 * (a + b);
}

// Segment class
class Segment2
{
protected:
  Vector2 a = Vector2::Null, b = Vector2::X; //!< End vertices of the segment.
public:
  //! Empty.
  Segment2() {}
  explicit Segment2(const Vector2&, const Vector2&);
  explicit Segment2(const Segment&);
  //! Empty.
  ~Segment2() {}

  // End vertices
  Vector2 Vertex(int) const;
  Vector2& Vertex(int);

  // Vertex along segment
  Vector2 VertexAt(const double&) const;
  Vector2 Center() const;

  Vector2 GetAxis() const;
  Box2 GetBox() const;

  Vector2 Orthogonal() const;

  Segment ToSegment(const double&) const;
  Segment ToSegment(const double&, const double&) const;

  void Translate(const Vector2&);
  void Rotate(const Matrix2&);
  void Scale(const double&);

  Segment2 Translated(const Vector2&) const;
  Segment2 Extended(const double&) const;

  double Length() const;
  friend std::ostream& operator<<(std::ostream&, const Segment2&);

  friend bool operator==(const Segment2&, const Segment2&);

  void Draw(QGraphicsScene&, const QPen & = QPen()) const;
  void DrawArrow(QGraphicsScene&, const double&, const QPen & = QPen(), const QBrush & = QBrush()) const;

  double R(const Vector2&) const;
  double R(const Vector2&, double&) const;
  Vector2 Normal(const Vector2&) const;

  bool Intersect(const Segment2&) const;
  bool IntersectOpen(const Segment2&) const;
  bool Intersection(const Segment2&, Vector2&) const;

  static constexpr const double epsilon = 1e-8; //!< Epsilon value for intersection test
  static Vector2 Intersect(const Vector2&, const Vector2&, const double, const double, const double = 0.0);
};

/*!
\brief Creates a planar segment given end vertices.
\param a, b End vertices of the segment.
*/
inline Segment2::Segment2(const Vector2& a, const Vector2& b) :a(a), b(b)
{
}

/*!
\brief Creates a planar segment given a three dimensional segment.
\param s %Segment.
*/
inline Segment2::Segment2(const Segment& s) :a(s.Vertex(0)), b(s.Vertex(1))
{
}

//! Return one of the end vertex of the axis.
inline Vector2 Segment2::Vertex(int i) const
{
  if (i == 0) return a;
  else return b;
}

//! Return one of the end vertex of the axis.
inline Vector2& Segment2::Vertex(int i)
{
  if (i == 0) return a;
  else return b;
}

//! Return axis length.
inline double Segment2::Length() const
{
  return Norm(b - a);
}

//! Returns the normalized axis vector. 
inline Vector2 Segment2::GetAxis() const
{
  return (b - a) / Norm(b - a);
}

/*!
\brief Compute a point on the segment.
\param t Interpolation parameter.
*/
inline Vector2 Segment2::VertexAt(const double& t) const
{
  return (1.0 - t) * a + t * b;
}

/*!
\brief Compute the center of the segment.
*/
inline Vector2 Segment2::Center() const
{
  return 0.5 * (a + b);
}

/*!
\brief Compare two segments.
\param a,b Segments.
*/
inline bool operator==(const Segment2& a, const Segment2& b)
{
  return (a.a == b.a) && (a.b == b.b);
}

class Line :protected Segment {
protected:
public:
  //! Empty
  Line() {}
  explicit Line(const Vector&, const Vector&);

  Vector Vertex(int) const;
  Vector& Vertex(int);

  double R(const Vector&) const;
  double R(const Line&) const;
};

/*!
\brief Creates a line.
\param a,b Two points on the line.
*/
inline Line::Line(const Vector& a, const Vector& b) :Segment(a, b)
{
}

//! Return one of the vertices of the line.
inline Vector Line::Vertex(int i) const
{
  if (i == 0) return a;
  else return b;
}

//! Return one of the vertices of the line.
inline Vector& Line::Vertex(int i)
{
  if (i == 0) return a;
  else return b;
}

class Circle2;

class Line2 :public Segment2
{
protected:
public:
  //! Empty
  Line2() {}
  explicit Line2(const Vector2&, const Vector2&);
  explicit Line2(const Segment2&);

  // Vertices
  Vector2 Vertex(int) const;
  Vector2& Vertex(int);

  // Inbtersection
  bool Intersection(const Line2&, Vector2&) const;
  bool Intersection(const Segment2&, Vector2&) const;

  Vector2 Symmetry(const Vector2&) const;
  Box2 Symmetry(const Box2&) const;
  Circle2 Symmetry(const Circle2&) const;

  // Distance
  double R(const Vector2&) const;

  bool IsLeftOrOn(const Vector2&, const double& = 0.0) const;
  bool IsRightOrOn(const Vector2&, const double& = 0.0) const;
public:
  static const Line2 X; //!< Horizontal.
  static const Line2 Y; //!< Vertical.
};

//! Return one of the vertices of the line.
inline Vector2 Line2::Vertex(int i) const
{
  if (i == 0) return a;
  else return b;
}

//! Return one of the vertices of the line.
inline Vector2& Line2::Vertex(int i)
{
  if (i == 0) return a;
  else return b;
}

class Polygon2;
class Polygons2;

class SegmentSet2
{
protected:
  QVector<Vector2> vertices; //!< Vertices.
  QVector<int> indices; //!< %Segment vertex indices.
public:
  explicit SegmentSet2();
  explicit SegmentSet2(const QVector<Vector2>&, const QVector<int>&);

  QVector<Segment2> GetSegments() const;
  Segment2 GetSegment(int) const;
  Vector2 Vertex(int) const;

  bool IsEmpty() const;
  int SegmentSize() const;
  Polygons2 GetPolygons() const;

  void Draw(QGraphicsScene&, const QPen & = QPen()) const;
};

/*!
\brief Return the k-th vertex.
\param k Index.
*/
inline Vector2 SegmentSet2::Vertex(int k) const
{
  return vertices.at(k);
}

/*!
\brief Return the k-th segment.
\param k Index.
*/
inline Segment2 SegmentSet2::GetSegment(int k) const
{
  return Segment2(vertices.at(indices.at(2 * k)), vertices.at(indices.at(2 * k + 1)));
}

/*!
\brief Check is the set of segments is empty.
*/
inline bool SegmentSet2::IsEmpty() const
{
  return vertices.size() == 0;
}

/*!
\brief Return the number of segments.
*/
inline int SegmentSet2::SegmentSize() const
{
  return indices.size() / 2;
}

