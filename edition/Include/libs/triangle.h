// Triangle

#pragma once

#include <QtCore/QVector>

#include "libs/frame.h"
#include "libs/sphere.h"

// Forward declare classes: Plane, Segment
class Plane;
class Segment;
class Segment2;

// Triangle
class Triangle
{
protected:
  Vector p[3]; //!< Array of vertices.
public:
  //! Empty.
  Triangle() {}
  explicit Triangle(const Vector&, const Vector&, const Vector&);

  //! Empty.
  ~Triangle() {}

  Vector operator[] (int) const;

  // Point in triangle
  Vector Vertex(const double&, const double&) const;

  // Intersection
  bool Intersect(const Ray&, double&) const;
  bool Intersect(const Ray&, double&, double&, double&) const;
  bool Intersect(const Sphere&) const;

  void Translate(const Vector&);
  void Rotate(const Matrix&);
  void Scale(const Vector&);
  void Transform(const Frame&);
  Triangle Transformed(const FrameScaled&) const;
  void Shrink(const double&);

  // Geometry
  Vector Normal() const;
  Vector AreaNormal() const;
  Vector Center() const;

  // Distance
  double R(const Vector&) const;
  Vector Normal(const Vector&) const;

  // Spheres
  Sphere Inscribed() const;
  Sphere Circumscribed() const;

  double Area() const;
  double Aspect() const;
  Box GetBox() const;
  int Intersect(const Triangle&, int&, Vector&, Vector&) const;

  int Intersect(const Plane&, Segment&) const;
  int Intersect(const double&, Segment&) const;

  // Stream
  friend std::ostream& operator<<(std::ostream&, const Triangle&);

  Vector BarycentricCoordinates(const Vector&) const;

  void Subdivide(int, QVector<Vector>&, QVector<int>&) const;
  double InscribedRadius() const;
  double CircumscribedRadius() const;

  Vector RandomInside(Random & = Random::R239) const;
protected:
  static constexpr const double epsilon = 1.0e-7; //!< Internal epsilon constant.
};

/*!
\brief Return the i-th vertex.
\param i Index.
*/
inline Vector Triangle::operator[] (int i) const
{
  return p[i];
}

//! Compute the barycenter of the triangle.
inline Vector Triangle::Center() const
{
  return (p[0] + p[1] + p[2]) / 3.0;
}

//! Compute the area of the triangle.
inline double Triangle::Area() const
{
  return 0.5 * Norm((p[0] - p[1]) / (p[2] - p[0]));
}

/*!
\brief Create a triangle.
\param a,b,c Vertices of the triangle.
*/
inline Triangle::Triangle(const Vector& a, const Vector& b, const Vector& c) :p{ a,b,c }
{
}

// Forward declare classes: Circle2
class Circle2;

class Triangle2 {
protected:
  Vector2 p[3]; //!< Array of vertices.
public:
  //! Empty.
  Triangle2() {}
  explicit Triangle2(const Vector2&, const Vector2&, const Vector2&);
  explicit Triangle2(const Triangle&);

  //! Empty.
  ~Triangle2() {}

  Vector2 operator[] (int) const;

  // Geometry
  Vector2 Center() const;
  Vector2 BaryCenter(const Vector&) const;
  Vector2 OrthoCenter() const;

  // Spheres
  Circle2 Inscribed() const;
  Circle2 Circumscribed() const;

  double R(const Vector2&) const;
  Vector2 Normal(const Vector2&) const;
  bool Intersect(const Circle2&) const;
  bool Intersect(const Box2&) const;
  bool Intersect(const Segment2&) const;
  bool Intersect(const Ray2&, double&, double&) const;

  double Area() const;
  double SignedArea() const;
  double Perimeter() const;

  Box2 GetBox() const;

  QVector<Vector2> Poisson(const double&, int, Random & = Random::R239) const;

  Vector BarycentricCoordinates(const Vector2&) const;


  double InscribedRadius() const;
  double CircumscribedRadius() const;
  double Aspect() const;

  bool Inside(const Box2&) const;
  bool Inside(const Circle2&) const;
  bool Inside(const Vector2&) const;
  bool Inside(const Vector2&, double&, double&) const;

  Vector2 RandomInside(Random & = Random::R239) const;

  void Rotate(const Matrix2&);
  void Translate(const Vector2&);
  void Scale(const double&);
  void Scale(const Vector2&);
  void Shrink(const double&);

  bool Overlap(const Triangle2&) const;
  bool Intersect(const Triangle2&) const;

  // Stream
  friend std::ostream& operator<<(std::ostream&, const Triangle2&);
  void Draw(QGraphicsScene&, const QPen & = QPen(), const QBrush & = QBrush()) const;

  // Equilateral
  static Triangle2 Equilateral(const Vector2&, const double& = 1.0, const double& = 0.0);
protected:
  bool IntersectEdge(const Triangle2&) const;
  bool IntersectVertex(const Triangle2&) const;
protected:
  static double epsilon; //!< Internal epsilon constant.
};

/*!
\brief Create a triangle given three points.
\param a, b, c Vertices.
*/
inline Triangle2::Triangle2(const Vector2& a, const Vector2& b, const Vector2& c) :p{ a,b,c }
{
}

/*!
\brief Return the i-th vertex.
*/
inline Vector2 Triangle2::operator[] (int i) const
{
  return p[i];
}

/*!
\brief Compute the circle circumscribing the triangle.

This function directly relies on the constructor of the class.

\sa Circle2::Circle2(const Vector2&, const Vector2&, const Vector2&);
*/
inline Circle2 Triangle2::Circumscribed() const
{
  return Circle2(p[0], p[1], p[2]);
}

