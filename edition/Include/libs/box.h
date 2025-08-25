// Box
#pragma once

// Forward declare seldom used pre-processing intensive Qt class
#include <QtWidgets/QGraphicsScene>

#include <QtGui/QPen>
#include <QtGui/QBrush>
#include <QtCore/QRect>
#include <QtCore/QSize>

#include "libs/random.h"
#include "libs/frame.h"
#include "libs/axis.h"

class Segment;
class Plane;
class Axis;

class Box
{
protected:
  Vector a = Vector(0.0), b = Vector(1.0); //!< Lower and upper vertexex.
public:
  //! Empty.
  Box() {}
  explicit Box(const double&);
  explicit Box(const Vector&);
  explicit Box(const Vector&, const Vector&);
  explicit Box(const Vector&, const double&);
  explicit Box(const Vector&, const double&, const double&, const double&);
  explicit Box(const double&, const double&, const double&);
  explicit Box(const Vector*, int);
  explicit Box(const QVector<Vector>&);
  explicit Box(const Box&, const Box&);
  explicit Box(const Box&, const Frame&);
  explicit Box(const Box&, const FrameScaled&);
  explicit Box(const Box&, const Matrix&);

  //! Empty.
  ~Box() {}

  // Access vertexes
  Vector& operator[] (int);
  Vector operator[] (int) const;

  // Comparison
  friend bool operator==(const Box&, const Box&);
  friend bool operator!=(const Box&, const Box&);

  // Acces to vertexes
  Vector Center() const;
  Vector Vertex(int) const;
  Vector Vertex(int, int, int, int, int, int) const;
  Segment Edge(int) const;
  Plane Face(int) const;

  Vector Size() const;
  Vector Diagonal() const;
  double Radius() const;

  double R(const Vector&) const;
  double Signed(const Vector&) const;
  Vector Normal(const Vector&) const;

  // Distance
  double R(const Box&) const;
  double RInfinity(const Box&) const;

  // Intersection with a ray 
  int Intersect(const Ray&) const;
  int Intersect(const Ray&, double&, double&) const;
  int Intersect(const Ray&, double&) const;

  int Intersect(const Ray&, double&, double&, Vector&, Vector&) const;
  int Intersect(const Ray&, double&, Vector&) const;
  int Intersect(const Vector&, const Vector&) const;

  // Intersection with another box
  bool Intersect(const Box&) const;

  bool Inside(const Box&) const;
  bool Inside(const Vector&) const;

  // Boolean
  Box Intersection(const Box&) const;
  int Difference(const Box&, Box*) const;

  bool Empty() const;
  double Volume() const;
  double Area() const;

  // Extend box to cube
  void SetCubic();
  void SetInscribedCubic();

  Box Cube() const;

  void Extend(const double&);
  Box Extended(const double&) const;
  void Extend(const Vector&);
  Box Extended(const Vector&) const;

  void SetParallelepipedic(const double&, int&, int&, int&);
  void SetParallelepipedic(int, int&, int&, int&);

  // Compute sub-box
  Box Sub(int) const;

  // Classification
  int Octant(const Vector&) const;

  // Minkowski sum
  friend Box operator+(const Box&, const Box&);
  friend Box operator*(const double&, const Box&);
  friend Box operator*(const Box&, const double&);

  // Translation, scale
  void Translate(const Vector&);
  Box Translated(const Vector&) const;
  void Scale(const double&);
  void Scale(const Vector&);
  Box Scaled(const Vector&) const;
  Box Offsetted(const Vector&) const;

  Box Centered() const;
  Box Cut(const Vector&, const Vector&) const;

  friend std::ostream& operator<<(std::ostream&, const Box&);

  int IntegerAxis() const;

  // Random
  Vector RandomInside(Random & = Random::R239) const;
  Vector RandomSurface(Random & = Random::R239) const;

  QVector<Vector> Poisson(const double&, int, Random & = Random::R239) const;
public:
  static constexpr const double epsilon = 1.0e-5; //!< Epsilon value used to check intersections and some round off errors and for ray intersection tests.
  static const Box Infinity; //!< Largest box.
  static const Box Null; //!< Null box, equivalent to: \code Box(Vector(0.0)); \endcode 
  static const Box Unit; //!< Unit box.
  static const int edge[24]; //!< Edge vertex indexes.
  static const Vector normal[6]; //!< Face normals.
};

//! Returns either end vertex of the box.
inline Vector& Box::operator[] (int i)
{
  if (i == 0) return a;
  else return b;
}

//! Overloaded.
inline Vector Box::operator[] (int i) const
{
  if (i == 0) return a;
  else return b;
}

/*!
\brief Computes the Minkowski sum of two boxes.
\param a, b Argument boxes.
*/
inline Box operator+(const Box& a, const Box& b)
{
  return Box(a.a + b.a, a.b + b.b);
}

//! Scales a box by a scalar factor.
inline Box operator*(const double& t, const Box& box)
{
  if (t < 0.0)
    return Box(t * box.b, t * box.a);
  else
    return Box(t * box.a, t * box.b);
}

//! Overloaded.
inline Box operator*(const Box& box, const double& t)
{
  if (t < 0.0)
    return Box(t * box.b, t * box.a);
  else
    return Box(t * box.a, t * box.b);
}

//! Returns the center of the box.
inline Vector Box::Center() const
{
  return 0.5 * (a + b);
}

/*!
\brief Returns the diagonal of the box.
*/
inline Vector Box::Diagonal() const
{
  return (b - a);
}

/*!
\brief Compute the size (width, length and height) of a box.
\sa Box::Diagonal()
*/
inline Vector Box::Size() const
{
  return b - a;
}

/*!
\brief Returns the radius of the box, i.e. the length of the half diagonal of the box.
*/
inline double Box::Radius() const
{
  return 0.5 * Norm(b - a);
}

/*!
\brief Returns the k-th vertex of the box.

\image html box-indexes.png

The returned vector is computed by analysing the first three bits of k as follows:
\code
Vector vertex=Vector((k&1)?b[0]:a[0],(k&2)?b[1]:a[1],(k&4)?b[2]:a[2]);
\endcode
*/
inline Vector Box::Vertex(int k) const
{
  return Vector((k & 1) ? b[0] : a[0], (k & 2) ? b[1] : a[1], (k & 4) ? b[2] : a[2]);
}

//! Compute the volume of a box.
inline double Box::Volume() const
{
  Vector side = b - a;
  return side[0] * side[1] * side[2];
}

/*!
\brief Compute the surface area of a box.
*/
inline double Box::Area() const
{
  Vector side = b - a;
  return 2.0 * (side[0] * side[1] + side[0] * side[2] + side[1] * side[2]);
}

/*!
\brief Computes the squared Euclidean distance between the box and a point.
\param p Point.
*/
inline double Box::R(const Vector& p) const
{
  double r = 0.0;

  for (int i = 0; i < 3; i++)
  {
    if (p[i] < a[i])
    {
      double s = p[i] - a[i];
      r += s * s;
    }
    else if (p[i] > b[i])
    {
      double s = p[i] - b[i];
      r += s * s;
    }
  }
  return r;
}

/*!
\brief Computes the normal vector between a point and a box.

Let <b>q</b> the projection of <b>p</b> onto the box, the normal vector is defined as <b>n</b>=<b>p</b>-<b>q</b>.
\param p Point.
*/
inline Vector Box::Normal(const Vector& p) const
{
  Vector n;

  for (int i = 0; i < 3; i++)
  {
    if (p[i] < a[i])
    {
      n[i] = p[i] - a[i];
    }
    else if (p[i] > b[i])
    {
      n[i] = p[i] - b[i];
    }
    else
    {
      n[i] = 0.0;
    }
  }
  return n;
}

/*!
\brief Check if the box intersects another box.

\param box Argument box.
*/
inline bool Box::Intersect(const Box& box) const
{
  if (((a[0] >= box.b[0]) || (a[1] >= box.b[1]) || (a[2] >= box.b[2]) || (b[0] <= box.a[0]) || (b[1] <= box.a[1]) || (b[2] <= box.a[2])))
    return false;
  else
    return true;
}

/*!
\brief Check if a box is empty.

The box is empty if one of the coordinates of its lower point is
greater than the coordinates of the opposite point.
*/
inline bool Box::Empty() const
{
  if ((a[0] > b[0]) || (a[1] > b[1]) || (a[2] > b[2]))
    return true;
  else
    return false;
}

/*!
\brief Check if an argument box is inside the box.
\param box The box.
*/
inline bool Box::Inside(const Box& box) const
{
  return ((a < box.a) && (b > box.b));
}

/*!
\brief Check if a point is inside the box.
\param p Point.
*/
inline bool Box::Inside(const Vector& p) const
{
  return ((a < p) && (b > p));
}

/*!
\brief Check if two boxes are (strictly) equal.
\param a, b Boxes.
*/
inline bool operator==(const Box& a, const Box& b)
{
  return (a.a == b.a) && (a.b == b.b);
}

/*!
\brief Check if two boxes are (strictly) different.
\param a, b Boxes.
*/
inline bool operator!=(const Box& a, const Box& b)
{
  return !(a == b);
}

class Line2;
class Segment2;

class Box2
{
protected:
  Vector2 a = Vector2(0.0), b = Vector2(1.0); //!< Lower and upper vertexes of the box.
public:
  //! Empty
  Box2() {}
  explicit Box2(const double&);
  explicit Box2(const double&, const double&);
  explicit Box2(const Vector2&);
  explicit Box2(const Vector2&, const Vector2&);
  explicit Box2(const Vector2&, const double&);
  explicit Box2(const Vector2&, const double&, const double&);
  explicit Box2(const Box&);
  explicit Box2(const Box2&, const Box2&);
  explicit Box2(const Box2&, const Matrix2&);
  explicit Box2(const Box2&, const Frame2&);
  explicit Box2(const QVector<Vector2>&);
  explicit Box2(const QSize&);

  Vector2& operator[] (int);
  Vector2 operator[] (int) const;

  Vector2 Size() const;
  double Width() const;
  double Height() const;
  Vector2 Diagonal() const;
  double Radius() const;

  double Area() const;
  double Perimeter() const;

  // Access to vertexes
  Vector2 Center() const;
  Vector2 Vertex(int) const;
  Vector2 Vertex(int, int, int, int) const;

  Vector2 VertexUV(const Vector2&) const;

  Box2 Sub(int) const;

  double R(const Vector2&) const;
  double Signed(const Vector2&) const;
  double R(const Box2&) const;

  void Translate(const Vector2&);
  void Scale(const Vector2&);
  void Scale(const double&);

  Box2 Translated(const Vector2&) const;
  Box2 Scaled(const double&) const;
  Box2 Scaled(const Vector2&) const;
  Box2 Scaled(const QSize&) const;
  Box2 ScaledTo(const double&) const;

  Box2 Rotated(const double&) const;
  Box2 Rotated(const Matrix2&) const;

  Box2 Centered() const;

  // Change shape
  void SetCubic();
  void SetInscribedCubic();

  Box2 Cube() const;

  void Extend(const double&);
  void Extend(const Vector2&);
  Box2 Extended(const double&) const;

  void SetParallelepipedic(int, int&, int&);
  void SetParallelepipedic(const double&, int&, int&);

  // Inside
  bool Inside(const Vector2&) const;
  bool Inside(const Vector2&, const double&) const;

  // Intersection with another box
  bool Intersect(const Box2&) const;
  double OverlapArea(const Box2&) const;

  // Boolean
  Box2 Intersection(const Box2&) const;

  Box ToBox(const double&, const double&) const;

  // Comparison
  friend bool operator==(const Box2&, const Box2&);
  friend bool operator!=(const Box2&, const Box2&);

  // Classification
  int Quadrant(const Vector2&) const;

  // Minkowski sum
  friend Box2 operator+(const Box2&, const Box2&);
  friend Box2 operator*(const double&, const Box2&);
  friend Box2 operator*(const Box2&, const double&);

  bool Intersect(const Ray2&, double&, double&) const;
  bool Intersect(const Ray2&) const;

  bool Intersect(const Line2&, double&, double&) const;
  bool Intersect(const Segment2&, double&, double&) const;
  bool Intersect(const Line2&) const;
  bool Intersect(const Segment2&) const;

  friend std::ostream& operator<<(std::ostream&, const Box2&);
  void Draw(QGraphicsScene&, const QPen & = QPen(), const QBrush & = QBrush()) const;

  // Random
  Vector2 RandomInside(Random & = Random::R239) const;
  Vector2 RandomOn(Random & = Random::R239) const;

  QRectF GetQtRect() const;

  QRect TileRange(const Box2&) const;
  Box2 Tile(int, int) const;
  Box2 Tile(const QRect&) const;

  Segment2 GetSegment(const double&, bool = true) const;

  QVector<Vector2> Poisson(const double&, int, Random & = Random::R239) const;
  QVector<Vector2> Poisson(const double&, int, const QVector<Vector2>&, bool = true, Random & = Random::R239) const;

  static Box2 MinMax(const Vector2&, const Vector2&);
public:
  static const double epsilon;
  static const Box2 Infinity;
  static const Box2 Null;
  static const Box2 Unit;
};

/*!
\brief Create a box.

Note that is possible to create a box using Vector as parameters
as the compiler will call the constructor Vector2::Vector2(const Vector&).

\param a,b Points.
*/
inline Box2::Box2(const Vector2& a, const Vector2& b)
{
  Box2::a = a;
  Box2::b = b;
}

/*!
\brief Create a box in the plane given a box.
\param box The box.
*/
inline Box2::Box2(const Box& box)
{
  a = box[0];
  b = box[1];
}

/*!
\brief Creates a box.

\param c Center.
\param r Radius.
*/
inline Box2::Box2(const Vector2& c, const double& r)
{
  a = c - Vector2(r);
  b = c + Vector2(r);
}

/*!
\brief Create a box given a center point and its width, length, and height.

\param c center.
\param x,y Width and height.
*/
inline Box2::Box2(const Vector2& c, const double& x, const double& y)
{
  Vector2 r = 0.5 * Vector2(x, y);
  a = c - r;
  b = c + r;
}

/*!
\brief Create a square box centered at the origin and of given half side length.

This is equivalent to:
\code
Box2 box(Vector2(0.0),2.0);  // Simplified constructor Box2(2.0);
\endcode
\param r Half side length.
*/
inline Box2::Box2(const double& r)
{
  a = -Vector2(r);
  b = Vector2(r);
}
/*!
\brief Create a box centered at the origin and of given dimensions.

This is equivalent to:
\code
double width, height;
Box2 box(Vector2(-width/2.0,-height/2.0),Vector2(width/2.0,height/2.0);  // Simplified constructor Box2(width,height);
\endcode
\param x,y %Size.
*/
inline Box2::Box2(const double& x, const double& y)
{
  a = -Vector2(x / 2.0, y / 2.0);
  b = Vector2(x / 2.0, y / 2.0);
}

/*!
\brief Create a box from a single vertex.
\param p Point.
*/
inline Box2::Box2(const Vector2& p) :a(p), b(p)
{
}

/*!
\brief Create a box from any two points.
\param a,b Points.
*/
inline Box2 Box2::MinMax(const Vector2& a, const Vector2& b)
{
  return Box2(Vector2::Min(a, b), Vector2::Max(a, b));
}

//! Returns either end vertex of the box.
inline Vector2& Box2::operator[] (int i)
{
  if (i == 0) return a;
  else return b;
}

//! Overloaded.
inline Vector2 Box2::operator[] (int i) const
{
  if (i == 0) return a;
  else return b;
}

/*!
\brief Compute the size (width and height) of a box.
*/
inline Vector2 Box2::Size() const
{
  return b - a;
}

/*!
\brief Compute the width of a box.

\sa Box2::Size()
*/
inline double Box2::Width() const
{
  return b[0] - a[0];
}

/*!
\brief Compute the height of a box.

\sa Box2::Size()
*/
inline double Box2::Height() const
{
  return b[1] - a[1];
}

/*!
\brief Returns the radius of the box, i.e. the length of the half diagonal of the box.
*/
inline double Box2::Radius() const
{
  return 0.5 * Norm(b - a);
}

/*!
\brief Compute the perimeter of the box.
*/
inline double Box2::Perimeter() const
{
  return 2.0 * (b[0] - a[0] + b[1] - a[1]);
}

/*!
\brief Computes the Minkowski sum of two boxes.
\param a, b Argument boxes.
*/
inline Box2 operator+(const Box2& a, const Box2& b)
{
  return Box2(a.a + b.a, a.b + b.b);
}

/*!
\brief Scales a box by a scalar factor.
\param t Scaling factor.
\param box The box.
*/
inline Box2 operator*(const double& t, const Box2& box)
{
  if (t < 0.0)
    return Box2(t * box.b, t * box.a);
  else
    return Box2(t * box.a, t * box.b);
}

//! Overloaded.
inline Box2 operator*(const Box2& box, const double& t)
{
  if (t < 0.0)
    return Box2(t * box.b, t * box.a);
  else
    return Box2(t * box.a, t * box.b);
}

/*!
\brief Returns the diagonal of the box.
*/
inline Vector2 Box2::Diagonal() const
{
  return (b - a);
}

/*!
\brief Returns the center of the box.
*/
inline Vector2 Box2::Center() const
{
  return 0.5 * (a + b);
}

//! Compute the surface area of a box.
inline double Box2::Area() const
{
  return Width() * Height();
}

/*!
\brief Returns the k-th vertex of the box.

The returned vector is computed by analysing the first two bits of k as follows:
\code
Vector2 vertex=Vector2((k&1)?b[0]:a[0],(k&2)?b[1]:a[1]);
\endcode
*/
inline Vector2 Box2::Vertex(int k) const
{
  return Vector2((k & 1) ? b[0] : a[0], (k & 2) ? b[1] : a[1]);
}

/*!
\brief Check if the box intersects another box.

\param box Argument box.
*/
inline bool Box2::Intersect(const Box2& box) const
{
  if (((a[0] >= box.b[0]) || (a[1] >= box.b[1]) || (b[0] <= box.a[0]) || (b[1] <= box.a[1])))
    return false;
  else
    return true;
}

/*!
\brief Computes the overlapping area of two boxes.

This is a convenience function for the following code:

\code
Box2 a,b;
double area=(a.Intersection(b)).Area();
\endcode
\sa Box2::Intersection(const Box2&) const

\param box
*/
inline double Box2::OverlapArea(const Box2& box) const
{
  if (!Intersect(box))
    return 0.0;
  return Box2(Vector2::Max(a, box[0]), Vector2::Min(b, box[1])).Area();
}

/*!
\brief Check if two boxes are (strictly) equal.
\param a, b Boxes.
*/
inline bool operator==(const Box2& a, const Box2& b)
{
  return (a.a == b.a) && (a.b == b.b);
}

/*!
\brief Check if two boxes are (strictly) different.
\param a, b Boxes.
*/
inline bool operator!=(const Box2& a, const Box2& b)
{
  return !(a == b);
}

/*!
\brief Compute the coordinates of a point inside the box given uv coordinates.
\param uv Coordinates.
*/
inline Vector2 Box2::VertexUV(const Vector2& uv) const
{
  return Vector2((1.0 - uv[0]) * a[0] + uv[0] * b[0], (1.0 - uv[1]) * a[1] + uv[1] * b[1]);
}


class BoxEmpty :public Box
{
public:
  //! Empty.
  BoxEmpty() {}
  explicit BoxEmpty(const Vector&, const Vector&);
  explicit BoxEmpty(const Vector&, const double&);
  explicit BoxEmpty(const Box&);

  //! Empty.
  ~BoxEmpty() {}

  // Distance
  double R(const Vector&) const;
  Vector Normal(const Vector&) const;
  double Signed(const Vector&) const;
};
