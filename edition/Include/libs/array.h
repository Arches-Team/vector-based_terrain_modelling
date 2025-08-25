// Fundamentals
#pragma once

#include "libs/box.h"
#include "libs/circle.h"

class Triangle2;

class Array2 :protected Box2
{
protected:
  int nx = 0, ny = 0; //!< Sizes.
  Vector2 celldiagonal = Vector2::Null; //!< Cell diagonal.
  Vector2 inversecelldiagonal = Vector2::Null; //!< Inverse cell diagonal.
public:
  Array2();
  explicit Array2(const Box2&, int, int);
  explicit Array2(const Box2&, int);

  Array2 DownSample(int) const;
  Array2 UpSample(int) const;

  using Box2::VertexUV;

  //! Empty
  ~Array2() {}

  int VertexSize() const;
  int GetSizeX() const;
  int GetSizeY() const;

  int VertexBorderSize() const;

  int CellSize() const;
  int CellSizeX() const;
  int CellSizeY() const;
  QSize GetQtSize() const;

  Box2 Cell(int) const;
  Box2 Cell(int, int) const;
  void HalfCell(int, int, bool, Triangle2&, Triangle2&) const;

  Vector2 CellCenter(int, int) const;

  void Scale(const Vector2&);
  void Scale(const double&);

  Box2 GetBox() const;

  Box2 UnitCell() const;

  Vector2 CellDiagonal() const;
  double CellArea() const;

  Vector2 ArrayVertex(int, int) const;
  Vector2 ArrayVertex(const QPoint&) const;
  QVector <Vector2> ArrayVertexes(const QVector<QPoint>&) const;

  Vector2 Size() const;

  QRect AreaInteger() const;

  QImage ImageGrid(const Box2&, const double&) const;
  void Draw(QGraphicsScene&, const QPen & = QPen(), const QBrush & = QBrush()) const;

  // Subdivision
  void Subdivide();

  // Empty
  bool IsEmpty() const;

  QPoint VertexBorder(int) const;
  int VertexBorderIndex(int, int) const;

  // Vertex queries
  void VertexInteger(const Vector2&, int&, int&) const;
  QPoint VertexInteger(const Vector2&) const;
  QPoint VertexInteger(const Vector2&, double&, double&) const;
  QRect VertexIntegerArea(const Box2&) const;
  QRect VertexIntegerArea(const Circle2&) const;
  QRect VertexIntegerArea(const QRect&) const;
  QRect VertexIntegerArea(int,int,int) const;

  // Cell queries
  void CellInteger(const Vector2&, int&, int&) const;
  void CellInteger(const Vector2&, int&, int&, double&, double&) const;
  QPoint CellInteger(const Vector2&) const;
  QPoint CellInteger(const Vector2&, double&, double&) const;

  QRect CellIntegerArea(const Box2&) const;
  QRect CellIntegerArea(const Circle2&) const;

  friend std::ostream& operator<<(std::ostream&, const Array2&);

  void Translate(const Vector2&);
  Vector2 Center() const;

  // Domain queries
  constexpr bool InsideVertexIndex(int, int) const;
  constexpr bool OutsideVertexIndex(int, int) const;
  constexpr bool InsideVertexIndex(const QPoint&) const;
  constexpr bool InsideVertexIndex(int, int, int) const;

  bool Inside(const Vector2&) const;

  constexpr bool BorderVertexIndex(int, int) const;
  constexpr bool BorderVertexIndex(const QPoint&) const;

  // Indexes for storing elements at vertices
  constexpr int VertexIndex(int, int) const;
  constexpr int VertexIndex(const QPoint&) const;

  QPoint Next(const QPoint&, int) const;

  QImage CreateEmptyImage() const;

  QString Statistics() const;

  void ClampVertexIndex(int&, int&) const;
protected:
  constexpr void InverseVertexIndex(int, int&, int&) const;
  constexpr QPoint InverseVertexIndex(int) const;

  // Indexes for storing elements at cells
  constexpr int CellIndex(int, int) const;
  constexpr int CellIndex(const QPoint&) const;

  constexpr bool InsideCellIndex(int, int) const;
  constexpr bool InsideCellIndex(const QPoint&) const;
  void InverseCellIndex(int, int&, int&) const;
public:
  static int NeighborCode(int, int);
  static QPoint CodeToDir(int);
protected:
  static const QPoint next[8]; //!< Array of points in the 1-ring neighborhood.
  static const double length[8]; //!< Length to the i-th neighbor.
  static const double inverselength[8]; //!< Inverse length.
};

/*!
\brief Downsample the array.
\param d Sampling factor.
*/
inline Array2 Array2::DownSample(int d) const
{
  return Array2(Box2(a, b), nx / d, ny / d);
}

/*!
\brief Upsample the array.
\param u Sampling factor.
*/
inline Array2 Array2::UpSample(int u) const
{
  return Array2(Box2(a, b), nx * u, ny * u);
}

/*!
\brief Return the size of the array.
*/
inline Vector2 Array2::Size() const
{
  return Box2::Size();
}

/*!
\brief Detect if the array is empty, i.e., any dimension equal to zero.
*/
inline bool Array2::IsEmpty() const
{
  return (nx <= 0) || (ny <= 0);
}

/*!
\brief Get the vertex size of the array for x axis.
*/
inline int Array2::GetSizeX() const
{
  return nx;
}

/*!
\brief Get the vertex size of the array for y axis.
*/
inline int Array2::GetSizeY() const
{
  return ny;
}

/*!
\brief Create an image with the same size as the array.

This is a convenience function for:
\code
Array2 a;
  QImage image(a.GetSizeX(), a.GetSizeY(), QImage::Format_ARGB32);

\endcode
*/
inline QImage Array2::CreateEmptyImage() const
{
  return QImage(nx, ny, QImage::Format_ARGB32);
}

/*!
\brief Get the cell size of the array for x axis.
*/
inline int Array2::CellSizeX() const
{
  return nx - 1;
}

/*!
\brief Get the cell size of the array for y axis.
*/
inline int Array2::CellSizeY() const
{
  return ny - 1;
}

/*!
\brief Return the size of the vertex array.
*/
inline int Array2::VertexSize() const
{
  return nx * ny;
}
/*!
\brief Return the number of vertices on the boundary of the rectangle.
*/
inline int Array2::VertexBorderSize() const
{
  return nx * 2 + ny * 2 - 4;
}

/*!
\brief Return the Qt size of the vertex array.
*/
inline QSize Array2::GetQtSize() const
{
  return QSize(nx, ny);
}

/*!
\brief Return the size of the cell array.
*/
inline int Array2::CellSize() const
{
  return (nx - 1) * (ny - 1);
}

/*!
\brief Compute the coordinates of a point on the grid.
\param i,j Integer coordinates.
*/
inline Vector2 Array2::ArrayVertex(int i, int j) const
{
  return Vector2(a[0] + i * celldiagonal[0], a[1] + j * celldiagonal[1]);
}

/*!
\brief Compute the coordinates of a point on the grid.
\param p Point.
*/
inline Vector2 Array2::ArrayVertex(const QPoint& p) const
{
  return Vector2(a[0] + p.x() * celldiagonal[0], a[1] + p.y() * celldiagonal[1]);
}

/*!
\brief Get the box of the array.
*/
inline Box2 Array2::GetBox() const
{
  return Box2(a, b);
}

/*!
\brief Compute the index of a given cell.
\param i,j Integer coordinates of the cell.
*/
inline constexpr int Array2::VertexIndex(int i, int j) const
{
  return i + nx * j;
}
/*!
\brief Compute the index of a given cell.
\param p Point.
*/
inline constexpr int Array2::VertexIndex(const QPoint& p) const
{
  return p.x() + nx * p.y();
}

/*!
\brief Compute the coordinates of a given cell.
\param c Index of the cell.
\param i,j Integer coordinates of the cell.
*/
inline constexpr void Array2::InverseVertexIndex(int c, int& i, int& j) const
{
  i = c % nx;
  j = c / nx;
}

/*!
\brief Compute the coordinates of a given cell.
\param c Index of the cell.
*/
inline constexpr QPoint Array2::InverseVertexIndex(int c) const
{
  return QPoint(c % nx, c / nx);
}

/*!
\brief Check if the indexes are within range.
\param i,j Integer coordinates of the vertex.
*/
inline constexpr bool Array2::InsideCellIndex(int i, int j) const
{
  return (i >= 0) && (i < nx - 1) && (j >= 0) && (j < ny - 1);
}

/*!
\brief Check if the indexes are within range.
\param p Point.
*/
inline constexpr bool Array2::InsideCellIndex(const QPoint& p) const
{
  return (p.x() >= 0) && (p.x() < nx - 1) && (p.y() >= 0) && (p.y() < ny - 1);
}

/*!
\brief Check if the indexes are within range.
\param i,j Integer coordinates of the vertex.
*/
inline constexpr bool Array2::InsideVertexIndex(int i, int j) const
{
  return (i >= 0) && (i < nx) && (j >= 0) && (j < ny);
}

/*!
\brief Check if the indexes are outside or on the border.
\param i,j Integer coordinates of the vertex.
*/
inline constexpr bool Array2::OutsideVertexIndex(int i, int j) const
{
  if (i <= 0 || j <= 0 || i >= nx - 1 || j >= ny - 1) return false;
  return true;
}

/*!
\brief Check if the indexes are within range.
\param p Point.
*/
inline constexpr bool Array2::InsideVertexIndex(const QPoint& p) const
{
  return (p.x() >= 0) && (p.x() < nx) && (p.y() >= 0) && (p.y() < ny);
}

/*!
\brief Check if the indexes are on the border.
\param i,j Integer coordinates of the vertex.
*/
inline constexpr bool Array2::BorderVertexIndex(int i, int j) const
{
  return (i == 0) || (i == nx - 1) || (j == 0) || (j == ny - 1);
}

/*!
\brief Check if the indexes are on the border.
\param p Point.
*/
inline constexpr bool Array2::BorderVertexIndex(const QPoint& p) const
{
  return (p.x() == 0) || (p.x() == nx - 1) || (p.y() == 0) || (p.y() == ny - 1);
}

/*!
\brief Check if a point is in the rectangular domain.
\param p Point.
*/
inline bool Array2::Inside(const Vector2& p) const
{
  return Box2::Inside(p);
}

/*!
\brief Check if the indexes are within k-range.
\param i,j Integer coordinates of the vertex.
\param k Thickness of the boundary around the domain
*/
inline constexpr bool Array2::InsideVertexIndex(int i, int j, int k) const
{
  return (i >= 0 + k) && (i < nx - k) && (j >= 0 + k) && (j < ny - k);
}

/*!
\brief Compute the coordinates of a given cell.
\param c Index of the cell.
\param i,j Integer coordinates of the cell.
*/
inline void Array2::InverseCellIndex(int c, int& i, int& j) const
{
  i = c % (nx - 1);
  j = c / (nx - 1);
}

/*!
\brief Compute the point on the grid given an input point.
\param p Point.
\param i,j Integer coordinates of the cell.
\param u,v Coordinates of the point in the corresponding cell.
*/
inline void Array2::CellInteger(const Vector2& p, int& i, int& j, double& u, double& v) const
{
  Vector2 q = p - a;

  /*
    Vector2 d = b - a;

    u = q[0] / d[0];
    v = q[1] / d[1];

    // Scale
    u *= (nx - 1);
    v *= (ny - 1);
  */
  u = q[0] * inversecelldiagonal[0];
  v = q[1] * inversecelldiagonal[1];

  // Integer coordinates
  i = int(u);
  j = int(v);

  // Local coordinates within cell
  u -= i;
  v -= j;
}

/*!
\brief Compute the point on the grid given an input point.
\param p Point.
\param u,v Coordinates of the point in the corresponding cell.
*/
inline QPoint Array2::CellInteger(const Vector2& p, double& u, double& v) const
{
  Vector2 q = p - a;

  /*
    Vector2 d = b - a;

    u = q[0] / d[0];
    v = q[1] / d[1];

    // Scale
    u *= (nx - 1);
    v *= (ny - 1);
  */
  u = q[0] * inversecelldiagonal[0];
  v = q[1] * inversecelldiagonal[1];

  // Integer coordinates
  int i = int(u);
  int j = int(v);

  // Local coordinates within cell
  u -= i;
  v -= j;

  return QPoint(i, j);
}

/*!
\brief Compute the point on the grid given an input point.
\param p Point.
\param u, v Coordinates of the point in the corresponding cell.
*/
inline QPoint Array2::VertexInteger(const Vector2& p, double& u, double& v) const
{
  Vector2 q = p - a;

  /*
    Vector2 d = b - a;

    u = q[0] / d[0];
    v = q[1] / d[1];

    // Scale
    u *= (nx - 1);
    v *= (ny - 1);
  */
  u = q[0] * inversecelldiagonal[0];
  v = q[1] * inversecelldiagonal[1];

  // Integer coordinates
  int i = int(u);
  int j = int(v);

  // Local coordinates within cell
  u -= i;
  v -= j;

  return QPoint(i, j);
}

/*!
\brief Clamp vertex indexes to the size of the array.
\param i,j %Vertex indexes
*/
inline void Array2::ClampVertexIndex(int& i, int& j) const
{
  if (i < 0) { i = 0; }
  if (i > nx - 1) { i = nx - 1; }
  if (j < 0) { j = 0; }
  if (j > nx - 1) { j = nx - 1; }
}

/*!
\brief Compute the index of a given cell.
\param i,j Integer coordinates of the cell.
*/
inline constexpr int Array2::CellIndex(int i, int j) const
{
  return i + (nx - 1) * j;
}
/*!
\brief Compute the index of a given cell.
\param p Point.
*/
inline constexpr int Array2::CellIndex(const QPoint& p) const
{
  return p.x() + (nx - 1) * p.y();
}

/*!
\brief Compute the coordinates of a vertex inside a cell.
\param p Point.
\param i,j Integer coordinates of the cell.
*/
inline void Array2::VertexInteger(const Vector2& p, int& i, int& j) const
{
  Vector2 q = p - a;
  /*
  Vector2 d = b - a;

  double u = q[0] / d[0];
  double v = q[1] / d[1];

  i = int(u * (nx - 1));
  j = int(v * (ny - 1));
  */
  i = int(q[0] * inversecelldiagonal[0]);
  j = int(q[1] * inversecelldiagonal[1]);
}

/*!
\brief Compute the coordinates of a vertex inside a cell.
\param p Point.
*/
inline QPoint Array2::VertexInteger(const Vector2& p) const
{
  int i, j;
  VertexInteger(p, i, j);
  return QPoint(i, j);
}

/*!
\brief Compute the coordinates of a vertex inside a cell.
\param p Point.
*/
inline QPoint Array2::CellInteger(const Vector2& p) const
{
  int i, j;
  CellInteger(p, i, j);
  return QPoint(i, j);
}

/*!
\brief Return the range of integer values for the domain.
*/
inline QRect Array2::AreaInteger() const
{
  return QRect(0, 0, nx - 1, ny - 1);
}

/*!
\brief Compute the coordinates of a vertex inside a cell.
\param p Point.
\param i,j Integer coordinates of the cell.
*/
inline void Array2::CellInteger(const Vector2& p, int& i, int& j) const
{
  Vector2 q = p - a;
  /*
  Vector2 d = b - a;
  double u = q[0] / d[0];
  double v = q[1] / d[1];
  i = int(u * (nx - 1));
  j = int(v * (ny - 1));
  */
  i = int(q[0] * inversecelldiagonal[0]);
  j = int(q[1] * inversecelldiagonal[1]);
}

/*!
\brief Compute the point next to another one.
\param p Point.
\param n Next neighbor, should be in [0,7].
*/
inline QPoint Array2::Next(const QPoint& p, int n) const
{
  return p + next[n];
}

class Array :protected Box
{
protected:
  int nx, ny, nz; //!< Sizes
  Vector celldiagonal; //!< Cell diagonal.
  Vector inversecelldiagonal; //!< Inverse cell diagonal.
public:
  //! Empty.
  Array();
  explicit Array(const Box&, int, int, int);

  //! Empty
  ~Array() {}

  Array Extract(int, int, int, int, int, int) const;

  int VertexSize() const;
  int VertexSize(int) const;

  Box Cell(int) const;
  Box Cell(int, int, int) const;
  Vector CellCenter(int, int, int) const;

  void Scale(const Vector&);
  void Scale(const double&);

  Box GetBox() const;

  Vector CellDiagonal() const;
  double CellVolume() const;

  Vector Vertex(int, int, int) const;

  int GetSizeX() const;
  int GetSizeY() const;
  int GetSizeZ() const;

  int CellSize() const;
  int CellSizeX() const;
  int CellSizeY() const;
  int CellSizeZ() const;

  // Empty
  bool IsEmpty() const;

  Box UnitCell() const;

  void Resize(const Box&, int, int, int);

  // Vertex queries
  void VertexInteger(const Vector&, int&, int&, int&) const;
  void VertexIntegerVolume(const Box&, int&, int&, int&, int&, int&, int&) const;

  // Cell queries
  void CellInteger(const Vector&, int&, int&, int&) const;
  void CellInteger(const Vector&, int&, int&, int&, double&, double&, double&) const;
  void CellIntegerVolume(const Box&, int&, int&, int&, int&, int&, int&) const;

  friend std::ostream& operator<<(std::ostream&, const Array&);
  void OutStream(QDataStream&) const;
  void InStream(QDataStream&);

  int Memory() const;
  static Array CubicCells(Box, int);

public:
  constexpr int VertexIndex(int, int, int) const;
  constexpr bool InsideVertexIndex(int, int, int) const;
  constexpr bool InsideCellIndex(int, int, int) const;
  void InverseVertexIndex(int, int&, int&, int&) const;

  constexpr int CellIndex(int, int, int) const;
  void InverseCellIndex(int, int&, int&, int&) const;
};

/*!
\brief Detect if the array is empty.
*/
inline bool Array::IsEmpty() const
{
  return (nx <= 0) || (ny <= 0) || (nz <= 0);
}

/*!
\brief Get the size of the array for every axis.
\param i %Axis.
*/
inline int Array::VertexSize(int i) const
{
  if (i == 0) return nx;
  else if (i == 1) return ny;
  else return nz;
}

/*!
\brief Return the size of the array.
*/
inline int Array::VertexSize() const
{
  return nx * ny * nz;
}

/*!
\brief Get the box of the array.
*/
inline Box Array::GetBox() const
{
  return Box(a, b);
}

/*!
\brief Compute the coordinates of a point on the grid.
\param i,j,k Integer coordinates.
*/
inline Vector Array::Vertex(int i, int j, int k) const
{
  return Vector(a[0] + i * celldiagonal[0], a[1] + j * celldiagonal[1], a[2] + k * celldiagonal[2]);
}

/*!
\brief Compute the index of a given vertex.
\param i,j,k Integer coordinates of the vertex.
*/
inline constexpr int Array::VertexIndex(int i, int j, int k) const
{
  return (i * ny + j) * nz + k;
}

/*!
\brief Check if the indexes are within range.
\param i,j,k Integer coordinates of the vertex.
*/
inline constexpr bool Array::InsideCellIndex(int i, int j, int k) const
{
  return (i >= 0) && (i < nx - 1) && (j >= 0) && (j < ny - 1) && (k >= 0) && (k < nz - 1);
}

/*!
\brief Check if the indexes are within range.
\param i,j,k Integer coordinates of the vertex.
*/
inline constexpr bool Array::InsideVertexIndex(int i, int j, int k) const
{
  return (i >= 0) && (i < nx) && (j >= 0) && (j < ny) && (k >= 0) && (k < nz);
}

/*!
\brief Compute the coordinates of a given vertex.
\param c Index of the vertex.
\param i,j,k Integer coordinates of the vertex.
*/
inline void Array::InverseVertexIndex(int c, int& i, int& j, int& k) const
{
  // Extract last coefficient
  k = c % nz;
  // At this step, both i and j coefficients are still wrapped together
  j = c / nz;
  // Avoid temporary variable by getting i first
  i = j / ny;
  // Get j coefficent
  j = j % ny;
}

/*!
\brief Compute the coordinates of a given cell.
\param c Index of the cell.
\param i,j,k Integer coordinates of the cell.
*/
inline void Array::InverseCellIndex(int c, int& i, int& j, int& k) const
{
  // Extract last coefficient
  k = c % (nz - 1);
  // At this step, both i and j coefficients are still wrapped together
  j = c / (nz - 1);
  // Avoid temporary variable by getting i first
  i = j / (ny - 1);
  // Get j coefficent
  j = j % (ny - 1);
}

/*!
\brief Compute the point on the grid given an input point.
\param p Point.
\param i,j,k Integer coordinates of the cell.
\param u, v,w Coordinates of the point in the corresponding cell.
*/
inline void Array::CellInteger(const Vector& p, int& i, int& j, int& k, double& u, double& v, double& w) const
{
  Vector q = p - a;
  Vector d = b - a;

  u = q[0] / d[0];
  v = q[1] / d[1];
  w = q[2] / d[2];

  // Scale
  u *= (nx - 1);
  v *= (ny - 1);
  w *= (nz - 1);

  // Integer coordinates
  i = int(u);
  j = int(v);
  k = int(w);

  // Local coordinates within cell
  u -= i;
  v -= j;
  w -= k;
}

/*!
\brief Compute the index of a given cell.
\param i,j,k Integer coordinates of the cell.
*/
inline constexpr int Array::CellIndex(int i, int j, int k) const
{
  return (i * (ny - 1) + j) * (nz - 1) + k;
}

/*!
\brief Compute the coordinates of a vertex inside a cell.
\param p Point.
\param i, j ,k Integer coordinates of the cell.
*/
inline void Array::VertexInteger(const Vector& p, int& i, int& j, int& k) const
{
  Vector q = p - a;

  i = int(q[0] * inversecelldiagonal[0]);
  j = int(q[1] * inversecelldiagonal[1]);
  k = int(q[2] * inversecelldiagonal[2]);
}

/*!
\brief Compute the coordinates of a vertex inside a cell.
\param p Point.
\param i,j,k Integer coordinates of the cell.
*/
inline void Array::CellInteger(const Vector& p, int& i, int& j, int& k) const
{
  Vector q = p - a;
  double u = q[0] * inversecelldiagonal[0];
  double v = q[1] * inversecelldiagonal[1];
  double w = q[2] * inversecelldiagonal[2];

  i = int(u);
  j = int(v);
  k = int(w);
}

/*!
\brief Get the size of the array for x axis.
*/
inline int Array::GetSizeX() const
{
  return nx;
}

/*!
\brief Get the size of the array for y axis.
*/
inline int Array::GetSizeY() const
{
  return ny;
}

/*!
\brief Get the size of the array for z axis.
*/
inline int Array::GetSizeZ() const
{
  return nz;
}

/*!
\brief Return the size of the voxel.
*/
inline int Array::CellSize() const
{
  return (nx - 1) * (ny - 1) * (nz - 1);
}

/*!
\brief Get the size of the voxel for x axis.
*/
inline int Array::CellSizeX() const
{
  return nx - 1;
}

/*!
\brief Get the size of the voxel for y axis.
*/
inline int Array::CellSizeY() const
{
  return ny - 1;
}

/*!
\brief Get the size of the voxel for z axis.
*/
inline int Array::CellSizeZ() const
{
  return nz - 1;
}

/*!
\brief Create an array with enforced cubic cells.
\param box The box.
\param n Maximum subdivision.

\sa Box::SetParallelepipedic
*/
inline Array Array::CubicCells(Box box, int n)
{
  int nx, ny, nz;
  box.SetParallelepipedic(n, nx, ny, nz);
  return Array(box, nx, ny, nz);
}


