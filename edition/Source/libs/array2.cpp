// Arrays

#include <QtGui/QPainter>

#include "libs/array.h"
#include "libs/segment.h"
#include "libs/triangle.h"

/*!
\class Array2 array.h
\brief A base two-dimensional array structure.

\image html array.png

An array is a two-dimensional rectangle virtually decomposed into cells.

\ingroup StructureGroup
*/

const QPoint Array2::next[8] = { QPoint(1, 0), QPoint(1, 1), QPoint(0, 1), QPoint(-1, 1), QPoint(-1, 0), QPoint(-1, -1), QPoint(0, -1), QPoint(1, -1) };
const double Array2::length[8] = { 1.0, sqrt(2.0), 1.0, sqrt(2.0), 1.0, sqrt(2.0), 1.0, sqrt(2.0) };
const double Array2::inverselength[8] = { 1.0, 1.0 / sqrt(2.0), 1.0, 1.0 / sqrt(2.0), 1.0, 1.0 / sqrt(2.0), 1.0, 1.0 / sqrt(2.0) };

/*!
\brief Empty array, with empty box.
*/
Array2::Array2() :Box2(Box2::Null), nx(0), ny(0), celldiagonal(Vector2::Null), inversecelldiagonal(Vector2::Null)
{
}

/*!
\brief Create the array structure.
\param box The box.
\param x,y Size of the array.
*/
Array2::Array2(const Box2& box, int x, int y) :Box2(box), nx(x), ny(y), celldiagonal((b[0] - a[0]) / (nx - 1), (b[1] - a[1]) / (ny - 1))
{
  inversecelldiagonal = celldiagonal.Inverse();
}

/*!
\brief Create the array structure.
\param box The box.
\param n Size of the array.
*/
Array2::Array2(const Box2& box, int n) :Array2(box, n, n)
{
}

/*!
\brief Scales the array structure.

\param s Scaling vector.
*/
void Array2::Scale(const Vector2& s)
{
  Box2::Scale(s);
  celldiagonal = Vector2((b[0] - a[0]) / (nx - 1), (b[1] - a[1]) / (ny - 1));
  inversecelldiagonal = celldiagonal.Inverse();
}

/*!
\brief Scales the array structure.

\param s Scaling factor.
*/
void Array2::Scale(const double& s)
{
  Box2::Scale(Vector2(s));
  celldiagonal = Vector2((b[0] - a[0]) / (nx - 1), (b[1] - a[1]) / (ny - 1));
  inversecelldiagonal = celldiagonal.Inverse();
}


/*!
\brief Compute the k-th point on the boundary of the domain.

\image html borders.png

\param k Integer index.
*/
QPoint Array2::VertexBorder(int k) const
{
  int x;
  int y;
  if (k < 2 * nx)
  {
    if (k < nx)
    {
      x = k;
      y = 0;
    }
    else
    {
      x = k - nx;
      y = ny - 1;
    }
  }
  else
  {
    k -= 2 * nx;
    if (k < ny - 2)
    {
      x = 0;
      y = k + 1;
    }
    else
    {
      x = nx - 1;
      y = k + 1 - (ny - 2);
    }
  }

  return QPoint(x, y);
}


/*!
\brief Compute the index of a point on the boundary of the domain.
\param i,j Coordinates.
*/
int Array2::VertexBorderIndex(int i, int j) const
{
  if (j == 0)
  {
    return i;
  }
  else if (j == ny - 1)
  {
    return i + nx;
  }
  else
  {
    // Already processed corners in previous if statements
    if (i == 0)
    {
      return nx * 2 + j - 1;
    }
    else
    {
      return nx * 2 + ny - 2 + j - 1;
    }
  }
}

/*!
\brief Return the geometry of the cell.
\param c Integer cell index.
*/
Box2 Array2::Cell(int c) const
{
  int i, j;
  InverseCellIndex(c, i, j);
  return Cell(i, j);
}

/*!
\brief Return the geometry of the cell.
\param i,j Integer coordinates of the cell.
*/
Box2 Array2::Cell(int i, int j) const
{
  Vector2 e = ArrayVertex(i, j);
  return Box2(e, e + celldiagonal);
}

/*!
\brief Return the geometry of the cell.
\param i, j Integer coordinates of the cell.
\param d Boolean, direct diagonal if true, opposite diagonal otherwize.
\param ta, tb Returned triangles.
*/
void Array2::HalfCell(int i, int j, bool d, Triangle2& ta, Triangle2& tb) const
{
  // Box of the cell
  Box2 b = Cell(i, j);

  Vector2 a10 = Vector2(b[1][0], b[0][1]);
  Vector2 a01 = Vector2(b[0][0], b[1][1]);

  if (d == true)
  {
    ta = Triangle2(b[0], a10, b[1]);
    tb = Triangle2(b[0], b[1], a01);
  }
  else
  {
    ta = Triangle2(a01, b[0], a01);
    tb = Triangle2(a01, a10, b[1]);
  }
}

/*!
\brief Return the center of the cell.
\param i,j Integer coordinates of the cell.
*/
Vector2 Array2::CellCenter(int i, int j) const
{
  return  a + celldiagonal.Scaled(Vector2(i, j)) + 0.5 * celldiagonal;
}

/*!
\brief Return the geometry of a generic cell.

This is a box, centered at origin, with the right size.
*/
Box2 Array2::UnitCell() const
{
  // Half diagonal scaled by grid size
  Vector2 h = 0.5 * celldiagonal;
  return Box2(-h, h);
}

/*!
\brief Overloaded.
\param s Stream.
\param a The array.
*/
std::ostream& operator<<(std::ostream& s, const Array2& a)
{
  s << "Array2(" << a.a << ',' << a.b << "," << a.nx << ',' << a.ny << ')';
  return s;
}


/*!
\brief Return the cell diagonal vector.
*/
Vector2 Array2::CellDiagonal() const
{
  return celldiagonal;
}

/*!
\brief Return the area of a cell.

This is convenience function, which is the same as the following code, but more efficient:
\code
Array2 array;
double a=array.UnitCell().Area();
\endcode
*/
double Array2::CellArea() const
{
  return Box2::Area() / ((nx - 1) * (ny - 1));
}

/*!
\brief Translate the array.

Centering the array can be coded as:
\code
Array2 a;
a.Translate(-a.Center());
\endcode

\param t Translation vector.
*/
void Array2::Translate(const Vector2& t)
{
  Box2::Translate(t);
}

/*!
\brief Return the center of the array.
*/
Vector2 Array2::Center() const
{
  return Box2::Center();
}

/*!
\brief Create a texture with a regular grid.
\param box The box.
\param a Distance between lines.
*/
QImage Array2::ImageGrid(const Box2& box, const double& a) const
{
  QImage image(nx, ny, QImage::Format_RGB32);
  double xa = fmod(box[0][0], a) * a;
  double ya = fmod(box[0][1], a) * a;

  QPainter paint(&image);
  // Grey color
  paint.setPen(QColor(128, 128, 128, 128));

  for (double x = xa; x < box[1][0]; x += a)
  {
    QPoint p = VertexInteger(Vector2(x, 0.0));
    paint.drawLine(p.x(), 0, p.x(), ny - 1);
  }
  for (double y = ya; y < box[1][1]; y += a)
  {
    QPoint p = VertexInteger(Vector2(0.0, y));
    paint.drawLine(0, p.y(), nx - 1, p.y());
  }
  return image;
}

/*!
\brief Draw the cells of the grid.
\param scene Graphics scene.
\param pen The pen.
\param brush The brush.
*/
void Array2::Draw(QGraphicsScene& scene, const QPen& pen, const QBrush& brush) const
{
  // Contour
  Box2::Draw(scene, pen, brush);

  for (int i = 1; i < nx - 1; i++)
  {
    Segment2(ArrayVertex(i, 0), ArrayVertex(i, ny - 1)).Draw(scene, pen);
  }
  for (int j = 1; j < ny - 1; j++)
  {
    Segment2(ArrayVertex(0, j), ArrayVertex(nx - 1, j)).Draw(scene, pen);
  }
}

/*!
\brief Compute the integer coordinates of the vertices embedding a cell with prescribed integer radius.

The coordinates are clamped to the size of the array.
\param i,j Cell.
\param r %Radius.
*/
QRect Array2::VertexIntegerArea(int i, int j, int r) const
{
  QRect area(i - r, j - r, 2 * r + 1, 2 * r + 1);

  // Limit to domain
  QRect mask(0, 0, nx - 1, ny - 1);
  return area.intersected(mask);
}

/*!
\brief Compute the integer coordinates of the vertices embedding a box.

The coordinates are clamped to the size of the array.
\param box %Box.
*/
QRect Array2::VertexIntegerArea(const Box2& box) const
{
  QPoint pa = VertexInteger(box[0]);
  QPoint pb = VertexInteger(box[1]);

  // Rectangle
  QRect area(pa.x(), pa.y(), pb.x() - pa.x() + 1, pb.y() - pa.y() + 1);

  // Limit to domain
  QRect mask(0, 0, nx - 1, ny - 1);
  return area.intersected(mask);
}

/*!
\brief Compute the integer coordinates of the vertices embedding a circle.

The coordinates are clamped to the size of the array.
\param circle %Circle.
*/
QRect Array2::VertexIntegerArea(const Circle2& circle) const
{
  return VertexIntegerArea(circle.GetBox());
}

/*!
\brief Compute the integer coordinates span of the vertices embedding a rectangle.

\param rect %Rectangle.
*/
QRect Array2::VertexIntegerArea(const QRect& rect) const
{
  QRect mask = AreaInteger();
  return rect.intersected(mask);
}

/*!
\brief Compute the integer coordinates of the cells of the array embedding an input box.

The coordinates are clamped to the size of the array.
\param box %The box.
*/
QRect Array2::CellIntegerArea(const Box2& box) const
{
  QPoint ca = CellInteger(box[0]);
  QPoint cb = CellInteger(box[1]);

  // Rectangle
  QRect area(ca.x(), ca.y(), cb.x() - ca.x() + 1, cb.y() - ca.y() + 1);

  // Limit to domain
  QRect mask(0, 0, nx - 1, ny - 1);

  return area.intersected(mask);
}

/*!
\brief Compute the integer coordinates of the cells of the array embedding an input circle.

The coordinates are clamped to the size of the array.
\param circle %The circle.
*/
QRect Array2::CellIntegerArea(const Circle2& circle) const
{
  return CellIntegerArea(circle.GetBox());
}

/*!
\brief Get statistics.
*/
QString Array2::Statistics() const
{
  QString statistics = QString("Array2 %1 x %2 \nResolution %3 x %4\n Cell %5 x %6\n").arg(Size()[0]).arg(Size()[1]).arg(nx).arg(ny).arg(CellDiagonal()[0]).arg(CellDiagonal()[1]);
  return statistics;
}

/*!
\brief Convert a direction to an 8 bit integer code.
\param a,b Neighbor direction vector.
*/
int Array2::NeighborCode(int a, int b)
{
  int c = 0;

  if (a == 1 && b == 0) 		c = 1;
  if (a == 1 && b == -1)		c = 2;
  if (a == 0 && b == -1)		c = 4;
  if (a == -1 && b == -1)		c = 8;
  if (a == -1 && b == 0)		c = 16;
  if (a == -1 && b == 1)		c = 32;
  if (a == 0 && b == 1)		c = 64;
  if (a == 1 && b == 1)		c = 128;

  return c;
  // static const int x[9] = { 32,64,128,16,0,1,8,4,2 };
  // return x[(b + 1) * 3 + a + 1];
}

/*!
\brief Convert an 8 bit integer code direction into a direction vector.
\param c Direction code.
*/
QPoint Array2::CodeToDir(int c)
{
  /*
  static const QPoint q[8 * 2] = {
      QPoint(1, 0),
      QPoint(1,-1),
      QPoint(0,-1),
      QPoint(-1,-1),
      QPoint(-1, 0),
      QPoint(-1, 1),
      QPoint(0, 1),
      QPoint(1, 1) };
  return q[c];
  */
  int a = 0;
  int b = 0;

  if (c == 1) { a = 1; b = 0; }
  if (c == 2) { a = 1; b = -1; }
  if (c == 4) { a = 0; b = -1; }
  if (c == 8) { a = -1; b = -1; }
  if (c == 16) { a = -1; b = 0; }
  if (c == 32) { a = -1; b = 1; }
  if (c == 64) { a = 0; b = 1; }
  if (c == 128) { a = 1; b = 1; }
  return QPoint(a, b);
}

/*!
\brief Change the resolution.

Increases the resolution by two.
*/
void Array2::Subdivide()
{
  // Change resolution
  nx = nx * 2 - 1;
  ny = ny * 2 - 1;

  // Change diagonals
  celldiagonal = Vector2((b[0] - a[0]) / (nx - 1), (b[1] - a[1]) / (ny - 1));
  inversecelldiagonal = celldiagonal.Inverse();
}

/*!
\brief Compute the coordinates of a set of points on the grid.

\sa ArrayVertex
\param s Set of point.
*/
QVector <Vector2> Array2::ArrayVertexes(const QVector<QPoint>& s) const
{
  QVector<Vector2> v(s.size());
  for (int i = 0; i < s.size(); i++)
  {
    v[i] = ArrayVertex(s.at(i));
  }
  return v;
}
