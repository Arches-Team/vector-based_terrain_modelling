// Arrays

#include "libs/array.h"

/*!
\class Array array.h
\brief A core three dimensional grid structure.

An array is a three-dimensional box virtually decomposed into cells.

\ingroup StructureGroup
*/

/*!
\brief Empty array, with empty box.
*/
Array::Array() :Box(Box::Null), nx(0), ny(0), nz(0), celldiagonal(Vector::Null), inversecelldiagonal(Vector::Null)
{
}

/*!
\brief Create the array structure.
\param box The box.
\param x,y,z Size of the array.
*/
Array::Array(const Box& box, int x, int y, int z) :Box(box), nx(x), ny(y), nz(z)
{
  celldiagonal = Vector((b[0] - a[0]) / (nx - 1), (b[1] - a[1]) / (ny - 1), (b[2] - a[2]) / (nz - 1));
  inversecelldiagonal = celldiagonal.Inverse();
}

/*!
\brief Scales the array structure.

\param s Scaling vector.
*/
void Array::Scale(const Vector& s)
{
  Box::Scale(s);
  celldiagonal = Vector((b[0] - a[0]) / (nx - 1), (b[1] - a[1]) / (ny - 1), (b[2] - a[2]) / (nz - 1));
  inversecelldiagonal = celldiagonal.Inverse();
}

/*!
\brief Scales the array structure.

\param s Scaling factor.
*/
void Array::Scale(const double& s)
{
  Box::Scale(s);
  celldiagonal = Vector((b[0] - a[0]) / (nx - 1), (b[1] - a[1]) / (ny - 1), (b[2] - a[2]) / (nz - 1));
  inversecelldiagonal = celldiagonal.Inverse();
}

/*!
\brief Return the geometry of the cell.
\param c Integer cell index.
*/
Box Array::Cell(int c) const
{
  int i, j, k;
  InverseCellIndex(c, i, j, k);
  return Cell(i, j, k);
}

/*!
\brief Return the cell diagonal.
*/
Vector Array::CellDiagonal() const
{
  return celldiagonal;
}

/*!
\brief Return the volume of a cell.

This is convenience function, which is the same as the following code, but more efficient:
\code
Array array;
double a=array.UnitCell().Volume();
\endcode
*/
double Array::CellVolume() const
{
  return celldiagonal[0] * celldiagonal[1] * celldiagonal[2];
}

/*!
\brief Return the geometry of the cell.
\param i,j,k Integer coordinates of the cell.
*/
Box Array::Cell(int i, int j, int k) const
{
  Vector e = a + celldiagonal.Scaled(Vector(i, j, k));
  return Box(e, e + celldiagonal);
}

/*!
\brief Return the center of the cell.
\param i,j,k Integer coordinates of the cell.
*/
Vector Array::CellCenter(int i, int j, int k) const
{
  return  a + celldiagonal.Scaled(Vector(i, j, k)) + 0.5 * celldiagonal;
}

/*!
\brief Return the geometry of a generic cell.

This is a box, centered at origin, with the right size.
*/
Box Array::UnitCell() const
{
  Vector h = 0.5 * celldiagonal;
  return Box(-h, h);
}

/*!
\brief Resize the voxel with the given parameters.
\param box The box.
\param x,y,z Sizes.
*/
void Array::Resize(const Box& box, int x, int y, int z)
{
  nx = x;
  ny = y;
  nz = z;

  Box::a = box[0];
  Box::b = box[1];
}

/*!
\brief Overloaded.
\param s Stream.
\param a The array.
*/
std::ostream& operator<<(std::ostream& s, const Array& a)
{
  s << "Array(" << a.a << ',' << a.b << "," << a.nx << ',' << a.ny << ',' << a.nz << ')';
  return s;
}

/*!
\brief Write data to a data stream.
\param s Stream.
*/
void Array::OutStream(QDataStream& s) const
{
  s << a[0] << a[1] << a[2];
  s << b[0] << b[1] << b[2];
  s << nx;
  s << ny;
  s << nz;
}

/*!
\brief Read data to a data stream.
\param s Stream.
*/
void Array::InStream(QDataStream& s)
{
  s >> a[0] >> a[1] >> a[2];
  s >> b[0] >> b[1] >> b[2];
  s >> nx;
  s >> ny;
  s >> nz;
}

/*!
\brief Compute the size of the array.
*/
int Array::Memory() const
{
  return sizeof(Array);
}

/*!
\brief Compute the integer coordinates of the vertices of a box embedding an input box.

The coordinates are clamped to the size of the array.
\param box %Box.
\param ax, ay, az, bx, by, bz Integer coordinates of the two points a and b.
*/
void Array::VertexIntegerVolume(const Box& box, int& ax, int& ay, int& az, int& bx, int& by, int& bz) const
{
  VertexInteger(box[0], ax, ay, az);
  VertexInteger(box[1], bx, by, bz);

  // Limit to domain
  ax = max(ax, 0);
  ay = max(ay, 0);
  az = max(az, 0);
  bx = min(bx, nx - 1);
  by = min(by, ny - 1);
  bz = min(bz, nz - 1);
}

/*!
\brief Compute the integer coordinates of the cells of the array embedding an input box.

The coordinates are clamped to the size of the array.
\param box %Box.
\param ax, ay, az, bx, by, bz Integer coordinates of the two points a and b.
*/
void Array::CellIntegerVolume(const Box& box, int& ax, int& ay, int& az, int& bx, int& by, int& bz) const
{
  CellInteger(box[0], ax, ay, az);
  CellInteger(box[1], bx, by, bz);

  // Limit to domain
  ax = max(ax, 0);
  ay = max(ay, 0);
  az = max(az, 0);
  bx = min(bx, nx - 1);
  by = min(by, ny - 1);
  bz = min(bz, nz - 1);
}

/*!
\brief Extract an sub-array.

Compute the coordinates of the box and set the subdivision.

\param xa,ya,za,xb,yb,zb Integer position of the array vertices.
*/
Array Array::Extract(int xa, int ya, int za, int xb, int yb, int zb) const
{
  return Array(Box(Vertex(xa, ya, za), Vertex(xb, yb, zb)), xb - xa + 1, yb - ya + 1, zb - za + 1);
}
