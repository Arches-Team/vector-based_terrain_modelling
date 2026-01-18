// Voxel

#pragma once

#include <QtCore/QVector>

#include "libs/array.h"

class Voxel :public Array
{
protected:
  QVector<int> voxel; //!< %Array of integer storing the data.
public:
  Voxel() { }
  explicit Voxel(const Array&);
  explicit Voxel(const Box&, int);
  explicit Voxel(const Box&, int, int, int);
  explicit Voxel(const Box&, int, int, int, const QVector<int>&);

  // Access to elements
  int At(int, int, int) const;
  int At(int) const;

  int AtExtended(int, int, int) const;

  int& operator()(int, int, int);
  int& operator()(int);

  int& operator[](int);

  bool Inside(const Vector&) const;
  bool Inside(int, int, int) const;

  bool IsVertex(int, int, int) const;

  double Signed(const Vector&) const;
  int OctantCells(const Vector&, int[24], int&, Vector&) const;

  QVector<Vector> GetCubes(Box&, bool = false) const;

  int Accessibility(int, int, int) const;

  unsigned int Memory() const;

  // Octant distance
  static double OctantSigned(const Vector&, const Vector&, int);
protected:
  static const int octantcellindex[8][24];
};

/*!
\brief Return the data in the cell of the voxel.
\param c Cell index.
*/
inline int Voxel::At(int c) const
{
  return voxel.at(c);
}

/*!
\brief Return the data in the cell of the voxel.

\sa AtExtended()

\param i,j,k Integer coordinates of the cell.
*/
inline int Voxel::At(int i, int j, int k) const
{
  return voxel.at(CellIndex(i, j, k));
}

/*!
\brief Return the data in the cell of the voxel.

This function handles queries outside of the voxel.
If the coordinates are outside, function returns 0.

\sa At()

\param i,j,k Integer coordinates of the cell.
*/
inline int Voxel::AtExtended(int i, int j, int k) const
{
  if (!InsideCellIndex(i, j, k))
  {
    return 0;
  }

  return voxel.at(CellIndex(i, j, k));
}

/*!
\brief Return the data in the cell of the voxel.
\param c Cell index.
*/
inline int& Voxel::operator()(int c)
{
  return voxel[c];
}

/*!
\brief Return the data in the cell of the voxel.
\param i,j,k Integer coordinates of the cell.
*/
inline int& Voxel::operator()(int i, int j, int k)
{
  return voxel[CellIndex(i, j, k)];
}

/*!
\brief Return the data in the cell of the voxel.
\param c Integer cell index.
*/
inline int& Voxel::operator[](int c)
{
  return voxel[c];
}


