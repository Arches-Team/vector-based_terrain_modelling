// Voxels

#include "libs/voxel.h"

/*!
\class Voxel array.h
\brief A base three-dimensional voxel.

\ingroup StructureGroup

*/

/*!<
Set of eight cells, sorted.
*/
const int Voxel::octantcellindex[8][24] =
{
  // X  Y  Z  X  Y  Z  X  Y  Z  X  Y  Z
    {0, 0, 0,-1, 0, 0, 0,-1, 0,-1,-1, 0,
     0, 0,-1,-1, 0,-1, 0,-1,-1,-1,-1,-1},
    {0, 0, 0, 1, 0, 0, 0,-1, 0, 1,-1, 0,
     0, 0,-1, 1, 0,-1, 0,-1,-1, 1,-1,-1},
    {0, 0, 0,-1, 0, 0, 0, 1, 0,-1, 1, 0,
     0, 0,-1,-1, 0,-1, 0, 1,-1,-1, 1,-1},
    {0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0,
     0, 0,-1, 1, 0,-1, 0, 1,-1, 1, 1,-1},
    {0, 0, 0,-1, 0, 0, 0,-1, 0,-1,-1, 0,
     0, 0, 1,-1, 0, 1, 0,-1, 1,-1,-1, 1},
    {0, 0, 0, 1, 0, 0, 0,-1, 0, 1,-1, 0,
     0, 0, 1, 1, 0, 1, 0,-1, 1, 1,-1, 1},
    {0, 0, 0,-1, 0, 0, 0, 1, 0,-1, 1, 0,
     0, 0, 1,-1, 0, 1, 0, 1, 1,-1, 1, 1},
    {0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0,
     0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1}
};

/*!
\brief Create the voxel structure.
\param a The array.
*/
Voxel::Voxel(const Array& a) :Array(a)
{
  voxel.resize(CellSize());
}

/*!
\brief Create the voxel structure.
\param box The cubic box.
\param n Discretization.
*/
Voxel::Voxel(const Box& box, int n) :Voxel(box, n, n, n)
{
}

/*!
\brief Create the voxel structure.
\param box The box.
\param x,y,z Discretization, i.e. number of cells.
*/
Voxel::Voxel(const Box& box, int x, int y, int z) :Array(box, x + 1, y + 1, z + 1)
{
  voxel.resize(CellSize());
}

/*!
\brief Create and initialize the voxel structure.
\param box The box.
\param x,y,z Size of the array, i.e. number of cells.
\param v Voxel data.
*/
Voxel::Voxel(const Box& box, int x, int y, int z, const QVector<int>& v) :Array(box, x + 1, y + 1, z + 1)
{
  voxel = v;
}

/*!
\brief Check if a point is inside a non-empty cell.

Note that if the point is not inside the box, then it is defined as outside.

\param p Point.
*/
bool Voxel::Inside(const Vector& p) const
{
  if (!Array::Inside(p)) return false;
  int i, j, k;
  CellInteger(p, i, j, k);

  return At(i, j, k) > 0;
}

/*!
\brief Compute the integer coordinates of the (at most) eight octant cells closest to p.

\param p Point.
\param cells Returned array of integer coordinates.
\param configuration Returned configuration.
\param q Returned center of the octant that can be used for computing the distance.
\return The number of cells in the voxel structure.
*/
int Voxel::OctantCells(const Vector& p, int cells[24], int& configuration, Vector& q) const
{
  int i, j, k;
  CellInteger(p, i, j, k);

  Box cell = Cell(i, j, k);
  Vector center = cell.Center();
  int octant = center.Octant(p);

  // Returned center of the octant
  q = cell.Vertex(octant);

  int nc = 0;
  configuration = 0;
  for (int u = 0; u < 8; u++)
  {
    int iu = i + octantcellindex[octant][u * 3 + 0];
    int ju = j + octantcellindex[octant][u * 3 + 1];
    int ku = k + octantcellindex[octant][u * 3 + 2];
    if (InsideCellIndex(iu, ju, ku))
    {
      cells[3 * nc] = iu;
      cells[3 * nc + 1] = ju;
      cells[3 * nc + 2] = ku;
      nc++;

      if (At(iu, ju, ku) > 0)
      {
        configuration |= 1 << u;
      }
    }
  }

  return nc;
}

/*!
\brief Check if a vertex is a corner of the surface of voxel.

This function handles queries outside of the voxel.
If the coordinates are outside, function returns false.

\param x, y, z Integer coordinates of the cell.
*/
bool Voxel::IsVertex(int x, int y, int z) const
{
  int c = 0;
  c += AtExtended(x, y, z) > 0 ? 1 : 0;
  c += AtExtended(x - 1, y, z) ? 2 : 0;
  c += AtExtended(x, y - 1, z) ? 4 : 0;
  c += AtExtended(x - 1, y - 1, z) ? 8 : 0;
  c += AtExtended(x, y, z - 1) ? 1 : 16;
  c += AtExtended(x - 1, y, z - 1) ? 1 : 32;
  c += AtExtended(x, y - 1, z - 1) ? 1 : 64;
  c += AtExtended(x - 1, y - 1, z - 1) ? 128 : 0;

  if ((c == 0) || (c == 255))
  {
    return false;
  }
  return true;
}


/*!
\brief Compute the <I>approximate</I> signed distance to the voxel.

\param p Point.
*/
double Voxel::Signed(const Vector& p) const
{
  const double u = 0.5 * celldiagonal[0];

  // Outside 
  if (!Array::Inside(p))
  {
    // Distance to box
    double d = Box::Signed(p);

    if (d > u)
    {
      return d;
    }
  }

  // Inside
  double d = u;
  int c[24];
  int configuration = 0;
  Vector q;
  int nc = OctantCells(p, c, configuration, q);

  // The following simple code does not work: it leaves 0 distance on internal (shared) faces.
  /*
  for (int i = 0; i < nc; i++)
  {
    if (At(c[i * 3], c[i * 3 + 1], c[i * 3 + 2]) > 0)
    {
      d = Math::Min(d, Cell(c[i * 3], c[i * 3 + 1], c[i * 3 + 2]).Signed(p));
    }
  }*/
  d = OctantSigned(p, q, configuration);
  if (d > 0.0)
  {
    d = Math::Min(d, u);
  }
  else
  {
    d = Math::Max(d, -u);
  }
  return d;
}

/*!
\brief Compute the amount of memory used by the voxel structure.
*/
unsigned int Voxel::Memory() const
{
  return int(voxel.size()) * sizeof(int) + sizeof(Array);
}

/*!
\brief Compute the integer accessibility of a voxel.

Uses a 1-neighborhood, thus checking 3<SUP>3</SUP>3-1=26 cells.

\param x, y, z Integer coordinates of the cell.
*/
int Voxel::Accessibility(int x, int y, int z) const
{
  // Range
  int ax = max(x - 1, 0);
  int bx = min(x + 2, nx - 1);
  int ay = max(y - 1, 0);
  int by = min(y + 2, ny - 1);
  int az = max(z - 1, 0);
  int bz = min(z + 2, nz - 1);

  int accessibility = 0;

  for (int i = ax; i < bx; i++)
  {
    for (int j = ay; j < by; j++)
    {
      for (int k = az; k < bz; k++)
      {
        if (At(x + i, y + j, z + k) > 0)
        {
          accessibility++;
        }
      }
    }
  }
  return accessibility;
}

/*!
\brief Test if a cell is inside or outside of the voxel.

This function handles queries outside of the voxel.
If the coordinates are outside, function returns false.

\param x, y, z Integer coordinates of the cell.
*/
bool Voxel::Inside(int x, int y, int z) const
{
  if (!InsideCellIndex(x, y, z))
  {
    return false;
  }
  return At(x, y, z) > 0;
}

/*!
\brief Computes the distance between a point a eight octant volumes.

The octant configuration is provided as an integer whose bits are set to 0 if octants are empty, 1 otherwise.

\sa Box::Octant(const Vector&)
\param p Point.
\param q Origin of the eight octant volumes.
\param c Octant configuration.
*/
double Voxel::OctantSigned(const Vector& p, const Vector& q, int c)
{
  Vector e = p - q;
  int o = q.Octant(p);

  // Handle relative octant position of p with respect to q (8 cases)
  if ((o & 1) == 0)
  {
    // Symmetry: x
    e[0] = -e[0];
    // Swap bits
    o = ((o & 0b01010101) << 1) | ((o & 0b10101010) >> 1);
  }
  if ((o & 2) == 0)
  {
    // Symmetry: y
    e[1] = -e[1];
    // Swap bits
    o = ((o & 0b00110011) << 2) | ((o & 0b11001100) >> 2);
  }
  if ((o & 4) == 0)
  {
    // Symmetry: z
    e[2] = -e[2];
    // Swap bits
    o = ((o & 0b00001111) << 4) | ((o & 0b11110000) >> 4);
  }

  // At this stage, e is the relative vector as if p had been located in the positive octant
  // The (signed) distance can be any of the following: distance to an axis, a plane, or the center point

  // Point p is inside a voxel
  if (c & 0b10000000)
  {
    double d = -Math::Infinity;

    // Center
    if (!(c & 0b00000001))
    {
      d = -Norm(e);
    }

    // Check edges
    if (!(c & 0b00000010))
    {
      d = Math::Min(d, Norm(Vector2(e[1], e[2]))); // Edge x
    }
    if (!(c & 0b00000100))
    {
      d = Math::Min(d, Norm(Vector2(e[0], e[2]))); // Edge y
    }
    if (!(c & 0b00010000))
    {
      d = Math::Min(d, Norm(Vector2(e[0], e[1]))); // Edge z
    }

    // Check planes
    if (!(c & 0b00100000))
    {
      d = Math::Max(d, -e[1]); // Plane xz
    }
    if (!(c & 0b01000000))
    {
      d = Math::Max(d, -e[0]); // Plane yz
    }
    if (!(c & 0b00001000))
    {
      d = Math::Max(d, -e[2]); // Plane xy
    }
    return d;
  }
  else
  {
    double d = Math::Infinity;

    // Center
    if (c & 0b00000001)
    {
      d = Norm(e);
    }

    // Check edges
    if (c & 0b00000010)
    {
      d = Math::Min(d, Norm(Vector2(e[1], e[2]))); // Edge x
    }
    if (c & 0b00000100)
    {
      d = Math::Min(d, Norm(Vector2(e[0], e[2]))); // Edge y
    }
    if (c & 0b00010000)
    {
      d = Math::Min(d, Norm(Vector2(e[0], e[1]))); // Edge z
    }
    // Check planes
    if (c & 0b00100000)
    {
      d = Math::Min(d, e[1]); // Plane xz
    }
    if (c & 0b01000000)
    {
      d = Math::Min(d, e[0]); // Plane yz
    }
    if (c & 0b00001000)
    {
      d = Math::Min(d, e[2]); // Plane xy
    }
    return d;
  }
}

/*!
\brief Create a set of cubes representing a voxel.
\param u Unit reference box of the voxel.
\param surface Only produce cubes that participate to the visible surface.
*/
QVector<Vector> Voxel::GetCubes(Box& u, bool surface) const
{
  u = UnitCell();

  QVector<Vector> cubes;

  if (surface == false)
  {
    // Parse all cells
    for (int i = 0; i < CellSize(); i++)
    {
      if (voxel.at(i) != 0)
      {
        cubes.append(Cell(i).Center());
      }
    }
  }
  // Discard cubes that do not have a frontier 
  else
  {
    // Parse all cells
    for (int i = 0; i < CellSizeX(); i++)
    {
      for (int j = 0; j < CellSizeY(); j++)
      {
        for (int k = 0; k < CellSizeZ(); k++)
        {
          // Compute configuration
          if (At(i, j, k) != 0)
          {
            bool s = false;
            if (i - 1 < 0)
            {
              s = true;
            }
            else if (At(i - 1, j, k) == 0)
            {
              s = true;
            }

            if (i + 1 > CellSizeX() - 1)
            {
              s = true;
            }
            else if (At(i + 1, j, k) == 0)
            {
              s = true;
            }
            if (j - 1 < 0)
            {
              s = true;
            }
            else if (At(i, j - 1, k) == 0)
            {
              s = true;
            }
            if (j + 1 > CellSizeY() - 1)
            {
              s = true;
            }
            else if (At(i, j + 1, k) == 0)
            {
              s = true;
            }
            if (k - 1 < 0)
            {
              s = true;
            }
            else if (At(i, j, k - 1) == 0)
            {
              s = true;
            }

            if (k + 1 > CellSizeZ() - 1)
            {
              s = true;
            }
            else if (At(i, j, k + 1) == 0)
            {
              s = true;
            }

            if (s == true)
            {
              cubes.append(Cell(i, j, k).Center());
            }
          }
        }
      }
    }
  }
  return cubes;
}
