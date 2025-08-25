
// Fields

#include "libs/scalarfield.h"
#include "libs/mesh.h"
#include "libs/meshcolor.h"

/*!
\brief Compute the polygonal mesh approximating the implicit surface.

This member function implements the marching cube tesselation algorithm, using
a fixed look-up table. This technique does suffer from the ambiguities in
the traditional marching cubes algorithms where some special cases result
in unexpected holes in the resulting surface.

<BR>The first part of the algorithm uses a table (edgeTable) which maps the
vertices under the isosurface to the intersecting edges. An 8 bit index
is formed where each bit corresponds to a vertex.

<P>Looking up the edge table returns a 12 bit number, each bit corresponding to an edge, 0 if the edge isn't cut by the
isosurface, 1 if the edge is cut by the isosurface. If none of the edges are cut the table returns a 0, this occurs when
cubeindex is 0 (all vertices below the isosurface) or 0xff (all vertices above the isosurface).

<P>The last part of the algorithm involves forming the correct facets from the positions that the isosurface intersects the edges
of the grid cell. Again a table is used which this time uses the same cubeindex but allows the vertex sequence to
be looked up for as many triangular facets are necessary to represent the isosurface within the grid cell. There at most 5
triangular facets necessary.

<P>Lets say vertex 0 and 3 are below the isosurface. cubeindex will then be 0000 1001==9.
The 9th entry into the egdeTable is 905 = 1001 0000 0101: edges 11, 8, 2, and 0 are
cut and so we work out the vertices of the intersection of the isosurface with those edges.
Next 9 in the triangle table is 0, 11, 2, 8, 11, 0. This corresponds to 2 triangular facets, one between the intersection of edge 0
11 and 2; the other between the intersections along edges 8, 11, and 0.

\param box %Box defining the region that will be polygonized.
\param n Discretization parameter.
\param g Returned geometry.
\param epsilon Epsilon value for computing vertices on straddling edges.
*/
void AnalyticScalarField::Polygonize(int n, Mesh& g, const Box& box, const double& epsilon) const
{
  QVector<Vector> vertex;
  QVector<Vector> normal;

  QVector<int> triangle;

  vertex.reserve(200000);
  normal.reserve(200000);
  triangle.reserve(200000);

  int nv = 0;
  int nt = 0;
  int nx, ny, nz;

  Box clipped = box;
  clipped.SetParallelepipedic(n, nx, ny, nz);

  // Number of processed cubes
  ncubes = (nx - 1) * (ny - 1) * (nz - 1);

  // Clamped integer values
  int nax = 0;
  int nbx = nx;
  int nay = 0;
  int nby = ny;
  int naz = 0;
  int nbz = nz;

  const int size = nx * ny;

  // Intensities
  double* a = new double[size];
  double* b = new double[size];

  // Vertex
  Vector2* u = new Vector2[size];

  // Edges
  int* eax = new int[size];
  int* eay = new int[size];
  int* ebx = new int[size];
  int* eby = new int[size];
  int* ez = new int[size];

  // Diagonal of a cell
  Vector d = Array(clipped, nx, ny, nz).CellDiagonal();

  double za = clipped[0][2];

  // Compute field inside lower Oxy plane
  for (int i = nax; i < nbx; i++)
  {
    for (int j = nay; j < nby; j++)
    {
      u[i * ny + j] = clipped[0] + Vector2(i * d[0], j * d[1]);
      a[i * ny + j] = Value(u[i * ny + j].ToVector(za));
    }
  }

  // Compute straddling edges inside lower Oxy plane
  for (int i = nax; i < nbx - 1; i++)
  {
    for (int j = nay; j < nby; j++)
    {
      // We need a xor b, which can be implemented a == !b 
      if (!((a[i * ny + j] < 0.0) == !(a[(i + 1) * ny + j] >= 0.0)))
      {
        vertex.append(Dichotomy(u[i * ny + j].ToVector(za), u[(i + 1) * ny + j].ToVector(za), a[i * ny + j], a[(i + 1) * ny + j], d[0], epsilon));
        normal.append(Normal(vertex.last()));
        eax[i * ny + j] = nv;
        nv++;
      }
    }
  }
  for (int i = nax; i < nbx; i++)
  {
    for (int j = nay; j < nby - 1; j++)
    {
      if (!((a[i * ny + j] < 0.0) == !(a[i * ny + (j + 1)] >= 0.0)))
      {
        vertex.append(Dichotomy(u[i * ny + j].ToVector(za), u[i * ny + (j + 1)].ToVector(za), a[i * ny + j], a[i * ny + (j + 1)], d[1], epsilon));
        normal.append(Normal(vertex.last()));
        eay[i * ny + j] = nv;
        nv++;
      }
    }
  }

  // Array for edge vertices
  int e[12];

  // For all layers
  for (int k = naz; k < nbz; k++)
  {
    double zb = za + d[2];

    for (int i = nax; i < nbx; i++)
    {
      for (int j = nay; j < nby; j++)
      {
        b[i * ny + j] = Value(u[i * ny + j].ToVector(zb));
      }
    }

    // Compute straddling edges inside lower Oxy plane
    for (int i = nax; i < nbx - 1; i++)
    {
      for (int j = nay; j < nby; j++)
      {
        //   if (((b[i*ny + j] < 0.0) && (b[(i + 1)*ny + j] >= 0.0)) || ((b[i*ny + j] >= 0.0) && (b[(i + 1)*ny + j] < 0.0)))
        if (!((b[i * ny + j] < 0.0) == !(b[(i + 1) * ny + j] >= 0.0)))
        {
          vertex.append(Dichotomy(u[i * ny + j].ToVector(zb), u[(i + 1) * ny + j].ToVector(zb), b[i * ny + j], b[(i + 1) * ny + j], d[0], epsilon));
          normal.append(Normal(vertex.last()));
          ebx[i * ny + j] = nv;
          nv++;
        }
      }
    }

    for (int i = nax; i < nbx; i++)
    {
      for (int j = nay; j < nby - 1; j++)
      {
        // if (((b[i*ny + j] < 0.0) && (b[i*ny + (j + 1)] >= 0.0)) || ((b[i*ny + j] >= 0.0) && (b[i*ny + (j + 1)] < 0.0)))
        if (!((b[i * ny + j] < 0.0) == !(b[i * ny + (j + 1)] >= 0.0)))
        {
          vertex.append(Dichotomy(u[i * ny + j].ToVector(zb), u[i * ny + (j + 1)].ToVector(zb), b[i * ny + j], b[i * ny + (j + 1)], d[1], epsilon));
          normal.append(Normal(vertex.last()));
          eby[i * ny + j] = nv;
          nv++;
        }
      }
    }

    // Create vertical straddling edges
    for (int i = nax; i < nbx; i++)
    {
      for (int j = nay; j < nby; j++)
      {
        // if ((a[i*ny + j] < 0.0) && (b[i*ny + j] >= 0.0) || (a[i*ny + j] >= 0.0) && (b[i*ny + j] < 0.0))
        if (!((a[i * ny + j] < 0.0) == !(b[i * ny + j] >= 0.0)))
        {
          vertex.append(Dichotomy(u[i * ny + j].ToVector(za), u[i * ny + j].ToVector(zb), a[i * ny + j], b[i * ny + j], d[2], epsilon));
          normal.append(Normal(vertex.last()));
          ez[i * ny + j] = nv;
          nv++;
        }
      }
    }
    // Create mesh
    for (int i = nax; i < nbx - 1; i++)
    {
      for (int j = nay; j < nby - 1; j++)
      {
        int cubeindex = 0;
        if (a[i * ny + j] < 0.0)       cubeindex |= 1;
        if (a[(i + 1) * ny + j] < 0.0)   cubeindex |= 2;
        if (a[i * ny + j + 1] < 0.0)     cubeindex |= 4;
        if (a[(i + 1) * ny + j + 1] < 0.0) cubeindex |= 8;
        if (b[i * ny + j] < 0.0)       cubeindex |= 16;
        if (b[(i + 1) * ny + j] < 0.0)   cubeindex |= 32;
        if (b[i * ny + j + 1] < 0.0)     cubeindex |= 64;
        if (b[(i + 1) * ny + j + 1] < 0.0) cubeindex |= 128;

        // Cube is straddling the surface
        if ((cubeindex != 255) && (cubeindex != 0))
        {
          e[0] = eax[i * ny + j];
          e[1] = eax[i * ny + (j + 1)];
          e[2] = ebx[i * ny + j];
          e[3] = ebx[i * ny + (j + 1)];
          e[4] = eay[i * ny + j];
          e[5] = eay[(i + 1) * ny + j];
          e[6] = eby[i * ny + j];
          e[7] = eby[(i + 1) * ny + j];
          e[8] = ez[i * ny + j];
          e[9] = ez[(i + 1) * ny + j];
          e[10] = ez[i * ny + (j + 1)];
          e[11] = ez[(i + 1) * ny + (j + 1)];

          for (int h = 0; ScalarField::TriangleTable[cubeindex][h] != -1; h += 3)
          {
            triangle.append(e[ScalarField::TriangleTable[cubeindex][h + 0]]);
            triangle.append(e[ScalarField::TriangleTable[cubeindex][h + 1]]);
            triangle.append(e[ScalarField::TriangleTable[cubeindex][h + 2]]);
            nt += 3;
          }
        }
      }
    }

    Math::Swap(a, b);

    za += d[2];
    Math::Swap(eax, ebx);
    Math::Swap(eay, eby);
  }

  delete[]a;
  delete[]b;
  delete[]u;

  delete[]eax;
  delete[]eay;
  delete[]ebx;
  delete[]eby;
  delete[]ez;

  QVector<int> normals = triangle;

  g = Mesh(vertex, normal, triangle, normals);
}

/*!
\brief Colorize the vertexes of a mesh according to the analytic scalar field.

\param m %Mesh.
\param c Colored mesh.
*/
void AnalyticScalarField::Colorize(const Mesh& m, MeshColor& c) const
{
  c = MeshColor(m);
  for (int i = 0; i < c.Vertexes(); i++)
  {
    c.SetColor(i, GetMaterial(c.Vertex(i), c.Normal(i)));
  }
}

#include <stack>

class Octree
{
public:
  Box box; //!< Cubic domain.
  int depth; //!< Depth.
  int n; //!< Number of vertexes, set to 2<SUP>d</SUP>+1 where d is the depth.

  Octree(const Box& b, int l) :box(b), depth(l)
  {
    n = (1 << l) + 1;
  }
  static inline int64_t OctreeVertexIndex(int x, int y, int z, int l)
  {
    return int64_t(x) | (int64_t(y) << (l + 1)) | (int64_t(z) << (2 * l + 2));
  }

  static inline int64_t OctreeEdgeIndex(int x, int y, int z, int t, int l)
  {
    return OctreeVertexIndex(x, y, z, l) | ((int64_t(t) << (3 * l + 3)));
  }

  double EdgeLength() const {
    return (box[1][0] - box[0][0]) / double((1 << depth)); // Same as n-1
  }
};

class OctreeCell
{
public:
  Box box;
  int celllevel;
  int depth;
  int x, y, z;
public:
  OctreeCell() :box(Box(1.0)), depth(0), celllevel(0), x(0), y(0), z(0) {}
  OctreeCell(const Box& b, int u) :box(b), depth(u), celllevel(u) { x = 0; y = 0; z = 0; }
  OctreeCell(const Box& b, int u, int l) :box(b), depth(u), celllevel(l) { x = 0; y = 0; z = 0; }
  OctreeCell(const Box& b, int u, int l, int x, int y, int z) :box(b), depth(u), celllevel(l), x(x), y(y), z(z) {}
  ~OctreeCell() {}

  const int origin[12 * 3] =
  {
  0,0,0, 0,1,0, 0,0,1, 0,1,1, // Along x-axis
   0,0,0, 1,0,0, 0,0,1, 1,0,1, // y-axis
  0,0,0, 1,0,0, 0,1,0, 1,1,0, //  z-axis
  };

  int64_t VertexCode(int v) const
  {
    return Octree::OctreeVertexIndex((v & 1) ? x + (1 << celllevel) : x, (v & 2) ? y + (1 << celllevel) : y, (v & 4) ? z + (1 << celllevel) : z, depth);
  }
  int64_t CenterCode() const
  {
    return Octree::OctreeVertexIndex(x + (1 << (celllevel - 1)), y + (1 << (celllevel - 1)), z + (1 << (celllevel - 1)), depth);
  }
  inline int64_t EdgeCode(int e) const
  {
    return Octree::OctreeEdgeIndex(x + origin[e * 3], y + origin[e * 3 + 1], z + origin[e * 3 + 2], e / 4, depth); // Orientation is e/4 : x for the first 4, then y, then z
  }
  inline OctreeCell Sub(int i) const
  {
    return OctreeCell(box.Sub(i), depth, celllevel - 1, x + ((i & 1) ? 1 << (celllevel - 1) : 0), y + ((i & 2) ? 1 << (celllevel - 1) : 0), z + ((i & 4) ? 1 << (celllevel - 1) : 0));
  }
};

/*!
\brief
*/
double OctreeGetSharedDoubleIndex(const AnalyticScalarField& sf, std::map<int64_t, double>& vertexmap, int64_t key, const Vector& p)
{
  // have we already computed the dual point?
  auto iterator = vertexmap.find(key);
  if (iterator != vertexmap.end())
  {
    // just return the dual point index
    return iterator->second;
  }
  else
  {
    double value = sf.Value(p);
    vertexmap[key] = value;
    return value;
  }
}

/*!
\brief Compute the polygonal mesh approximating the implicit surface using an octree decomposition of space.

\param box Box defining the region that will be polygonized.
\param level Octree depth, the corresponding grid subdivision is 2<SUP>n</SUP>+1.
\param epsilon Epsilon value for computing vertices on straddling edges.
*/
Mesh AnalyticScalarField::PolygonizeOctree(const Box& box, int level, const double& epsilon) const
{
  QVector<Vector> vertex;
  QVector<int> indexes;

  Box clipped = box;
  clipped.SetCubic();
  Octree octree(clipped, level);

  const double length = octree.EdgeLength();

  std::stack<OctreeCell> cells;

  std::map<int64_t, int64_t> medges;
  std::map<int64_t, double> mvertexes;

  cells.push(OctreeCell(clipped, level));

  while (!cells.empty())
  {
    OctreeCell c = cells.top();
    cells.pop();

    // Terminal cell: generate mesh
    if (c.celllevel == 0)
    {
      if (Norm(c.box[0] - Vector(1, 2, 4)) < 0.1)
      {
        std::cout << "debug";

      }
      double x[8];
      // For all vertexes
      for (int i = 0; i < 8; i++)
      {
        // Compute or get field values
        x[i] = OctreeGetSharedDoubleIndex(*this, mvertexes, c.VertexCode(i), c.box.Vertex(i));
      }
      //  Determine the index into the edge table which tells us which vertices are inside of the surface
      int cubeindex = 0;
      if (x[0] < 0.0) cubeindex |= 1;
      if (x[1] < 0.0) cubeindex |= 2;
      if (x[2] < 0.0) cubeindex |= 4;
      if (x[3] < 0.0) cubeindex |= 8;
      if (x[4] < 0.0) cubeindex |= 16;
      if (x[5] < 0.0) cubeindex |= 32;
      if (x[6] < 0.0) cubeindex |= 64;
      if (x[7] < 0.0) cubeindex |= 128;

      // Cube is entirely in/out of the surface
      if (ScalarField::edgeTable[cubeindex] == 0)
        continue;

      // Vertex index on edges
      int e[12];

      // Process 12 edges
      for (int i = 0; i < 12; i++)
      {
        // Find the vertices where the surface intersects the cube 
        if (ScalarField::edgeTable[cubeindex] & (1 << i))
        {
          int64_t key = c.EdgeCode(i);
          auto iterator = medges.find(key);
          if (iterator != medges.end())
          {
            // just return the dual point index
            e[i] = iterator->second;
          }
          else
          {
            e[i] = vertex.size();

            // Insert vertex ID into map 
            medges[key] = e[i];

            // Compute vertex
            int a = Box::edge[i * 2];
            int b = Box::edge[i * 2 + 1];
            vertex.append(Dichotomy(c.box.Vertex(a), c.box.Vertex(b), x[a], x[b], length, epsilon));
          }
        }
      }

      // Create the triangles
      for (int i = 0; ScalarField::TriangleTable[cubeindex][i] != -1; i += 3)
      {
        indexes.append(e[ScalarField::TriangleTable[cubeindex][i + 0]]);
        indexes.append(e[ScalarField::TriangleTable[cubeindex][i + 1]]);
        indexes.append(e[ScalarField::TriangleTable[cubeindex][i + 2]]);
      }
    }
    // Not a terminal cell: check Lipschitz condition and recurse it cell is not empty
    else
    {
      // Compute or get field value at center
      double v = OctreeGetSharedDoubleIndex(*this, mvertexes, c.CenterCode(), c.box.Center());

      // On pourrait amÃ©liorer en ne pushant la valeur au centre qu'apres le test de non vide
      double k = K(c.box);
      if (fabs(v) > k * c.box.Radius())
      {
        continue;
      }

      // Recurse  
      for (int i = 0; i < 8; i++)
      {
        OctreeCell oci = c.Sub(i);
        cells.push(c.Sub(i));
      }
    }
  }
  // Set normals
  QVector<Vector> normals;
  QVector<int> normalindexes = indexes;
  normals.reserve(vertex.size());

  // Compute normal for all vertexes
  for (int i = 0; i < vertex.size(); i++)
  {
    normals.append(Normal(vertex.at(i)));
  }
  return Mesh(vertex, normals, indexes, normalindexes);
}


void AnalyticScalarField::PolygonizeOctree(int n, Mesh& g, const Box& box, const double& epsilon) const
{
  g = PolygonizeOctree(box, int(log2(n)), epsilon);
  std::cout << "nb subdivision : " << int(log2(n)) << ", n : " << n << ", nb_cell : " << (1 << int(log2(n))) << std::endl;
  return;
}