// Mesh
#include <QtCore/QFile>
#include <QtCore/QTextStream>

#include "libs/mesh.h"
#include "libs/curvepoint.h"
#include "libs/torus.h"
#include "libs/cone.h"
#include "libs/voxel.h"
#include "libs/rectangle.h"
#include "libs/octahedron.h"
#include "libs/pyramid.h"
#include "libs/sphereset.h"
#include "libs/dodecahedron.h"
#include "libs/icosahedron.h"
#include "libs/curve.h"
#include "libs/curveset.h"
#include "libs/cuboctahedron.h"
#include "libs/icosidodecahedron.h"
#include "libs/hexagon.h"
#include "libs/cuboid.h"
#include "libs/array.h"
#include "libs/quadricsurface.h"
#include "libs/arrayinteger.h"

#include <QtCore/QRegularExpression>

/*!
\class Mesh mesh.h

\brief Core triangle mesh class.
\ingroup ExtendedKernelGroup
*/

/*!
\brief Initialize the mesh to empty.
*/
Mesh::Mesh()
{
}

/*!
\brief Create a mesh from a set of vertices.

The set defines the 3n vertexes of n triangles. The array of triangle indexes is generated internally.
\param vertices List of vertices.
*/
Mesh::Mesh(const QVector<Vector>& vertices) :vertices(vertices)
{
  normals.resize(vertices.size());
  varray.resize(vertices.size());

  // Set vertex indexes
  for (int i = 0; i < varray.size(); i++)
  {
    varray[i] = i;
  }

  // Compute normals
  for (int i = 0; i < varray.size() / 3; i++)
  {
    Vector nt = GetTriangle(i).Normal();
    normals[3 * i] = nt;
    normals[3 * i + 1] = nt;
    normals[3 * i + 2] = nt;
  }

  // Set normal indexes
  narray = varray;
}

/*!
\brief Create a mesh from a set of vertices and a set of triangle indexes.

Indices must have a size multiple of three. Normals are computed to define flat triangles.

\param vertices List of vertices.
\param indices List of vertex indexes.
\param flat Boolean, set to true if face should be flat, set to false to generate smooth normals.
*/
Mesh::Mesh(const QVector<Vector>& vertices, const QVector<int>& indices, bool flat) :vertices(vertices), varray(indices)
{
  int nfaces = varray.size() / 3;
  narray.resize(indices.size());

  if (flat == true)
  {
    // Number of normals is equal to the number of faces
    normals.resize(nfaces);

    for (int i = 0; i < nfaces; i++)
    {
      Vector nt = GetTriangle(i).Normal();
      normals[i] = nt;
      narray[3 * i] = i;
      narray[3 * i + 1] = i;
      narray[3 * i + 2] = i;
    }
  }
  else
  {
    // Number of normals is equal to the number of vertexes
    normals.fill(Vector::Null, vertices.size());

    for (int i = 0; i < nfaces; i++)
    {
      Vector nt = GetTriangle(i).AreaNormal();
      normals[varray[3 * i]] += nt;
      normals[varray[3 * i + 1]] += nt;
      normals[varray[3 * i + 2]] += nt;
    }

    for (int i = 0; i < normals.size(); i++)
    {
      Normalize(normals[i]);
    }
    // The normal indexes are the same
    narray = indices;
  }
}

/*!
\brief Initialize the structure from a two dimensional mesh.

\param mesh Planar mesh.
*/
Mesh::Mesh(const Mesh2& mesh) :varray(mesh.Indexes())
{
  // The vertex index array was directly copied from the argument mesh
  vertices.reserve(mesh.VertexSize());

  for (int i = 0; i < mesh.VertexSize(); i++)
  {
    vertices.append(mesh.Vertex(i).ToVector(0.0));
  }

  // Normals
  normals.reserve(mesh.TriangleSize());
  narray.reserve(mesh.TriangleSize() * 3);

  for (int i = 0; i < mesh.TriangleSize(); i++)
  {
    normals.append(Vector::Z);

    narray.append(i);
    narray.append(i);
    narray.append(i);
  }
}

/*!
\brief Create the mesh.

\param vertices, normals Array of vertices and normals.
\param va, na Array of vertex and normal indexes.
*/
Mesh::Mesh(const QVector<Vector>& vertices, const QVector<Vector>& normals, const QVector<int>& va, const QVector<int>& na) :vertices(vertices), normals(normals), varray(va), narray(na)
{
}

/*!
\brief Reserve memory for arrays.
\param nv,nn,nvi,nvn Number of vertices, normals, vertex indexes and vertex normals.
*/
void Mesh::Reserve(int nv, int nn, int nvi, int nvn)
{
  vertices.reserve(nv);
  normals.reserve(nn);
  varray.reserve(nvi);
  narray.reserve(nvn);
}

/*!
\brief Empty
*/
Mesh::~Mesh()
{
}

/*!
\brief Generate a quadrangle.
\param quad %Quadrangle, should be flat.
*/
Mesh::Mesh(const Quadrangle& quad)
{
  vertices.resize(4);
  normals.resize(1);

  varray.resize(6);
  narray.resize(6);

  normals[0] = quad.Normal();

  for (int i = 0; i < 4; i++)
  {
    vertices[i] = quad.Vertex(i);
  }

  // Vertex indexes
  varray[0] = 0;
  varray[1] = 1;
  varray[2] = 2;

  varray[3] = 0;
  varray[4] = 2;
  varray[5] = 3;

  // Normal indexes
  for (int i = 0; i < 6; i++)
  {
    narray[i] = 0;
  }
}

/*!
\brief Generate a grid geometry over an input quadrangle.
\param quad %Quadrangle.
\param nx, ny Number of points in the subdivision.
*/
Mesh::Mesh(const Quadrangle& quad, int nx, int ny)
{
  vertices.resize(nx * ny);
  normals.resize(nx * ny);

  varray.resize(3 * (nx - 1) * (ny - 1) * 2);
  narray.resize(3 * (nx - 1) * (ny - 1) * 2);

  // Vertices and normals
  for (int j = 0; j < nx; j++)
  {
    double u = Math::Unit(j, nx);
    for (int i = 0; i < ny; i++)
    {
      double v = Math::Unit(i, ny);

      vertices[j * nx + i] = quad.Vertex(u, v);
      normals[j * nx + i] = quad.Normal(u, v);
    }
  }

  // Triangle indexes
  int n = 0;
  for (int j = 0; j < nx - 1; j++)
  {
    for (int i = 0; i < ny - 1; i++)
    {
      narray[n + 0] = varray[n + 0] = j * nx + i;
      narray[n + 1] = varray[n + 1] = j * nx + i + 1;
      narray[n + 2] = varray[n + 2] = (j + 1) * nx + i;

      n += 3;

      narray[n + 0] = varray[n + 0] = j * nx + i + 1;
      narray[n + 1] = varray[n + 1] = (j + 1) * nx + i + 1;
      narray[n + 2] = varray[n + 2] = (j + 1) * nx + i;

      n += 3;
    }
  }
}

/*!
\brief Generate a grid geometry over an input rectangle.
\param rect %Rectangle.
\param nx, ny Number of points in the subdivision.
*/
Mesh::Mesh(const Rectangles& rect, int nx, int ny)
{
  vertices.resize(nx * ny);
  normals.resize(nx * ny);

  varray.resize(3 * (nx - 1) * (ny - 1) * 2);
  narray.resize(3 * (nx - 1) * (ny - 1) * 2);

  // Vertices and normals
  for (int j = 0; j < nx; j++)
  {
    double u = Math::Unit(j, nx);
    for (int i = 0; i < ny; i++)
    {
      double v = Math::Unit(i, ny);

      vertices[j * nx + i] = rect.Vertex(u, v);
      normals[j * nx + i] = rect.Normal();
    }
  }

  // Triangle indexes
  int n = 0;
  for (int j = 0; j < nx - 1; j++)
  {
    for (int i = 0; i < ny - 1; i++)
    {
      narray[n + 0] = varray[n + 0] = j * nx + i;
      narray[n + 1] = varray[n + 1] = j * nx + i + 1;
      narray[n + 2] = varray[n + 2] = (j + 1) * nx + i;

      n += 3;

      narray[n + 0] = varray[n + 0] = j * nx + i + 1;
      narray[n + 1] = varray[n + 1] = (j + 1) * nx + i + 1;
      narray[n + 2] = varray[n + 2] = (j + 1) * nx + i;

      n += 3;
    }
  }
}

/*!
\brief Creates a capsule.

The mesh has 2n(n-1)+2 vertices, 2n(n-1)+2 normals and 4n(n-1) triangles.
\param c Capsule.

\param n Number of points on circle, the number of slices on a half sphere will be E(n/4).
\param m Number of slices on the cylinder part (should be 2 at least).
*/
Mesh::Mesh(const Capsule& c, int n, int m)
{
  const Vector a = c.Vertex(0);
  const Vector b = c.Vertex(1);
  const double r = c.Radius();

  Vector z = c.GetAxis();

  // Frame
  Vector x, y;
  z.Orthonormal(x, y);

  // Perimeter
  int p = n;

  // Sphere slices
  n = int(n / 4);

  // Size: two half spheres and k circles defining the cylinder

  // Number of points (and normals, overall mesh is not optimized)
  int s = p * (2 * n + m) + 2;

  vertices.resize(s);
  normals.resize(s);

  // Create set of vertices
  double dt = Math::TwoPi / p;
  double df = 0.5 * Math::Pi / n;

  int k = 0;
  // Phi
  double f = -Math::HalfPi;

  // First half sphere, does not include equator
  for (int j = 0; j < n; j++)
  {
    f += df;
    // Theta
    double t = 0.0;
    for (int i = 0; i < p; i++)
    {
      Vector u = cos(t) * cos(f) * x + sin(t) * cos(f) * y + sin(f) * z;
      normals[k] = u;
      vertices[k] = a + r * u;
      k++;
      t += dt;
    }
  }
  // Set of circles
  for (int j = 0; j < m; j++)
  {
    Vector cj = a + (b - a) * double(j) / double(m - 1);
    // Theta
    double t = 0.0;
    for (int i = 0; i < p; i++)
    {
      Vector u = cos(t) * x + sin(t) * y;
      normals[k] = u;
      vertices[k] = cj + r * u;
      k++;
      t += dt;
    }
  }

  f = 0.0;
  // Second half sphere, does not include equator
  for (int j = 0; j < n; j++)
  {
    f += df;
    // Theta
    double t = 0.0;
    for (int i = 0; i < p; i++)
    {
      Vector u = cos(t) * cos(f) * x + sin(t) * cos(f) * y + sin(f) * z;
      normals[k] = u;
      vertices[k] = b + r * u;
      k++;
      t += dt;
    }
  }

  // First apex
  normals[s - 2] = z;
  vertices[s - 2] = b + r * z;

  // Seconde apex
  normals[s - 1] = -z;
  vertices[s - 1] = a - r * z;

  // Reserve space for the smooth triangle array
  varray.reserve(p * (n - 1 + m - 2) + 2 * p);
  narray.reserve(p * (n - 1 + m - 2) + 2 * p);

  // Apexes

  for (int i = 0; i < p; i++)
  {
    AddSmoothTriangle(s - 1, s - 1, (i + 1) % p, (i + 1) % p, i, i);
  }

  for (int i = 0; i < p; i++)
  {
    AddSmoothTriangle(s - 2, s - 2, s - 2 - p + i, s - 2 - p + i, s - 2 - p + (i + 1) % p, s - 2 - p + (i + 1) % p);
  }

  // Cylinder part
  for (int j = 1; j < (2 * n + m); j++)
  {
    for (int i = 0; i < p; i++)
    {
      AddSmoothQuadrangle((j - 1) * p + i, (j - 1) * p + i, (j - 1) * p + (i + 1) % p, (j - 1) * p + (i + 1) % p, j * p + (i + 1) % p, j * p + (i + 1) % p, j * p + i, j * p + i);
    }
  }
}

/*!
\brief Compute the histogram of the aspect ratio of the triangles of the mesh.
\param n Number of classes.
*/
QVector<int> Mesh::AspectRatios(int n) const
{
  QVector<int> a(n, 0);
  for (int i = 0; i < varray.size() / 3; i++)
  {
    double r = GetTriangle(i).Aspect();

    int k = int(Math::Clamp(r, 0.0, 1.0 - 0.000001) * n);
    a[k]++;
  }

  return a;
}

/*!
\brief Compute the average edge length.
*/
double Mesh::AverageEdgeLength() const
{
  double length = 0.0;

  // For all edges of the triangle, sum lengths
  for (int i = 0; i < varray.size(); i += 3)
  {
    const Vector& a = vertices[varray[i + 0]];
    const Vector& b = vertices[varray[i + 1]];
    const Vector& c = vertices[varray[i + 2]];

    // Edge lengths
    length += Norm(b - a);
    length += Norm(c - b);
    length += Norm(a - c);
  }

  length /= 3.0 * varray.size();
  return length;
}

/*!
\brief Compute the range of edge lengths.
\param a, b Range.
*/
void Mesh::EdgeLengthRange(double& a, double& b) const
{
  a = b = 0.0;

  if (varray.size() == 0)
    return;

  a = b = Norm(vertices[varray[1]] - vertices[varray[0]]);

  // For all edges of the triangle
  for (int i = 0; i < varray.size(); i += 3)
  {
    const Vector& va = vertices[varray[i + 0]];
    const Vector& vb = vertices[varray[i + 1]];
    const Vector& vc = vertices[varray[i + 2]];

    // Edge lengths
    Math::SetMinMax(Norm(vb - va), a, b);
    Math::SetMinMax(Norm(vc - vb), a, b);
    Math::SetMinMax(Norm(va - vc), a, b);
  }
}

/*!
\brief Compute the next index of an edge in a triangle.
\param i Index of the starting vertex of an edge.
*/
inline int Mesh::NextIndex(int i) const
{
  return ((i % 3) == 2) ? i - 2 : i + 1;
}

/*!
\brief Compute the previous index of an edge in a triangle.
\param i Index of the starting vertex of an edge.
*/
inline int Mesh::PrevIndex(int i) const
{
  return ((i % 3) == 0) ? i + 2 : i - 1;
}

/*!
\brief Compute the base triangle index given an edge index.
\param i Edge index.
*/
inline int Mesh::BaseIndex(int i) const
{
  return i - (i % 3);
}


bool Mesh::FindEdge(int a, int b, int& k, int s) const
{
  std::cout << "Mesh::FindEdge(" << a << "," << b << ", at " << s << ")" << std::endl;
  for (k = s; k < varray.size(); k++)
  {
    if ((varray[k] == a) && (varray[NextIndex(k)] == b))
    {
      std::cout << " Found at " << k << "," << NextIndex(k) << std::endl;
      // Found opposite half edges 
      return true;
    }
  }
  std::cout << " Not found " << std::endl;
  return false;
}

/*!
\brief Check if a mesh is manifold.
*/
bool Mesh::Manifold() const
{
  int rsize = varray.size();

  int truc[20];
  // For all triangles
  for (int i = 0; i < rsize; i++)
  {
    int ia = varray[i];
    int ib = varray[NextIndex(i)];
    int n = 0;
    for (int j = 0; j < rsize; j++)
    {
      if ((varray[j] == ib) && (varray[NextIndex(j)] == ia))
      {
        // Found opposite half edges 
        truc[n] = j;
        n++;
      }
    }
    if (n != 1)
    {
      //std::cout << "n=" << n << std::endl;
      for (int k = 0; k < n; k++)
      {
        //std::cout << "  truc=" << truc[k] << std::endl;
      }
      return false;
    }
  }
  return true;
}

/*!
\brief Run a simple edge collapse algorithm over a triangle mesh.

This algorithm runs in O(n<SUP>2</SUP>) where n denotes the number
of triangles in the mesh.

\param epsilon Threshold length for collapsing an edge.
*/
void Mesh::EdgeCollapse(const double& epsilon)
{
  std::cout << "Mesh::EdgeCollapse" << std::endl;
  // std::cout << "      Manifold: " << Manifold()<<std::endl;
  const double e = epsilon * epsilon;
  int count = 0;

  // For all triangles
  for (int i = 0; i < varray.size(); i++)
  {
    //std::cout << "i ="<<i << std::endl;
    //std::cout << "   "<<varray.size() << std::endl;
    int ia = varray[i];
    int ib = varray[NextIndex(i)];
    // Edge
    Vector ab = vertices[ia] - vertices[ib];

    // Size of edge lower than threshold : collapse edge
    if (ab * ab < e)
    {
      std::cout << "Count: " << count++ << std::endl;

      EdgeCollapse(i);
    }
  }
  // Should be removed when finally works
  narray = varray;
}

void Mesh::EdgeCollapse(int i)
{
  // std::cout << *this;
  std::cout << "Mesh::EdgeCollapse(int=" << i << " )" << std::endl;
  std::cout << "   which means i= " << i << " and next(i)= " << NextIndex(i) << std::endl;

  int ia = varray[i];
  int ib = varray[NextIndex(i)];
  std::cout << "   which means ia= " << ia << " and ib= " << ib << std::endl;

  int k = 0;

  // Find the (possible) other face in the remaining faces that share the same edge
  bool vs = FindEdge(ib, ia, k, 0);

  if (vs == false)
  {
    std::cout << "NONNNNN" << std::endl;

    return;
  }

  //Math::Sort(ia, ib);

  // Take the minimum of indexes to avoid putting new vertex at the end of the array (which will be removed later)
  int ic = min(ia, ib);
  int id = max(ia, ib);

  // Collapse to center and change normal to average
  vertices[ic] = 0.5 * (vertices[ia] + vertices[ib]);
  normals[ic] = Normalized(normals[ia] + normals[ib]);

  // Copy last vertex to position d
  vertices[id] = vertices.last();
  normals[id] = normals.last();

  std::cout << "(i=" << i << ", k=" << k << ") ( size=" << varray.size() << ")  ";
  std::cout << "ic=" << ic << " id=" << id << "  size=" << vertices.size() << std::endl;
  vertices.removeLast();
  normals.removeLast();

  // Update vertex indexing since deleted vertex d has collapsed to c, and last vertex was moved to d
  // // d->c & n->d
  // Also update normals
  for (int j = 0; j < varray.size(); j++)
  {
    if (varray[j] == id)
    {
      varray[j] = ic;
      narray[j] = ic;
    }
    if (varray[j] == vertices.size())
    {
      varray[j] = id;
      narray[j] = id;
    }
  }

  // std::cout << *this << std::endl;


  std::cout << "Base of i " << BaseIndex(i) << " and of k " << BaseIndex(k) << std::endl;

  Sort(i, k);
  std::cout << " Sorted i " << BaseIndex(i) << " and of k " << BaseIndex(k) << std::endl;

  /*
  i = BaseIndex(i);
  k = BaseIndex(k);
 // Move two triangles
 varray[i] = varray[l];
  varray[i+1] = varray[l+1];
  varray[i+2] = varray[l+2];

  varray[k] = varray[l+3];
  varray[k + 1] = varray[l + 4];
  varray[k + 2] = varray[l + 5];


  for (int ii = 0; ii < 6; ii++)
  {
    varray.removeLast();
    narray.removeLast();
  }
 // varray.resize(varray.size()-6);
 // narray.resize(narray.size()-6);
  */


  if (i < varray.size() - 6)
  {
    // Last two elements can be moved and deleted
    if (k < varray.size() - 6)
    {
      std::cout << "Delete type 1" << std::endl;
      // Warning: swap order to preserve order
      varray[NextIndex(i)] = varray.takeLast();
      varray[i] = varray.takeLast();
      varray[PrevIndex(i)] = varray.takeLast();

      varray[NextIndex(k)] = varray.takeLast();
      varray[k] = varray.takeLast();
      varray[PrevIndex(k)] = varray.takeLast();

      narray[NextIndex(i)] = narray.takeLast();
      narray[i] = narray.takeLast();
      narray[PrevIndex(i)] = narray.takeLast();

      narray[NextIndex(k)] = narray.takeLast();
      narray[k] = narray.takeLast();
      narray[PrevIndex(k)] = narray.takeLast();
    }
    else
    {
      // Move last, remove the other
      if (k < varray.size() - 3)
      {
        std::cout << "Delete type 2 A" << std::endl;
        // Warning: swap order to preserve order
        varray[NextIndex(i)] = varray.takeLast();
        varray[i] = varray.takeLast();
        varray[PrevIndex(i)] = varray.takeLast();

        varray.removeLast();
        varray.removeLast();
        varray.removeLast();

        narray[NextIndex(i)] = narray.takeLast();
        narray[i] = narray.takeLast();
        narray[PrevIndex(i)] = narray.takeLast();

        narray.removeLast();
        narray.removeLast();
        narray.removeLast();

      }
      else
      {
        std::cout << "Delete type 2 B" << std::endl;
        varray.removeLast();
        varray.removeLast();
        varray.removeLast();

        // Warning: swap order to preserve order
        varray[NextIndex(i)] = varray.takeLast();
        varray[i] = varray.takeLast();
        varray[PrevIndex(i)] = varray.takeLast();

        narray.removeLast();
        narray.removeLast();
        narray.removeLast();

        // Warning: swap order to preserve order
        narray[NextIndex(i)] = narray.takeLast();
        narray[i] = narray.takeLast();
        narray[PrevIndex(i)] = narray.takeLast();
      }
    }
  }
  // Elements to be removed are the last two
  else
  {
    std::cout << "Delete type 3" << std::endl;
    varray.removeLast();
    varray.removeLast();
    varray.removeLast();
    varray.removeLast();
    varray.removeLast();
    varray.removeLast();


    narray.removeLast();
    narray.removeLast();
    narray.removeLast();
    narray.removeLast();
    narray.removeLast();
    narray.removeLast();
  }


  // std::cout << *this << std::endl;
   /*
   if (Manifold() == false)
   {
     std::cout << "      Manifold? NO ! " << std::endl;
     getchar();
   }
   */
}

/*!
\brief Smooth the normals of the mesh.

This function weights the normals of the faces by their corresponding area.
\sa Triangle::AreaNormal()
*/
void Mesh::SmoothNormals()
{
  // Set size
  normals.resize(vertices.size());

  // Initialize 
  normals.fill(Vector(0.0), vertices.size());

  narray = varray;

  // Sum normals
  for (int i = 0; i < varray.size(); i += 3)
  {
    Vector tn = Triangle(vertices[varray.at(i)], vertices[varray.at(i + 1)], vertices[varray.at(i + 2)]).AreaNormal();
    normals[narray[i + 0]] += tn;
    normals[narray[i + 1]] += tn;
    normals[narray[i + 2]] += tn;
  }

  // Normalize 
  for (int i = 0; i < normals.size(); i++)
  {
    Normalize(normals[i]);
  }
}

/*!
\brief Check if the i<sup>th</sup> triangle in the mesh is a flat or smooth triangle.

This is performed by comparing the indexes of the normals, not the normals themselves.
Therefore, should the normal indexes be different but the normal vector equal, the
function will return that the triangle is a smooth one anyway.
\param i Index of the triangle.
*/
int Mesh::IsSmooth(int i) const
{
  return !((narray[i * 3] == narray[i * 3 + 1]) && (narray[i * 3] == narray[i * 3 + 2]));
}

/*!
\brief Add a smooth triangle to the geometry.
\param a, b, c Indexes of the vertices.
\param na, nb, nc Indexes of the normals.
*/
void Mesh::AddSmoothTriangle(int a, int na, int b, int nb, int c, int nc)
{
  varray.append(a);
  narray.append(na);
  varray.append(b);
  narray.append(nb);
  varray.append(c);
  narray.append(nc);
}
/*!
\brief Add a triangle to the geometry.
\param a, b, c Indexes of the vertices.
\param n Index of the normal.
*/
void Mesh::AddTriangle(int a, int b, int c, int n)
{
  varray.append(a);
  narray.append(n);
  varray.append(b);
  narray.append(n);
  varray.append(c);
  narray.append(n);
}

/*!
\brief Add a grid structure.
\param nx, ny Grid size.
*/
void Mesh::AddArray(int nx, int ny)
{
  for (int i = 0; i < nx - 1; i++)
  {
    for (int j = 0; j < ny - 1; j++)
    {
      AddQuadrangle(i * ny + j, i * ny + j + 1, (i + 1) * ny + j + 1, (i + 1) * ny + j);
    }
  }

}

/*!
\brief Add a smmoth quadrangle to the geometry.

Creates two smooth triangles abc and acd.

\param a, b, c, d  Indexes of the vertexes.
\param na, nb, nc, nd Indexes of the normals for the vertexes.
*/
void Mesh::AddSmoothQuadrangle(int a, int na, int b, int nb, int c, int nc, int d, int nd)
{
  // First triangle
  AddSmoothTriangle(a, na, b, nb, c, nc);

  // Second triangle
  AddSmoothTriangle(a, na, c, nc, d, nd);
}

/*!
\brief Add a quadrangle to the geometry.

\param a, b, c, d  Indexes of the vertices and normals.
*/
void Mesh::AddQuadrangle(int a, int b, int c, int d)
{
  AddSmoothQuadrangle(a, a, b, b, c, c, d, d);
}

/*!
\brief Add a quadrangle to the geometry.

\param a, b, c, d  Indexes of the vertexes.
\param n  Index of the (shared) normal.
*/
void Mesh::AddQuadrangle(int a, int b, int c, int d, int n)
{
  AddSmoothQuadrangle(a, n, b, n, c, n, d, n);
}

/*!
\brief Add a pentagon to the geometry.

It creates three triangles abc, acd and ade.
\param a, b, c, d, e  Indexes of the vertices.
\param n Index of the normal shared by all vertices.
*/
void Mesh::AddPentagon(int a, int b, int c, int d, int e, int n)
{
  // First triangle
  AddTriangle(a, b, c, n);

  // Second triangle
  AddTriangle(a, c, d, n);

  // Third triangle
  AddTriangle(a, d, e, n);
}

/*!
\brief Transforms the mesh given a transformation.
\param frame Transformation.
*/
void Mesh::Transform(const Frame& frame)
{
  // Apply transformation to vertices
  for (int i = 0; i < vertices.size(); i++)
  {
    vertices[i] = frame.Transform(vertices[i]);
  }

  // Compute normal, no need to normalize
  for (int i = 0; i < normals.size(); i++)
  {
    normals[i] = frame.TransformDirection(normals[i]);
  }
}

/*!
\brief Transforms the mesh given a transformation.
\param frame Transformation.
*/
void Mesh::Transform(const FrameScaled& frame)
{
  // Apply transformation to vertices
  for (int i = 0; i < vertices.size(); i++)
  {
    vertices[i] = frame.Transform(vertices[i]);
  }

  // Compute normal
  for (int i = 0; i < normals.size(); i++)
  {
    normals[i] = Normalized(frame.TransformDirection(normals[i]));
  }
}

/*!
\brief Transforms the mesh given a transformation.
\param frame Transformation.
*/
Mesh Mesh::Transformed(const FrameScaled& frame) const
{
  QVector<Vector> vi = vertices;
  QVector<Vector> ni = normals;

  // Apply transformation to vertices
  for (int i = 0; i < vertices.size(); i++)
  {
    vi[i] = frame.Transform(vertices[i]);
  }

  // Compute normal
  for (int i = 0; i < normals.size(); i++)
  {
    ni[i] = Normalized(frame.TransformDirection(normals[i]));
  }
  return Mesh(vi, ni, varray, narray);
}

/*!
\brief Inverse transform.
\sa Mesh::Transform()
\param frame Transformation.
*/
void Mesh::InverseTransform(const FrameScaled& frame)
{
  // Apply transformation to vertices
  for (int i = 0; i < vertices.size(); i++)
  {
    vertices[i] = frame.InverseTransform(vertices[i]);
  }

  // Compute normal
  for (int i = 0; i < normals.size(); i++)
  {
    normals[i] = Normalized(frame.InverseTransformDirection(normals[i]));
  }
}

/*!
\brief Compute the bounding box of the object.
*/
Box Mesh::GetBox() const
{
  if (vertices.size() == 0)
  {
    return Box::Null;
  }
  return Box(vertices);
}

/*!
\brief Compute the bounding sphere of the object.
*/
Sphere Mesh::GetSphere() const
{
  if (vertices.size() == 0)
  {
    return Sphere::Null;
  }
  return Sphere(vertices);
}

/*!
\brief Compute the intersection between a ray and the mesh.

Very computationally intensive function: compute all the ray-triangle intersections.

\param ray The ray.
\param depth Intersection depth.
\param index Index of the intersected triangle.
*/
bool Mesh::Intersect(const Ray& ray, double& depth, int& index) const
{
  depth = 0.0;
  index = -1;

  bool check = false;

  for (int i = 0; i < varray.size() / 3; i++)
  {
    // Temporary depth
    double t;
    // Edge
    if (GetTriangle(i).Intersect(ray, t))
    {
      if (t > 0.0)
      {
        if (check == false)
        {
          check = true;
          depth = t;
          index = i;
        }
        else
        {
          if (t < depth)
          {
            depth = t;
            index = i;
          }
        }
      }
    }
  }

  return check;
}

/*!
\brief Compute the number of intersections between a ray and the mesh.

Very computationally intensive function: compute all the ray-triangle intersections.

\param ray The ray.
*/
int Mesh::Intersections(const Ray& ray) const
{
  const int n = varray.size() / 3;
  int c = 0;

  for (int i = 0; i < n; i++)
  {
    // Temporary depth
    double t;
    // Edge
    if (GetTriangle(i).Intersect(ray, t))
    {
      if (t > 0.0)
      {
        c++;
      }
    }
  }

  return c;
}

/*!
\brief Check if a point is inside or outside of the object.

This is a brute force algorithm, that computes all the intersections between the ray and the triangles of the object.

The mesh should be closed, i.e. without holes.
\param p Point.
*/
bool Mesh::Inside(const Vector& p) const
{
  // Ray used for intersecting the triangles
  const Ray ray(p, Vector(1.0, 0.0, 0.0));
  bool check = false;

  for (int i = 0; i < varray.size() / 3; i++)
  {
    // Temporary depth
    double t;
    // Edge
    if (GetTriangle(i).Intersect(ray, t))
    {
      if (t >= 0.0) check = !check;
    }
  }

  return check;
}


/*!
\brief Creates the mesh of an octaedron.
\param octa %Octahedron.
*/
Mesh::Mesh(const Octahedron& octa) :vertices(6), normals(8)
{
  // Vertices
  for (int i = 0; i < 6; i++)
  {
    vertices[i] = octa.Vertex(i);
  }

  // Normals
  for (int i = 0; i < 8; i++)
  {
    normals[i] = octa.Normal(i);
  }

  for (int i = 0; i < 8; i++)
  {
    AddSmoothTriangle(Octahedron::VertexIndex(i, 0), i, Octahedron::VertexIndex(i, 1), i, Octahedron::VertexIndex(i, 2), i);
  }
}


/*!
\brief Creates the mesh of an octaedron.

\param pyramid %Pyramid.
*/
Mesh::Mesh(const Pyramid& pyramid) :vertices(5), normals(6)
{
  // Vertices
  for (int i = 0; i < 5; i++)
  {
    vertices[i] = pyramid.Vertex(i);
  }

  // Normals
  // Faces #0 and #1 share the same normal, so skip #0
  for (int i = 1; i < 6; i++)
  {
    normals[i - 1] = pyramid.Normal(i);
  }

  // Triangles
  AddTriangle(Pyramid::VertexIndex(0, 0), Pyramid::VertexIndex(0, 1), Pyramid::VertexIndex(0, 2), 0);
  for (int i = 1; i < 6; i++)
  {
    AddTriangle(Pyramid::VertexIndex(i, 0), Pyramid::VertexIndex(i, 1), Pyramid::VertexIndex(i, 2), i - 1);
  }
}

/*!
\brief Creates an axis aligned box.

The object has 8 vertices, 6 normals and 12 triangles.
\param box The box.
*/
Mesh::Mesh(const Box& box) :vertices(8), normals(6)
{
  // Vertices
  for (int i = 0; i < 8; i++)
  {
    vertices[i] = box.Vertex(i);
  }

  // Normals
  normals[0] = Vector(-1, 0, 0);
  normals[1] = Vector(1, 0, 0);
  normals[2] = Vector(0, -1, 0);
  normals[3] = Vector(0, 1, 0);
  normals[4] = Vector(0, 0, -1);
  normals[5] = Vector(0, 0, 1);

  // Reserve space for the triangle array
  varray.reserve(12 * 3);
  narray.reserve(12 * 3);

  AddTriangle(0, 2, 1, 4);
  AddTriangle(1, 2, 3, 4);

  AddTriangle(4, 5, 6, 5);
  AddTriangle(5, 7, 6, 5);

  AddTriangle(0, 4, 2, 0);
  AddTriangle(4, 6, 2, 0);

  AddTriangle(1, 3, 5, 1);
  AddTriangle(3, 7, 5, 1);

  AddTriangle(0, 1, 5, 2);
  AddTriangle(0, 5, 4, 2);

  AddTriangle(3, 2, 7, 3);
  AddTriangle(6, 7, 2, 3);
}

/*!
\brief Creates a triangle.

\param t The triangle.
*/
Mesh::Mesh(const Triangle& t) :vertices(3), normals(1)
{
  // Vertices
  for (int i = 0; i < 3; i++)
  {
    vertices[i] = t[i];
  }

  // Normals
  normals[0] = t.Normal();

  // Reserve space for the triangle array
  varray.reserve(1 * 3);
  narray.reserve(1 * 3);

  AddTriangle(0, 1, 2, 0);
}

/*!
\brief Extrude a contour along a point curve.
\param curve %Curve.
\param profile %Contour.
*/
Mesh Mesh::Extrusion(const PointCurve& curve, const PointCurve2& profile)
{
  Mesh mesh;

  mesh.Reserve(curve.Size() * profile.Size(), curve.Size() * profile.Size(), (curve.Size() - 1 * profile.Size() - 1) * 2 * 3, (curve.Size() - 1 * profile.Size() - 1) * 2 * 3);

  for (int i = 0; i < curve.Size(); i++)
  {
    PointCurve p = profile.Transform(Frame(Matrix::Frenet(curve.Tangent(i)), curve.At(i)));
    for (int j = 0; j < profile.Size(); j++)
    {
      mesh.vertices.append(p[j]);
    }

    Vector ti = curve.Tangent(i);
    for (int j = 0; j < profile.Size(); j++)
    {
      Vector tj = p.Tangent(j);
      mesh.normals.append(Normalized(tj / ti));
    }
  }

  int nc = profile.Size();

  for (int i = 0; i < curve.Size() - 1; i++)
  {
    for (int j = 0; j < nc - 1; j++)
    {
      mesh.AddQuadrangle(i * nc + j, i * nc + j + 1, i * nc + j + nc + 1, i * nc + j + nc);
    }
  }
  return mesh;
}

/*!
\brief Extrude a set of contour along a point curve.

The curve and the set of contours should have the same size.

\param curve% Curve.
\param profile %Contour.
*/
Mesh Mesh::Extrusion(const PointCurve& curve, const QVector<PointCurve2>& profile)
{
  Mesh mesh;

  int nc = profile.at(0).Size();

  mesh.Reserve(curve.Size() * nc, curve.Size() * nc, (curve.Size() - 1 * nc - 1) * 2 * 3, (curve.Size() - 1 * nc - 1) * 2 * 3);

  for (int i = 0; i < curve.Size(); i++)
  {
    PointCurve p = profile.at(i).Transform(Frame(Matrix::Frenet(curve.Tangent(i)), curve.At(i)));
    for (int j = 0; j < nc; j++)
    {
      mesh.vertices.append(p[j]);
    }

    Vector ti = curve.Tangent(i);
    for (int j = 0; j < nc; j++)
    {
      Vector tj = p.Tangent(j);
      mesh.normals.append(Normalized(tj / ti));
    }
  }

  for (int i = 1; i < curve.Size() - 2; i++)
  {
    for (int j = 0; j < nc - 1; j++)
    {
      mesh.AddQuadrangle(i * nc + j, i * nc + j + 1, i * nc + j + nc + 1, i * nc + j + nc);
    }
  }
  return mesh;
}



/*!
\brief Import a mesh from an .obj file.
\param filename File name.
*/
void Mesh::Load(const QString& filename)
{
  vertices.clear();
  normals.clear();
  varray.clear();
  narray.clear();

  QFile data(filename);

  if (!data.open(QFile::ReadOnly))
    return;
  QTextStream in(&data);

  // Set of regular expressions : Vertex, Normal, Triangle
  QRegularExpression rexv("v\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)");
  QRegularExpression rexn("vn\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)");
  QRegularExpression rext("f\\s*(\\d*)/\\d*/(\\d*)\\s*(\\d*)/\\d*/(\\d*)\\s*(\\d*)/\\d*/(\\d*)");
  while (!in.atEnd())
  {
    QString line = in.readLine();
    QRegularExpressionMatch match = rexv.match(line);
    QRegularExpressionMatch matchN = rexn.match(line);
    QRegularExpressionMatch matchT = rext.match(line);
    if (match.hasMatch())//rexv.indexIn(line, 0) > -1)
    {
      Vector q = Vector(match.captured(1).toDouble(), match.captured(2).toDouble(), match.captured(3).toDouble()); vertices.append(q);
    }
    else if (matchN.hasMatch())//rexn.indexIn(line, 0) > -1)
    {
      Vector q = Vector(matchN.captured(1).toDouble(), matchN.captured(2).toDouble(), matchN.captured(3).toDouble());  normals.append(q);
    }
    else if (matchT.hasMatch())//rext.indexIn(line, 0) > -1)
    {
      varray.append(matchT.captured(1).toInt() - 1);
      varray.append(matchT.captured(3).toInt() - 1);
      varray.append(matchT.captured(5).toInt() - 1);
      narray.append(matchT.captured(2).toInt() - 1);
      narray.append(matchT.captured(4).toInt() - 1);
      narray.append(matchT.captured(6).toInt() - 1);
    }
  }
  data.close();
}

/*!
\brief Creates a sphere, using polar coordinates.

The object has (2 n) (n - 1) + 2 vertices, (2 n) (n - 1) + 2 normals and 4n(n-1) triangles.

\param sphere %Sphere.
\param n Discretization parameter.
*/
Mesh::Mesh(const Sphere& sphere, int n)
{
  Vector c = sphere.Center();
  double r = sphere.Radius();

  // Perimeter
  int p = 2 * n;

  // Size
  int s = (2 * n) * (n - 1) + 2;

  vertices.resize(s);
  normals.resize(s);

  // Create set of vertices
  double dt = Math::Pi / n;
  double df = Math::Pi / n;

  int k = 0;

  // Phi
  double f = -Math::HalfPi;

  for (int j = 1; j < n; j++)
  {
    f += df;
    // Theta
    double t = 0.0;
    for (int i = 0; i < 2 * n; i++)
    {
      Vector u = Vector(cos(t) * cos(f), sin(t) * cos(f), sin(f));
      normals[k] = u;
      vertices[k] = c + r * u;
      k++;
      t += dt;
    }
  }
  // North pole
  normals[s - 2] = Vector(0.0, 0.0, 1.0);
  vertices[s - 2] = c + Vector(0.0, 0.0, r);

  // South
  normals[s - 1] = Vector(0.0, 0.0, -1.0);
  vertices[s - 1] = c - Vector(0.0, 0.0, r);

  // Reserve space for the smooth triangle array
  varray.reserve(4 * n * (n - 1) * 3);
  narray.reserve(4 * n * (n - 1) * 3);

  // South cap
  for (int i = 0; i < 2 * n; i++)
    AddSmoothTriangle(s - 1, s - 1, (i + 1) % p, (i + 1) % p, i, i);

  // North cap
  for (int i = 0; i < 2 * n; i++)
    AddSmoothTriangle(s - 2, s - 2, 2 * n * (n - 2) + i, 2 * n * (n - 2) + i, 2 * n * (n - 2) + (i + 1) % p, 2 * n * (n - 2) + (i + 1) % p);

  // Sphere
  for (int j = 1; j < n - 1; j++)
  {
    for (int i = 0; i < 2 * n; i++)
    {
      AddSmoothQuadrangle((j - 1) * p + i, (j - 1) * p + i, (j - 1) * p + (i + 1) % p, (j - 1) * p + (i + 1) % p, j * p + (i + 1) % p, j * p + (i + 1) % p, j * p + i, j * p + i);
    }
  }
}


/*!
\brief Creates a set of spheres.

\param set %Set of spheres.
\param n Discretization parameter.
*/
Mesh::Mesh(const SphereSet& set, int n) :Mesh()
{
  for (int i = 0; i < set.Size(); i++)
  {
    Sphere s = set.GetSphere(i);
    Append(Mesh(s, n));
  }
}

/*!
\brief Creates a cylinder.

The cylinder with caps has 2n+2 vertices, n+2 normals and 4n triangles.

\param cylinder %Cylinder.
\param ca, cb Booleans defining whether the cylinder should have caps.
\param n Discretization.
*/
Mesh::Mesh(const Cylinder& cylinder, int n, bool ca, bool cb)
{
  const Vector a = cylinder.Vertex(0);
  const Vector b = cylinder.Vertex(1);
  const double r = cylinder.Radius();

  Vector z = Normalized(b - a);
  Vector x, y;
  z.Orthonormal(x, y);

  // Set dimension of array of vertices and normals
  vertices.resize(2 * n + (ca ? 1 : 0) + (cb ? 1 : 0));
  normals.resize(n + (ca ? 1 : 0) + (cb ? 1 : 0));

  // Vertices and normals
  for (int i = 0; i < n; i++)
  {
    double t = Math::Angle(i, n);
    Vector u = x * cos(t) + y * sin(t);

    vertices[i] = a + r * u;
    vertices[n + i] = b + r * u;
    normals[i] = u;
  }

  // Indexes for apex vertices and normals
  int iav = 2 * n;
  int ibv = 2 * n;
  int ian = n;
  int ibn = n;
  // caps
  if (ca)
  {
    vertices[iav] = a;
    normals[ian] = -z;
  }

  if (cb)
  {
    if (ca)
    {
      ibv = iav + 1;
      ibn = ian + 1;
    }
    vertices[ibv] = b;
    normals[ibn] = z;
  }

  // Reserve space for the triangle and smooth triangle arrays
  varray.reserve(4 * n * 3);
  narray.reserve(4 * n * 3);

  // Caps
  if (ca)
  {
    for (int i = 0; i < n; i++)
    {
      AddTriangle(iav, (i + 1) % n, i, ian);
    }
  }
  if (cb)
  {
    for (int i = 0; i < n; i++)
    {
      AddTriangle(ibv, n + i, n + (i + 1) % n, ibn);
    }
  }
  // Side
  for (int i = 0; i < n; i++)
  {
    AddSmoothTriangle(i, i, (i + 1) % n, (i + 1) % n, n + i, i);
    AddSmoothTriangle((i + 1) % n, (i + 1) % n, n + (i + 1) % n, (i + 1) % n, n + i, i);
  }
}

/*!
\brief Creates a disc.

The disc with has n+1 vertices.

\param disc %Disc.
\param n Discretization.
*/
Mesh::Mesh(const Disc& disc, int n)
{
  Vector z = Normalized(disc.Axis());
  Vector x, y;
  z.Orthonormal(x, y);

  const Vector c = disc.Center();
  const double r = disc.Radius();

  // Set dimension of array of vertices and normals
  vertices.resize(n + 1);

  // Vertices 
  for (int i = 0; i < n; i++)
  {
    double t = Math::Angle(i, n);
    Vector u = x * cos(t) + y * sin(t);

    vertices[i] = c + r * u;
  }
  vertices[n] = c;

  normals.append(z);

  // Reserve space for the triangle array
  varray.reserve(n * 3);
  narray.reserve(n * 3);

  for (int i = 0; i < n; i++)
  {
    AddTriangle(n, (i + 1) % n, i, 0);
  }
}



/*!
\brief This function creates the mesh of a dodecahedron.

The object has 20 vertices, 12 normals and 36 triangles
\param dodecahedron %Dodecahedron.
*/
Mesh::Mesh(const Dodecahedron& dodecahedron) :vertices(20), normals(12)
{
  // Vertices
  for (int i = 0; i < 20; i++)
  {
    vertices[i] = dodecahedron.Vertex(i);
  }
  // Normals
  for (int i = 0; i < 12; i++)
  {
    normals[i] = dodecahedron.Normal(i);
  }

  // Triangles
  for (int i = 0; i < 12; i++)
  {
    for (int j = 1; j < 4; j++)
    {
      AddTriangle(Dodecahedron::VertexIndex(i, 0), Dodecahedron::VertexIndex(i, j), Dodecahedron::VertexIndex(i, j + 1), i);
    }
  }
}



/*!
\brief Create the mesh of an icosahedron.
\param icosahedron %Icosahedron.
*/
Mesh::Mesh(const Icosahedron& icosahedron) :vertices(12), normals(20)
{
  for (int i = 0; i < 12; i++)
  {
    vertices[i] = icosahedron.Vertex(i);
  }

  for (int i = 0; i < 20; i++)
  {
    normals[i] = icosahedron.Normal(i);
  }

  for (int i = 0; i < 20; i++)
  {
    AddTriangle(Icosahedron::VertexIndex(i, 0), Icosahedron::VertexIndex(i, 1), Icosahedron::VertexIndex(i, 2), i);
  }
}

/*!
\brief Creates the mesh of a cone.

A truncated cone has 2n+2 vertices, n+2 normals and 4n triangles.

A real cone has n+2 vertices, n+2 normals and 2n triangles.

\param cone The cone.
\param n Discretization parameter.
*/
Mesh::Mesh(const Cone& cone, int n)
{
  double ra = cone.Radius(0);
  double rb = cone.Radius(1);
  Vector a = cone.Vertex(0);
  Vector b = cone.Vertex(1);

  Vector z = Normalized(b - a);
  Vector x, y;
  z.Orthonormal(x, y);

  if ((ra != 0.0) && (rb != 0.0))
  {
    // Set dimension of array of vertices and normals
    vertices.resize(2 * n + 2);
    normals.resize(n + 2);

    // Vertices and normals
    for (int i = 0; i < n; i++)
    {
      double t = Math::Angle(i, n);
      Vector u = x * cos(t) + y * sin(t);

      vertices[i] = a + ra * u;
      vertices[n + i] = b + rb * u;
      normals[i] = Normalized((ra - rb) * z + Norm(b - a) * u);
    }
    vertices[2 * n] = a;
    normals[n] = -z;
    vertices[2 * n + 1] = b;
    normals[n + 1] = z;

    // Space for the triangle and smooth triangle arrays
    varray.reserve(4 * n * 3);
    narray.reserve((n + 2) * 3);

    // Caps
    for (int i = 0; i < n; i++)
    {
      AddTriangle(2 * n, (i + 1) % n, i, n);
    }

    for (int i = 0; i < n; i++)
    {
      AddTriangle(2 * n + 1, n + i, n + (i + 1) % n, n + 1);
    }
    // Side
    for (int i = 0; i < n; i++)
    {
      AddSmoothTriangle(i, i, (i + 1) % n, (i + 1) % n, n + i, i);
      AddSmoothTriangle((i + 1) % n, (i + 1) % n, n + (i + 1) % n, (i + 1) % n, n + i, i);
    }

  }
  else
  {
    // Swap so that vertex b should be the apex
    if (ra == 0.0)
    {
      Swap(a, b);
      Math::Swap(ra, rb);
    }
    // Set dimension of array of vertices and normals
    vertices.resize(n + 2);
    normals.resize(n + 2);

    // Vertices and normals
    for (int i = 0; i < n; i++)
    {
      double t = Math::Angle(i, n);
      Vector u = x * cos(t) + y * sin(t);

      vertices[i] = a + ra * u;
      normals[i] = Normalized(ra * z + Norm(b - a) * u);
    }
    vertices[n] = a;
    normals[n] = -z;

    vertices[n + 1] = b;
    normals[n + 1] = z;

    // Space for the triangle and smooth triangle arrays
    varray.reserve(2 * n * 3);
    narray.reserve((n + 1) * 3);

    // Caps
    for (int i = 0; i < n; i++)
    {
      AddTriangle(n, (i + 1) % n, i, n);
    }

    // Side
    for (int i = 0; i < n; i++)
    {
      AddSmoothTriangle(i, i, (i + 1) % n, (i + 1) % n, n + 1, n + 1);
    }

  }
}

/*!
\brief Apply a planar symmetry to the object.
\param plane Plane.
*/
void Mesh::Symmetry(const Plane& plane)
{
  // Vertexes
  for (int i = 0; i < vertices.size(); i++)
  {
    vertices[i] = plane.Symmetry(vertices[i]);
  }

  // Normals
  for (int i = 0; i < normals.size(); i++)
  {
    normals[i] = plane.Reflect(normals[i]);
  }

  // Index should be swapped
  for (int i = 0; i < varray.size(); i += 3)
  {
    int t = varray[i];
    varray[i] = varray[i + 1];
    varray[i + 1] = t;
  }
}

/*!
\brief Rotate the mesh.
\param r Rotation matrix.
*/
void Mesh::Rotate(const Matrix& r)
{
  // Vertexes
  for (int i = 0; i < vertices.size(); i++)
  {
    vertices[i] = r * vertices[i];
  }
  // Normals
  for (int i = 0; i < normals.size(); i++)
  {
    normals[i] = r * normals[i];
  }
}

/*!
\brief Translate the mesh.
\param t Translation vector.
*/
void Mesh::Translate(const Vector& t)
{
  // Vertexes
  for (int i = 0; i < vertices.size(); i++)
  {
    vertices[i] += t;
  }
}

/*!
\brief Scale the mesh.
\param s Scaling factor.
*/
void Mesh::Scale(const double& s)
{
  // Vertexes
  for (int i = 0; i < vertices.size(); i++)
  {
    vertices[i] *= s;
  }

  if (s < 0.0)
  {
    // Normals
    for (int i = 0; i < normals.size(); i++)
    {
      normals[i] = -normals[i];
    }
  }
}

/*!
\brief Scale the mesh.
\param s Scaling vector, terms should be >0.
*/
void Mesh::Scale(const Vector& s)
{
  // Vertexes
  for (int i = 0; i < vertices.size(); i++)
  {
    vertices[i] *= s;
  }
}

/*!
\brief Save the mesh in .obj format, with vertices and normals.
\param url Filename.
\param meshName %Mesh name in .obj file.
*/
void Mesh::SaveObj(const QString& url, const QString& meshName) const
{
  QFile data(url);
  if (!data.open(QFile::WriteOnly))
    return;
  QTextStream out(&data);
  out << "# From Core::Mesh::Save " << "\n";
  out << "# " << vertices.size() << " Vertexes, " << varray.size() / 3 << " Triangles" << "\n";
  out << "g " << meshName << "\n";
  for (int i = 0; i < vertices.size(); i++)
    out << "v " << vertices.at(i)[0] << " " << vertices.at(i)[1] << " " << vertices.at(i)[2] << QString('\n');
  for (int i = 0; i < normals.size(); i++)
    out << "vn " << normals.at(i)[0] << " " << normals.at(i)[1] << " " << normals.at(i)[2] << QString('\n');
  for (int i = 0; i < varray.size(); i += 3)
  {
    out << "f " << varray.at(i) + 1 << "//" << narray.at(i) + 1 << " "
      << varray.at(i + 1) + 1 << "//" << narray.at(i + 1) + 1 << " "
      << varray.at(i + 2) + 1 << "//" << narray.at(i + 2) + 1 << " "
      << "\n";
  }
  out.flush();
  data.close();
}



/*!
\brief Create a swept-sphere along a quadric curve.
\param c %Quadric curve
\param r Radius.
\param n Number of points on circle.
\param m Number of slices on the curve (should be 2 at least).
*/
Mesh::Mesh(const QuadricCurve& c, const double& r, int n, int m)
{
  // Perimeter
  int p = n;

  // Sphere slices
  n = int(n / 4);

  // Size: two half spheres and k circles defining the cylinder

  // Number of points (and normals, overall mesh is not optimized)
  int s = p * (2 * n + m) + 2;

  vertices.resize(s);
  normals.resize(s);

  // Create set of vertices
  const double dt = Math::TwoPi / p;
  const double df = 0.5 * Math::Pi / n;

  int k = 0;

  // Phi
  double f = -Math::HalfPi;

  // Frame
  Frame frame = c.GetFrame(0.0);

  // Apex
  normals[s - 1] = -frame.R().C(0);
  vertices[s - 1] = frame.T() - r * frame.R().C(0);

  // First half sphere, does not include equator
  for (int j = 0; j < n; j++)
  {
    f += df;
    // Theta
    double t = 0.0;
    for (int i = 0; i < p; i++)
    {
      normals[k] = frame.SphereNormal(t, f);
      vertices[k] = frame.T() + r * normals[k];
      k++;
      t += dt;
    }
  }

  // Set of circles
  for (int j = 0; j < m; j++)
  {
    frame = c.GetFrame(double(j) / (m - 1));
    Vector a = frame.T();
    Vector x = frame.R().C(1);
    Vector y = frame.R().C(2);
    Vector z = frame.R().C(0);
    // Theta
    double t = 0.0;
    for (int i = 0; i < p; i++)
    {
      Vector u = cos(t) * x + sin(t) * y;
      normals[k] = u;
      vertices[k] = a + r * u;
      k++;
      t += dt;
    }
  }

  f = 0.0;

  // Frame
  frame = c.GetFrame(1.0);
  Vector a = frame.T();
  Vector x = frame.R().C(1);
  Vector y = frame.R().C(2);
  Vector z = frame.R().C(0);

  // Apex
  normals[s - 2] = z;
  vertices[s - 2] = a + r * z;

  // Second half sphere, does not include equator
  for (int j = 0; j < n; j++)
  {
    f += df;

    // Theta
    double t = 0.0;
    for (int i = 0; i < p; i++)
    {
      Vector u = cos(t) * cos(f) * x + sin(t) * cos(f) * y + sin(f) * z;
      normals[k] = u;
      vertices[k] = a + r * u;
      k++;
      t += dt;
    }
  }

  // Reserve space for the smooth triangle array
  varray.reserve(p * (n - 1 + m - 2) + 2 * p);
  narray.reserve(p * (n - 1 + m - 2) + 2 * p);

  // Apexes
  for (int i = 0; i < p; i++)
    AddSmoothTriangle(s - 1, s - 1, (i + 1) % p, (i + 1) % p, i, i);
  for (int i = 0; i < p; i++)
    AddSmoothTriangle(s - 2, s - 2, s - 2 - p + i, s - 2 - p + i, s - 2 - p + (i + 1) % p, s - 2 - p + (i + 1) % p);

  // Cylinder part
  for (int j = 1; j < (2 * n + m); j++)
  {
    for (int i = 0; i < p; i++)
      AddSmoothQuadrangle((j - 1) * p + i, (j - 1) * p + i, (j - 1) * p + (i + 1) % p, (j - 1) * p + (i + 1) % p, j * p + (i + 1) % p, j * p + (i + 1) % p, j * p + i, j * p + i);
  }
}

/*!
\brief Create a swept-sphere along a cubic curve.
\param c %Cubic curve
\param r Radius.
\param n Number of points on circle.
\param m Number of slices on the curve (should be 2 at least).
*/
Mesh::Mesh(const CubicCurve& c, const double& r, int n, int m)
{
  // Perimeter
  int p = n;

  // Sphere slices
  n = int(n / 4);

  // Size: two half spheres and k circles defining the cylinder

  // Number of points (and normals, overall mesh is not optimized)
  int s = p * (2 * n + m) + 2;

  vertices.resize(s);
  normals.resize(s);

  // Create set of vertices
  const double dt = Math::TwoPi / p;
  const double df = 0.5 * Math::Pi / n;

  int k = 0;

  // Phi
  double f = -Math::HalfPi;

  // Frame
  Frame frame = c.GetFrame(0.0);

  // Apex
  normals[s - 1] = -frame.R().C(0);
  vertices[s - 1] = frame.T() - r * frame.R().C(0);

  // First half sphere, does not include equator
  for (int j = 0; j < n; j++)
  {
    f += df;
    // Theta
    double t = 0.0;
    for (int i = 0; i < p; i++)
    {
      normals[k] = frame.SphereNormal(t, f);
      vertices[k] = frame.T() + r * normals[k];
      k++;
      t += dt;
    }
  }

  // Set of circles
  for (int j = 0; j < m; j++)
  {
    frame = c.GetFrame(double(j) / (m - 1));
    Vector a = frame.T();
    Vector x = frame.R().C(1);
    Vector y = frame.R().C(2);
    Vector z = frame.R().C(0);
    // Theta
    double t = 0.0;
    for (int i = 0; i < p; i++)
    {
      Vector u = cos(t) * x + sin(t) * y;
      normals[k] = u;
      vertices[k] = a + r * u;
      k++;
      t += dt;
    }
  }

  f = 0.0;

  // Frame
  frame = c.GetFrame(1.0);
  Vector a = frame.T();
  Vector x = frame.R().C(1);
  Vector y = frame.R().C(2);
  Vector z = frame.R().C(0);

  // Apex
  normals[s - 2] = z;
  vertices[s - 2] = a + r * z;

  // Second half sphere, does not include equator
  for (int j = 0; j < n; j++)
  {
    f += df;

    // Theta
    double t = 0.0;
    for (int i = 0; i < p; i++)
    {
      Vector u = cos(t) * cos(f) * x + sin(t) * cos(f) * y + sin(f) * z;
      normals[k] = u;
      vertices[k] = a + r * u;
      k++;
      t += dt;
    }
  }

  // Reserve space for the smooth triangle array
  varray.reserve(p * (n - 1 + m - 2) + 2 * p);
  narray.reserve(p * (n - 1 + m - 2) + 2 * p);

  // Apexes
  for (int i = 0; i < p; i++)
    AddSmoothTriangle(s - 1, s - 1, (i + 1) % p, (i + 1) % p, i, i);
  for (int i = 0; i < p; i++)
    AddSmoothTriangle(s - 2, s - 2, s - 2 - p + i, s - 2 - p + i, s - 2 - p + (i + 1) % p, s - 2 - p + (i + 1) % p);

  // Cylinder part
  for (int j = 1; j < (2 * n + m); j++)
  {
    for (int i = 0; i < p; i++)
      AddSmoothQuadrangle((j - 1) * p + i, (j - 1) * p + i, (j - 1) * p + (i + 1) % p, (j - 1) * p + (i + 1) % p, j * p + (i + 1) % p, j * p + (i + 1) % p, j * p + i, j * p + i);
  }
}


/*!
\brief Create a swept-sphere along a cubic curve.
\param c Set of cubic curves
\param r Radius.
\param n Number of points on circle.
\param m Number of slices on the curve (should be 2 at least).
*/
Mesh::Mesh(const CubicCurveSet& c, const double& r, int n, int m)
{
  // Perimeter
  int p = n;

  // Sphere slices
  n = int(n / 4);

  // Size: two half spheres and k circles defining the cylinder

  // Number of points (and normals, overall mesh is not optimized)
  int s = p * (2 * n + m) + 2;

  vertices.resize(s);
  normals.resize(s);

  // Create set of vertices
  double dt = Math::TwoPi / p;
  double df = 0.5 * Math::Pi / n;

  int k = 0;
  // Phi
  double f = -Math::HalfPi;

  Frame frame;
  Vector a;
  Vector x;
  Vector y;
  Vector z;

  // Frame
  frame = c.GetFrame(0.0);
  a = frame.T();
  x = frame.R().C(1);
  y = frame.R().C(2);
  z = frame.R().C(0);

  // Apex
  normals[s - 1] = -z;
  vertices[s - 1] = a - r * z;

  // First half sphere, does not include equator
  for (int j = 0; j < n; j++)
  {
    f += df;
    // Theta
    double t = 0.0;
    for (int i = 0; i < p; i++)
    {
      Vector u = cos(t) * cos(f) * x + sin(t) * cos(f) * y + sin(f) * z;
      normals[k] = u;
      vertices[k] = a + r * u;
      k++;
      t += dt;
    }
  }

  // Set of circles
  for (int j = 0; j < m; j++)
  {
    frame = c.GetFrame(double(j) * c.Size() / (m - 1));
    a = frame.T();
    x = frame.R().C(1);
    y = frame.R().C(2);
    z = frame.R().C(0);
    // Theta
    double t = 0.0;
    for (int i = 0; i < p; i++)
    {
      Vector u = cos(t) * x + sin(t) * y;
      normals[k] = u;
      vertices[k] = a + r * u;
      k++;
      t += dt;
    }
  }

  f = 0.0;
  // Frame
  frame = c.GetFrame(double(c.Size()));
  a = frame.T();
  x = frame.R().C(1);
  y = frame.R().C(2);
  z = frame.R().C(0);

  // Apex
  normals[s - 2] = z;
  vertices[s - 2] = a + r * z;

  // Second half shpere, does not include equator
  for (int j = 0; j < n; j++)
  {
    f += df;
    // Theta
    double t = 0.0;
    for (int i = 0; i < p; i++)
    {
      Vector u = cos(t) * cos(f) * x + sin(t) * cos(f) * y + sin(f) * z;
      normals[k] = u;
      vertices[k] = a + r * u;
      k++;
      t += dt;
    }
  }

  // Space for the smooth triangle array
  varray.reserve(p * (n - 1 + m - 2) + 2 * p);
  narray.reserve(p * (n - 1 + m - 2) + 2 * p);

  // Apexes
  for (int i = 0; i < p; i++)
  {
    AddSmoothTriangle(s - 1, s - 1, (i + 1) % p, (i + 1) % p, i, i);
  }

  for (int i = 0; i < p; i++)
  {
    AddSmoothTriangle(s - 2, s - 2, s - 2 - p + i, s - 2 - p + i, s - 2 - p + (i + 1) % p, s - 2 - p + (i + 1) % p);
  }

  // Cylinder part
  for (int j = 1; j < (2 * n + m); j++)
  {
    for (int i = 0; i < p; i++)
    {
      AddSmoothQuadrangle((j - 1) * p + i, (j - 1) * p + i, (j - 1) * p + (i + 1) % p, (j - 1) * p + (i + 1) % p, j * p + (i + 1) % p, j * p + (i + 1) % p, j * p + i, j * p + i);
    }
  }
}

/*!
\brief Create a swept-sphere along a set of quadric curves.
\param c Set of quadric curves
\param r Radius.
\param n Number of points on circle.
\param m Number of slices on the curve (should be 2 at least).
*/
Mesh::Mesh(const QuadricCurveSet& c, const double& r, int n, int m)
{
  // Perimeter
  int p = n;

  // Sphere slices
  n = int(n / 4);

  // Size: two half spheres and k circles defining the cylinder

  // Number of points (and normals, overall mesh is not optimized)
  int s = p * (2 * n + m) + 2;

  vertices.resize(s);
  normals.resize(s);

  // Create set of vertices
  const double dt = Math::TwoPi / p;
  const double df = 0.5 * Math::Pi / n;

  int k = 0;
  // Phi
  double f = -Math::HalfPi;

  Frame frame;
  Vector a;
  Vector z;

  // Frame
  frame = c.GetFrame(0.0);
  a = frame.T();
  z = frame.R().C(0);

  // Apex
  normals[s - 1] = -z;
  vertices[s - 1] = a - r * z;

  // First half sphere, does not include equator
  for (int j = 0; j < n; j++)
  {
    f += df;
    // Theta
    double t = 0.0;
    for (int i = 0; i < p; i++)
    {
      Vector u = frame.SphereNormal(t, f);
      normals[k] = u;
      vertices[k] = a + r * u;
      k++;
      t += dt;
    }
  }

  // Set of circles
  for (int j = 0; j < m; j++)
  {
    frame = c.GetFrame(double(j) * c.Size() / (m - 1));
    a = frame.T();
    // Theta
    double t = 0.0;
    for (int i = 0; i < p; i++)
    {
      Vector u = frame.CircleNormal(t);
      //Vector u = cos(t) * x + sin(t) * y;
      normals[k] = u;
      vertices[k] = a + r * u;
      k++;
      t += dt;
    }
  }

  f = 0.0;
  // Frame
  frame = c.GetFrame(double(c.Size()));
  a = frame.T();
  z = frame.R().C(0);

  // Apex
  normals[s - 2] = z;
  vertices[s - 2] = a + r * z;

  // Second half shpere, does not include equator
  for (int j = 0; j < n; j++)
  {
    f += df;
    // Theta
    double t = 0.0;
    for (int i = 0; i < p; i++)
    {
      Vector u = frame.SphereNormal(t, f);
      normals[k] = u;
      vertices[k] = a + r * u;
      k++;
      t += dt;
    }
  }

  // Space for the smooth triangle array
  varray.reserve(p * (n - 1 + m - 2) + 2 * p);
  narray.reserve(p * (n - 1 + m - 2) + 2 * p);

  // Apexes
  for (int i = 0; i < p; i++)
  {
    AddSmoothTriangle(s - 1, s - 1, (i + 1) % p, (i + 1) % p, i, i);
  }

  for (int i = 0; i < p; i++)
  {
    AddSmoothTriangle(s - 2, s - 2, s - 2 - p + i, s - 2 - p + i, s - 2 - p + (i + 1) % p, s - 2 - p + (i + 1) % p);
  }

  // Cylinder part
  for (int j = 1; j < (2 * n + m); j++)
  {
    for (int i = 0; i < p; i++)
    {
      AddQuadrangle((j - 1) * p + i, (j - 1) * p + (i + 1) % p, j * p + (i + 1) % p, j * p + i);
    }
  }
}

/*!
\brief Create a generalized cylinder with star-shaped contours.
\param ca, cb Star-shaped contours, they should have the same size.
\param a, b Centers of the two contours.
*/
Mesh::Mesh(const QVector<Vector>& ca, const QVector<Vector>& cb, const Vector& a, const Vector& b)
{
  int n = ca.size();

  // Vertices
  vertices += ca; // [0,n-1]
  vertices += cb; // [n,2n-1]
  vertices.append(a);  // 2n
  vertices.append(b);  // 2n+1

  // Normals
  for (int i = 0; i < n; i++)
  {
    normals.append(Triangle(ca[i], ca[(i + 1) % n], cb[i]).Normal());
  }
  normals.append(Triangle(a, ca[1], ca[0]).Normal()); // n
  normals.append(Triangle(b, cb[0], cb[1]).Normal()); // n+1

  // Indexes
  // Caps
  for (int i = 0; i < n; i++)
  {
    AddTriangle(2 * n, (i + 1) % n, i, n);
  }

  for (int i = 0; i < n; i++)
  {
    AddTriangle(2 * n + 1, n + i, n + (i + 1) % n, n + 1);
  }
  // Side
  for (int i = 0; i < n; i++)
  {
    AddTriangle(i, (i + 1) % n, n + i, i);
    AddTriangle((i + 1) % n, n + (i + 1) % n, n + i, i);
  }
}

/*!
\brief Create the mesh of an torus.
\param t Torus.
\param n,p Discretization of the generating circles.
\param ia, ib Starting and ending indexes of the sections of the torus
*/
Mesh::Mesh(const Torus& t, int n, int p, int ia, int ib)
{
  if (ia == -1 || ib == -1)
  {
    ia = 0;
    ib = n;
  }

  const Vector c = t.Center();
  const Vector a = t.Axis();
  const double r = t.Radius();
  const double s = t.Small();

  // Number of cylinder segments, used for index j
  int js = (ib - ia) > n ? n : ib - ia - 1;

  // Size 
  int size = (ib - ia) * p;

  vertices.resize(size);
  normals.resize(size);

  int k = 0;

  // Frame
  Vector z = Normalized(a);
  Vector x, y;
  z.Orthonormal(x, y);

  for (int j = ia; j < ib; j++)
  {
    // Theta
    double theta = Math::Angle(j, n);

    for (int i = 0; i < p; i++)
    {
      double phi = Math::Angle(i, p);

      Vector u = cos(theta) * x + sin(theta) * y;
      Vector v = cos(phi) * u + sin(phi) * z;
      normals[k] = v;
      vertices[k] = c + r * u + s * v;
      k++;
    }
  }

  // Reserve space for the smooth triangle array
  varray.reserve(3 * 2 * size);
  narray.reserve(3 * 2 * size);

  // Sections of the torus
  for (int j = 0; j < js; j++)
  {
    for (int i = 0; i < p; i++)
    {
      AddSmoothQuadrangle(j * p + i, j * p + i, ((j + 1) % n) * p + i, ((j + 1) % n) * p + i, ((j + 1) % n) * p + (i + 1) % p, ((j + 1) % n) * p + (i + 1) % p, j * p + (i + 1) % p, j * p + (i + 1) % p);
    }
  }
}


/*!
\brief Create the mesh of a cuboctahedron.
\param c %Cuboctahedron.
*/
Mesh::Mesh(const Cuboctahedron& c)
{
  vertices.resize(12);
  normals.resize(14);

  for (int i = 0; i < 12; i++)
  {
    vertices[i] = c.Vertex(i);
  }

  for (int i = 0; i < 14; i++)
  {
    normals[i] = c.Normal(i);
  }

  // Triangles : 8
  for (int i = 0; i < 8; i++)
  {
    AddTriangle(Cuboctahedron::VertexIndex(i, 0), Cuboctahedron::VertexIndex(i, 1), Cuboctahedron::VertexIndex(i, 2), i);
  }
  // Squares : 6
  for (int i = 8; i < 20; i++)
  {
    AddTriangle(Cuboctahedron::VertexIndex(i, 0), Cuboctahedron::VertexIndex(i, 1), Cuboctahedron::VertexIndex(i, 2), 8 + (i - 8) / 2);
  }
}


/*!
\brief Create the mesh of an icosidodecahedron.
\param icosi Icosidodecahedron.
*/
Mesh::Mesh(const Icosidodecahedron& icosi)
{
  vertices.resize(30);
  normals.resize(32);

  for (int i = 0; i < 30; i++)
  {
    vertices[i] = icosi.Vertex(i);
  }

  for (int i = 0; i < 32; i++)
  {
    normals[i] = icosi.Normal(i);
  }

  // Pentagons : 12
  for (int i = 0; i < 12; i++)
  {
    AddPentagon(Icosidodecahedron::VertexPentagonIndex(i, 0), Icosidodecahedron::VertexPentagonIndex(i, 1), Icosidodecahedron::VertexPentagonIndex(i, 2), Icosidodecahedron::VertexPentagonIndex(i, 3), Icosidodecahedron::VertexPentagonIndex(i, 4), i);
  }

  // Triangles : 20
  for (int i = 0; i < 20; i++)
  {
    AddTriangle(Icosidodecahedron::VertexTriangleIndex(i, 0), Icosidodecahedron::VertexTriangleIndex(i, 1), Icosidodecahedron::VertexTriangleIndex(i, 2), 12 + i);
  }
}

/*!
\brief Scales (shrinks) all the triangles of the model.

It should be possible to look into the geometry as holes
will appear at the edges and vertices of the model.

\param e Erosion.
*/
Mesh Mesh::ShrinkedTriangles(const double& e) const
{
  QVector<Vector> v;
  QVector<Vector> n;
  QVector<int> vi;
  QVector<int> ni;

  v.reserve(varray.size());
  n.reserve(varray.size() / 3);
  vi.reserve(varray.size());
  ni.reserve(varray.size());

  for (int i = 0; i < varray.size(); i += 3)
  {
    Vector a = vertices.at(varray.at(i));
    Vector b = vertices.at(varray.at(i + 1));
    Vector c = vertices.at(varray.at(i + 2));
    Triangle t(a, b, c);
    t.Shrink(e);

    v.append(t[0]);
    v.append(t[1]);
    v.append(t[2]);
    n.append(t.Normal());
    vi.append(i);
    vi.append(i + 1);
    vi.append(i + 2);
    ni.append(i / 3);
    ni.append(i / 3);
    ni.append(i / 3);
  }
  return Mesh(v, n, vi, ni);
}


/*!
\brief Create the mesh of a hexagonal prism.
\param hexa Hexagonal prism.
*/
Mesh::Mesh(const Hexagonal& hexa)
{
  vertices.resize(14);
  normals.resize(8);

  for (int i = 0; i < 12; i++)
  {
    vertices[i] = hexa.Vertex(i);
  }
  vertices[12] = hexa.Center(0);
  vertices[13] = hexa.Center(1);

  for (int i = 0; i < 6; i++)
  {
    normals[i] = hexa.Normal(i);
  }
  normals[6] = -Vector::Z;
  normals[7] = Vector::Z;

  // Quadrangles : 6
  for (int i = 0; i < 6; i++)
  {
    AddQuadrangle(i, (i + 1) % 6, 6 + (i + 1) % 6, 6 + i, i);
  }

  // Triangles : 6+6
  for (int i = 0; i < 6; i++)
  {
    AddTriangle(12, (i + 1) % 6, i, 6);
    AddTriangle(13, 6 + i, 6 + (i + 1) % 6, 7);
  }
}

/*!
\brief Returns the geometry as a set of triangles.
*/
QVector<Triangle> Mesh::GetTriangles() const
{
  QVector<Triangle> triangles;

  for (int i = 0; i < varray.size(); i += 3)
  {
    triangles.append(Triangle(vertices.at(varray.at(i + 0)), vertices.at(varray.at(i + 1)), vertices.at(varray.at(i + 2))));
  }

  return triangles;
}

/*!
\brief Overloaded.

\param s Stream.
\param mesh The mesh.
*/
std::ostream& operator<<(std::ostream& s, const Mesh& mesh)
{
  s << "Mesh(" << mesh.vertices.size() << ',' << mesh.normals.size() << ',' << mesh.varray.size() << ',' << mesh.narray.size() << ')' << std::endl;
  for (int i = 0; i < mesh.vertices.size(); i++)
  {
    s << mesh.vertices.at(i) << std::endl;
  }
  for (int i = 0; i < mesh.varray.size(); i += 3)
  {
    s << mesh.varray.at(i) << ' ' << mesh.varray.at(i + 1) << ' ' << mesh.varray.at(i + 2) << std::endl;
  }
  for (int i = 0; i < mesh.narray.size(); i += 3)
  {
    s << mesh.narray.at(i) << ' ' << mesh.narray.at(i + 1) << ' ' << mesh.narray.at(i + 2) << std::endl;
  }
  return s;
}
/*!
\brief Compute the squared distance to the mesh.

This function computes the squared distance to the surface, not the squared distance to the solid object.

\param p Point.
*/
double Mesh::R(const Vector& p) const
{
  const int n = varray.size() / 3;

  double d = Math::Infinity;
  for (int i = 0; i < n; i++)
  {
    d = Math::Min(d, GetTriangle(i).R(p));
  }

  return d;
}

/*!
\brief Compute the signed distance to the mesh.

The mesh should be closed.

\param p Point.
*/
double Mesh::Signed(const Vector& p) const
{
  const int n = varray.size() / 3;

  // Squared distance to the surface
  double d = Math::Infinity;
  for (int i = 0; i < n; i++)
  {
    d = Math::Min(d, GetTriangle(i).R(p));
  }

  // Distance
  d = sqrt(d);

  // Exterior check
  int c = Intersections(Ray(p, Vector::X));

  if (c % 2 == 1)
  {
    d = -d;
  }
  return d;

}

/*!
\brief Merge two meshes into a single one.
\param mesh Argument mesh.
*/
void Mesh::Append(const Mesh& mesh)
{
  // Shifts
  const int v = vertices.size();
  const int n = normals.size();

  // Vertexes
  vertices.append(mesh.vertices);

  // Normals
  normals.append(mesh.normals);

  // Include shifting for vertex and normal indexes
  for (int i = 0; i < mesh.varray.size(); i++)
  {
    varray.append(mesh.varray.at(i) + v);
  }
  for (int i = 0; i < mesh.narray.size(); i++)
  {
    narray.append(mesh.narray.at(i) + n);
  }
}

/*!
\brief Creates an arrow.
The object has 3n+2 vertices, 2n+2 normals and 6n triangles.

Note that this shape is more efficient in terms of storage than a cylinder and a cone.
\param a, b End vertices.
\param r, R %Cylinder and cone radii.
\param l Length of the cylinder, the size of the cone will be derived from the distance <b>b</b>-<b>a</b> minus this length.
\param n Discretization.
*/
Mesh Mesh::Arrow(const Vector& a, const Vector& b, const double& r, const double& R, const double& l, int n)
{
  Vector z = Normalized(b - a);

  // Frame
  Vector x, y;
  z.Orthonormal(x, y);

  Mesh mesh;

  // Set dimension of array of vertices and normals
  mesh.vertices.resize(3 * n + 2);
  mesh.normals.resize(2 * n + 2);

  // Vertices and normals
  for (int i = 0; i < n; i++)
  {
    double t = Math::Angle(i, n);
    Vector u = x * cos(t) + y * sin(t);

    mesh.vertices[i] = a + r * u;
    mesh.vertices[n + i] = a + l * z + r * u;
    mesh.vertices[2 * n + i] = a + l * z + R * u;
    mesh.normals[i] = u;
    mesh.normals[n + i] = Normalized(R * z + (Norm(b - a) - l) * u);
  }

  // End vertices
  mesh.vertices[3 * n] = a;
  mesh.vertices[3 * n + 1] = b;

  // End normals
  mesh.normals[2 * n] = -z;
  mesh.normals[2 * n + 1] = z;

  // Reserve space for the triangle and smooth triangle arrays
  mesh.varray.reserve(6 * n * 3);
  mesh.narray.reserve(6 * n * 3);

  // Lower cap of cylinder
  for (int i = 0; i < n; i++)
  {
    mesh.AddTriangle(3 * n, (i + 1) % n, i, 2 * n);
  }

  // Cone side
  for (int i = 0; i < n; i++)
  {
    mesh.AddSmoothTriangle(2 * n + i, n + i, 2 * n + (i + 1) % n, n + (i + 1) % n, 3 * n + 1, 2 * n + 1);
  }

  // Cylinder side
  for (int i = 0; i < n; i++)
  {
    mesh.AddSmoothTriangle(i, i, (i + 1) % n, (i + 1) % n, n + i, i);
    mesh.AddSmoothTriangle((i + 1) % n, (i + 1) % n, n + (i + 1) % n, (i + 1) % n, n + i, i);
  }

  // Mid cylider closing cap
  for (int i = 0; i < n; i++)
  {
    mesh.AddTriangle(n + i, n + (i + 1) % n, 2 * n + i, 2 * n);
    mesh.AddTriangle(n + (i + 1) % n, 2 * n + (i + 1) % n, 2 * n + i, 2 * n);
  }
  return mesh;
}

/*!
\brief Creates the mesh of a cuboid.
\param cuboid %The cuboid.
*/
Mesh::Mesh(const Cuboid& cuboid) :vertices(8), normals(12)
{
  // Vertices
  for (int i = 0; i < 8; i++)
  {
    vertices[i] = cuboid.Vertex(i);
  }

  // Vertices
  for (int i = 0; i < 6; i++)
  {
    Quadrangle q = cuboid.GetQuadrangle(i);
    normals[2 * i] = q.Normal();
    normals[2 * i + 1] = q.Normal();
  }

  // Reserve space for the triangle array
  varray.reserve(12 * 3);
  narray.reserve(12 * 3);

  AddTriangle(0, 2, 1, 8);
  AddTriangle(1, 2, 3, 9);

  AddTriangle(4, 5, 6, 10);
  AddTriangle(5, 7, 6, 11);

  AddTriangle(0, 4, 2, 0);
  AddTriangle(4, 6, 2, 1);

  AddTriangle(1, 3, 5, 2);
  AddTriangle(3, 7, 5, 3);

  AddTriangle(0, 1, 5, 4);
  AddTriangle(0, 5, 4, 5);

  AddTriangle(3, 2, 7, 6);
  AddTriangle(6, 7, 2, 7);
}


/*!
\brief Generate a grid geometry and polygonize the quadric surface.
\param s %The quadric surface.
\param box %The box.
\param nx, ny Number of points in the subdivision.
*/
Mesh::Mesh(const QuadricSurface& s, const Box2& box, int nx, int ny)
{
  Array2 a(box, nx, ny);
  vertices.resize(nx * ny);
  normals.resize(nx * ny);

  varray.resize(3 * (nx - 1) * (ny - 1) * 2);
  narray.resize(3 * (nx - 1) * (ny - 1) * 2);

  // Vertices and normals
  for (int j = 0; j < nx; j++)
  {
    for (int i = 0; i < ny; i++)
    {
      int ij = a.VertexIndex(i, j);

      Vector2 p = a.ArrayVertex(i, j);

      vertices[ij] = p.ToVector(s.Value(p));
      Vector2 g = s.Gradient(Vector2(p));
      normals[ij] = Normalized(Vector(-g[0], -g[1], 1.0));
    }
  }

  // Triangle indexes
  int n = 0;
  for (int j = 0; j < nx - 1; j++)
  {
    for (int i = 0; i < ny - 1; i++)
    {
      narray[n + 0] = varray[n + 0] = a.VertexIndex(i, j);
      narray[n + 1] = varray[n + 1] = a.VertexIndex(i + 1, j);
      narray[n + 2] = varray[n + 2] = a.VertexIndex(i, j + 1);

      n += 3;

      narray[n + 0] = varray[n + 0] = a.VertexIndex(i + 1, j);
      narray[n + 1] = varray[n + 1] = a.VertexIndex(i + 1, j + 1);
      narray[n + 2] = varray[n + 2] = a.VertexIndex(i, j + 1);

      n += 3;
    }
  }

}


/*!
\brief Compute the surface of the voxel.
\param voxel The voxel.
*/
Mesh::Mesh(const Voxel& voxel)
{
  const int nvx = voxel.GetSizeX();
  const int nvy = voxel.GetSizeY();
  const int nvz = voxel.GetSizeZ();

  int nv = 0;

  // Define normals
  normals = { Vector::X, Vector::Y, Vector::Z, -Vector::X, -Vector::Y, -Vector::Z };

  // Sweeping horizontal grid of vertex indexes
  Array2 sweep(Box2(voxel.GetBox()), nvx, nvy);

  Array2I* pa = new Array2I(sweep);
  Array2I* pb = new Array2I(sweep);

  // Initialize plane a
  for (int i = 0; i < nvx; i++)
  {
    for (int j = 0; j < nvy; j++)
    {
      if (voxel.IsVertex(i, j, 0))
      {
        vertices.append(voxel.Vertex(i, j, 0));
        (*pa)(i, j) = nv;
        nv++;
      }
      else
      {
        (*pa)(i, j) = -1;
      }
    }
  }
  // Compute xy quadrangles in plane a
  for (int i = 0; i < nvx - 1; i++)
  {
    for (int j = 0; j < nvy - 1; j++)
    {
      if (((*pa)(i, j) != -1) && ((*pa)(i + 1, j) != -1) && ((*pa)(i, +1) != -1) && ((*pa)(i + 1, j + 1) != -1))
      {
        // Add face -Z]
        AddQuadrangle((*pa)(i, j), (*pa)(i + 1, j), (*pa)(i + 1, +1), (*pa)(i, j + 1), 5);
      }
    }
  }

  for (int k = 1; k < nvz; k++)
  {
    // Compute vertex indexes for plane b
    for (int i = 0; i < nvx; i++)
    {
      for (int j = 0; j < nvy; j++)
      {
        if (voxel.IsVertex(i, j, k))
        {
          vertices.append(voxel.Vertex(i, j, k));
          (*pb)(i, j) = nv;
          nv++;
        }
        else
        {
          (*pb)(i, j) = -1;
        }
      }
    }

    // Compute xy quadrangles in plane b
    for (int i = 0; i < nvx - 1; i++)
    {
      for (int j = 0; j < nvy - 1; j++)
      {
        if (((*pb)(i, j) != -1) && ((*pb)(i + 1, j) != -1) && ((*pb)(i, j + 1) != -1) && ((*pb)(i + 1, +1) != -1))
        {
          // Add face +/-Z
          if (voxel.At(i, j, k - 1) > 0)
          {
            AddQuadrangle((*pb)(i, j), (*pb)(i + 1, j), (*pb)(i + 1, +1), (*pb)(i, j + 1), 5);
          }
          else
          {
            AddQuadrangle((*pb)(i, j), (*pb)(i + 1, j), (*pb)(i + 1, +1), (*pb)(i, j + 1), 2);
          }
        }
      }
    }

    // Compute xz quadrangles 
    for (int i = 0; i < nvx - 1; i++)
    {
      for (int j = 0; j < nvy; j++)
      {
        if (((*pa)(i, j) != -1) && ((*pa)(i, +1) != -1) && ((*pb)(i, j) != -1) && ((*pb)(i, j + 1) != -1))
        {
          // Add face +/-X
          if (1)
          {
            AddQuadrangle((*pa)(i, j), (*pa)(i, +1), (*pb)(i, j + 1), (*pb)(i, j), 3);
          }
          else
          {
            AddQuadrangle((*pa)(i, j), (*pa)(i, +1), (*pb)(i, j + 1), (*pb)(i, j), 0);
          }
        }
      }
    }
    // Compute yz quadrangles 
    for (int j = 0; j < nvy - 1; j++)
    {
      for (int i = 0; i < nvx; i++)
      {
        if (((*pa)(i, j) != -1) && ((*pb)(i + 1, j) != -1) && ((*pb)(i, j) != -1) && ((*pb)(i + 1, j) != -1))
        {
          // Add face +/-Y
          if (1)
          {
            AddQuadrangle((*pa)(i, j), (*pb)(i + 1, j), (*pb)(i + 1, j), (*pb)(i, j), 4);
          }
          else
          {
            AddQuadrangle((*pa)(i, j), (*pb)(i + 1, j), (*pb)(i + 1, j), (*pb)(i, j), 1);
          }
        }
      }
    }

    std::swap(pa, pb);
  }

  delete pa;
  delete pb;
}
