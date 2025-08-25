#pragma once


#include <QtCore/QVector>

#include "libs/triangle.h"
#include "libs/quadrangle.h"
#include "libs/capsule.h"
#include "libs/circle.h"
#include "libs/polygon.h"
#include "libs/plane.h"

class PointCurve;
class PointCurve2;
class QuadricCurve;
class QuadricCurveSet;
class CubicCurve;
class CubicCurveSet;
class CubicCurve2Set;
class Circle2;
class Ellipse2;

class Torus;
class Icosidodecahedron;
class Icosahedron;
class Dodecahedron;
class Cuboctahedron;
class Octahedron;
class Pyramid;
class Sphere;
class Hexagonal;
class Cone;
class Rectangles;
class SphereSet;
class Cuboid;
class Voxel;

class QuadricSurface;

class Mesh2;

class Mesh
{
protected:
  QVector<Vector> vertices; //!< Vertices.
  QVector<Vector> normals;  //!< Normals.
  QVector<int> varray;      //!< Vertex indexes.
  QVector<int> narray;      //!< Normal indexes.
public:
  explicit Mesh();
  explicit Mesh(const Mesh2&);
  explicit Mesh(const QVector<Vector>&);
  explicit Mesh(const QVector<Vector>&, const QVector<int>&, bool = true);
  explicit Mesh(const QVector<Vector>&, const QVector<Vector>&, const QVector<int>&, const QVector<int>&);
  ~Mesh();

  // Convex hull
  static Mesh Hull(const QVector<Vector>&);

  void Reserve(int, int, int, int);

  Triangle GetTriangle(int) const;
  const Vector Vertex(int) const;
  Vector& Vertex(int);
  Vector Vertex(int, int) const;

  Vector Normal(int) const;

  // Number of elements
  int Triangles() const;
  int Vertexes() const;
  int Normals() const;

  bool Manifold() const;

  void Append(const Mesh&);

  QVector<Vector> GetVertices() const;
  QVector<Vector> GetNormals() const;
  QVector<int> VertexIndexes() const;
  QVector<int> NormalIndexes() const;

  QVector<Triangle> GetTriangles() const;

  int VertexIndex(int, int) const;
  int VertexIndex(int) const;
  int NormalIndex(int, int) const;
  int NormalIndex(int) const;

  Vector operator[](int) const;

  double R(const Vector&) const;
  double Signed(const Vector&) const;

  bool Intersect(const Ray&, double&, int&) const;
  int Intersections(const Ray&) const;

  bool Inside(const Vector&) const;

  Box GetBox() const;
  Sphere GetSphere() const;

  double AverageEdgeLength() const;
  void EdgeLengthRange(double&, double&) const;
  QVector<int> AspectRatios(int = 4) const;

  void EdgeCollapse(const double&);
  void EdgeCollapse(int);

  void Transform(const Frame&);
  void Transform(const FrameScaled&);
  Mesh Transformed(const FrameScaled&) const;
  void InverseTransform(const FrameScaled&);

  void Symmetry(const Plane&);
  void Rotate(const Matrix&);
  void Translate(const Vector&);
  void Scale(const Vector&);
  void Scale(const double&);

  void SmoothNormals();
  int IsSmooth(int) const;

  bool FindEdge(int, int, int&, int = 0) const;

  // Constructors from core classes
  explicit Mesh(const Capsule&, int, int = 2);
  explicit Mesh(const Quadrangle&, int, int);
  explicit Mesh(const Quadrangle&);
  explicit Mesh(const Octahedron&);
  explicit Mesh(const Pyramid&);
  explicit Mesh(const Box&);
  explicit Mesh(const Triangle&);
  explicit Mesh(const Sphere&, int);
  explicit Mesh(const SphereSet&, int);
  explicit Mesh(const Cylinder&, int, bool = true, bool = true);
  explicit Mesh(const Disc&, int);
  explicit Mesh(const Dodecahedron&);
  explicit Mesh(const Icosahedron&);
  explicit Mesh(const QVector<Vector>&, const QVector<Vector>&, const Vector&, const Vector&);
  explicit Mesh(const Torus&, int, int, int = -1, int = -1);
  explicit Mesh(const Cuboctahedron&);
  explicit Mesh(const Icosidodecahedron&);
  explicit Mesh(const Hexagonal&);
  explicit Mesh(const Cone&, int);
  explicit Mesh(const Rectangles&, int, int);
  explicit Mesh(const Cuboid&);
  explicit Mesh(const QuadricSurface&, const Box2&, int, int);
  explicit Mesh(const Voxel&);

  // Swept spheres
  explicit Mesh(const QuadricCurve&, const double&, int, int);
  explicit Mesh(const CubicCurve&, const double&, int, int);
  explicit Mesh(const CubicCurveSet&, const double&, int, int);
  explicit Mesh(const QuadricCurveSet&, const double&, int, int = 2);

  // Extrusion
  static Mesh Extrusion(const PointCurve&, const PointCurve2&);
  static Mesh Extrusion(const PointCurve&, const QVector<PointCurve2>&);
  static Mesh ExtrusionRotation(const PointCurve2&, const double&, const double&, int, const Vector & = Vector::Null);
  static Mesh ExtrusionRotation(const PointCurve2&, int, const Vector & = Vector::Null);

  static Mesh Extrusion4new(Vector, Vector, const QVector<PointCurve>&);
  static Mesh Extrusion3new(const PointCurve&, const QVector<PointCurve>&, QVector<int>, QVector<double>);
  static Mesh Extrusion10new(const PointCurve&, const QVector<PointCurve>&);

  static Mesh extrusion5new(Vector, Vector, const QVector<PointCurve>&, const QVector<PointCurve>&);
  static Mesh extrusion8new(QVector<Vector>, double, double, double);

  static Mesh Arrow(const Vector&, const Vector&, const double&, const double&, const double&, int);

  Mesh ShrinkedTriangles(const double&) const;

  void Load(const QString&);
  void SaveObj(const QString&, const QString & = QString("mesh")) const;

  friend std::ostream& operator<<(std::ostream&, const Mesh&);
protected:
  void AddTriangle(int, int, int, int);
  void AddSmoothTriangle(int, int, int, int, int, int);
  void AddSmoothQuadrangle(int, int, int, int, int, int, int, int);
  void AddQuadrangle(int, int, int, int);
  void AddQuadrangle(int, int, int, int, int);
  void AddPentagon(int, int, int, int, int, int);

  void AddArray(int, int);
protected:
  int NextIndex(int) const;
  int PrevIndex(int) const;
  int BaseIndex(int) const;
  friend class MeshStack;
};

/*!
\brief Return the set of vertex indexes.
*/
inline QVector<int> Mesh::VertexIndexes() const
{
  return varray;
}

/*!
\brief Return the set of normal indexes.
*/
inline QVector<int> Mesh::NormalIndexes() const
{
  return narray;
}

/*!
\brief Get the vertex index of a given triangle.
\param t Triangle index.
\param i Vertex index.
*/
inline int Mesh::VertexIndex(int t, int i) const
{
  return varray.at(t * 3 + i);
}

/*!
\brief Get the normal index of a given triangle.
\param t Triangle index.
\param i Normal index.
*/
inline int Mesh::NormalIndex(int t, int i) const
{
  return narray.at(t * 3 + i);
}

/*!
\brief Get the vertex index.
\param i Index.
*/
inline int Mesh::VertexIndex(int i) const
{
  return varray.at(i);
}

/*!
\brief Get the normal index.
\param i Index.
*/
inline int Mesh::NormalIndex(int i) const
{
  return narray.at(i);
}

/*!
\brief Get a triangle.
\param i Index.
\return The triangle.
*/
inline Triangle Mesh::GetTriangle(int i) const
{
  return Triangle(vertices.at(varray.at(i * 3 + 0)), vertices.at(varray.at(i * 3 + 1)), vertices.at(varray.at(i * 3 + 2)));
}

/*!
\brief Get a vertex.
\param i The index of the vertex.
\return The vertex position.
*/
inline const Vector Mesh::Vertex(int i) const
{
  return vertices[i];
}

/*!
\brief Acces to a vertex.
\param i The index of the vertex.
*/
inline Vector& Mesh::Vertex(int i)
{
  return vertices[i];
}

/*!
\brief Get a vertex from a specific triangle.
\param t The number of the triangle wich contain the wanted vertex.
\param v The triangle vertex: 0, 1, or 2.
\return The wanted vertex (as a 3D Vector).
*/
inline Vector Mesh::Vertex(int t, int v) const
{
  return vertices[varray[t * 3 + v]];
}

/*!
\brief Get the number of vertices in the geometry.
\return The number of vertices in the geometry, in other words the size of vertices.
*/
inline int Mesh::Vertexes() const
{
  return varray.size();
}

/*!
\brief Get the number of normals in the geometry.
\return The number of normals in the geometry, in other words the size of vertices.
*/
inline int Mesh::Normals() const
{
  return narray.size();
}

/*!
\brief Get a normal.
\param i Index of the wanted normal.
\return The normal.
*/
inline Vector Mesh::Normal(int i) const
{
  return normals[i];
}

/*!
\brief Get the number of triangles.
*/
inline int Mesh::Triangles() const
{
  return varray.size() / 3;
}

/*!
\brief Get a vertex.
\param i The index of the wanted vertex.
\return The wanted vertex (as a 3D Vector).
\see vertex(int i) const
*/
inline Vector Mesh::operator[](int i) const
{
  return vertices[i];
}

/*!
\brief Get the array of vertices.
*/
inline QVector<Vector> Mesh::GetVertices() const
{
  return vertices;
}

/*!
\brief Get the array of normals.
*/
inline QVector<Vector> Mesh::GetNormals() const
{
  return normals;
}

class MeshStack :public Mesh
{
protected:
  FrameScaled frame; //!< Current transformation.
  QVector<FrameScaled> stackframe; //!< Stack of frames.
public:
  MeshStack();
  MeshStack(const Mesh&);
  //! Empty
  ~MeshStack() {}

  // Stack management
  void Push();
  void Pop();
  void PopAll();

  void Transform(const FrameScaled&);
  void Translate(const Vector&);
  void Rotate(const Vector&);
  void Scale(const Vector&);

  // Add simple primitives
  void AddTriangle(const Vector&, const Vector&, const Vector&);
  void AddQuadrangle(const Vector&, const Vector&, const Vector&, const Vector&);

  // Add a complete geometry
  void AddGeometry(const Mesh&);
  friend MeshStack& operator<<(MeshStack&, const Mesh&);

  FrameScaled GetFrame() const;
protected:
  void Reserve(int, int, int);
};

class Mesh2
{
protected:
  QVector<Vector2> vertices; //!< Vertices.
  QVector<int> indices; //!< %Triangle vertex indices.
public:
  explicit Mesh2();
  explicit Mesh2(const QVector<Vector2>&, const QVector<int>&);

  // Constructors from core classes
  explicit Mesh2(const Quadrangle2&, int, int);
  explicit Mesh2(const Box2&, int, int);
  explicit Mesh2(const Circle2&, int);
  explicit Mesh2(const Ellipse2&, int);
  explicit Mesh2(const Polygon2&);

  ~Mesh2();

  Triangle2 GetTriangle(int) const;
  Vector2 Vertex(int) const;
  Vector2 Vertex(int, int) const;

  Box2 GetBox() const;

  int VertexSize() const;
  int TriangleSize() const;
  int IndexSize() const;

  int index(int) const;
  int index(int, int) const;

  QVector<int> Indexes() const;
  QVector<int>& Indexes();
  QVector<Vector2> Vertices() const;

  QVector<int> Valences() const;

  void Transform(const Frame2&);
  void Translate(const Vector2&);
  void Scale(const Vector2&);
  void Rotate(const Matrix2&);

  void Draw(QGraphicsScene&, const QPen & = QPen(), bool = true) const;

  double Aspect() const;

  Mesh2 Intersect(const Box2&) const;

  static Mesh2 Delaunay(const QVector<Vector2>&);
  static Mesh2 Alpha(const QVector<Vector2>&,const double&);
  static Mesh2 DelaunayN4(const QVector<Vector2>&);
public:
  static int NextIndex(int);
};

/*!
\brief Return the set of indexes.
*/
inline QVector<int> Mesh2::Indexes() const
{
  return indices;
}

/*!
\brief Return the set of indexes.
*/
inline QVector<int>& Mesh2::Indexes()
{
  return indices;
}

/*!
\brief Get the array of vertices.
*/
inline QVector<Vector2> Mesh2::Vertices() const
{
  return vertices;
}

/*!
\brief Get a triangle.
\param i Index.
\return The triangle.
*/
inline Triangle2 Mesh2::GetTriangle(int i) const
{
  return Triangle2(vertices.at(indices.at(i * 3 + 0)), vertices.at(indices.at(i * 3 + 1)), vertices.at(indices.at(i * 3 + 2)));
}

/*!
\brief Get a vertex.
\param i The index of the wanted vertex.
\return The wanted vertex (as a 3D Vector).
*/
inline Vector2 Mesh2::Vertex(int i) const
{
  return vertices.at(i);
}

/*!
\brief Get a vertex from a specific triangle.
\param t The number of the triangle that contains the wanted vertex.
\param v The triangle vertex: 0, 1, or 2.
\return The wanted vertex.
*/
inline Vector2 Mesh2::Vertex(int t, int v) const
{
  return vertices.at(indices.at(t * 3 + v));
}

/*!
\brief Get the vertex index.
\param i The number in indices of the wanted vertex.
\return The index of the wanted vertex.
*/
inline int Mesh2::index(int i) const
{
  return indices.at(i);
}

/*!
\brief Get a vertex/normal index according to triangle information.
\param t The triangle number.
\param v The triangle vertex index.
\return The index of the wanted vertex.
*/
inline int Mesh2::index(int t, int v) const
{
  return indices.at(t * 3 + v);
}

/*!
\brief Get the number of triangles.

This is the size of the array of indexes divided by three.

\return The effective number of triangle in the geometry.
*/
inline int Mesh2::TriangleSize() const
{
  return indices.size() / 3;
}

/*!
\brief Get the size of the index array.

\sa Mesh2::TriangleSize()
*/
inline int Mesh2::IndexSize() const
{
  return indices.size();
}

/*!
\brief Get the number of vertices in the geometry.
\return The number of vertices in the geometry, in other words the size of vertices.
*/
inline int Mesh2::VertexSize() const
{
  return vertices.size();
}

class MeshTopo : protected Mesh
{
protected:
  QVector<int> tarray; //!< Neighbouring half edge, from which it is possible to derive the neighboring triangle.
public:
  explicit MeshTopo();
  explicit MeshTopo(const Mesh&);
  ~MeshTopo();

  void EdgeCollapse(const double&);
protected:

  void EdgeCollapse(int);
  std::vector<int> SharedVertex(int) const;
};

/*!
\brief Compute the next index of an edge in a triangle.
\param i Index of the starting vertex of an edge.
*/
inline int Mesh2::NextIndex(int i)
{
  return ((i % 3) == 2) ? i - 2 : i + 1;
}
