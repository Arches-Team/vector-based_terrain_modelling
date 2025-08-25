#include <QtCore/QFile>

#include "libs/mayageometry.h"
#include "libs/mesh.h"
#include "libs/cylinder.h"

/*!
\mainpage MayaLib
<center>
\image html maya-small.png
<P>Latest stable version : 2024.01.12
</center>
<h3>License</h3>
Copyright &copy; by
<a href="mailto:eric.galin@liris.cnrs.fr">Eric Galin</a>, <a href="mailto:adrien.peytavie@liris.cnrs.fr">Adrien Peytavie</a>, <a href="mailto:eric.guerin@liris.cnrs.fr">Eric Guérin</a>.
No representations are made about the suitability of this software
for any purpose. It is provided as is without express or implied
warranty.
<h3>Features overview</h3>
This library implements core classes for handling meshes.
The overall architecture of the system is as follows:
\image html architecture.png

<h2>Fundamentals</h2>

The Maya library is built upon a set of core classes, provided by the Core library. These include the following:
-# The MayaGeometry classes, a set of classes for storing the geometry of objects.
-# The MayaInstance classes for handling objects stored on the disc and identified by their name.
-# The MayaGpu classes for handling objects stored on the graphics card.

<h2>Graphical User Interface</h2>

The Maya library includes support for Qt development with QtGLWidget.
This class provides a simple user interface for visualizing geometric objects defined using MayaGeometry or MayaInstance elements.

*/

/*!
\defgroup MayaCore Maya geometry core classes

\brief Maya core geometric classes include several core classes such as MayaGeometry, MayaInstance,
MayaInstanceSet and many others that are usefull for handling meshes.
*/

/*!
\defgroup MayaGpuGroup Maya GPU classes

\brief Maya GPU classes implement Vertex Buffer Objects, shading and intancing on the GPU.
*/

/*!
\defgroup MayaInstanceGroup Maya instance classes

\brief Maya instance classes implement a framework for working with named instances.
MayaInstance, MayaInstanceSet and MayaInstanceAll classes are used in the MayaResources
which manages atlases of assets referenced by their names.
*/

/*!
\class MayaGeometry maya.h
\brief A simple triangle mesh representation.

This class has been designed for the manipulation of geometric instances.

Example of a MayaGeometry that stores a tetrahedron:
\code
QVector<Vector> vertex;
QVector<Vector> normal;
QVector<int> index;
vertex.append(Vector(0,0,0));
vertex.append(Vector(1,0,0));
vertex.append(Vector(0,1,0));
vertex.append(Vector(0,0,1));
normal.append(Vector(-1,0,0));
normal.append(Vector(0,-1,0));
normal.append(Vector(0,0,-1));
normal.append(Normalized(Vector(1,1,1)));
// Set the indexes defining the oblique triangle
index.append(1); // Index of first vertex
index.append(3); // Oblique normal
index.append(2); // Second
index.append(3); // Same oblique normal
index.append(3); // Third
index.append(3); // Again the same normal
//  and so on for the other triangles
MayaGeometry object("tetrahedron",vertex,normal,index);
\endcode
Note that in this case, the vertices and normals are uniquely defined
and correctly shared.
Another possibility is to create the tetrahedron by adding triangles
directly to the MayaGeometry structure. In that case however, the
object is a triangle soup: points (and normals) are replicated.
\code
MayaGeometry object("tetrahedron");
object.AddTriangle(Vector(1,0,0),Vector(0,1,0),Vector(0,0,1)); // ... and so on for the other triangles
\endcode
\ingroup MayaCore
*/

/*!
\brief Creates an geometry from a mesh.
\param name Name of the object.
\param mesh The mesh.
\param mo Material.
*/
MayaGeometry::MayaGeometry(const QString& name, const Mesh& mesh, const MayaMaterial& mo) :name(name), mat(mo)
{
  MayaGeometry::vertex = mesh.GetVertices();
  MayaGeometry::normal = mesh.GetNormals();

  indexes.resize(mesh.Triangles() * 3);

  for (int i = 0; i < mesh.Triangles(); i++)
  {
    indexes[i * 3 + 0] = MayaIndexVertexData(mesh.VertexIndex(i, 0), mesh.NormalIndex(i, 0));
    indexes[i * 3 + 1] = MayaIndexVertexData(mesh.VertexIndex(i, 1), mesh.NormalIndex(i, 1));
    indexes[i * 3 + 2] = MayaIndexVertexData(mesh.VertexIndex(i, 2), mesh.NormalIndex(i, 2));
  }
}

/*!
\brief Creates an geometry from a mesh.
\param name Name of the object.
\param mesh The mesh.
\param mo Material.
*/
MayaGeometry::MayaGeometry(const QString& name, const MeshColor& mesh, const MayaMaterial& mo) :name(name), mat(mo)
{
  MayaGeometry::vertex = mesh.GetVertices();
  MayaGeometry::normal = mesh.GetNormals();

  MayaGeometry::color.resize(mesh.GetColors().size());
  for (int i = 0; i < mesh.GetColors().size(); i++)
  {
    Color c = mesh.GetColor(i);
    MayaGeometry::color[i] = Vector(c[0], c[1], c[2]);
  }

  indexes.resize(mesh.Triangles() * 3);

  for (int i = 0; i < mesh.Triangles(); i++)
  {
    indexes[i * 3 + 0] = MayaIndexVertexData(mesh.VertexIndex(i, 0), mesh.NormalIndex(i, 0), mesh.ColorIndex(i, 0));
    indexes[i * 3 + 1] = MayaIndexVertexData(mesh.VertexIndex(i, 1), mesh.NormalIndex(i, 1), mesh.ColorIndex(i, 1));
    indexes[i * 3 + 2] = MayaIndexVertexData(mesh.VertexIndex(i, 2), mesh.NormalIndex(i, 2), mesh.ColorIndex(i, 2));
  }
}

/*!
\brief Creates an geometry from a mesh.
\param name Name of the object.
\param mesh The mesh.
\param mo Material.
*/
MayaGeometry::MayaGeometry(const QString& name, const Mesh2& mesh, const MayaMaterial& mo) :name(name), mat(mo)
{
  vertex.resize(mesh.VertexSize());
  indexes.resize(mesh.TriangleSize() * 3);
  normal.fill(Vector::Z, mesh.VertexSize());

  for (int i = 0; i < mesh.VertexSize(); i++)
  {
    vertex[i] = mesh.Vertex(i).ToVector();
  }

  for (int i = 0; i < mesh.TriangleSize(); i++)
  {
    indexes[i * 3 + 0] = MayaIndexVertexData(mesh.index(i, 0), mesh.index(i, 0));
    indexes[i * 3 + 1] = MayaIndexVertexData(mesh.index(i, 1), mesh.index(i, 1));
    indexes[i * 3 + 2] = MayaIndexVertexData(mesh.index(i, 2), mesh.index(i, 2));
  }

  // MayaGeometry::mat.materialNode = None;
}

/*!
\brief Creates a geometry from a set of triangles.
\param name Name of the object.
\param tris the triangles.
\param mo Material.
*/
MayaGeometry::MayaGeometry(const QString& name, const QVector<Triangle>& tris, const MayaMaterial& mo) :name(name), mat(mo)
{
  for (int i = 0; i < tris.size(); i++)
  {
    AddTriangle(tris[i][0], tris[i][1], tris[i][2]);
  }
  // MayaGeometry::mat.materialNode = None;
}

/*!
\brief Creates an empty mesh structure.
\param name Name of the object.
\param mo Material.
*/
MayaGeometry::MayaGeometry(const QString& name, const MayaMaterial& mo) :name(name), mat(mo)
{
}

/*!
\brief Empty.
*/
MayaGeometry::~MayaGeometry()
{
}

/*!
\brief Creates a mesh given an array of vertices and normals, and an array of integers
defining the (possibly smooth) triangles of the mesh.

The array of integers stores the index of the vertices interlaced with the indexes of the normals.
\param name Name of the object.
\param vertex Array of vertices.
\param normal Array of normals.
\param index Indexes.
\param mo Material.
*/
MayaGeometry::MayaGeometry(const QString& name, const QVector<Vector>& vertex, const QVector<Vector>& normal, const QVector<MayaIndexVertexData>& index, const MayaMaterial& mo) :name(name), mat(mo)
{
  MayaGeometry::vertex = vertex;
  MayaGeometry::normal = normal;
  MayaGeometry::indexes = index;

  // MayaGeometry::mat = mo;
  // MayaGeometry::mat.materialNode = None;
}

/*!
\brief Create a mesh.

This function is the same as the constructor with QVector parameters, but uses C++ arrays instead.

\param name Name of the object.
\param vertex Array of vertices.
\param nv Number of vertices.
\param normal Array of normals.
\param nn Number of normals.
\param index Indexes defining the (possibly smooth) triangles of the mesh.
\param ni Number of indexes.
\param mo Material.
*/
MayaGeometry::MayaGeometry(const QString& name, const Vector* vertex, int nv, const Vector* normal, int nn, const int* index, int ni, const MayaMaterial& mo) :name(name), mat(mo)
{

  // Set size
  MayaGeometry::vertex.resize(nv);
  MayaGeometry::normal.resize(nn);
  MayaGeometry::indexes.resize(ni / 2);

  // Copy vertices
  for (int i = 0; i < nv; i++)
  {
    MayaGeometry::vertex[i] = vertex[i];
  }
  // Copy 
  for (int i = 0; i < nn; i++)
  {
    MayaGeometry::normal[i] = normal[i];
  }
  // Copy 
  for (int i = 0; i < ni; i += 2)
  {
    MayaIndexVertexData tmpindex;
    tmpindex.setIV(index[i + 0]);
    tmpindex.setIN(index[i + 1]);

    MayaGeometry::indexes[i / 2] = tmpindex;
  }
  //  MayaGeometry::mat.materialNode = None;
}

/*!
\brief Creates a mesh.
\param name Name of the object.
\param vertex Array of vertices.
\param normal Array of normals.
\param color  Array of colors.
\param index Array of integers defining the (possibly smooth) triangles of the mesh.
\param mo Material.
*/
MayaGeometry::MayaGeometry(const QString& name, const QVector<Vector>& vertex, const QVector<Vector>& normal, const QVector<Vector>& color, const QVector<MayaIndexVertexData>& index, const MayaMaterial& mo) :name(name), mat(mo)
{
  MayaGeometry::vertex.clear();
  MayaGeometry::normal.clear();
  MayaGeometry::color.clear();
  MayaGeometry::UVmap.clear();
  MayaGeometry::indexes.clear();

  MayaGeometry::vertex = vertex;
  MayaGeometry::normal = normal;
  MayaGeometry::color = color;
  MayaGeometry::indexes = index;

  MayaGeometry::mat.materialNode = VertexColor;
}

/*!
\brief Creates a mesh given an array of vertices, normals and colors, and an array of integers defining the (possibly smooth) triangles of the mesh.
\param name Name of the object.
\param vertex Array of vertices.
\param normal Array of normals.
\param UVmap  Array of UVmaps.
\param index Indexes.
\param mo Material.
*/
MayaGeometry::MayaGeometry(const QString& name, const QVector<Vector>& vertex, const QVector<Vector>& normal, const QVector<Vector2>& UVmap, const QVector<MayaIndexVertexData>& index, const MayaMaterial& mo) :name(name), mat(mo)
{
  MayaGeometry::vertex = vertex;
  MayaGeometry::normal = normal;
  MayaGeometry::UVmap = UVmap;
  MayaGeometry::indexes = index;

  MayaGeometry::mat.materialNode = UVMapping;
}

/*!
\brief Get the i-th triangle from the mesh as a flat triangle.

Note that the normals at the vertices of the triangle in the mesh structure is not taken into account.
\param i Index.
*/
Triangle MayaGeometry::GetTriangle(int i) const
{
  return Triangle(vertex[indexes[i * 3 + 0].getIV()], vertex[indexes[i * 3 + 1].getIV()], vertex[indexes[i * 3 + 2].getIV()]);
}

/*!
\brief Returns the geometry as a set of triangles.
*/
QVector<Triangle> MayaGeometry::GetTriangles() const
{
  QVector<Triangle> triangles;

  for (int i = 0; i < indexes.size(); i += 3)
  {
    triangles.append(Triangle(vertex[indexes[i + 0].getIV()], vertex[indexes[i + 1].getIV()], vertex[indexes[i + 2].getIV()]));
  }

  return triangles;
}

/*!
\brief Adds a triangle to the mesh.
\param a,b,c Vertices.
*/
void MayaGeometry::AddTriangle(const Vector& a, const Vector& b, const Vector& c)
{
  // Append triangle
  AddTriangle(vertex.size(), vertex.size() + 1, vertex.size() + 2, normal.size());

  // Append vertices 
  vertex.append(a);
  vertex.append(b);
  vertex.append(c);

  // Append normal
  normal.append(Normalized((b - a) / (c - a)));
}

/*!
\brief Adds a smooth triangle to the mesh.
\param a, b, c Vertices
\param na, nb, nc Normals
*/
void MayaGeometry::AddSmoothTriangle(const Vector& a, const Vector& b, const Vector& c, const Vector& na, const Vector& nb, const Vector& nc)
{
  // Append triangle
  AddSmoothTriangle(vertex.size(), normal.size(), vertex.size() + 1, normal.size() + 1, vertex.size() + 2, normal.size() + 2);

  // Append vertices 
  vertex.append(a);
  vertex.append(b);
  vertex.append(c);

  // Append normals
  normal.append(na);
  normal.append(nb);
  normal.append(nc);
}

/*!
\brief Add a smooth triangle to the geometry.
\param a, b, c Index of the vertices.
\param na, nb, nc Index of the normals.
*/
void MayaGeometry::AddSmoothTriangle(int a, int na, int b, int nb, int c, int nc)
{
  MayaIndexVertexData mi1(a, na);
  indexes.append(mi1);
  MayaIndexVertexData mi2(b, nb);
  indexes.append(mi2);
  MayaIndexVertexData mi3(c, nc);
  indexes.append(mi3);
}


/*!
\brief Add a triangle to the geometry.
\param a, b, c Index of the vertices.
\param n Index of the normal shared by all vertices.
*/
void MayaGeometry::AddTriangle(int a, int b, int c, int n)
{
  MayaIndexVertexData mi1(a, n);
  indexes.append(mi1);
  MayaIndexVertexData mi2(b, n);
  indexes.append(mi2);
  MayaIndexVertexData mi3(c, n);
  indexes.append(mi3);
}


/*!
\brief Add a quadrangle to the geometry.

It creates two triangles abc and acd.
\param a, b, c, d  Index of the vertices.
\param n Index of the normal shared by all vertices.
*/
void MayaGeometry::AddQuadrangle(int a, int b, int c, int d, int n)
{
  // First triangle
  AddTriangle(a, b, c, n);

  // Second triangle
  AddTriangle(a, c, d, n);
}

/*!
\brief Merges an argument mesh with the existing mesh.

Note that this function may be slow as it may require
a time consuming copy of several arrays.

\param mesh Geometry that will be merged.
*/
void MayaGeometry::Merge(const MayaGeometry& mesh)
{
  // Set size
  vertex.resize(vertex.size() + mesh.vertex.size());
  normal.resize(normal.size() + mesh.normal.size());
  color.resize(color.size() + mesh.color.size());
  UVmap.resize(UVmap.size() + mesh.UVmap.size());

  indexes.resize(indexes.size() + mesh.indexes.size());

  // Merge vertices, normals, colors, UVCoord 
  for (int i = 0; i < mesh.vertex.size(); i++)
  {
    vertex[vertex.size() - mesh.vertex.size() + i] = mesh.vertex[i];
  }
  for (int i = 0; i < mesh.normal.size(); i++)
  {
    normal[normal.size() - mesh.normal.size() + i] = mesh.normal[i];
  }
  for (int i = 0; i < mesh.color.size(); i++)
  {
    color[color.size() - mesh.color.size() + i] = mesh.color[i];
  }
  for (int i = 0; i < mesh.UVmap.size(); i++)
  {
    UVmap[UVmap.size() - mesh.UVmap.size() + i] = mesh.UVmap[i];
  }

  // Indexes should be shifted
  for (int i = 0; i < mesh.indexes.size(); i++)
  {
    MayaIndexVertexData mi = mesh.indexes[i];
    mi.setIV(mi.getIV() + vertex.size() - mesh.vertex.size());
    mi.setIN(mi.getIN() + normal.size() - mesh.normal.size());
    mi.setIC(mi.getIC() != -1 ? (mi.getIC() + color.size() - mesh.color.size()) : -1);
    mi.setIUV(mi.getIUV() != -1 ? (mi.getIUV() + UVmap.size() - mesh.UVmap.size()) : -1);

    indexes[indexes.size() - mesh.indexes.size() + i] = mi;
  }
}

/*!
\brief Compute the geometry information to Maya index lists.
\param te Edge indexes
\param e Vertex indexes
\param n Normal indexes
*/
void MayaGeometry::ConvertToVertexEdgeFaceMaya(QVector<int>& te, QVector<int>& col, QVector<int>& e, QVector<Vector>& n) const
{
  te.clear();
  e.clear();
  n.clear();
  col.clear();

  // Get Edge
  for (int i = 0; i < indexes.size(); i = i + 3)
  {
    int iv1 = indexes[i + 0].getIV();
    int iv2 = indexes[i + 1].getIV();
    int iv3 = indexes[i + 2].getIV();

    int ie1, ie2, ie3;
    {ie1 = e.size(); e.append(iv1); e.append(iv2); e.append(0); }

    {ie2 = e.size(); e.append(iv2); e.append(iv3); e.append(0); }

    {ie3 = e.size(); e.append(iv3); e.append(iv1); e.append(0); }

    te.append(ie1 / 3);
    te.append(ie2 / 3);
    te.append(ie3 / 3);

    col.append(iv1);
    col.append(iv2);
    col.append(iv3);
  }

  // Transform Triangles
  for (int i = 0; i < te.size(); i = i + 3)
  {
    int ie1 = te[i + 0];
    int ie2 = te[i + 1];
    int ie3 = te[i + 2];

    if (!(e[ie1 + 1] == e[ie2 + 0] || e[ie1 + 1] == e[ie2 + 1])) te[i + 0] = -(ie1 + 1);
    if (!(e[ie2 + 1] == e[ie3 + 0] || e[ie2 + 1] == e[ie3 + 1])) te[i + 1] = -(ie2 + 1);
    if (!(e[ie3 + 1] == e[ie1 + 0] || e[ie3 + 1] == e[ie1 + 1])) te[i + 2] = -(ie3 + 1);
  }

  // Get Normal
  for (int i = 0; i < indexes.size(); i = i + 3)
  {
    int in1 = indexes[i + 0].getIN();
    int in2 = indexes[i + 1].getIN();
    int in3 = indexes[i + 2].getIN();

    n.append(normal[in1]);
    n.append(normal[in2]);
    n.append(normal[in3]);
  }
}

/*!
\brief Transforms the geometry given a transformation operator.
\param frame Transformation.
*/
MayaGeometry& MayaGeometry::Transform(const FrameScaled& frame)
{
  // Apply transformation to vertices
  for (int i = 0; i < vertex.size(); i++)
  {
    vertex[i] = frame.Transform(vertex[i]);
  }

  // Compute normal
  for (int i = 0; i < normal.size(); i++)
  {
    normal[i] = Normalized(frame.TransformDirection(normal[i]));
  }

  return *this;
}

/*!
\brief Inverse transforms the geometry given a transformation operator.
\param frame Transformation.
*/
MayaGeometry& MayaGeometry::InverseTransform(const FrameScaled& frame)
{
  // Apply transformation to vertices
  for (int i = 0; i < vertex.size(); i++)
  {
    vertex[i] = frame.InverseTransform(vertex[i]);
  }

  // Compute normal
  for (int i = 0; i < normal.size(); i++)
  {
    normal[i] = Normalized(frame.InverseTransformDirection(normal[i]));
  }

  return *this;
}

/*!
\brief Compute the bounding box of the object.
*/
Box MayaGeometry::GetBox() const
{
  int n = vertex.size();
  if (n > 0)
  {
    // Starting box
    Box box = Box(vertex[0]);

    // Extend
    for (int i = 1; i < n; i++)
    {
      box.Extend(vertex[i]);
    }

    return box;
  }
  else
  {
    return Box::Null;
  }
}

/*!
\brief Clear every data and set a neutral grey material
*/
void MayaGeometry::Clear()
{
  name = "NoName";
  vertex.clear();
  normal.clear();
  color.clear();
  UVmap.clear();
  indexes.clear();

  mat.shaderNode = ShaderPhong;
  mat.materialNode = None;
  mat.ambient = Color(0.5, 0.5, 0.5);
  mat.diffuse = Color(0.5, 0.5, 0.5);
  mat.specular = Color(0.5, 0.5, 0.5);
  mat.shininess = 10;

  mat.texture.clear();
}

/*!
\brief Compute the statistics of the object.

This member computes the number of vertices and triangles of the object.
Since vertices may not be shared, those numbers may be larger than the real
number of vertices of the geometric object.

\return The statistics.
*/
MayaStatistics MayaGeometry::GetStatistics() const
{
  return MayaStatistics(1, vertex.size(), indexes.size() / 3, 0, 0, 0);
}


/*!
\brief Apply a planar mapping onto the object.

Delete ancient UVmaps.
\param box Domain.
*/
void MayaGeometry::generatePlannarZ_Mapping(Box2 box)
{
  UVmap.clear();
  Vector2 size = box.Size();

  for (int i = 0; i < vertex.size(); i++)
  {
    UVmap << Vector2((vertex[i][0] - box[0][0]) / size[0], 1.0 - (vertex[i][1] - box[0][1]) / size[1]);
  }

  for (int i = 0; i < indexes.size(); i++)
  {
    indexes[i].setIUV(indexes[i].getIV());
  }
}

/*!
\brief Apply a planar mapping onto the object using a given normal direction

Delete old UVmaps.
\param dir Normal direction.
\param step Size in spatial coordinates of the tiling texture
*/
void MayaGeometry::generatePlannarMapping(const Vector& dir, const double& step)
{
  UVmap.clear();
  Box box = GetBox();
  //Vector2 size = GetBox().Size();
  Vector x, y;
  dir.Orthonormal(x, y);
  x /= step;
  y /= step;
  for (int i = 0; i < vertex.size(); i++)
  {
    UVmap << Vector2(vertex[i] * x, vertex[i] * y);
  }

  for (int i = 0; i < indexes.size(); i++)
  {
    indexes[i].setIUV(indexes[i].getIV());
  }
}


/*!
\brief Get the text information.
\param spaces Spacing for indenting the text.
\param html Flag to specify syntax highlighting.
*/
QString MayaGeometry::GetText(int spaces, bool html) const
{
  return QString(spaces, ' ') + (html == true ? QString("<span style=\"color:#248\">") : QString("")) + QString("MayaGeometry") + (html == true ? QString("</span>") : QString("")) + QString("( ") + GetName() + QString(" )");
}


/*
\brief Create .cpp with static data from .obj
\param className corresponding class name
*/
void MayaGeometry::WriteCpp(const QString& className)
{
  for (int i = 0; i < vertex.size(); i++)
  {
    vertex[i][0] *= -1;
  }

  QString writeString;

  // Vertices
  writeString += QString("QList<Vector>" + className + "VL=QList<Vector>()\n");
  for (int i = 0; i < vertex.size(); i++)
  {
    writeString += "<<" + vertex[i].ToString(3) + "\n";
  }
  writeString += QString(";\n");
  writeString += QString("const QVector<Vector> " + className + "::vertex=QVector<Vector>::fromList(" + className + "VL);\n");
  writeString += QString("\n");

  // Normals
  writeString += QString("QList<Vector>" + className + "NL=QList<Vector>()\n");
  for (int i = 0; i < normal.size(); i++)
  {
    Vector tempNormal = normal.at(i);
    tempNormal[0] *= -1.0;
    writeString += "<<" + tempNormal.ToString(3) + "\n";
  }
  writeString += QString(";\n");
  writeString += QString("const QVector<Vector> " + className + "::normal=QVector<Vector>::fromList(" + className + "NL);\n");
  writeString += QString("\n");

  // Vertex Index
  writeString += QString("QList<MayaIndexVertexData>" + className + "IL=QList<MayaIndexVertexData>()\n");

  for (int i = 0; i < indexes.size(); i++)
  {
    writeString += QString("<<MayaIndexVertexData(%1, %2)\n").arg(QString::number(indexes[i].getIV()), QString::number(indexes[i].getIN()));
  }
  writeString += QString(";\n");
  writeString += QString("const QVector<MayaIndexVertexData> " + className + "::index=QVector<MayaIndexVertexData>::fromList(" + className + "IL);\n");

  QFile file("../" + className + ".cpp");
  file.open(QIODevice::WriteOnly | QIODevice::Text);
  file.write(writeString.toLatin1().data());
  file.close();
}

