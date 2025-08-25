#include "libs/meshcolor.h"

/*!
\class MeshColor meshcolor.h

\brief Core triangle mesh class with color support.
\ingroup ExtendedKernelGroup
*/

/*!
\brief Create an empty mesh.
*/
MeshColor::MeshColor()
{
}

/*!
\brief Constructor from a Mesh with color array and indices.
\param m Base mesh.
\param cols Color array.
\param carr Color indexes.
*/
MeshColor::MeshColor(const Mesh& m, const QVector<Color>& cols, const QVector<int>& carr) : Mesh(m), colors(cols), carray(carr)
{
}

/*!
\brief Constructor from a Mesh with color array for every vertex.
\param m Base mesh.
\param cols Color array.
*/
MeshColor::MeshColor(const Mesh& m, const QVector<Color>& cols) : Mesh(m), colors(cols),carray(varray)
{
  //carray = varray;
}

/*!
\brief Constructor from a Mesh.
\param m the base mesh
*/
MeshColor::MeshColor(const Mesh& m) : Mesh(m)
{
  colors.resize(vertices.size());
  colors.fill(Color::Black);
  carray = varray;
}

/*!
\brief Empty.
*/
MeshColor::~MeshColor()
{
}

/*!
\brief Merge two meshes into a single one.
\param mesh Argument mesh.
*/
void MeshColor::Merge(const MeshColor& mesh)
{
  // Base class
  Mesh::Append(mesh);

  // Shifts
  const int c = colors.size();

  // Colors
  colors.append(mesh.colors);

  // Include shifting for colors
  for (int i = 0; i < mesh.varray.size(); i++)
  {
    carray.append(mesh.carray.at(i) + c);
  }
}
