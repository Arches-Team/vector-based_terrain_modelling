// Heightfield

#include "libs/heightfield.h"

/*!
\brief Create the geometry of the border of the terrain.

\param z The height of the base of the block.
*/
Mesh HeightField::CreateMeshSide(const double& z) const
{
  int n = 2 * nx + 2 * ny - 4;

  QVector<Vector> vertex;
  vertex.reserve(2 * n);
  QVector<Vector> normal;
  normal.resize(4);

  QVector<int> indexes;
  QVector<int> normals;
  indexes.reserve(6 * 2 * n);

  // Base contour
  for (int i = 0; i < nx - 1; i++)
  {
    Vector p = Vertex(i, 0);
    Vector q = p;
    p[2] = z;
    vertex.append(p);
    vertex.append(q);
  }

  for (int i = 0; i < ny - 1; i++)
  {
    Vector p = Vertex(nx - 1, i);
    Vector q = p;
    p[2] = z;
    vertex.append(p);
    vertex.append(q);
  }

  for (int i = nx - 1; i > 0; i--)
  {
    Vector p = Vertex(i, ny - 1);
    Vector q = p;
    p[2] = z;
    vertex.append(p);
    vertex.append(q);
  }

  for (int i = ny - 1; i > 0; i--)
  {
    Vector p = Vertex(0, i);
    Vector q = p;
    p[2] = z;
    vertex.append(p);
    vertex.append(q);
  }

  // Compute normals
  normal[0] = Vector(0.0, 1.0, 0.0);
  normal[1] = Vector(-1.0, 0.0, 0.0);
  normal[2] = Vector(0.0, -1.0, 0.0);
  normal[3] = Vector(1.0, 0.0, 0.0);

  // Compute indexes : loop over all vertices
  for (int i = 0; i < 2 * n; i += 2)
  {
    int side = 3;
    if ((i / 2 >= 0) && (i / 2 < nx - 1)) { side = 0; }
    else if ((i / 2 >= nx - 1) && (i / 2 < nx + ny - 2)) { side = 1; }
    else if ((i / 2 >= nx + ny - 2) && (i / 2 < 2 * nx + ny - 3)) { side = 2; }

    //Triangle 1
    indexes.append(i);
    indexes.append(i + 1);
    indexes.append((i + 3) % (2 * n));

    //Triangle 2
    indexes.append(i);
    indexes.append((i + 3) % (2 * n));
    indexes.append((i + 2) % (2 * n));

    normals.append(side);
    normals.append(side);
    normals.append(side);

    normals.append(side);
    normals.append(side);
    normals.append(side);
  }

  return Mesh(vertex, normal, indexes, normals);
}


/*!
\brief Create the geometry of the heightfield.

\param surface Generate the surface elevation of the heightfield.
\param side Generate the side of the heightfield.

\param z The height of the base of the block.
*/
Mesh HeightField::CreateMesh(bool surface, bool side, const double& z) const
{
  Mesh mesh;

  if (surface == true)
  {
    QVector<Vector> vertex;
    vertex.resize(nx * ny);
    QVector<Vector> normal;
    normal.resize(nx * ny);

    QVector<int> indexes;

    indexes.reserve(6 * (nx - 1) * (ny - 1));

    // Height of vertices
    for (int i = 0; i < nx; i++)
    {
      for (int j = 0; j < ny; j++)
      {
        vertex[i * ny + j] = Vertex(i, j);
      }
    }

    // Compute normals
    for (int i = 0; i < nx; i++)
    {
      for (int j = 0; j < ny; j++)
      {
        normal[i * ny + j] = Normal(i, j);
      }
    }

    // Compute indexes
    for (int i = 0; i < nx - 1; i++)
    {
      for (int j = 0; j < ny - 1; j++)
      {
        //Triangle 1
        indexes.append(i * ny + j);
        indexes.append((i + 1) * ny + j);
        indexes.append((i + 1) * ny + j + 1);

        //Triangle 2
        indexes.append(i * ny + j);
        indexes.append((i + 1) * ny + j + 1);
        indexes.append(i * ny + j + 1);
      }
    }
    QVector<int> normals = indexes;
    mesh = Mesh(vertex, normal, indexes, normals);
  }

  if (side == true)
  {
    Mesh sidemesh = CreateMeshSide(z);
    mesh.Append(sidemesh);
  }
  return mesh;
}

/*!
\brief Create a stack representation of the model.

Stacks are defined as scaled boxes.

\param cube Reference cube.
\param frames Set of instances.
\param z Base elevation.
*/
void HeightField::CreateCubes(Box& cube, QVector<FrameScaled>& frames, const double& z) const
{
  // Could use Box::Unit as well
  cube = Box(Vector(0.0, 0.0, 0.0), Vector(1.0, 1.0, 1.0));

  frames.reserve(nx * ny);

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      Vector2 p = ArrayVertex(i, j);
      double scale = at(i, j) - z;
      frames.append(FrameScaled(Matrix::Identity, p.ToVector(z), celldiagonal.ToVector(scale)));
    }
  }
}
