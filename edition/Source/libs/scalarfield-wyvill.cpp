
#include "libs/scalarfield.h"
#include "libs/mesh.h"
#include "libs/ivector.h"

#include <unordered_map>
#include <queue>
#include <iostream>

// Points
//     5 ----- 6
//    /|     / |
//   1 ----- 2 |	z
//   | 4 ----|-7	| y
//   |/      |/		|/
//   0 ----- 3		0 --- x

// Faces
//     . ----- .
//    /|  3   /|
//   . ----5 .4|
//   |1. 0---|-.	
//   |/   2  |/		
//   . ----- .		


//               Wyvill   Polygonize   Factor
// Bear : 1000 :   3.81    43:92       11.5
//         700 :   1:98    15:44        7.8
//         500 :   0:69     3:03        4.4

class Cube
{
public:
  int x, y, z;	//!< Integer coordinates (from the first cube)
  double values[8]; //!< Vertex values
  double intersections[12]; // < -1.0 if no intersection, between 0.0 and 1.0 else (distance from the first vertice index)

  // Links between intersections
  //int links[12][4]; // -1 if no link with this intersection, the number of the second point if some link with another intersection
  // Maximum four link per intersection

  int intersectionsCoordinates[12];	// True intersection coordinates in the mesh vertexes list

  // Constructor
  Cube(int x, int y, int z);
  void setValue(int vertex, double value);
  double getValue(int vertex) const;


protected:

  // Cube vertices, defined between two vertexes
  static const int vertices[12][2];

  // Vertexes (points) on each face
  static const int facesP[6][4];

  // Vertices on each face
  static const int facesV[6][4];

public:

  void setIntersection(int vertice, double value)
  {
    //std::cout << "Intersection vertice nÂ° " << vertice << " : " << value << std::endl;
    intersections[vertice] = value;
  }

  // Compute intersections by linear interpolation between the values
  double computeIntersections(int vertice)
  {
    // Intersection
    if (values[vertices[vertice][0]] * values[vertices[vertice][1]] <= 0.0)
    {
      // Linear interpolation between the values
      if (values[vertices[vertice][0]] <= 0.0)
      {
        intersections[vertice] = -values[vertices[vertice][0]] / (values[vertices[vertice][1]] - values[vertices[vertice][0]]);
      }
      else
      {
        intersections[vertice] = 1.0 + values[vertices[vertice][1]] / (values[vertices[vertice][0]] - values[vertices[vertice][1]]);
      }
    }
    //std::cout << "Intersection vertice nÂ° " << vertice << " : " << intersections[vertice] << std::endl;
    return intersections[vertice];
  }

  double getIntersection(int vertice)
  {
    return intersections[vertice];
  }

  void setIntersectionCoordinates(int vertice, int coordinate)
  {
    //std::cout << "True intersection vertice nÂ° " << vertice << " : " << coordinate << std::endl;
    intersectionsCoordinates[vertice] = coordinate;
  }

  int getIntersectionCoordinates(int vertice)
  {
    return intersectionsCoordinates[vertice];
  }

  // Return the two vertexes corresponding to the vertice
  void getVerticesPoints(int vertice, int& a, int& b)
  {
    a = vertices[vertice][0];
    b = vertices[vertice][1];
  }

  // Return the values of the two vertexes corresponding to the vertice
  void getVerticesValues(int vertice, double& a, double& b)
  {
    a = values[vertices[vertice][0]];
    b = values[vertices[vertice][1]];
  }

  // Compute the triangles
  QVector<int> computePolygones()
  {
    // Version 2: using TrianglesTable
    QVector<int> triangles;

    int cubeIndex = 0;
    if (values[0] < 0.0) cubeIndex |= 1;
    if (values[1] < 0.0) cubeIndex |= 4;
    if (values[2] < 0.0) cubeIndex |= 8;
    if (values[3] < 0.0) cubeIndex |= 2;
    if (values[4] < 0.0) cubeIndex |= 16;
    if (values[5] < 0.0) cubeIndex |= 64;
    if (values[6] < 0.0) cubeIndex |= 128;
    if (values[7] < 0.0) cubeIndex |= 32;

    if (cubeIndex != 255 && cubeIndex != 0)
    {
      for (int h = 0; ScalarField::TriangleTable[cubeIndex][h] != -1; h++)
      {
        switch (ScalarField::TriangleTable[cubeIndex][h])
        {
        case 0: triangles << 1;
          break;
        case 1: triangles << 3;
          break;
        case 2: triangles << 9;
          break;
        case 3: triangles << 10;
          break;
        case 4: triangles << 0;
          break;
        case 5: triangles << 5;
          break;
        case 6: triangles << 8;
          break;
        case 7: triangles << 11;
          break;
        case 8: triangles << 2;
          break;
        case 9: triangles << 7;
          break;
        case 10: triangles << 4;
          break;
        case 11: triangles << 6;
          break;
        default: break;
        }
      }
    }

    return triangles;

  }


  // Compute the triangles
  int computePolygones(int pppp[20])
  {
    // Version 2: using TrianglesTable
    int nbp = 0;

    int cubeIndex = 0;
    if (values[0] < 0.0) cubeIndex |= 1;
    if (values[1] < 0.0) cubeIndex |= 4;
    if (values[2] < 0.0) cubeIndex |= 8;
    if (values[3] < 0.0) cubeIndex |= 2;
    if (values[4] < 0.0) cubeIndex |= 16;
    if (values[5] < 0.0) cubeIndex |= 64;
    if (values[6] < 0.0) cubeIndex |= 128;
    if (values[7] < 0.0) cubeIndex |= 32;

    if (cubeIndex != 255 && cubeIndex != 0)
    {
      for (int h = 0; ScalarField::TriangleTable[cubeIndex][h] != -1; h++)
      {
        switch (ScalarField::TriangleTable[cubeIndex][h])
        {
        case 0: pppp[nbp] = 1;
          break;
        case 1: pppp[nbp] = 3;
          break;
        case 2: pppp[nbp] = 9;
          break;
        case 3: pppp[nbp] = 10;
          break;
        case 4: pppp[nbp] = 0;
          break;
        case 5: pppp[nbp] = 5;
          break;
        case 6: pppp[nbp] = 8;
          break;
        case 7: pppp[nbp] = 11;
          break;
        case 8: pppp[nbp] = 2;
          break;
        case 9: pppp[nbp] = 7;
          break;
        case 10: pppp[nbp] = 4;
          break;
        case 11: pppp[nbp] = 6;
          break;
        default: break;
        }
        nbp++;
      }
    }

    return nbp;
  }
};

Cube::Cube(int x, int y, int z) : x(x), y(y), z(z)
{
  // Initialisations
  for (int i = 0; i < 8; i++)
  {
    values[i] = 0.0;
  }
  for (int i = 0; i < 12; i++)
  {
    intersections[i] = -1.0;
    //links[i][0] = -1;
    //links[i][1] = -1;
    //links[i][2] = -1;
    //links[i][3] = -1;
    intersectionsCoordinates[i] = -1;
  }
}

inline void Cube::setValue(int vertex, double value)
{
  //std::cout << "Value vertex nÂ° " << sommet << " : " << value << std::endl;
  values[vertex] = value;
}

inline double Cube::getValue(int vertex) const
{
  return values[vertex];
}

const int Cube::vertices[12][2] = {
    {0, 1}, //0
    {0, 3}, //1
    {0, 4}, //2
    {1, 2}, //3
    {1, 5}, //4
    {2, 3}, //5
    {2, 6}, //6
    {3, 7}, //7
    {4, 5}, //8
    {4, 7}, //9
    {5, 6}, //10
    {6, 7}  //11
};

const int Cube::facesP[6][4] = {
    {0, 1, 2, 3},
    {4, 5, 1, 0},
    {4, 0, 3, 7},
    {1, 5, 6, 2},
    {3, 2, 6, 7},
    {7, 6, 5, 4}
};

const int Cube::facesV[6][4] = {
    {0, 3, 5, 1},
    {8, 4, 0, 2},
    {2, 1, 7, 9},
    {4, 10, 6, 3},
    {5, 6, 11, 7},
    {11, 10, 8, 9}
};

// Compute the id of the cube of the indicated face
int getCubeId(int face, int x, int y, int z)
{
  int cubeId = 0;
  switch (face)
  {
  case 0: cubeId = (1024 * x + (y - 1)) * 1024 + z;
    break;
  case 1: cubeId = (1024 * (x - 1) + y) * 1024 + z;
    break;
  case 2: cubeId = (1024 * x + y) * 1024 + (z - 1);
    break;
  case 3: cubeId = (1024 * x + y) * 1024 + (z + 1);
    break;
  case 4: cubeId = (1024 * (x + 1) + y) * 1024 + z;
    break;
  case 5: cubeId = (1024 * x + (y + 1)) * 1024 + z;
    break;
  default: break;
  }
  return cubeId;
}

// Compute the true coordinates a vertex
inline Vector computeVertexPosition(Vec3I vertex, Vec3I pos, Vector start, Vector size)
{
  return start + Vector(((double)pos[0] + (double)vertex[0]) * size[0], ((double)pos[1] + (double)vertex[1]) * size[1], ((double)pos[2] + (double)vertex[2]) * size[2]);
}

// Compute the true coordinate of an intersection
inline Vector computeIntersectionPosition(Vec3I vertex1, Vec3I vertex2, double intersection, Vec3I pos, Vector start, Vector size)
{
  return computeVertexPosition(vertex1, pos, start, size) + intersection * Vector((double)(vertex2 - vertex1)[0] * size[0], (double)(vertex2 - vertex1)[1] * size[1], (double)(vertex2 - vertex1)[2] * size[2]);
}

bool boxIntersection(Box box, Vec3I pos, Vector start, Vector size)
{
  Vector a = start + Vector((double)pos[0] * size[0], (double)pos[1] * size[1], (double)pos[2] * size[2]);
  Vector b = a + size;
  return box.Intersect(Box(a, b));
}

/*!
\brief Implementation of the original continuation polygonization method.
\param box Starting cube straddling the surface.
\param mesh %Mesh.
\param dichotomy Boolean for computing vertexes using bisection.
\param computeNormals Boolean for computing normals.

The original method is described in: Data structure forsoft objects. Geoff Wyvill, Craig McPheeters, Brian Wyvill. <i>The Visual Computer</i>, <b>2</b>, 227â€“234, 1986.
*/
void AnalyticScalarField::PolygonizeLucie(const Box& box, Mesh& mesh, bool dichotomy, bool computeNormals) const
{
  std::unordered_map<int, Cube> cubes; // Contains the created cubes
  std::queue<Vec3I> unvisited; // Contains the coordinates of the cube to visit next
  std::queue<int> parent; // Contains the face of the parent cube of each cube to visit

  // Mesh vertexes, triangles and normals
  QVector<Vector> vertex;
  QVector<int> triangles;
  QVector<Vector> normals;

  vertex.reserve(200000);
  triangles.reserve(400000);
  normals.reserve(200000);

  Vector size = box.Size();

  Vector start = (box)[0];

  Box treeBox = GetBox();

  unvisited.push(Vec3I(0, 0, 0)); // First cube
  parent.push(-1); // This cube has no parent

  // Vertexes of each face
  int faceVertex[6][4] =
  {
    {0, 1, 2, 3},
    {0, 1, 4, 5},
    {0, 3, 4, 7},
    {1, 2, 5, 6},
    {2, 3, 6, 7},
    {4, 5, 6, 7}
  };

  // Vertexes of the cube at the indicated face
  int parentVertex[6][4] =
  {
    {4, 5, 6, 7},
    {3, 2, 7, 6},
    {1, 2, 5, 6},
    {0, 3, 4, 7},
    {1, 0, 5, 4},
    {0, 1, 2, 3}
  };

  // Faces vertice intersects
  int verticeFaces[12][2] =
  {
    {0, 1},
    {0, 2},
    {1, 2},
    {0, 3},
    {1, 3},
    {0, 4},
    {3, 4},
    {2, 4},
    {1, 5},
    {2, 5},
    {3, 5},
    {4, 5}
  };

  // Same as in cube
  // The vertices numbers ov each face
  int faceVertices[6][4] =
  {
    {0, 3, 5, 1},
    {8, 4, 0, 2},
    {2, 1, 7, 9},
    {4, 10, 6, 3},
    {5, 6, 11, 7},
    {11, 10, 8, 9}
  };

  // Coordinates of the vertices inside the cube
  Vec3I vertices[12][2] = {
      {Vec3I(0, 0, 0), Vec3I(0, 0, 1)},
      {Vec3I(0, 0, 0), Vec3I(1, 0, 0)},
      {Vec3I(0, 0, 0), Vec3I(0, 1, 0)},
      {Vec3I(0, 0, 1), Vec3I(1, 0, 1)},
      {Vec3I(0, 0, 1), Vec3I(0, 1, 1)},
      {Vec3I(1, 0, 1), Vec3I(1, 0, 0)},
      {Vec3I(1, 0, 1), Vec3I(1, 1, 1)},
      {Vec3I(1, 0, 0), Vec3I(1, 1, 0)},
      {Vec3I(0, 1, 0), Vec3I(0, 1, 1)},
      {Vec3I(0, 1, 0), Vec3I(1, 1, 0)},
      {Vec3I(0, 1, 1), Vec3I(1, 1, 1)},
      {Vec3I(1, 1, 1), Vec3I(1, 1, 0)}
  };

  // Coordinates of each vertex of the cube
  Vec3I cubeVertex[8] = {
    Vec3I(0, 0, 0),
    Vec3I(0, 0, 1),
    Vec3I(1, 0, 1),
    Vec3I(1, 0, 0),
    Vec3I(0, 1, 0),
    Vec3I(0, 1, 1),
    Vec3I(1, 1, 1),
    Vec3I(1, 1, 0)
  };

  // Corresponding parent vertices
  int parentVertices[6][12] =
  {
    {8, 9, -1, 10, -1, 11, -1, -1, -1, -1, -1, -1},
    {5, -1, 7, -1, 6, -1, -1, -1, 11, -1, -1, -1},
    {-1, 3, 4, -1, -1, -1, -1, 6, -1, 10, -1, -1},
    {-1, -1, -1, 1, 2, -1, 7, -1, -1, -1, 9, -1},
    {-1, -1, -1, -1, -1, 0, 4, 2, -1, -1, -1, 8},
    {-1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 3, 5},
  };

  int nbCubes = 0;	// Number of computed cubes

  // Continue while there are still cubes to check
  while (!unvisited.empty())
  {
    // Cubes position (from the first cube)
    Vec3I pos = unvisited.front();
    int x = pos[0];
    int y = pos[1];
    int z = pos[2];

    //std::cout << "Unvisited : cube nÂ° " << x << " ; " << y << " ; " << z << std::endl;

    // Cube ID in the map
    int id = (1024 * x + y) * 1024 + z;

    // Check if the cube is alreeady treated or not and if the cube is inside de tree box
    if (cubes.count(id) == 0 && boxIntersection(treeBox, pos, start, size))
    {
      //std::cout << "Nouveau cube ! Cube nÂ°" << nbCubes << std::endl;
      Cube cube(x, y, z);
      nbCubes++;

      // Check of the already computed vertex
      int parentFace = parent.front();
      bool toCompute[8] = { true, true, true, true, true, true, true, true }; // Vertexes to compute
      int intersectionsToCompute[12] = { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }; // Intersections to compute
      bool parentToCompute[6] = { true, true, true, true, true, true }; // Cubes to add to the list

      // Eliminate the first cube
      if (parentFace != -1)
      {
        //std::cout << "Parent face : " << parentFace << std::endl;

        // If it's the parent cube, vertexes and vertices of this face are already calculated, no need to check
        for (int j = 0; j < 4; j++)
        {
          int sommet = faceVertex[parentFace][j];
          toCompute[sommet] = false;
          cube.setValue(sommet, cubes.at(getCubeId(parentFace, x, y, z)).getValue(parentVertex[parentFace][j]));
          intersectionsToCompute[faceVertices[parentFace][j]] = parentFace;
        }
        parentToCompute[parentFace] = false;

        // Check each neighboor cube
        for (int i = 0; i < 6; i++)
        {
          if (i != parentFace)
          {
            int cubeId = getCubeId(i, x, y, z);

            // Check if this neightboor is already computed
            if (cubes.count(cubeId) > 0)
            {
              //std::cout << "Un autre parent a ete trouve sur la face : " << i << std::endl;

              // If it's already computed, put the values of the vertices and intersections
              for (int j = 0; j < 4; j++)
              {
                int sommet = faceVertex[i][j];
                if (toCompute[sommet])
                {
                  toCompute[sommet] = false;
                  cube.setValue(sommet, cubes.at(cubeId).getValue(parentVertex[i][j]));
                }
                intersectionsToCompute[faceVertices[i][j]] = i;
              }
              parentToCompute[i] = false;
            }
          }
        }
      }

      // Compute the vertexex not already computed
      for (int i = 0; i < 8; i++)
      {
        if (toCompute[i])
        {
          //std::cout << "Le sommet " << i << " doit etre calcule." << std::endl;
          //std::cout << "Sa position relle est " << computeVertexPosition(cubeVertex[i], pos, start, size) << std::endl;

          cube.setValue(i, Value(computeVertexPosition(cubeVertex[i], pos, start, size)));
        }
      }

      // Compute the intersections
      int intersectedFaces[6] = { false, false, false, false, false, false }; // Face with a new intersection
      for (int i = 0; i < 12; i++)
      {
        // New intersection
        if (intersectionsToCompute[i] == -1)
        {
          //std::cout << "Intersection a calculer en " << i  << std::endl;
          double intersection = -1.0;
          Vector intersectionPosition;
          if (dichotomy)
          {
            // Dichotomy computing
            double va, vb;
            cube.getVerticesValues(i, va, vb);
            if (va * vb <= 0.0)
            {
              int sa, sb;
              cube.getVerticesPoints(i, sa, sb);
              Vector a = computeVertexPosition(cubeVertex[sa], pos, start, size);
              if (i == 0 || i == 8)
              {
                Vector b = a + Vector(0.0, 0.0, size[2]);
                double length = size[2];
                intersectionPosition = Dichotomy(a, b, va, vb, length);
                intersection = (intersectionPosition[2] - a[2]) / length;
              }
              else if (i == 5 || i == 11)
              {
                Vector b = a - Vector(0.0, 0.0, size[2]);
                double length = size[2];
                intersectionPosition = Dichotomy(a, b, va, vb, length);
                intersection = (a[2] - intersectionPosition[2]) / length;
              }
              else if (i == 1 || i == 3 || i == 9 || i == 10)
              {
                Vector b = a + Vector(size[0], 0.0, 0.0);
                double length = size[0];
                intersectionPosition = Dichotomy(a, b, va, vb, length);
                intersection = (intersectionPosition[0] - a[0]) / length;
              }
              else
              {
                Vector b = a + Vector(0.0, size[1], 0.0);
                double length = size[1];
                intersectionPosition = Dichotomy(a, b, va, vb, length);
                intersection = (intersectionPosition[1] - a[1]) / length;
              }

              cube.setIntersection(i, intersection);
            }
          }
          else
          {
            // Linear computing
            intersection = cube.computeIntersections(i);
            intersectionPosition = computeIntersectionPosition(vertices[i][0], vertices[i][1], intersection, pos, start, size);
          }

          // If there is an intersection
          if (intersection >= 0.0)
          {
            // Check intersected faces to add next cubes to the list
            intersectedFaces[verticeFaces[i][0]] = true;
            intersectedFaces[verticeFaces[i][1]] = true;

            //std::cout << "C'est une nouvelle intersection, de nÂ° " << vertex.size() << " et de position reelle " << intersectionPosition << std::endl;

            // Add the new vertex
            vertex << intersectionPosition;
            Vector origin = computeVertexPosition(Vec3I(0, 0, 0), pos, start, size);
            if (computeNormals)
            {
              normals << Normal(intersectionPosition);
            }
            // Indicates the coordinates (in the vertex list) of this intersection
            cube.setIntersectionCoordinates(i, vertex.size() - 1);
          }
        }
        else // Already computed intersection
        {
          //std::cout << "Intersection deja calculee en " << i << std::endl;

          // Get the intersection from the parent cube indicated
          int cubeId = getCubeId(intersectionsToCompute[i], x, y, z);
          double intersection = cubes.at(cubeId).getIntersection(parentVertices[intersectionsToCompute[i]][i]);
          if (intersection >= 0.0)
          {
            intersectedFaces[verticeFaces[i][0]] = true;
            intersectedFaces[verticeFaces[i][1]] = true;
            cube.setIntersection(i, intersection);
            cube.setIntersectionCoordinates(i, cubes.at(cubeId).getIntersectionCoordinates(parentVertices[intersectionsToCompute[i]][i]));
          }
        }
      }

      // Add the next cubes to visit
      for (int i = 0; i < 6; i++)
      {
        if (parentToCompute[i] && intersectedFaces[i])
        {
          //std::cout << "Cube de face " << i << " ajoute a la liste." << std::endl;
          switch (i)
          {
          case 0:
          {
            unvisited.push(pos - Vec3I(0, 1, 0));
            parent.push(5);
          }
          break;
          case 1:
          {
            unvisited.push(pos - Vec3I(1, 0, 0));
            parent.push(4);
          }
          break;
          case 2:
          {
            unvisited.push(pos - Vec3I(0, 0, 1));
            parent.push(3);
          }
          break;
          case 3:
          {
            unvisited.push(pos + Vec3I(0, 0, 1));
            parent.push(2);
          }
          break;
          case 4:
          {
            unvisited.push(pos + Vec3I(1, 0, 0));
            parent.push(1);
          }
          break;
          case 5:
          {
            unvisited.push(pos + Vec3I(0, 1, 0));
            parent.push(0);
          }
          break;
          default: break;
          }
        }
      }

      // Polygonize the cube
      /*
      QVector<int> polygones = cube.computePolygones();

      //std::cout << "Polygones du cube : ";

      // Add the triangles coordinates

      for (int i = 0; i < polygones.size(); i++)
      {
        //std::cout << polygones[i] << "(" << cube->getIntersectionCoordinates(polygones[i]) << ") ; ";
        triangles << cube.getIntersectionCoordinates(polygones[i]);
      }
      */
      int polygos[20];
      int npolygos = cube.computePolygones(polygos);

      //std::cout << "Polygones du cube : ";

      // Add the triangles coordinates
      for (int i = 0; i < npolygos; i++)
      {
        //std::cout << polygones[i] << "(" << cube->getIntersectionCoordinates(polygones[i]) << ") ; ";
        triangles << cube.getIntersectionCoordinates(polygos[i]);
      }

      // Put the computed cube in the map
      cubes.insert({ id, cube });
    }

    // Clear the cube and the parent
    unvisited.pop();
    parent.pop();
  }

  if (computeNormals)
  {
    mesh = Mesh(vertex, normals, triangles, triangles);
  }
  else
  {
    mesh = Mesh(vertex, triangles);
  }
  // Number of processed cubes
  ncubes = nbCubes;
}

/*!
\brief Wyvill's marching cubes algorithm

Find a cube straddling the surface, and apply Wyvill's marching algorithm.

\param size Size of the cubes that will be straddling the surface.
\param mesh %Mesh.
\param dichotomy Boolean for computing vertexes using bisection.
\param computeNormals Boolean for computing normals.
*/
void AnalyticScalarField::PolygonizeLucie(const double& size, Mesh& mesh, bool dichotomy, bool computeNormals) const
{
  Vector p;
  if (GetSample(p, GetBox()))
  {
    PolygonizeLucie(Box(p, size / 2.0), mesh, dichotomy, computeNormals);
  }
}
