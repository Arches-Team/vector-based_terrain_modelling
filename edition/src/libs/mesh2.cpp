// Mesh

#include "libs/mesh.h"
#include "libs/array.h"
#include "libs/ellipse.h"
#include "libs/polygon.h"
#include "libs/hexagon.h"

/*!
\class Mesh2 mesh.h

\brief Planar triangle mesh.
\ingroup PlanarGroup, ExtendedKernelGroup
*/

/*!
\brief Create an empty mesh.
*/
Mesh2::Mesh2()
{
}

/*!
\brief Create a mesh from a list of vertices and a indexes.

Indices should be a multiple of three; the number of triangles is derived
from the size of indices.

\sa Mesh2::TriangleSize

\param v Set of vertices.
\param i Indexes that represent the triangles.
*/
Mesh2::Mesh2(const QVector<Vector2>& v, const QVector<int>& i) :vertices(v), indices(i)
{
}

/*!
\brief Empty
*/
Mesh2::~Mesh2()
{
}

/*!
\brief Generate a grid geometry over an input quadrangle.
\param quad Quadrangle.
\param nx, ny Number of points in the subdivision.
*/
Mesh2::Mesh2(const Quadrangle2& quad, int nx, int ny)
{
  vertices.resize(nx * ny);
  indices.resize(3 * (nx - 1) * (ny - 1) * 2);

  // Vertices
  for (int j = 0; j < nx; j++)
  {
    double u = Math::Unit(j, nx);
    for (int i = 0; i < ny; i++)
    {
      double v = Math::Unit(i, ny);

      vertices[j * nx + i] = quad.Vertex(u, v);
    }
  }

  // Triangles
  int n = 0;
  for (int j = 0; j < nx - 1; j++)
  {
    for (int i = 0; i < ny - 1; i++)
    {
      indices[n + 0] = j * nx + i;
      indices[n + 1] = j * nx + i + 1;
      indices[n + 2] = (j + 1) * nx + i;
      n += 3;

      indices[n + 0] = j * nx + i + 1;
      indices[n + 1] = (j + 1) * nx + i + 1;
      indices[n + 2] = (j + 1) * nx + i;
      n += 3;
    }
  }
}

/*!
\brief Generate a grid geometry.
\param box The box.
\param x, y Number of points in the subdivision.
*/
Mesh2::Mesh2(const Box2& box, int x, int y)
{
  Array2 a(box, x, y);

  vertices.resize(x * y);
  indices.resize(3 * (x - 1) * (y - 1) * 2);

  // Vertices
  for (int j = 0; j < x; j++)
  {
    for (int i = 0; i < y; i++)
    {
      vertices[j * x + i] = a.ArrayVertex(i, j);
    }
  }

  // Triangles
  int n = 0;
  for (int j = 0; j < x - 1; j++)
  {
    for (int i = 0; i < y - 1; i++)
    {
      indices[n + 0] = j * x + i;
      indices[n + 1] = j * x + i + 1;
      indices[n + 2] = (j + 1) * x + i;
      n += 3;
      indices[n + 0] = j * x + i + 1;
      indices[n + 1] = (j + 1) * x + i + 1;
      indices[n + 2] = (j + 1) * x + i;
      n += 3;
    }
  }
}


static inline int A(int k, int i, int j, int n, int p)
{
  if (i == 0) return 0;
  else
  {
    if (j < i)
    {
      return 1 + (k % p) * (n * (n + 1) / 2) + i * (i - 1) / 2 + j;
    }
    else
    {
      return 1 + ((k + 1) % p) * (n * (n + 1) / 2) + i * (i - 1) / 2 + 0;
    }
  }
}


/*!
\brief Create an n-adic subdivision of a circle.

Creates an hexagon an applies an n-adic subdivision.

\param circle The circle.
\param n Subdivision level.
*/
Mesh2::Mesh2(const Circle2& circle, int n)
{
  Polygon2 hexagon(Hexagon2(circle.Center(), circle.Radius()));
  int s = 6;

  // Starting point.
  vertices.append(circle.Center().ToVector());

  // Append vertices
  for (int k = 0; k < 6; k++)
  {
    Vector2 e01 = (hexagon.Vertex(k) - circle.Center()) / n;
    Vector2 e12 = (hexagon.Vertex((k + 1) % 6) - hexagon.Vertex(k)) / n;
    for (int j = 1; j <= n; j++)
    {
      for (int i = 0; i < j; i++)
      {
        // Scale the position only for vertices strictly inside the triangle
        if (i > 0)
        {
          vertices.append((circle.Center() + circle.Radius() * Normalized(e01 * j + e12 * i) * (double(j) / double(n))).ToVector());
        }
        // directly compute the coordinates from the frame.
        else
        {
          vertices.append((circle.Center() + e01 * j + e12 * i).ToVector());
        }

      }
    }
  }
  // Indexes
  for (int k = 0; k < 6; k++)
  {
    // Counter clockwise oriented triangle connecting the center : i=1
    indices.append(0);
    indices.append(A(k, 1, 0, n, s));
    indices.append(A(k + 1, 1, 0, n, s));

    // Other counter clockwise oriented triangles
    for (int i = 2; i <= n; i++)
    {
      for (int j = 0; j < i; j++)
      {
        indices.append(A(k, i, j, n, s));
        indices.append(A(k, i, j + 1, n, s));
        indices.append(A(k, i - 1, j, n, s));
      }
    }

    // Other clockwise oriented triangles
    for (int i = 2; i <= n; i++)
    {
      for (int j = 1; j < i; j++)
      {
        indices.append(A(k, i, j, n, s));
        indices.append(A(k, i - 1, j, n, s));
        indices.append(A(k, i - 1, j - 1, n, s));
      }
    }
  }
}

/*!
\brief Create an n-adic subdivision of an ellipse.

Create a circle and apply scaling and translation.

\sa Mesh2::Mesh2(const Circle&, int)

This constructor is the same as the following:
\code
Ellipse2 ellipse(Vector2(2.0,3.0),4.0,2.0);
Mesh2 mesh(Circle(1.0),n);
mesh.Scale(ellipse.A(),ellipse.B());
mesh.Translate(ellipse.Center());
\endcode

\param ellipse The ellipse.
\param n Subdivision level.
*/
Mesh2::Mesh2(const Ellipse2& ellipse, int n) :Mesh2(Circle2(1.0), n)
{
  Scale(Vector2(ellipse.A(), ellipse.B()));
  Vector2 u = ellipse.Axis();
  if (u != Vector::X)
  {
    Rotate(Matrix2(u, u.Orthogonal()));
  }
  Translate(ellipse.Center());
}

/*!
\brief Create triangle mesh from a polygon.

Perform a decomposition of the polygon into triangles using the ear-clipping algorithm.

\sa Polygon::EarClip()
\param polygon The polygon.
*/
Mesh2::Mesh2(const Polygon2& polygon)
{
  indices = polygon.EarClip();
  vertices = polygon.Vertices();
}

Mesh2 Mesh2::Alpha(const QVector<Vector2>& points, const double& radius)
{
  Mesh2 mesh = Delaunay(points);
  int n = 0;

  while (n<mesh.TriangleSize())
  {
    // Compute radius
    Triangle2 t= mesh.GetTriangle(n);
    double r = t.CircumscribedRadius();

    if (r > radius)
    {
      mesh.indices[3*n + 0] = mesh.indices[mesh.indices.size() - 3 + 0];
      mesh.indices[3*n + 1] = mesh.indices[mesh.indices.size() - 3 + 1];
      mesh.indices[3*n + 2] = mesh.indices[mesh.indices.size() - 3 + 2];

      mesh.indices.pop_back();
      mesh.indices.pop_back();
      mesh.indices.pop_back();
    }
    else
    {
      n++;
    }
  }


  for (int i = 0; i < mesh.TriangleSize(); i++)
  {
    Triangle2 t = mesh.GetTriangle(i);
    double r = t.CircumscribedRadius();
    if (r>radius)
      std::cout << r << std::endl;
  }
  return mesh;
}

/*!
\brief Compute the Delaunay triangulation of a set of points.

\param p Point set.

Bowyer-Watson algorithm, O(n log n) complexity, C++ implementation of http://paulbourke.net/papers/triangulate.
*/
Mesh2 Mesh2::Delaunay(const QVector<Vector2>& p)
{
  if (p.size() < 3)
    return Mesh2();

  double eps = 1e-4;
  std::vector<int> ltr; // triangles = 3 indexes

  // Init Data
  // Bounding box
  Box2 box(p);
  /*
  
  double xmin = points[0][0];
  double xmax = xmin;
  double ymin = points[0][1];
  double ymax = ymin;
  for (int i = 0; i < points.size(); i++)
  {
    xmin = Math::Min(xmin, points[i][0]);
    xmax = Math::Max(xmax, points[i][0]);
    ymin = Math::Min(ymin, points[i][1]);
    ymax = Math::Max(ymax, points[i][1]);
  }
  */
  /*
  double xmin = box[0][0];
  double xmax = box[1][0];
  double ymin = box[0][1];
  double ymax = box[1][1];
  

  const double dx = xmax - xmin;
  const double dy = ymax - ymin;


  const double dmax = Math::Max(dx, dy);

  const double midx = (xmin + xmax) / (2.0);
  const double midy = (ymin + ymax) / (2.0);
  */
  const double d = NormInfinity(box.Diagonal());
  const Vector2 c = box.Center();

  // Initialize
  /*
  Vector2 p0 = Vector2(midx - 20 * dmax, midy - dmax);
  Vector2 p1 = Vector2(midx, midy + 20 * dmax);
  Vector2 p2 = Vector2(midx + 20 * dmax, midy - dmax);
  */
  // GUNDY CLEAN THIS OLD CODE DELETE

  Vector2 p0 = c + Vector2(-20.0 * d, -d);
  Vector2 p1 = c+Vector2(0.0, 20 * d);
  Vector2 p2 = c+Vector2(20 * d, -d);

  QVector<Vector2> points = p;

  points.append(p0); ltr.emplace_back(points.size() - 1);
  points.append(p1); ltr.emplace_back(points.size() - 1);
  points.append(p2); ltr.emplace_back(points.size() - 1);

  // Triangulation
  for (int i = 0; i < points.size() - 3; i++)
  {
    std::vector<QPair<int, int>> edges;
    std::vector<int> tmps;

    for (int j = 0; j < ltr.size(); j += 3)
    {
      // Check if the point is inside the triangle circumcircle.
      Circle2 circle = Triangle2(points[ltr[j + 0]], points[ltr[j + 1]], points[ltr[j + 2]]).Circumscribed();
      const double dist = SquaredNorm(circle.Center() - points[i]);
      if ((dist - circle.Radius() * circle.Radius()) <= eps)
      {
        edges.push_back(qMakePair(ltr[j + 0], ltr[j + 1]));
        edges.push_back(qMakePair(ltr[j + 0], ltr[j + 2]));
        edges.push_back(qMakePair(ltr[j + 1], ltr[j + 2]));
      }
      else {
        tmps.push_back(ltr[j + 0]);
        tmps.push_back(ltr[j + 1]);
        tmps.push_back(ltr[j + 2]);
      }
    }

    // Delete duplicate edges. 
    std::vector<bool> remove(edges.size(), false);
    for (auto it1 = edges.begin(); it1 != edges.end(); ++it1)
    {
      for (auto it2 = edges.begin(); it2 != edges.end(); ++it2)
      {
        if (it1 == it2) {
          continue;
        }
        if (*it1 == *it2) {
          remove[std::distance(edges.begin(), it1)] = true;
          remove[std::distance(edges.begin(), it2)] = true;
        }
      }
    }

    edges.erase(
      std::remove_if(edges.begin(), edges.end(),
        [&](auto const& e) { return remove[&e - &edges[0]]; }),
      edges.end());

    // Update triangulation.
    for (int j = 0; j < edges.size(); j++)
    {
      tmps.push_back(edges[j].first);
      tmps.push_back(edges[j].second);
      tmps.push_back(i);
    }
    ltr = tmps;
  }

  // Remove original super triangle.
  points.removeLast(); points.removeLast(); points.removeLast();
  for (int j = 0; j < ltr.size(); j += 3) {

    if ((ltr[j + 0] > points.size() - 1) || (ltr[j + 1] > points.size() - 1) || (ltr[j + 2] > points.size() - 1))
    {
      ltr.erase(ltr.begin() + j, ltr.begin() + j + 3);
      j -= 3;
    }
  }

  // Define Mesh2
  return Mesh2(points, QVector<int>(ltr.begin(), ltr.end()));
}

/*!
\brief Compute the Delaunay triangulation of a set of points in the plane.

This function return the set of index of vertices

This function is only provided as a
reference: it uses a brute force approach with an O(n<SUP>4</SUP>) complexity.

Basically it computes all the possible triangles and checks whether they are
valid by testing if some points lie inside the circumscribed circle.
\param p Set of input points.

\sa Mesh2::Delaunay
*/
Mesh2 Mesh2::DelaunayN4(const QVector<Vector2>& p)
{
  QVector<int> t;

  const int n = p.size();

  for (int i = 0; i < n - 2; i++)
  {
    for (int j = i + 1; j < n - 1; j++)
    {
      for (int k = j + 1; k < n; k++)
      {
        if (WhichSide(p[i], p[j], p[k]) == 0.0) continue; // points are aligned

        const Circle2 c(p.at(i), p.at(j), p.at(k));
        bool reject = false;

        for (int w = 0; w < n; w++)
        {
          if (w == i || w == j || w == k)
            continue;

          if (c.Inside(p[w]))
          {
            reject = true;
            break;
          }
        }
        if (reject == false)
        {
          t.append(i);
          t.append(j);
          t.append(k);
        }
      }
    }
  }
  return Mesh2(p, t);
}

/*!
\brief Compute the boudning box.

Simply eqsuivalent as:
\code
Mesh2 mesh;
Box2 box=Box2(mesh.Vertices());
\endcode
*/
Box2 Mesh2::GetBox() const
{
  return Box2(vertices);
}

/*!
\brief Draw the mesh.
\param graphics Graphics scene.
\param pen Pen used for the edges and vertexes.
\param vertexes Boolean defining whether vertexes should be drawn as discs.
*/
void Mesh2::Draw(QGraphicsScene& graphics, const QPen& pen, bool vertexes) const
{
  Box2 box = GetBox();
  box.SetCubic();
  double r = box.Radius() / 100.0;

  for (int i = 0; i < indices.size() / 3; i++)
  {
    // Triangle
    Triangle2 t = GetTriangle(i);
    t.Draw(graphics, pen);
  }

  if (vertexes == true)
  {
    for (int i = 0; i < vertices.size(); i++)
    {
      QBrush brush(QColor(255, 255, 255));
      Circle2(vertices.at(i), r).Draw(graphics, pen, brush);
    }
  }

}

/*!
\brief Transform a mesh.
\param f Transform.
*/
void Mesh2::Transform(const Frame2& f)
{
  for (int i = 0; i < vertices.size(); i++)
  {
    vertices[i] = f.Transform(vertices[i]);
  }
}

/*!
\brief Translate a mesh.
\param t Translation vector.
*/
void Mesh2::Translate(const Vector2& t)
{
  for (int i = 0; i < vertices.size(); i++)
  {
    vertices[i] += t;
  }
}

/*!
\brief Translate a mesh.
\param s Scaling vector.
*/
void Mesh2::Scale(const Vector2& s)
{
  for (int i = 0; i < vertices.size(); i++)
  {
    // Scale operator
    vertices[i] *= s;
  }
}

/*!
\brief Rotate a mesh.
\param r Rotation matrix.
*/
void Mesh2::Rotate(const Matrix2& r)
{
  for (int i = 0; i < vertices.size(); i++)
  {
    vertices[i] = r * vertices[i];
  }
}

/*!
\brief Compute the valence of the vertices.
\return A vector with the computed valences.
*/
QVector<int> Mesh2::Valences() const
{
  QVector<int> v(vertices.size(), 0);
  for (int i = 0; i < indices.size(); i++)
  {
    v[indices[i]]++;
  }
  return v;
}

/*!
\brief Compute the average aspect ratio of the triangles of the mesh.
*/
double Mesh2::Aspect() const
{
  const int nt = indices.size();

  double a = 0.0;
  for (int i = 0; i < nt; i += 3)
  {
    a += Triangle2(vertices.at(indices.at(i + 0)), vertices.at(indices.at(i + 1)), vertices.at(indices.at(i + 2))).Aspect();
  }
  return a / (nt/3);
}

/*!
\brief Keep triangles that intersect the a box.
\param box The box.
*/
Mesh2 Mesh2::Intersect(const Box2& box) const
{
  QVector<Vector2> v = vertices;
  
  QVector<int> t;
  
  t.reserve(vertices.size());

  const int nt = indices.size();
  for (int i = 0; i < nt; i += 3)
  {
    if (Triangle2(vertices.at(indices.at(i + 0)), vertices.at(indices.at(i + 1)), vertices.at(indices.at(i + 2))).Intersect(box))
    {
      t.append(indices.at(i + 0));
      t.append(indices.at(i + 1));
      t.append(indices.at(i + 2));
    }
  }

  return Mesh2(v, t);
}
