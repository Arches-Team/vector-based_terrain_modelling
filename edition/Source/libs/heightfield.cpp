// Heightfield
#include "libs/heightfield.h"
#include "libs/vectorfield.h"

#include "libs/triangle.h"
#include "libs/gputerrainsimu.h"

RandomFast FlowStruct::fast;

/*!
\mainpage LibHeightField
<center>

<P>Latest stable version : 2024.01.12
</center>
<h3>License</h3>
Copyright &copy; by
<a href="mailto:eric.galin@liris.cnrs.fr">Eric Galin</a>.
No representations are made about the suitability of this software
for any purpose. It is provided as is without express or implied
warranty.

<h3>Fundamental classes in heightfield</h3>
The heightfield library implements several classes for representing heightfields. The list includes two kinds
of classes : \ref HeightField for elevation maps based on regular grids, which inherits from ScalarField2, and \ref AnalyticHeightField which implements several functions for
procedurally defined elevation functions representing terrains.
*/

/*!
\defgroup HeightFieldGroup Heightfield structures and algorithms

\brief Heightfield classes implement several algorithms for processing heightfields.
*/


/*!
\class HeightField heightfield.h
\brief A simple height field.

\ingroup HeightFieldGroup
*/

const double HeightField::flat = 1.0e-8;

/*!
\brief Create a flat heightfield.
\param box Rectangle domain of the terrain.
\param nx, ny Samples.
\param v Constant elevation.
*/
HeightField::HeightField(const Box2& box, int nx, int ny, const double& v) :ScalarField2(box, nx, ny, v)
{
}

/*!
\brief Create a heightfield.
\param box Rectangle domain of the terrain.
\param nx, ny Samples.
\param v Set of elevation values.
*/
HeightField::HeightField(const Box2& box, int nx, int ny, const QVector<double>& v) : ScalarField2(box, nx, ny, v)
{
}

/*!
\brief Create a heightfield from a scalar field.

This constructor provides implicit conversion.

\param s Scalar field.
*/
HeightField::HeightField(const ScalarField2& s) :ScalarField2(s)
{
}

/*!
\brief Create a heightfield from an image.
\param box Rectangle domain of the terrain.
\param image Elevation image.
\param a, b Minimum and maximum elevation range.
\param grayscale Boolean set to false if the image is provided in color.
*/
HeightField::HeightField(const Box2& box, const QImage& image, const double& a, const double& b, bool grayscale) :ScalarField2(box, image, a, b, grayscale)
{
}

/*!
\brief Refine the terrain.

This function subdivides the samples and uses a smooth interpolation.
A random elevation may be added to the samples.
\param e Amplitude of the random elevation added to the samples.
\param r %Random number generator.
*/
void HeightField::Subdivide(const double& e, Random& r)
{
  // Extend grid and smooth values
  ScalarField2::Subdivide();

  // Add some random values
  for (int i = 1; i < nx; i += 2)
  {
    for (int j = 1; j < ny; j += 2)
    {
      field[VertexIndex(i, j)] += r.Uniform(-e, e);
    }
  }
}

/*!
\brief Compute the normal for a given position on the terrain.

Note that this function may be expensive to compute.

\param p Point.
\param triangular Boolean, use triangle normals if set to true, bilinear interpolation of normals at vertices otherwize.
*/
Vector HeightField::Normal(const Vector2& p, bool triangular) const
{
  double u, v;
  int i, j;
  CellInteger(p, i, j, u, v);

  // Test position
  if (!InsideCellIndex(i, j))
    return Vector::Z;

  if (triangular)
  {
    if (u > v)
    {
      return Triangle(Vertex(i, j), Vertex(i + 1, j), Vertex(i + 1, j + 1)).Normal();
    }
    else
    {
      return Triangle(Vertex(i, j), Vertex(i + 1, j + 1), Vertex(i, j + 1)).Normal();
    }
  }
  else
  {
    return Normalized(Vector::Bilinear(Normal(i, j), Normal(i + 1, j), Normal(i + 1, j + 1), Normal(i, j + 1), u, v));
  }
}

/*!
\brief Compute the normal at a given sample.

This function uses the weighted sum (area) of the normals of the
triangles sharing the point on the grid. The returned vector is normalized.

\param i,j Integer coordinates of the sample.
*/
Vector HeightField::Normal(int i, int j) const
{
  Vector n;
  if (i == 0)
  {
    if (j == 0)
    {
      // Corner: 0/1
      n = Triangle(Vertex(i, j), Vertex(i + 1, j), Vertex(i + 1, j + 1)).AreaNormal() + Triangle(Vertex(i, j), Vertex(i + 1, j + 1), Vertex(i, j + 1)).AreaNormal();
    }
    else if (j == ny - 1)
    {
      // Corner: 5
      n = Triangle(Vertex(i, j), Vertex(i, j - 1), Vertex(i + 1, j)).AreaNormal();
    }
    else
    {
      // Edge: 0/1/5
      n = Triangle(Vertex(i, j), Vertex(i + 1, j), Vertex(i + 1, j + 1)).AreaNormal() + Triangle(Vertex(i, j), Vertex(i + 1, j + 1), Vertex(i, j + 1)).AreaNormal()
        + Triangle(Vertex(i, j), Vertex(i, j - 1), Vertex(i + 1, j)).AreaNormal();
    }
  }
  else if (i == nx - 1)
  {
    if (j == 0)
    {
      // Corner: 2
      n = Triangle(Vertex(i, j), Vertex(i, j + 1), Vertex(i - 1, j)).AreaNormal();

    }
    else if (j == ny - 1)
    {
      // Corner: 3/4
      n = Triangle(Vertex(i, j), Vertex(i - 1, j - 1), Vertex(i, j - 1)).AreaNormal() + Triangle(Vertex(i, j), Vertex(i - 1, j), Vertex(i - 1, j - 1)).AreaNormal();
    }
    else
    {
      // Edge: 2/3/4
      n = Triangle(Vertex(i, j), Vertex(i, j + 1), Vertex(i - 1, j)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i - 1, j), Vertex(i - 1, j - 1)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i - 1, j - 1), Vertex(i, j - 1)).AreaNormal();
    }
  }
  else
  {
    if (j == 0)
    {
      // Edge: 0/1/2
      n = Triangle(Vertex(i, j), Vertex(i + 1, j), Vertex(i + 1, j + 1)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i + 1, j + 1), Vertex(i, j + 1)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i, j + 1), Vertex(i - 1, j)).AreaNormal();
    }
    else if (j == ny - 1)
    {
      // Edge: 3/4/5
      n = Triangle(Vertex(i, j), Vertex(i - 1, j), Vertex(i - 1, j - 1)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i - 1, j - 1), Vertex(i, j - 1)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i, j - 1), Vertex(i + 1, j)).AreaNormal();
    }
    else
    {
      // Face: 0/1/2/3/4/5
      n = Triangle(Vertex(i, j), Vertex(i + 1, j), Vertex(i + 1, j + 1)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i + 1, j + 1), Vertex(i, j + 1)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i, j + 1), Vertex(i - 1, j)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i - 1, j), Vertex(i - 1, j - 1)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i - 1, j - 1), Vertex(i, j - 1)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i, j - 1), Vertex(i + 1, j)).AreaNormal();
    }
  }
  return Normalized(n);
}

/*!
\brief Compute the vertex corresponding to a given sample.
\param i, j Integer coordinates of the sample.
*/
Vector HeightField::Vertex(int i, int j) const
{
  return Vector(a[0] + i * celldiagonal[0], a[1] + j * celldiagonal[1], at(i, j));
}

/*!
\brief Compute the coordinates of a set of points on the grid.

\sa Array2::ArrayVertexes
\param p Set of point.
*/
QVector<Vector> HeightField::ArrayVertexes(const QVector<QPoint>& p) const
{
  QVector<Vector> v(p.size());
  for (int i = 0; i < p.size(); i++)
  {
    const QPoint pi = p.at(i);
    v[i] = ArrayVertex(pi).ToVector(at(pi));
  }
  return v;
}

/*!
\brief Compute the vertex position on the terrain.
\param p Point.
\param triangular Boolean, use triangular interpolation if set to true, bilinear interpolation otherwize.
*/
Vector HeightField::Vertex(const Vector2& p, bool triangular) const
{
  double u, v;
  int i, j;
  CellInteger(p, i, j, u, v);

  double z = 0.0;
  // Escape
  if (!InsideCellIndex(i, j))
  {
  }
  else
  {
    if (triangular)
    {
      if (u > v)
      {
        z = (1.0 - u) * at(i, j) + (u - v) * at(i + 1, j) + v * at(i + 1, j + 1);
      }
      else
      {
        z = (1.0 - v) * at(i, j) + u * at(i + 1, j + 1) + (v - u) * at(i, j + 1);
      }
    }
    else
    {
      z = Math::Bilinear(at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1), u, v);
    }
  }
  return Vector(p[0], p[1], z);
}

#include <QtCore/QStack>

#include <algorithm>

class NeighborPoint
{
public:
  QPoint p; //!< Point.
  double height; //!< Height difference.
  double slope;  //!< Slope.
  double weight;  //!< Weighted slope.
public:
  NeighborPoint(const QPoint& p, double h, double s, double w) :p(p), height(h), slope(s), weight(w) {}
};

/*!
\brief Compute the stream length.

The stream length is computed, for all cells, as the maximum distance to a spring, which are the peaks.
*/
ScalarField2 HeightField::StreamLength(bool steepest) const
{
  ScalarField2 streamlength(Box2(a, b), nx, ny, 0.0);

  QVector<ScalarPoint2> queue = GetScalarPoints();
  std::sort(queue.begin(), queue.end());

  for (int i = queue.size() - 1; i >= 0; i--)
  {
    const QPoint p = queue.at(i).Point();

    FlowStruct flow;
    int n = CheckFlowSlope(p, flow);
    if (n > 0)
    {
      if (steepest) {
        QPoint q = flow.Steepest();
        double d = Norm(Vector2(q.x() - p.x(), q.y() - p.y()));
        streamlength(q) = Math::Max(streamlength(q), streamlength(p) + d);
      }
      else {
        const double l = streamlength(p);
        for (int j = 0; j < n; j++)
        {
          streamlength(flow.q[j]) = Math::Max(streamlength(flow.q[j]), l + Array2::length[flow.Next(j)]);
        }
      }
    }
  }

  return streamlength;
}



/*!
\brief Return the maximum source elevation of the flows converging to each point in the height field
\param steepest If true, uses D8 routing algorithm. Otherwise, MFD
*/
ScalarField2 HeightField::StreamSourceElev(bool steepest) const
{
  ScalarField2 sourceelev(*this);

  QVector<ScalarPoint2> QEE = GetScalarPoints();
  std::sort(QEE.begin(), QEE.end());

  for (int i = QEE.size() - 1; i >= 0; i--)
  {
    const QPoint p = QEE.at(i).Point();

    FlowStruct flow;
    int n = CheckFlowSlope(p, flow);
    if (n > 0)
    {
      const double h = sourceelev(p);
      if (steepest) {
        sourceelev(flow.Steepest()) = Math::Max(sourceelev(flow.Steepest()), h);
      }
      else {
        for (int j = 0; j < n; j++)
        {
          sourceelev(flow.q[j]) = Math::Max(sourceelev(flow.q[j]), h);
        }
      }
    }
  }

  return sourceelev;
}



/*!
\brief Compute the stream area field of the terrain.

The stream area is computed as a fraction of cells, independtly of the surface area of a cell.
*/
ScalarField2 HeightField::StreamArea() const
{
  ScalarField2 stream(Box2(a, b), nx, ny, 1.0);

  QVector<ScalarPoint2> queue = GetScalarPoints();
  std::sort(queue.begin(), queue.end());

  for (int i = queue.size() - 1; i >= 0; i--)
  {
    QPoint p = queue.at(i).Point();

    FlowStruct flow;
    int n = CheckFlowSlope(p, flow);
    if (n > 0)
    {
      const double sp = stream(p);
      for (int j = 0; j < n; j++)
      {
        stream(flow.q[j]) += sp * flow.sn[j];
      }
    }
  }

  return stream;
}


/*!
\brief Compute the stream area field of the terrain.

The stream area is computed as a fraction of cells, independtly of the surface area of a cell.
\param power Power.
*/
ScalarField2 HeightField::StreamAreaWeighted(const double& power) const
{
  ScalarField2 stream(Box2(a, b), nx, ny, 1.0);

  QVector<ScalarPoint2> queue = GetScalarPoints();
  std::sort(queue.begin(), queue.end());

  for (int i = queue.size() - 1; i >= 0; i--)
  {
    QPoint p = queue.at(i).Point();

    FlowStruct flow;
    int n = CheckFlowSlopeWeighted(p, flow, power);
    if (n > 0)
    {
      const double sp = stream(p);
      for (int j = 0; j < n; j++)
      {
        stream(flow.q[j]) += sp * flow.sn[j];
      }
    }
  }

  return stream;
}

/*!
\brief Compute
*/
Vector2 AverageFlowDirection(const QPoint& p, QPoint q[8], double s[8], int n)
{
  Vector2 d = Vector2::Null;
  for (int i = 0; i < n; i++)
  {
    QPoint pq = q[i] - p;
    Vector2 pq2 = Vector2(pq.x(), pq.y());
    d += s[i] * Normalized(pq2);
  }
  return Normalized(d);
}


/*!
\brief Compute the flow field of the terrain.

The direction if the average downhill flow.

This function differs slightly from HeightField::StreamArea(). Pits have a null flow vector,
although the stream area is a local maximum.

The stream area is computed as a fraction of cells, independtly of the surface area of a cell.

\param stream Returned stream area;
*/
VectorField2 HeightField::FlowField(ScalarField2& stream) const
{
  stream = ScalarField2(Box2(a, b), nx, ny, 1.0);
  VectorField2 direction(Box2(a, b), nx, ny);

  QVector<ScalarPoint2> queue = GetScalarPoints();
  std::sort(queue.begin(), queue.end());

  for (int i = queue.size() - 1; i >= 0; i--)
  {
    QPoint p = queue.at(i).Point();

    FlowStruct flow;

    int n = CheckFlowSlope(p, flow);
    // Flow gets out of cell
    if (n > 0)
    {
      const double sp = stream(p);
      for (int j = 0; j < n; j++)
      {
        direction(flow.q[j]) = AverageFlowDirection(p, flow.q, flow.s, n);
        stream(flow.q[j]) += sp * flow.sn[j];
      }
    }
  }

  return direction;
}

/*!
\brief Compute the stream area field of the terrain.
\sa HeightField::StreamArea()
*/
ScalarField2 HeightField::StreamAreaSteepest() const
{
  ScalarField2 stream(Box2(a, b), nx, ny, 1.0);

  QVector<ScalarPoint2> queue = GetScalarPoints();
  std::sort(queue.begin(), queue.end());

  for (int i = queue.size() - 1; i >= 0; i--)
  {
    QPoint p = queue.at(i).Point();

    FlowStruct flow;
    int n = CheckFlowSlope(p, flow);
    if (n > 0)
    {
      const double sp = stream(p);
      stream(flow.Steepest()) += sp;
    }
  }

  return stream;
}

/*!
\brief Compute the stream area field of the terrain using a limited flow propagation.

\param r Flow propagation range, this function is the same as HeightField::StreamArea() if range is infinitely large.
\sa HeightField::StreamArea()
*/
ScalarField2 HeightField::StreamAreaLimited(const double& r) const
{
  ScalarField2 water(GetArray(), 1.0);
  ScalarField2 flow(GetArray(), 1.0);

  int n = int(r / Array2::UnitCell().Size()[0]);

  for (int k = 0; k < n; k++)
  {
    ScalarField2 deltawater(GetArray(), 0.0);
    for (int i = 0; i < nx; i++)
    {
      for (int j = 0; j < ny; j++)
      {
        double w = water(i, j);
        if (w <= 0.0) continue;

        FlowStruct fl;
        int num = CheckFlowSlope(QPoint(i, j), fl);
        if (num > 0)
        {
          for (int ii = 0; ii < num; ii++)
          {
            deltawater(fl.q[ii]) += w * fl.sn[ii];
            flow(fl.q[ii]) += w * fl.sn[ii];
          }
          deltawater(i, j) -= w;
        }
      }
    }
    water += deltawater;
  }
  return flow;
}

/*!
\brief Compute the maximum slope field of the terrain.

Simply compute the norm of the gradient.

This function may be used to perform hill-slope erosion simulation.

\param a Boolean, use average slope if set to true, and slope otherwize.
\sa (), Slope(), AverageSlope(), GradientNorm()
*/
ScalarField2 HeightField::Slope(bool a) const
{
  if (a == true)
  {
    return AverageSlope();
  }
  else
  {
    return GradientNorm();
  }
}

/*!
\brief Compute the average slope field of the terrain.
*/
ScalarField2 HeightField::AverageSlope() const
{
  ScalarField2 slope(Box2(a, b), nx, ny);
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      slope(i, j) = AverageSlope(i, j);
    }
  }
  return slope;
}


/*!
\brief Compute the aspect field (face orientation) of the terrain.
*/
ScalarField2 HeightField::Aspect() const
{
  ScalarField2 aspect(Box2(a, b), nx, ny);
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      aspect(i, j) = Aspect(i, j);
    }
  }
  return aspect;
}


/*!
\brief Compute the slope at a given integer point on the terrain.

This is a convenience function corresponding to the following piece of code:
\code
double s=Norm(heightfield.Gradient(i, j));
\endcode
\param i, j Integer coordinates
*/
double HeightField::Slope(int i, int j) const
{
  return Norm(Gradient(i, j));
}

/*!
\brief Compute the slope between two points.

\param pa, pb Points.
*/
double HeightField::Slope(const QPoint& pa, const QPoint& pb) const
{
  double h = fabs(at(pa) - at(pb));
  QPoint pab = pb - pa;
  double l = Norm(Vector2(pab.x(), pab.y()).Scaled(celldiagonal));
  return h / l;
}

/*!
\brief Compute the average slope at a given integer point on the terrain.

Simply average the slope in 8 directions.
\param i, j Integer coordinates
*/
double HeightField::AverageSlope(int i, int j) const
{
  double s;

  double e = Norm(celldiagonal);

  if (i == 0)
  {
    if (j == 0)
    {
      // Corner
      s = fabs(at(i, j) - at(i + 1, j)) * inversecelldiagonal[0] + fabs(at(i, j) - at(i + 1, j + 1)) / e + fabs(at(i, j) - at(i, j + 1)) * inversecelldiagonal[1];
      s /= 3.0;
    }
    else if (j == ny - 1)
    {
      // Corner
      s = fabs(at(i, j) - at(i + 1, j)) * inversecelldiagonal[0] + fabs(at(i, j) - at(i + 1, j - 1)) / e + fabs(at(i, j) - at(i, j - 1)) * inversecelldiagonal[1];
      s /= 3.0;
    }
    else
    {
      // Edge
      s = fabs(at(i, j) - at(i, j - 1)) * inversecelldiagonal[1] + fabs(at(i, j) - at(i + 1, j - 1)) / e + fabs(at(i, j) - at(i + 1, j)) * inversecelldiagonal[0] + fabs(at(i, j) - at(i + 1, j + 1)) / e + fabs(at(i, j) - at(i, j + 1)) * inversecelldiagonal[1];
      s /= 5.0;
    }
  }
  else if (i == nx - 1)
  {
    if (j == 0)
    {
      // Corner
      s = fabs(at(i, j) - at(i - 1, j)) * inversecelldiagonal[0] + fabs(at(i, j) - at(i - 1, j + 1)) / e + fabs(at(i, j) - at(i, j + 1)) * inversecelldiagonal[1];
      s /= 3.0;
    }
    else if (j == ny - 1)
    {
      // Corner
      s = fabs(at(i, j) - at(i - 1, j)) * inversecelldiagonal[0] + fabs(at(i, j) - at(i - 1, j - 1)) / e + fabs(at(i, j) - at(i, j - 1)) * inversecelldiagonal[1];
      s /= 3.0;
    }
    else
    {
      // Edge
      s = fabs(at(i, j) - at(i, j - 1)) * inversecelldiagonal[1] + fabs(at(i, j) - at(i - 1, j - 1)) / e + fabs(at(i, j) - at(i - 1, j)) * inversecelldiagonal[0] + fabs(at(i, j) - at(i - 1, j + 1)) / e + fabs(at(i, j) - at(i, j + 1)) * inversecelldiagonal[1];
      s /= 5.0;
    }
  }
  else
  {
    if (j == 0)
    {
      // Edge
      s = fabs(at(i, j) - at(i - 1, j)) * inversecelldiagonal[0] + fabs(at(i, j) - at(i - 1, j + 1)) / e + fabs(at(i, j) - at(i, j + 1)) * inversecelldiagonal[1] + fabs(at(i, j) - at(i + 1, j + 1)) / e + fabs(at(i, j) - at(i + 1, j)) * inversecelldiagonal[0];
      s /= 5.0;
    }
    else if (j == ny - 1)
    {
      // Edge
      s = fabs(at(i, j) - at(i - 1, j)) * inversecelldiagonal[0] + fabs(at(i, j) - at(i - 1, j - 1)) / e + fabs(at(i, j) - at(i, j - 1)) * inversecelldiagonal[1] + fabs(at(i, j) - at(i + 1, j - 1)) / e + fabs(at(i, j) - at(i + 1, j)) * inversecelldiagonal[0];
      s /= 5.0;
    }
    else
    {
      // Vertex
      s = fabs(at(i, j) - at(i + 1, j)) * inversecelldiagonal[0] + fabs(at(i, j) - at(i + 1, j + 1)) / e + fabs(at(i, j) - at(i, j + 1)) * inversecelldiagonal[1] + fabs(at(i, j) - at(i - 1, j + 1)) / e + fabs(at(i, j) - at(i - 1, j)) * inversecelldiagonal[0] + fabs(at(i, j) - at(i - 1, j - 1)) / e + fabs(at(i, j) - at(i, j - 1)) * inversecelldiagonal[1] + fabs(at(i, j) - at(i + 1, j - 1)) / e;
      s /= 8.0;
    }
  }
  return s;
}



/*!
\brief Compute the aspect (face orientation) at a given integer point on the terrain.

\param i, j Integer coordinates
*/
double HeightField::Aspect(int i, int j) const
{
  Vector2 g = Gradient(i, j);
  return Math::Pi + Math::ArcTan(g[0], g[1]);
}


/*!
\brief Compute the wetness index field of the terrain.

The algorithm works as follows: first compute the stream area
for every terrain sample, and then compute the wetness index.

Note that the wetness index is defined as ln ( A / s ) where
s is the slope, however the definition does not work well
for small slopes, therefore the implementation defines the wetness
index as ln ( A / ( e + s ) ), where e is 0.001.

\param c Coefficent for slope factor, if set to 0 the function will compute the logarithm of the drainage area.
\param a Boolean, use average slope if set to true, and slope otherwize.
\sa Slope(), AverageSlope()
*/
ScalarField2 HeightField::WetnessIndex(const double& c, bool a) const
{
  // Compute wetness index field using stream area
  ScalarField2 wetness = StreamArea();
  ScalarField2 slope = Slope(a);

  // Modify with slope
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      wetness(i, j) = log(wetness(i, j) / (c + slope(i, j)));
    }
  }
  return wetness;
}

/*!
\brief Compute the stream power field of the terrain.

The algorithm works as follows: first compute the stream area
for every terrain sample, and then compute the stream power index.

\param a Boolean, use average slope if set to true, and slope otherwize.
*/
ScalarField2 HeightField::StreamPower(bool a) const
{
  // Initialize wetness index field using stream area
  ScalarField2 power = StreamArea();
  ScalarField2 slope = Slope(a);

  // Modify with slope
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      power(i, j) = sqrt(power(i, j)) * slope(i, j);
    }
  }
  return power;
}

/*!
\brief Compute the stream power field of the terrain.

The algorithm works as follows: first compute the stream area
for every terrain sample, and then compute the stream power index.

\param m,n Exponent for the stream area and slope, n should be twice as m.
*/
ScalarField2 HeightField::StreamPower(double m, double n) const
{
  // Initialize stream power field using stream area
  ScalarField2 power = StreamArea();
  ScalarField2 slope = Slope();

  // Modify with slope
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      power(i, j) = pow(power(i, j), m) * pow(slope(i, j), n);
    }
  }
  return power;
}

/*!
\brief Get the bounding box of the heightfield.

Note that although this function has the same name as Array2::GetBox(),
it computes the minimum and maximum elevation of the terrain (computationally intensive).

\sa Array2::GetBox()
*/
Box HeightField::GetBox() const
{
  double za, zb;
  GetRange(za, zb);
  return Array2::GetBox().ToBox(za, zb);
}

/*!
\brief Compute the signed distance do the heightfield.

This function uses ScalarField2::BiCubicValue to compute the elevation.

The volume delimited by the heightfield is embedded an infinite vertical box.

\param p Point.
*/
double HeightField::Signed(const Vector& p) const
{
  double d = p[2] - BiCubicValue(Vector2(p));
  double b = Box2::Signed(p);
  return Math::Max(d, b);
}

/*!
\brief Compute the intersection between a ray and the surface border of the heightfield.

\param ray The ray.
\param t Distance along the ray.
\param box %Box where intersection will be computed.
\param p Returned intersection point.
\param n Returned normal.
*/
bool HeightField::IntersectBox(const Ray& ray, double& t, Vector& p, Vector& n, const Box& box) const
{
  double ta, tb;
  Vector na, nb;

  // Check the intersection with the bounding box
  if (!box.Intersect(ray, ta, tb, na, nb)) return false;

  t = ta + 0.0001;
  n = na;
  if (ta < 0.0)
  {
    t = tb - 0.0001;
    n = nb;
    return false;
  }

  p = ray(t);
  double z = GetHeight(p);
  if (z > p[2])
  {
    return true;
  }

  return false;
}


/*!
\brief Compute the intersection between a ray and the surface of the heightfield.

The algorithm uses a ray marching approach, therefore this function may require
many iterations and the resulting intersection may not be accurate.

\param ray The ray.
\param t Returned distance along the ray.
\param box %Box where intersection will be computed.
\param q Returned intersection point.
\param k Lipschitz constant.
\param length Maximum distance along the ray.
\param epsilon Minimum stepping distance.
*/
bool HeightField::Intersect(const Ray& ray, double& t, Vector& q, const Box& box, const double& k, const double& length, const double& epsilon) const
{
  double ta, tb;

  // Check the intersection with the bounding box
  if (!box.Intersect(ray, ta, tb)) return false;

  if (ta<-1.0e8 || tb>+1.0e8) return false;

  tb = Math::Min(tb, length);

  t = Math::Max(ta + epsilon, 0.0);

  // Ray marching
  while (t < tb)
  {
    // Point along the ray
    Vector p = ray(t);

    // Heightfield elevation
    double z = GetHeight(p);
    double h = p[2] - z;
    if (h < epsilon)
    {
      q = Vector(p[0], p[1], z);
      return true;
    }
    else
    {
      t += Math::Max(h / k, epsilon);
    }
  }
  return false;
}

/*!
\brief Compute the intersection between a ray and the surface of the heightfield.

The algorithm uses a sphere tracing approach, therefore this function may require
many iterations and the resulting intersection may not be accurate.

\param ray The ray.
\param tv Returned distance along the ray.
\param box %Box where intersection will be computed.
\param q Returned intersection point.
\param k Lipschitz constant.
*/
bool HeightField::Intersect(const Ray& ray, QVector<double>& tv, QVector<Vector>& q, const Box& box, const double& k) const
{
  double ta, tb;

  // Check the intersection with the bounding box
  if (!box.Intersect(ray, ta, tb)) return false;

  if (ta<-1000000000 || tb>+1000000000) return false;

  double t = Math::Max(ta + 0.0001, 0.0);

  // Ray marching
  bool findIntersection = false;
  while (t < tb)
  {
    // Point along the ray
    Vector p = ray(t);

    // Heightfield elevation
    double z = GetHeight(p);
    double h = abs(p[2] - z);
    if (h < 0.0001)
    {
      tv.append(t);
      q.append(Vector(p[0], p[1], z));
      findIntersection = true;
      t += 1.;
    }

    t += Math::Max(h / k, 0.0001);
  }
  return findIntersection;
}

/*!
\brief Shade a vector of the terrain.
\param p Terrain coordinates.
\param lights Set of lights.
\param intensity Set of light intensities.
\param sum Sum.
*/
double HeightField::Shade(const Vector& p, const QVector<Vector>& lights, const QVector<double>& intensity, const double& sum) const
{
  Vector q = Vector2(p).ToVector(Value(p) + 0.001);
  Vector n = Normal(p);

  double s = 0.0;
  Box box = GetBox();

  // Maximum distance along the ray
  double length = Norm(box.Diagonal());

  double k = K();

  for (int i = 0; i < lights.size(); i++)
  {
    Ray ray(q, Normalized(lights.at(i) - q));

    // There will be an intersection 
    if (ray.Direction() * n < 0.0) continue;

    Vector qq;
    double tq;
    if (!Intersect(ray, tq, qq, box, k, length))
    {
      double t = n * ray.Direction();
      t = 0.5 * (1 + t);

      s += t * intensity.at(i);
    }
  }
  s = s / sum;
  return s;
}

/*!
\brief Draw a heightfield.
\param scene Graphics scene.
\param pen The pen.
\param brush The brush.
*/
#pragma warning(push)
#pragma warning(disable: 4100)  
void HeightField::Draw(QGraphicsScene& scene, const QPen& pen, const QBrush& brush) const
{
  double a, b;
  GetRange(a, b);

  QPen invisible = pen;
  invisible.setStyle(Qt::NoPen);

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      double z = Linear::Step(at(VertexIndex(i, j)), a, b);

      QBrush fill(Color(z, z, z).GetQt());
      Box2 box = Box2(ArrayVertex(i, j) - celldiagonal, ArrayVertex(i, j) + celldiagonal);
      box.Draw(scene, invisible, fill);
    }
  }
}
#pragma warning(pop)

/*!
\brief Uniform scaling.

\param s Scaling factor.
*/
void HeightField::Scale(const double& s)
{
  Scale(Vector(s));
}


/*!
\brief Scale the height field.

The domain is scaled using the x and y components of the scaling vector,
whereas the heights of the heightfield are scaled using the z component.
\param s Scaling factor.
*/
void HeightField::Scale(const Vector& s)
{
  // Box
  ScalarField2::Scale(Vector2(s));

  // Heights
  for (int i = 0; i < field.size(); i++)
  {
    field[i] *= s[2];
  }
}

/*!
\brief Scale the heightfield so that the compact support should fit within unit box.

Centers the terrain and scales it.
*/
void HeightField::Unity()
{
  Translate(-Box2::Center());

  // Get largest dimension
  double s = NormInfinity(0.5 * Box2::Diagonal());

  Scale(Vector(1.0 / s));
}

/*!
\brief Compute the cells that are uphill of input point.

\param x, y Integer coordinates of the point.
*/
QVector<QPoint> HeightField::Up(int x, int y) const
{
  QPoint p = QPoint(x, y);
  // Elevation
  double z = at(p);

  QVector<QPoint> up;

  for (int i = 0; i < 8; i++)
  {
    QPoint q = p + next[i];
    if (InsideVertexIndex(q))
    {
      if (at(q) > z + 0.0001)
      {
        up.append(q);
      }
    }
  }
  return up;
}

/*!
\brief Compute the cells that are downhill of input point.

\param x, y Integer coordinates of the point.
*/
QVector<QPoint> HeightField::Down(int x, int y) const
{
  QPoint p = QPoint(x, y);
  // Elevation
  double z = at(p);

  QVector<QPoint> down;

  for (int i = 0; i < 8; i++)
  {
    QPoint q = p + next[i];
    if (InsideVertexIndex(q))
    {
      if (at(q) < z)
      {
        down.append(q);
      }
    }
  }
  return down;
}

/*!
\brief Compute the gradient distance between two heightfields.
\param h Second heightfield.
*/
ScalarField2 HeightField::GradientDistance(const HeightField& h) const
{
  ScalarField2 error(Box2(a, b), nx, ny);

  for (int j = 0; j < ny; j++)
  {
    for (int i = 0; i < nx; i++)
    {
      error(i, j) = 1.0 - (Normal(i, j) * h.Normal(i, j));
    }
  }
  return error;
}

/*!
\brief Compute the water elevation in every cell according to a uniform sea level.

\return The amount of water in every cell.
*/
ScalarField2 HeightField::Sea(const double& height) const
{
  ScalarField2 water(GetArray());
  for (int i = 0; i < nx * ny; i++)
  {
    double z = at(i);
    // Terrain below water
    if (z < height)
    {
      // Set water height
      water[i] = height - z;
    }
  }
  return water;
}

/*!
\brief Compute the elevation along a segment.

\param a,b %Segment.
\param n Number of samples.
*/
QVector<double> HeightField::Cross(const Vector2& a, const Vector2& b, int n) const
{
  QVector<double> cut;

  cut.reserve(n);

  cut.append(GetHeight(a));

  for (int k = 1; k < n - 1; k++)
  {
    Vector2 p = Vector2::Lerp(a, b, Math::Unit(k, n));
    cut.append(GetHeight(p));
  }

  cut.append(GetHeight(b));

  return cut;
}

/*!
\brief Compute the flow topology.
\param steep Flag, set to true if we rely on the steepest slope, false for multiple convergent and divergent flows.
*/
Array2I HeightField::Flow(bool steep) const
{
  Array2I topology(Box2(a, b), nx, ny);

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      FlowStruct flow;
      QPoint p(i, j);
      int n = CheckFlowSlope(p, flow);

      if (steep)
      {
        if (n > 0)
        {
          topology(i, j) = 1 << flow.i[flow.steepest];
        }
      }
      else
      {
        topology(i, j) = flow.Topology();
      }
    }
  }

  return topology;
}

/*!
\brief Compute the geomorphon index map.

Geomorphon classify cells into ten categories: 0: flat, 1: peak, 2: ridge, 3: shoulder, 4:hollow=convex, 5: slope, 6: spur=concave, 7: footslope, 8: valley, and 9: pit.

\param t Slope threshold.
*/
Array2I HeightField::Geomorphons(const double& t) const
{
  Array2I geo(GetArray());

  int shift[8] = { 1, 3, 9, 27, 81, 243, 729, 2187 };
  // Geomorphon index lookup table
  static const int index[9][9] =
  {
    { 0, 0, 0, 7, 7, 8, 8, 8, 9},
    { 0, 0, 7, 7, 7, 8, 8, 8,-1},
    { 0, 3, 5, 5, 4, 4, 8,-1,-1},
    { 3, 3, 5, 5, 5, 4,-1,-1,-1},
    { 3, 3, 6, 5, 5,-1,-1,-1,-1},
    { 2, 2, 6, 6,-1,-1,-1,-1,-1},
    { 2, 2, 2,-1,-1,-1,-1,-1,-1},
    { 2, 2,-1,-1,-1,-1,-1,-1,-1},
    { 1,-1,-1,-1,-1,-1,-1,-1,-1},
  };

  const double height = t * celldiagonal[0];

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      QPoint p(i, j);
      int type = 0;

      double zp = at(p);

      // Number of pluses and minuses
      int np = 0;
      int nm = 0;
      for (int k = 0; k < 8; k++)
      {
        QPoint q = Next(p, k);
        if (!InsideVertexIndex(q)) continue;
        double zq = at(q);
        double z = zp - zq;
        int s = 1 + Math::IntegerSign(z, height);

        // Unique type in [0,6561[
        type += shift[k] * s;
        if (s == 0) { nm++; }
        if (s == 2) { np++; }
      }
      // Geomorphon type
      geo(p) = index[nm][np];
    }
  }
  return geo;
}

/*!
\brief Compute the geomorphon index map.

Geomorphon classify cells into ten categories: flat, peak, ridge, shoulder, hollow, slope, spur, footslope, valley, and pit.

\param maxDist maximum ray traversal distance.
\param flatTangent threshold for considering a ray angle as flat.
*/
Array2I HeightField::GeomorphonsTangent(double maxDist, double flatTangent) const
{
  Array2I geo(GetArray());

  int shift[8] = { 1, 3, 9, 27, 81, 243, 729, 2187 };
  // Geomorphon index lookup table
  // FL-0 -  flat; PK 1 - peak; RI 2- ridge; SH 3 - shoulder; HL -4- hollow=convex; SL -5 slope; SP - 6 spur=concave; FS - 7 - footslope; VL - 8 -valley; PT 9- pit.
  static const int index[9][9] =
  {
    { 0, 0, 0, 7, 7, 8, 8, 8, 9},
    { 0, 0, 7, 7, 7, 8, 8, 8,-1},
    { 0, 3, 5, 5, 4, 4, 8,-1,-1},
    { 3, 3, 5, 5, 5, 4,-1,-1,-1},
    { 3, 3, 6, 5, 5,-1,-1,-1,-1},
    { 2, 2, 6, 6,-1,-1,-1,-1,-1},
    { 2, 2, 2,-1,-1,-1,-1,-1,-1},
    { 2, 2,-1,-1,-1,-1,-1,-1,-1},
    { 1,-1,-1,-1,-1,-1,-1,-1,-1},
  };

  const double dt = celldiagonal[0];

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      double zp = at(i, j);

      // Number of pluses and minuses
      int np = 0;
      int nm = 0;
      int type = 0;
      for (int k = 0; k < 8; k++)
      {
        double t = dt;
        Vector2 dir = Vector2(next[k].x(), next[k].y());
        Vector2 p0 = ArrayVertex(i, j);
        double slope = -1e12;

        // traverse terrain
        Vector2 p = p0 + t * dir;
        while ((maxDist < 0 || t < maxDist) && Inside(p)) {
          double dh = Value(p) - zp;
          double s = dh / t;
          slope = std::max(slope, s);
          t += dt;
          p = p0 + t * dir;
        }

        int s = 1 + Math::IntegerSign(slope, flatTangent);

        // Unique type in [0,6561[
        type += shift[k] * s;
        if (s == 0) { nm++; }
        if (s == 2) { np++; }
      }

      // Geomorphon type
      geo(i, j) = index[nm][np];
    }
  }
  return geo;
}

/*!
\brief Compute the aggregated geomorphon index map.

Aggregation classifies into only five categories: flat (flat), ridges (ridge, peak, spur, shoulder), slope (slope), valley (valley, pit, hollow) and foot slope (foot slope).

\param t Slope threshold.
*/
Array2I HeightField::GeomorphonsAggregated(const double& t) const
{
  Array2I geo(GetArray());

  // Aggregated geomorphon index lookup table
  static const int index[9][9] =
  {
    { 0, 0, 0, 4, 4, 3, 3, 3, 3},
    { 0, 0, 4, 4, 4, 3, 3, 3,-1},
    { 0, 1, 2, 2, 3, 3, 3,-1,-1},
    { 1, 1, 2, 2, 2, 3,-1,-1,-1},
    { 1, 1, 1, 2, 2,-1,-1,-1,-1},
    { 1, 1, 1, 1,-1,-1,-1,-1,-1},
    { 1, 1, 1,-1,-1,-1,-1,-1,-1},
    { 1, 1,-1,-1,-1,-1,-1,-1,-1},
    { 1,-1,-1,-1,-1,-1,-1,-1,-1},
  };

  const double height = t * celldiagonal[0];

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      QPoint p(i, j);

      double zp = at(p);

      // Number of pluses and minuses
      int np = 0;
      int nm = 0;
      for (int k = 0; k < 8; k++)
      {
        QPoint q = Next(p, k);
        if (!InsideVertexIndex(q)) continue;
        double zq = at(q);
        double z = zp - zq;
        int s = 1 + Math::IntegerSign(z, height);
        // Unique type in [0,6561[
        if (s == 0) { nm++; }
        if (s == 2) { np++; }
      }
      // Geomorphon type
      geo(p) = index[nm][np];
    }
  }
  return geo;
}

/*!
\brief Carve the terrain using a cubic curve, whose elevation define the prescribed elevation along the curve.
\param c %Curve.
\param r %Radius.
\param b %Blending radius.
*/
void HeightField::Carve(const CubicCurve& c, const double& r, const double& b)
{
  // Flat curve
  CubicCurve2 xy(c);

  // Elevation
  Cubic z = c[2];

  // Bounding rectangle
  Box2 box = xy.GetBox().Extended(r + b);

  // Area
  QRect area = VertexIntegerArea(box);

  double rb = (r + b);
  // Compute thickness
  for (int y = area.y(); y <= area.y() + area.height(); y++)
  {
    for (int x = area.x(); x <= area.x() + area.width(); x++)
    {
      // Squared distance to curve
      double u;
      double d = xy.R(ArrayVertex(x, y), u);

      if (d < rb * rb)
      {
        double zu = z(u);
        double a = Cubic::SmoothCompact(sqrt(d), r, rb);

        field[VertexIndex(x, y)] = (1.0 - a) * zu + a * field[VertexIndex(x, y)];
      }
    }
  }
}

/*!
\brief Carve the terrain using a segment, whose elevation define the prescribed elevation along the curve.
\param c %Segment.
\param r %Radius.
\param b %Blending radius.
*/
void HeightField::Carve(const Segment& c, const double& r, const double& b)
{
  // Flat curve
  Segment2 xy(c);

  // Elevation
  Linear z(c.Vertex(1)[2] - c.Vertex(0)[2], c.Vertex(0)[2]);

  // Bounding rectangle
  Box2 box = xy.GetBox().Extended(r + b);

  // Area
  QRect area = VertexIntegerArea(box);

  double rb = (r + b);
  // Compute thickness
  for (int y = area.y(); y <= area.y() + area.height(); y++)
  {
    for (int x = area.x(); x <= area.x() + area.width(); x++)
    {
      // Squared distance to curve
      double u;
      double d = xy.R(ArrayVertex(x, y), u);

      if (d < rb * rb)
      {
        double zu = z(u);
        double a = Cubic::SmoothCompact(sqrt(d), r, rb);

        field[VertexIndex(x, y)] = (1.0 - a) * zu + a * field[VertexIndex(x, y)];
      }
    }
  }
}

/*!
\brief Carve the terrain using a quadric curve, whose elevation define the prescribed elevation along the curve.
\param c %Curve.
\param r %Radius.
\param b %Blending radius.
*/
void HeightField::Carve(const QuadricCurve& c, const double& r, const double& b)
{
  // Flat curve
  QuadricCurve2 xy(c);

  // Elevation
  Quadric z = c[2];

  // Bounding rectangle
  Box2 box = xy.GetBox().Extended(r + b);

  // Area
  QRect area = VertexIntegerArea(box);

  double rb = (r + b);
  // Compute thickness
  for (int y = area.y(); y <= area.y() + area.height(); y++)
  {
    for (int x = area.x(); x <= area.x() + area.width(); x++)
    {
      // Squared distance to curve
      double u;
      double d = xy.R(ArrayVertex(x, y), u);

      if (d < rb * rb)
      {
        double zu = z(u);
        double a = Cubic::SmoothCompact(sqrt(d), r, rb);

        field[VertexIndex(x, y)] = (1.0 - a) * zu + a * field[VertexIndex(x, y)];
      }
    }
  }
}

/*!
\brief Carve the terrain using a piecewise cubic curve, whose elevation define the prescribed elevation along the curve.
\param c %Curve.
\param r %Radius.
\param b %Blending radius.
*/
void HeightField::Carve(const CubicCurveSet& c, const double& r, const double& b)
{
  for (int i = 0; i < c.Size(); i++)
  {
    Carve(c(i), r, b);
  }
}

/*!
\brief Carve the terrain using a piecewise quadric curve, whose elevation define the prescribed elevation along the curve.
\param c %Curve.
\param r %Radius.
\param b %Blending radius.
*/
void HeightField::Carve(const QuadricCurveSet& c, const double& r, const double& b)
{
  for (int i = 0; i < c.Size(); i++)
  {
    Carve(c(i), r, b);
  }
}

/*!
\brief Smooth breaching with a given input series.
\param ni Set of blurring iterations.
\author Hugo Schott
*/
void HeightField::SmoothBreachMultiScale(const QVector<int>& ni)
{
  HeightField temp_hf;
  ScalarField2 delta;

  int n = ni.size();

  for (int i = 0; i < n; i++)
  {
    temp_hf = *this;
    temp_hf.CompleteBreach();
    delta = *this;
    delta -= temp_hf;
    //delta.Blur(ni.at(i));
    double d_min, d_max_1, d_max_2;
    delta.GetRange(d_min, d_max_1);
    if (ni.at(i)) delta = delta.GaussianBlur(ni.at(i));
    delta.GetRange(d_min, d_max_2);
    delta.Normalize();
    double coeff = 1. / std::pow(double(n - i), 0.2);
    std::cout << ni.at(i) << " " << d_max_1 << " " << d_max_2 << " " << coeff << std::endl;
    (*this) -= coeff * d_max_1 * delta;
  }
}

/*!
\brief Smooth breaching with a geometric series.
\param r Radius, breaching will be performed at the corresponding resolution.
\param a Attenuation.
\author Hugo Schott
*/
void HeightField::SmoothBreachGeometric(double r, double a)
{
  QVector <int> smooth_iter = { int(r * inversecelldiagonal[0]) };

  double last_iter = smooth_iter[0];
  int k = 0;
  while (last_iter >= 0 && k < 200)
  {
    last_iter *= a;
    if (int(last_iter) != smooth_iter.back()) smooth_iter.push_back(int(last_iter));
    k++;
  }

  SmoothBreachMultiScale(smooth_iter);
}

/*!
\brief Smooth breaching process.
\param r Radius, breaching will be performed at the corresponding resolution with linear attenuation.
\author Hugo Schott
*/
void HeightField::SmoothBreachLinear(double r)
{
  QVector <int> smooth_iter;
  const int n = r * inversecelldiagonal[0];

  for (int i = n; i > -1; i--) {
    smooth_iter.push_back(i);
  }

  SmoothBreachMultiScale(smooth_iter);
}


/*!
\brief Preparametered multi-scale medium scale erosion. Similar capacity limited erosion/deposition to Mei [2007], but with drainage-like flow routing instead of shallow water.
\param radii List of erosion radius for each scale.
\param iters List of number of iterations to apply at each scale.
*/
void HeightField::LowErosionMultiScale(const QVector<double>& radii, const  QVector<double>& iters)
{
  GPUHydraulicErosionGrid_Test gpu_erosion;
  gpu_erosion.Init(*this);

  int nb_scale = std::min(radii.size(), iters.size());
  for (int i = 0; i < nb_scale; i++)
  {
    gpu_erosion.SetUniforms(1.0f   /*erosion speed*/, 0.6f  /*deposition speed*/, 2.6f  /*rain*/, 0.81f /*evaporation*/);
    gpu_erosion.SetGlobalSpeed(0.1f);
    gpu_erosion.SetMaxWater(10.0f);
    gpu_erosion.SetFlowRate(0.32f);
    gpu_erosion.SetSmoothSteps(int(radii[i] * inversecelldiagonal[0]));

    gpu_erosion.Step(iters[i]);
  }
  gpu_erosion.GetData(*this);
}

/*!
\brief Preparametered multi-scale medium scale erosion. Similar capacity limited erosion/deposition to Mei [2007], but with drainage-like flow routing instead of shallow water.
\param max_radius Largest erosion radius.
\param radius_attenuation Attenuation coeff to compute next radii.
\param max_iter Number of iteration for the largest radius.
\param iter_attenuation Attenuation coeff to compute next number of iteration.
Stops when radius is one pixel or below.
*/
void HeightField::LowErosionGeometric(double max_radius, double radius_attenuation, double max_iter, double iter_attenuation)
{
  QVector<double> radii = { max_radius };
  QVector<double> iters = { max_iter };

  double cellSize = celldiagonal[0];

  double last_radii = max_radius;
  double last_iter = max_iter;
  int k = 0;
  while (last_radii >= cellSize && k < 200)
  {

    last_radii *= radius_attenuation;
    radii.push_back(last_radii);

    last_iter *= iter_attenuation;
    iters.push_back(last_iter);
    k++;
  }

  LowErosionMultiScale(radii, iters);
}

/*!
\brief Code the size and range parameters of the heightfield into a string.

output / encode a string : x=40800 y=40800 z=+127 +3687 in meters or c=25645-z= in centimeters for cell size and grey level

*/
QString HeightField::EncodeSize() const
{
  double a, b;
  GetRange(a, b);
  Vector2 d = Box2::Diagonal();
  QString s = QString("x=%1 y=%2 z=%3 %4").arg(d[0], 0, 'g', 1).arg(d[1], 0, 'g', 1).arg(a, 0, 'g', 1).arg(b, 0, 'g', 1);
  return s;
}

/*!
\brief Compute the box of the heightfield from a string.
\param s String.
*/
Box HeightField::Size(const QString& s)
{
  // Size
  double x = 0;
  double y = 0;
  // Range
  double a = 0.0;
  double b = 0.0;
  return Box2(x, y).ToBox(a, b);
}
