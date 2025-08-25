#pragma once

#include <QtCore/QString>

#include "libs/scalarfield.h"
#include "libs/arrayinteger.h"
#include "libs/mesh.h"
#include "libs/random.h"
#include "libs/curve.h"
#include "libs/curveset.h"
#include "libs/arrayinteger.h"
#include "libs/smooth.h"

class FlowStruct
{
public:
  int mask = 0; //!< Integer mask of downstream flow topology.
public:
  int i[8] = { -1 }; //!< Indexes to neighboring cells.
  QPoint q[8] = { QPoint(0,0) }; //!< Flow points.
  double h[8] = { 0.0 }; //!< Heights.
  double s[8] = { 0.0 }; //!< Slopes.
  double sn[8] = { 0.0 }; //!< Normalized slopes.
  int num = 0; //TODO ADD NUMBER OF DIRECTIONS IN IT
protected:
  int steepest = -1; //!< Steepest slope index, i.e., index in the arrays storing the steepest slope data.
public:
  explicit FlowStruct() {}
  ~FlowStruct() {}

  void Add(const QPoint&, int, const double&, const double&);

  QPoint Steepest() const;
  double SteepestSlope() const;
  int Topology() const;
  int Next(int) const;
  QPoint NextPoint(int) const;

  int SelectSlope(int) const;

  friend class HeightField;
protected:
  static RandomFast fast;
};


/*!
\brief Select direction according to slope distribution.
\param n Number of directions.
*/
inline int FlowStruct::SelectSlope(int n) const
{
  return fast.Pick(s, n);
}

/*!
\brief Get k-th downward point.
\param k Index.
*/
inline QPoint FlowStruct::NextPoint(int k) const
{
  return q[k];
}

/*!
\brief Steepest downward point.
*/
inline QPoint FlowStruct::Steepest() const
{
  return q[steepest];
}

/*!
\brief Steepest downward point.
*/
inline double FlowStruct::SteepestSlope() const
{
  return s[steepest];
}

/*!
\brief %Topology mask.
*/
inline int FlowStruct::Topology() const
{
  return mask;
}

/*!
\brief Get the index of the k-th downward point.
\param k Index.
*/
inline int FlowStruct::Next(int k) const
{
  return i[k];
}


class HeightField :public ScalarField2
{
protected:
public:
  //! Empty
  HeightField() {}

  HeightField(const ScalarField2&);
  explicit HeightField(const Box2&, const QImage&, const double& = 0.0, const double& = 256.0 * 256.0 - 1.0, bool = true);
  explicit HeightField(const Box2&, int, int, const double& = 0.0);
  explicit HeightField(const Box2&, int, int, const QVector<double>&);

  //! Empty
  ~HeightField() {}

  double GetHeight(const Vector2&) const;

  // Vertices
  Vector Vertex(const Vector2&, bool = true) const;
  Vector Vertex(int, int) const;
  QVector<Vector> ArrayVertexes(const QVector<QPoint>&) const;

  // Normals and slope
  Vector Normal(const Vector2&, bool = true) const;
  Vector Normal(int, int) const;
  double Slope(int, int) const;
  double Slope(const QPoint&, const QPoint&) const;
  double AverageSlope(int, int) const;
  double Aspect(int, int) const;

  double K() const;

  void Subdivide(const double&, Random & = Random::R239);

  // Stream area
  ScalarField2 StreamArea() const;
  ScalarField2 StreamAreaWeighted(const double& = 2.0) const;
  ScalarField2 StreamAreaSteepest() const;
  ScalarField2 StreamAreaLimited(const double&) const;
  ScalarField2 StreamLength(bool = false) const;
  ScalarField2 StreamSourceElev(bool = false) const;

  VectorField2 FlowField(ScalarField2&) const;

  // Slope
  ScalarField2 AverageSlope() const;
  ScalarField2 Slope(bool = false) const;
  ScalarField2 Aspect() const;

  // Wetness index, stream power and others derived from previously defined functions
  ScalarField2 WetnessIndex(const double& = 1.0e-3, bool = false) const;
  ScalarField2 StreamPower(bool = false) const;
  ScalarField2 StreamPower(double, double) const;

  // Relief metrics
  ScalarField2 ReliefSteepness(const double&) const;
  ScalarField2 Jut(const double&) const;
  ScalarField2 Rut(const double&) const;
  ScalarField2 Peakedness(const double&) const;
  ScalarField2 TopographicPositionIndex(const double& r) const;
  ScalarField2 LocalVariance(int = 3) const;
  ScalarField2 SurfaceRoughness(int = 3) const;
  ScalarField2 AreaRatio(int = 3) const;
  ScalarField2 TerrainRuggednessIndex(int = 3, bool = false) const;
  ScalarField2 HillslopeAsymmetry(double, double = 0, double = 45) const;

  ScalarField2 Accessibility(const double&, int = 16) const;
  ScalarField2 Light(const Vector&, bool = false) const;


  int FillDepressions();
  int FillDepressions(const double&);
  int FillDepressions(ScalarField2&) const;
  int FillDepressions(const double&, ScalarField2&) const;
  void CompleteBreach();

  // Ray intersection
  bool Intersect(const Ray&, double&, Vector&, const Box&, const double&, const double& = 1.0e8, const double& = 1.0e-4) const;
  bool Intersect(const Ray&, QVector<double>&, QVector<Vector>&, const Box&, const double&) const;
  bool IntersectBox(const Ray&, double&, Vector&, Vector&, const Box&) const;

  double Shade(const Vector&, const QVector<Vector>&, const QVector<double>&, const double&) const;
  ScalarField2 Sun(const double&, int = 5, int = 1, int = 0, int = 12) const;
  ScalarField2 SelfShadow(const Vector&) const;

  Array2I Flow(bool = false) const;

  Array2I Geomorphons(const double& = 0.02) const;
  Array2I GeomorphonsTangent(double, double = 0.02) const;
  Array2I GeomorphonsAggregated(const double& = 0.02) const;

  // Box
  Box GetBox() const;

  double Signed(const Vector&) const;

  void Scale(const Vector&);
  void Scale(const double&);
  void Unity();

  void Draw(QGraphicsScene&, const QPen & = QPen(), const QBrush & = QBrush()) const;

  QVector<QPoint> Up(int, int) const;
  QVector<QPoint> Down(int, int) const;

  // Flow
  int CheckFlowSlope(const QPoint&, FlowStruct&) const;
  int CheckFlowSlopeWeighted(const QPoint&, FlowStruct&, const double& = 2.0) const;
  int CheckFlowDirectionsAngle(const QPoint&, const double&, FlowStruct&) const;

  Mesh CreateMesh(bool = true, bool = false, const double& = 0.0) const;
  Mesh CreateMeshSide(const double& = 0.0) const;

  void CreateCubes(Box&, QVector<FrameScaled>&, const double& = -10.0) const;

  void Carve(const CubicCurve&, const double&, const double&);
  void Carve(const QuadricCurve&, const double&, const double&);
  void Carve(const CubicCurveSet&, const double&, const double&);
  void Carve(const QuadricCurveSet&, const double&, const double&);
  void Carve(const Segment&, const double&, const double&);

  void Carve(const SmoothDisc2&, const double&);

  ScalarField2 GradientDistance(const HeightField&) const;

  ScalarField2 Sea(const double& = 0.0) const;

  QVector<double> Cross(const Vector2&, const Vector2&, int) const;

  void SmoothBreachMultiScale(const QVector<int>&);
  void SmoothBreachGeometric(double, double);
  void SmoothBreachLinear(double);

  void LowErosionMultiScale(const QVector<double>&, const QVector<double>&);
  void LowErosionGeometric(double, double, double, double);

  QString EncodeSize() const;
  static Box Size(const QString&);
public:
  static const double flat; //!< Small negative epsilon value used in breaching and flow algorithms.
};


/*!
\brief Compute the height at a given position.
\param p Point.
*/
inline double HeightField::GetHeight(const Vector2& p) const
{
  return Value(p);
}

class HeightFieldNext8 :public Array2
{
protected:
  QVector<unsigned char> flow; //!< Array of directions
public:
  HeightFieldNext8(const HeightField&, bool);
  //! Empty
  HeightFieldNext8() {}
  //! Empty
  ~HeightFieldNext8() {}

  unsigned char At(int, int) const;
  unsigned char At(const QPoint&) const;

  QPoint FlowTo(const QPoint&) const;
  QPoint FlowTo(const QPoint&, int) const;
  int N(int, int) const;
  int N(const QPoint&) const;
protected:
  unsigned char Flow(const HeightField&, const QPoint&, bool) const;
};

/*!
\brief Return the flowing directions, compacted into one byte.
\param i,j Integer coordinates of the sample.
*/
inline unsigned char HeightFieldNext8::At(int i, int j) const
{
  return flow.at(VertexIndex(i, j));
}

/*!
\brief Return the flowing directions, compacted into one byte.
\param p Point.
*/
inline unsigned char HeightFieldNext8::At(const QPoint& p) const
{
  return flow.at(VertexIndex(p));
}

class HeightFieldErosion :public HeightField
{
protected:
public:
  explicit HeightFieldErosion(const HeightField&);
  ~HeightFieldErosion() {}

  void HillSlope(double, double = 1.0);
  void DebrisSlope(double, double = 1.0);
  void StreamPowerErosion(double = 0.0005, double = 0.01, double = 2., double = 0.8, double = 100.);
  void StreamPowerErosion(double, const ScalarField2&, double = 2., double = 0.8, double = 100.);
  void HillSlope(double, double, double, double);
  void DebrisSlope(double, double, double, double);
};


class TerrainDX :public Array2
{
protected:
  HeightField terrain; //!< Reference to elevation.
  Array2I dirx; //!< Flow Direction.
  ScalarField2 accumulaton; //!< Drainage area. 
public:
  explicit TerrainDX(const HeightField&);
  //! Empty.
  virtual ~TerrainDX() {}

  // Get Local Type
  bool isRiver(int, int) const;
  bool isSource(int, int) const;
  bool isRiverBorder(int, int) const;

  // Getters Local
  QVector<int> getRiverInput(int i, int j) const;
  QVector<int> getFlowInput(int i, int j) const;
  int getRiverOutput(int i, int j) const;

  double getMeanFlow(int, int) const;
};

/*!
\brief
*/
inline bool TerrainDX::isRiver(int i, int j) const
{
  if (!accumulaton.InsideVertexIndex(i, j)) return false;
  return accumulaton.Value(i, j) > 200;
}



