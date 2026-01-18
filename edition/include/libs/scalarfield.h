// Fundamentals
#pragma once

#include <functional>

#include <QtGui/QImage>

#include "libs/array.h"
#include "libs/tetrahedra.h"
#include "libs/palette.h"

class Color;

class Polygons2;
class SegmentSet2;
class Segment2;

class Mesh;
class MeshColor;
class Voxel;

class Sphere;
class Segment;

class ScalarField2;
class ScalarField;

class Rectangles;
class Rectangles2;

class Quadrangle;

class AnalyticScalarField
{
protected:
  bool sign = true; //!< Sign convention, used for normal computation.
  static int ncubes; //!< Number of cubes processed by polygonization algorithms.
public:
  AnalyticScalarField(bool = true);
  virtual double Value(const Vector&) const;
  virtual Vector Gradient(const Vector&) const;
  Matrix Hessian(const Vector&) const;

  // Normal
  virtual Vector Normal(const Vector&) const;

  // Inside 
  bool Inside(const double&) const;

  // Color
  virtual Color GetMaterial(const Vector&, const Vector & = Vector::Null) const;

  virtual Box GetBox() const;


  // Curvature
  void Curvature(const Vector&, double&, double&) const;

  // Projection on the surface
  Vector Cast(const Vector&, const double&, const double& = 1e-6, int n = 100) const;

  // Dichotomy
  Vector Dichotomy(Vector, Vector, double, double, double, const double& = 1.0e-4) const;

  // Sampling
  virtual bool GetSample(Vector&, const Box&, Random & = Random::R239) const;
  virtual int GetSamples(QVector<Vector>&, const Box&, int, Random & = Random::R239) const;
  virtual QVector<Vector> Poisson(double, int, Random & = Random::R239) const;

  virtual ScalarField2 Sample(const Rectangles&, int, int) const;
  virtual ScalarField2 Sample(const Quadrangle&, int) const;
  virtual ScalarField Sample(const Box&, int) const;

  // Ray intersection
  int Roots(const Ray&, const double&, const double&, const double&, double*, int = 1, const double& = 1.0e-4) const;
  int Roots(const Ray&, const double&, const double&, const double&, const double&, const double&, double*, int = 1, const double& = 1.0e-4) const;

  //! Compute the Lipschitz constant of the field function
  virtual double K() const;
  virtual double K(const Box&) const;
  virtual double K(const Sphere&) const;
  virtual double K(const Segment&) const;

  virtual void Polygonize(int, Mesh&, const Box&, const double& = 1e-4) const;
  virtual void PolygonizeOctree(int, Mesh&, const Box&, const double& = 1e-4) const;
  virtual void PolygonizeLucie(const Box&, Mesh&, bool, bool) const;
  virtual void PolygonizeLucie(const double&, Mesh&, bool, bool) const;
  virtual void Dual(int, Mesh&, const Box&) const;

  virtual bool SphereTrace(const Ray&, const double&, const double&, const double&, double&, int&, const double& = 1e-4) const;

  virtual void Polygonize(const Box&, QVector<Triangle>&, const double& = 1e-4, bool = false) const;
  virtual void Polygonize(const Tetrahedra&, QVector<Triangle>&, const double& = 1e-4) const;
  virtual Mesh PolygonizeOctree(const Box&, int, const double& = 1e-4) const;

  virtual void Voxelize(const Box&, int, Box&, QVector<Vector>&) const;
  virtual void Voxelize(int, Voxel&, const Box&) const;

  // Volume
  virtual double Volume(int = 5) const;
  virtual double Volume(const Box&, int) const;
  virtual double Volume(const Box&, int, double&) const;

  // Center of gravity
  virtual Sphere Center(int) const;
  virtual Sphere Center(const Box&, int) const;

  virtual void Colorize(const Mesh&, MeshColor&) const;

  static int Cubes();
protected:
  bool Find(Vector&, bool, const Box&, int = 10000, Random & = Random::R239) const;
protected:
  static const double Epsilon; //!< Epsilon value for partial derivatives
};

/*!
\brief Constructor.

\param s Sign, by default set sign convention to true, i.e., negative values inside and positive outside.
*/
inline AnalyticScalarField::AnalyticScalarField(bool s) :sign(s)
{
}

/*!
\brief Check if the value is considered as inside or outside.

Result depend on the sign definition of the implicit surface.
\param v %Value.
*/
inline bool AnalyticScalarField::Inside(const double& v) const
{
  return (sign == true) ? (v > 0.0) : (v < 0.0);
}

class AnalyticScalarField2
{
protected:
  bool sign = true; //!< Sign convention, used for normal computation.
  Box2 box = Box2::Infinity; //!< Domain, set as infinite for base class.
public:
  AnalyticScalarField2(bool = true);
  virtual double Value(const Vector2&) const;
  virtual Vector2 Gradient(const Vector2&) const;
  virtual Matrix2 Hessian(const Vector2&) const;
  virtual ScalarField2 Sample(const Array2&) const;

  // Inside 
  bool Inside(const double&) const;

  // Curvature
  void Curvature(const Vector2&, double&, double&) const;

  virtual Box2 GetBox() const;
protected:
  Matrix Local(const Vector2&) const;
  static const double epsilon; //!< \htmlonly\epsilon;\endhtmlonly value for partial derivatives
};

/*!
\brief Constructor.

\param s Sign, by default set sign convention to true, i.e., negative values inside and positive outside.
*/
inline AnalyticScalarField2::AnalyticScalarField2(bool s) :sign(s)
{
}

/*!
\brief Return the bounding box.
*/
inline Box2 AnalyticScalarField2::GetBox() const
{
  return box;
}

/*!
\brief Check if the value is considered as inside or outside.

Result depend on the sign definition of the implicit surface.
\param v %Value.
*/
inline bool AnalyticScalarField2::Inside(const double& v) const
{
  return (sign == true) ? (v > 0.0) : (v < 0.0);
}

// Forward declare 
class VectorField2;

class ScalarPoint2
{
protected:
  QPoint p; //!< Point.
  double z; //!< Elevation
public:
  //! Empty.
  ScalarPoint2() :p(0, 0), z(0.0) {}
  explicit ScalarPoint2(const QPoint&, const double&);
  friend bool operator<(const ScalarPoint2& a, const ScalarPoint2& b) { return a.z < b.z; }

  // Access to members
  QPoint Point() const;
  double Scalar() const;
};

/*!
\class ScalarPoint2
\brief Internal class for sorting points.
*/

/*!
\brief Create a scalar point.
\param p Point.
\param s Elevation.
*/
inline ScalarPoint2::ScalarPoint2(const QPoint& p, const double& s) :p(p), z(s)
{
}

/*!
\brief Get point.
*/
inline QPoint ScalarPoint2::Point() const
{
  return p;
}

/*!
\brief Get scalar value.
*/
inline double ScalarPoint2::Scalar() const
{
  return z;
}

class ScalarPointSlope2 :public ScalarPoint2
{
protected:
  double s = 0.0; //!< Slope.
  double w = 0.0; //!< Weight.
public:
  //! Empty.
  ScalarPointSlope2() {}
  ScalarPointSlope2(const QPoint&, const double&, const double&, const double&);
  friend bool operator<(const ScalarPointSlope2& a, const ScalarPointSlope2& b) { return a.z < b.z; }

  // Access to members
  double Slope() const;
  double Weight() const;
};

/*!
\class ScalarPointSlope2
\brief Internal class for sorting points.
*/

/*!
\brief Create a scalar point.
\param p Point.
\param z Elevation.
\param s Slope.
\param w Weight.
*/
inline ScalarPointSlope2::ScalarPointSlope2(const QPoint& p, const double& z, const double& s, const double& w) :ScalarPoint2(p, z), s(s), w(w)
{
}

/*!
\brief Get slope.
*/
inline double ScalarPointSlope2::Slope() const
{
  return s;
}

/*!
\brief Get weight.
*/
inline double ScalarPointSlope2::Weight() const
{
  return w;
}

class FloatArray {
protected:
  QVector<float> a; //!< %Array.
public:
  FloatArray() {}
  FloatArray(int);
  float operator()(int) const;
  float& operator[](int);
};

inline FloatArray::FloatArray(int n)
{
  a.resize(n);
}

/*!
\brief Access to element.
*/
inline float FloatArray::operator()(int i) const
{
  return a.at(i);
}

/*!
\brief Access to element.
*/
inline float& FloatArray::operator[](int i)
{
  return a[i];
}

class Mesh2;
class PointCurve2;
class Histogram;

class ScalarField2 :public Array2
{
protected:
  QVector<double> field; //!< Field samples.
public:
  //! Empty.
  ScalarField2() {  }
  ScalarField2(const Array2&, const double& = 0.0);
  explicit ScalarField2(const Box2&, int, int, const double& = 0.0);
  explicit ScalarField2(const Box2&, int, int, const QVector<double>&);
  explicit ScalarField2(const Box2&, const QImage&, const double& = 0.0, const double& = 256.0 * 256.0 - 1.0, bool = true);
  explicit ScalarField2(const QImage&, const double& = 1.0, const double& = 1.0, bool = true);
  explicit ScalarField2(const QString&, double = 0.0, double = 1.0);

  //! Empty
  virtual ~ScalarField2();

  Array2 GetArray() const;

  void GetRange(double&, double&) const;
  double Integral() const;
  double Sum() const;
  double Average() const;
  double ScaledNorm() const;

  void Symmetry(bool = true, bool = false);
  void Rotate();

  Histogram GetHistogram(int) const;
  Histogram GetHistogram(int, double, double) const;
  QVector<std::pair<double, int>> OldHistogram(int) const;
  QVector<std::pair<double, int>> CumulativeHistogram(int) const;
  QVector<std::pair<double, double>> CumulativeNormedHistogram(int) const;

  virtual Vector2 Gradient(int, int) const;
  virtual Vector2 GradientSmooth(int, int) const;
  virtual Matrix2 Hessian(int, int) const;
  virtual double Value(int, int) const;
  virtual double K() const;

  void CutEpsilon(const double& = 1e-6);

  // Access to elements
  double at(int, int) const;
  double at(const QPoint&) const;
  double& operator()(int, int);
  double& operator()(const QPoint&);

  double at(int) const;
  double& operator[](int);

  double StandardDeviation() const;
  double StandardDeviation(const double&) const;

  // Create image from the scalar field
  QImage CreateImage(bool = true) const;

  QImage CreateImage(const AnalyticPalette&, bool = false, bool = false) const;
  QImage CreateImage(const double&, const double&, const AnalyticPalette&, bool = false) const;

  QImage CreateImage(const Palette&) const;
  QImage CreateImage(const double&, const double&, const GenericPalette&, bool = false) const;
  QImage CreateImage(const GenericPalette&) const;
  QImage CreateImage(const double&, const double&, bool = true) const;

  // PGM format
  bool ExportPGM(const QString&) const;
  bool ExportPGM(const QString&, const double&, const double&) const;

  ScalarField2 SetResolution(int, int, bool = false) const;
  ScalarField2 Sample(const Box2&, int, int, bool = false) const;
  ScalarField2 Sample(const Segment2&, const double&, int, int) const;
  ScalarField2 Sample(const Segment2&, const double&) const;

  ScalarField2 Crop(const QPoint&, const QPoint&) const;

  ScalarField2 DownSample(int) const;
  ScalarField2 UpSample(int, bool = false) const;

  virtual double Value(const Vector2&) const;
  virtual double Triangular(const Vector2&) const;
  virtual double BiCubic(const Vector2&) const;
  virtual double BiCubicValue(const Vector2&) const;

  ScalarField2 SummedAreaTable() const;

  void Fill(const double&);
  void Subdivide();

  // Local editing
  void Smooth(const Vector2&, double);
  void Flatten(const Vector2&, const double&, const double& = 0.25);
  void Multiply(const Vector2&, const double&, const double& = 1.0);
  void Level(const Vector2&, const double&, const double&);
  void Gaussian(const Vector2&, const double&, const double&);
  void Add(const Vector2&, const double&);

  ScalarField2 MedianFilter() const;

  // Smoothing
  void Smooth();
  void Smooth(int);
  void Blur();
  void Blur(int);
  void SmoothSmall();
  void MedianFilter();
  ScalarField2 GaussianBlur(int) const;
  ScalarField2 GaussianBlur(const double&) const;

  void Translate(const Vector2&);
  void Scale(const Vector2&);
  void Scale(const double&);

  // Functions
  void Normalize();
  void Unitize();

  void Clamp(const double&, const double&);

  ScalarField2& operator*=(const ScalarField2&);
  ScalarField2& operator*= (const double&);
  ScalarField2& operator/= (const double&);

  ScalarField2& operator-=(const ScalarField2&);
  ScalarField2& operator-=(const double&);

  ScalarField2& operator+=(const ScalarField2&);
  ScalarField2& operator+=(const double&);

  ScalarField2 operator-() const;

  friend ScalarField2 operator-(const double&, const ScalarField2&);
  friend ScalarField2 operator-(const ScalarField2&, const double&);
  friend ScalarField2 operator-(const ScalarField2&, const ScalarField2&);

  friend ScalarField2 operator+(const ScalarField2&, const ScalarField2&);
  friend ScalarField2 operator+(const ScalarField2&, const double&);

  friend ScalarField2 operator*(const ScalarField2&, const double&);
  friend ScalarField2 operator*(const double&, const ScalarField2&);
  friend ScalarField2 operator*(const ScalarField2&, const ScalarField2&);

  friend ScalarField2 operator/(const ScalarField2&, const double&);

  void Pow(const double&);
  void Atan();
  void Step(const double&, const double&);
  void SetRange(const double&, const double&);
  void Binarize(const double&);
  void Invert();
  void Negate();

  void SmoothStep(const double&, const double&, bool = false);
  void Abs();

  void Lerp(const ScalarField2&, const ScalarField2&, const double&);
  void Lerp(const ScalarField2&, const ScalarField2&, const ScalarField2&);

  void AdaptiveBilateralFiltering(double, double);
  void AdaptiveBilateralFiltering(int, double);

  // Morphological operators (for binary fields)
  ScalarField2 Convolution(double[], int);

  void MorphDilate();
  void MorphErode();
  void MorphRemoveEnds();
  void MorphHitAndMiss(double[], int);
  void MorphThin(double[], int);
  ScalarField2 MorphSkeleton(int);
  ScalarField2 MorphSkeletonConnected(int);

  // Gradient
  VectorField2 Gradient() const;
  VectorField2 GradientSmooth() const;
  ScalarField2 GradientNorm() const;
  ScalarField2 GradientSmoothNorm() const;

  ScalarField2 Ln() const;

  ScalarField2 Sqrt() const;
  ScalarField2 Cbrt() const;
  ScalarField2 Qurt() const;
  void Sqrted();

  void SelfOp(const std::function<double(double)>&);
  ScalarField2 Op(const std::function<double(double)>&) const;
  ScalarField2 Op(const std::function<double(const ScalarField2&, int, int)>&) const;

  // Laplacian
  virtual ScalarField2 Laplacian() const;
  virtual double Laplacian(int, int) const;

  QVector<ScalarPoint2> GetScalarPoints() const;
  FloatArray GetAsFloats() const;

  Vector2 Dichotomy(Vector2, Vector2, double, double, double, const double&, const double&) const;

  int Memory() const;

  ScalarField2 DistanceTransform(const double&) const;

  void SignalError(const ScalarField2&, double&, double&) const;
  void SignalError(const ScalarField2&, double&, double&, double&) const;

  ScalarField2 CurvatureProfile() const;
  ScalarField2 CurvatureContour() const;
  ScalarField2 CurvatureTangential() const;
  ScalarField2 CurvatureGaussian() const;
  ScalarField2 CurvatureMean() const;
  ScalarField2 CurvatureMin() const;
  ScalarField2 CurvatureMax() const;

  ScalarField2 MaxFilter(int) const;
  ScalarField2 MinFilter(int) const;

  static ScalarField2 LoadFromR32(const Box2&, const QString&);
  static ScalarField2 Load(Box2, const QString&, double = 0.0, double = 1.0);

  Mesh2 Polygonize() const;
  SegmentSet2 LineSegments(const double&, bool = false) const;

  friend std::ostream& operator<<(std::ostream&, const ScalarField2&);

  void Debug() const;

public:
  // Provide range-based for loops
  auto begin() { return field.begin(); }
  auto end() { return field.end(); }
  auto cbegin() const { return field.begin(); }
  auto cend() const { return field.end(); }
  auto begin() const { return field.begin(); }
  auto end() const { return field.end(); }



protected:
  void SetBorder(const double& = 0.0);
  double SmoothPoint(int, int) const;
  Matrix Local(int, int) const;

  void CopyEdges(ScalarField2&) const;
};


/*!
\brief Return the field value at a given array vertex.
\param i,j Integer coordinates of the vertex.
*/
inline double ScalarField2::at(int i, int j) const
{
  return field.at(VertexIndex(i, j));
}

/*!
\brief Return the field value at a given array vertex.
\param q Point.
*/
inline double ScalarField2::at(const QPoint& q) const
{
  return field.at(VertexIndex(q.x(), q.y()));
}

/*!
\brief Return the field value at a given array vertex.
\param q Point.
*/
inline double& ScalarField2::operator()(const QPoint& q)
{
  return field[VertexIndex(q.x(), q.y())];
}

/*!
\brief Return the field value at a given array vertex.
\param i,j Integer coordinates of the vertex.
*/
inline double& ScalarField2::operator()(int i, int j)
{
  return field[VertexIndex(i, j)];
}

/*!
\brief Return the data in the field.
\param c Index.
*/
inline double ScalarField2::at(int c) const
{
  return field.at(c);
}

/*!
\brief Return the data in the field.
\param c Index.
*/
inline double& ScalarField2::operator[](int c)
{
  return field[c];
}

class ScalarField :public Array
{
protected:
  bool sign; //!< Sign convention, used for normal computation.
  QVector<double> field; //!< Field samples.
public:
  //! Empty.
  ScalarField(bool = true);
  ScalarField(const Array&, const double& = 0.0, bool = true);
  ScalarField(const Box&, int, int, int, const double& = 0.0, bool = true);
  ScalarField(const Box&, int, int, int, const QVector<double>&, bool = true);

  //! Empty
  ~ScalarField() {}

  void GetRange(double&, double&) const;
  double Integral() const;
  double Sum() const;

  virtual Vector Gradient(int, int, int) const;
  virtual Matrix Hessian(int, int, int) const;
  virtual double Value(int, int, int) const;
  virtual double K() const;

  // Normal
  Vector Normal(const Vector&) const;

  void CutEpsilon(const double& = 1e-6);

  // Access to elements
  double at(int, int, int) const;
  double& operator()(int, int, int);

  double at(int) const;
  double& operator[](int);

  ScalarField& operator+= (const ScalarField&);
  ScalarField& operator*= (const double&);

  virtual double Value(const Vector&) const;
  virtual double BicubicValue(const Vector&) const;
  Vector Gradient(const Vector&) const;
  virtual Matrix Hessian(const Vector&) const;

  void Translate(const Vector&);
  void Scale(const Vector&);

  int Memory() const;

  void Polygonize(Mesh&) const;
  void Dual(Mesh&) const;

  int DualCellCode(int, int, int) const;

  void Save(const QString&) const;
  void Read(const QString&);

  void DualVertex(int, int, int, int, Vector&) const;
public:
  static int TriangleTable[256][16]; //!< Two dimensionnal array storing the straddling edges for every marching cubes configuration.
  static int edgeTable[256];    //!< Array storing straddling edges for every marching cubes configuration.
protected:
  void Blend(const double&, double[4]) const;
protected:
  static constexpr double Epsilon = 1.0e-6; //!< Epsilon value for partial derivatives
};

/*!
\brief Constructor.

\param s Sign, by default, set sign convention to true, i.e., negative values inside and positive outside.
*/
inline ScalarField::ScalarField(bool s) :sign(s)
{
}

/*!
\brief Return the field value at a given array vertex.
\param i,j,k Integer coordinates of the vertex.
*/
inline double ScalarField::at(int i, int j, int k) const
{
  return field.at(VertexIndex(i, j, k));
}

/*!
\brief Return the field value at a given array vertex.
\param i,j,k Integer coordinates of the vertex.
*/
inline double& ScalarField::operator()(int i, int j, int k)
{
  return field[VertexIndex(i, j, k)];
}

/*!
\brief Return the data in the field.
\param c Index.
*/
inline double ScalarField::at(int c) const
{
  return field.at(c);
}

/*!
\brief Return the data in the field.
\param c Index.
*/
inline double& ScalarField::operator[](int c)
{
  return field[c];
}

/*!
\brief Subtract a constant to the values the scalar field.

Provided out of consistency, simply calls ScalarField2::operator+=(const double&) with the negated constant.
\param c Constant.
*/
inline ScalarField2& ScalarField2::operator-=(const double& c)
{
  (*this) += -c;
  return *this;
}

/*!
\brief Scale the values of a scalar field.

\param s Scaling factor.
*/
inline ScalarField2& ScalarField2::operator/=(const double& s)
{
  *this *= 1.0 / s;
  return *this;
}

/*!
\brief Negates the values of the scalar field.
*/
inline ScalarField2 ScalarField2::operator-() const
{
  ScalarField2 r = *this;
  for (int i = 0; i < field.size(); i++)
  {
    r[i] = -r[i];
  }
  return r;
}

/*!
\brief Subtraction.

\param a Scalar field.
\param x Real.
*/
inline ScalarField2 operator-(const ScalarField2& a, const double& x)
{
  ScalarField2 r = a;
  r -= x;
  return r;
}

/*!
\brief Subtraction.

\param x Real.
\param a Scalar field.
*/
inline ScalarField2 operator-(const double& x, const ScalarField2& a)
{
  ScalarField2 r = -a;
  r += x;
  return r;
}

/*!
\brief Subtraction.

\param a,b Scalar fields.
*/
inline ScalarField2 operator-(const ScalarField2& a, const ScalarField2& b)
{
  ScalarField2 r = a;
  r -= b;
  return r;
}

/*!
\brief Addition.

\param a,b Scalar fields.
*/
inline ScalarField2 operator+(const ScalarField2& a, const ScalarField2& b)
{
  ScalarField2 r = a;
  r += b;
  return r;
}

/*!
\brief Addition.

\param a Scalar field.
\param x Real.
*/
inline ScalarField2 operator+(const ScalarField2& a, const double& x)
{
  ScalarField2 r = a;
  r += x;
  return r;
}

/*!
\brief Multiplication.

\param a Scalar field.
\param x Real.
*/
inline ScalarField2 operator*(const ScalarField2& a, const double& x)
{
  ScalarField2 r = a;
  r *= x;
  return r;
}

/*!
\brief Division.

Rely on multiplication by inverse.

\param a Scalar field.
\param x Real.
*/
inline ScalarField2 operator/(const ScalarField2& a, const double& x)
{
  ScalarField2 r = a;
  r *= 1.0 / x;
  return r;
}

/*!
\brief Multiplication.

\param a Scalar field.
\param x Real.
*/
inline ScalarField2 operator*(const double& x, const ScalarField2& a)
{
  ScalarField2 r = a;
  r *= x;
  return r;
}

/*!
\brief Multiplication.

\param a Scalar field.
\param b Scalar field.
*/
inline ScalarField2 operator*(const ScalarField2& a, const ScalarField2& b)
{
  ScalarField2 r = a;
  r *= b;
  return r;
}
