// Matrix 

#pragma once

// Include Vector class
#include "libs/evector.h"
#include "libs/random.h"
#include "libs/cubic.h"

static RandomFast random239(239);

class Matrix2 {
protected:
  double r[4] = { 1.0,0.0,0.0,1.0 };  //!< The array storing the coefficients of the matrix.
public:
  //! Empty.
  Matrix2() {}
  explicit Matrix2(const double&);
  explicit Matrix2(const Vector2&);
  explicit Matrix2(const Vector2&, const Vector2&);
  explicit Matrix2(const double&, const double&, const double&, const double&);
  explicit Matrix2(const double&, const double&, const double&);

  static const Matrix2 Null;
  static const Matrix2 Identity;

  // Functions to access class components  
  constexpr double& operator[] (int);
  constexpr double operator[] (int) const;

  constexpr double& operator() (int, int);
  constexpr double operator() (int, int) const;

  // Access to column vectors
  Vector2 C(int) const;

  // Unary operators
  Matrix2 operator+ () const { return *this; }
  Matrix2 operator- () const;

  // Assignment operators
  Matrix2& operator+= (const Matrix2&);
  Matrix2& operator-= (const Matrix2&);
  Matrix2& operator*= (const Matrix2&);
  Matrix2& operator*= (double);
  Matrix2& operator/= (double);

  Matrix2 Scaled(const Matrix2&) const;

  // Binary operators
  friend Matrix2 operator+ (const Matrix2&, const Matrix2&);
  friend Matrix2 operator- (const Matrix2&, const Matrix2&);
  friend Matrix2 operator* (const Matrix2&, const Matrix2&);

  Vector2 operator*(const Vector2&) const;

  friend Matrix2 operator* (const Matrix2&, const double&);
  friend Matrix2 operator* (const double&, const Matrix2&);

  friend Matrix2 operator/ (const Matrix2&, const double&);
  friend Matrix2 operator/ (const double&, const Matrix2&);

  static Matrix2 Rotation(const double&);
  static Matrix2 Rotation(RandomFast & = random239);
  static Matrix2 Tensor2(const double&);
  static Matrix2 Tensor2(const double&, const double&);

  friend std::ostream& operator<<(std::ostream&, const Matrix2&);

  friend QString ToGLSL(const Matrix2&);

  // Transpose
  Matrix2 T() const;

  // Determinant
  double Determinant() const;

  Matrix2 Inverse() const;
  double Trace() const;

  double FrobeniusNorm() const;
  Quadric Characteristic() const;
  double SpectralNorm() const;
  double InfinityNorm() const;

  // Eigensolver, matrix must be symmetric
  void EigenSolveSymmetric(double[2], Vector2[2]) const;
  Vector2 Eigen() const;

  static const Matrix2 RotationHalfPi; //!< Half-pi rotation.
  static const Matrix2 RotationTwoThirdsPi; //!< Two thirds pi rotation.
private:
  static double epsilon; //!< Epsilon value used to check angles of a rotation matrix, and in the implementation of some algebra methods.
};

//! Direct access to the array of the matrix.
inline constexpr double& Matrix2::operator[] (int i)
{
  return r[i];
}

//! Overloaded.
inline constexpr double Matrix2::operator[] (int i) const
{
  return r[i];
}

/*!
\brief Get element (i,j) of the matrix.
\param i Row.
\param j Column.
*/
inline constexpr double& Matrix2::operator() (int i, int j)
{
  return r[i + j + j];
}

//! Overloaded.
inline constexpr double Matrix2::operator() (int i, int j) const
{
  return r[i + j + j];
}

/*!
\brief Creates a matrix with a set of double values.
*/
inline Matrix2::Matrix2(const double& a00, const double& a01, const double& a10, const double& a11)
{
  r[0] = a00;
  r[1] = a01;

  r[2] = a10;
  r[3] = a11;
}

/*!
\brief Create a symmetric matrix.
\param a,b Diagonal terms.
\param c Last term.
*/
inline Matrix2::Matrix2(const double& a, const double& b, const double& c) :Matrix2(a, c, c, b)
{
}

/*!
\brief Creates a matrix with the same diagonal value.
\param a The diagonal value.
*/
inline Matrix2::Matrix2(const double& a)
{
  r[1] = r[2] = 0.0;
  r[0] = r[3] = a;
}

/*!
\brief Create a diagonal matrix with diagonal terms set to the vector entries.

Note that this is a scaling matrix.

\param a %Vector of diagonal values
*/
inline Matrix2::Matrix2(const Vector2& a)
{
  r[1] = r[2] = 0.0;
  r[0] = a[0];
  r[3] = a[1];
}

/*!
\brief Get the i-th column vector of the matrix.
\param i Index.
*/
inline Vector2 Matrix2::C(int i) const
{
  return Vector2(r[i * 2 + 0], r[i * 2 + 1]);
}

/*!
\brief Creates a column vector matrix.

\param a,b Column vectors.
*/
inline Matrix2::Matrix2(const Vector2& a, const Vector2& b)
{
  r[0] = a[0];
  r[1] = a[1];

  r[2] = b[0];
  r[3] = b[1];
}

/*!
\brief Transpose the matrix.

Note that if the matrix is a rotation matrix, then the transpose is equal to its inverse.

\sa Matrix::T()
*/
inline Matrix2 Matrix2::T() const
{
  return Matrix2(r[0], r[2], r[1], r[3]);
}

/*!
\brief Computes the determinant of the matrix.
*/
inline double Matrix2::Determinant() const
{
  return r[0] * r[3] - r[1] * r[2];
}

//! Returns the opposite of a matrix -A.
inline Matrix2 Matrix2::operator-() const
{
  return Matrix2(-r[0], -r[1], -r[2], -r[3]);
}

/*!
\brief Compute the trace (sum of diagonal terms) of a matrix.
*/
inline double Matrix2::Trace() const
{
  return r[0] + r[3];
}

/*!
\brief Right multiply by a vector.
*/
inline Vector2 Matrix2::operator*(const Vector2& v) const
{
  return Vector2(v[0] * r[0] + v[1] * r[2], v[0] * r[1] + v[1] * r[3]);
}

//! Overloaded. 
inline Matrix2 operator-(const Matrix2& u, const Matrix2& v)
{
  return Matrix2(u[0] - v[0], u[1] - v[1], u[2] - v[2], u[3] - v[3]);
}

//! Overloaded. 
inline Matrix2 operator+(const Matrix2& u, const Matrix2& v)
{
  return Matrix2(u[0] + v[0], u[1] + v[1], u[2] + v[2], u[3] + v[3]);
}

//! Right multiply by a double. 
inline Matrix2 operator*(const Matrix2& u, const double& a)
{
  return Matrix2(a * u[0], a * u[1], a * u[2], a * u[3]);
}

//! Left multiply by a double. 
inline Matrix2 operator*(const double& a, const Matrix2& u)
{
  return Matrix2(a * u[0], a * u[1], a * u[2], a * u[3]);
}

//! Right divide by a double. 
inline Matrix2 operator/(const Matrix2& u, const double& a)
{
  const double ia = 1.0 / a;
  return u * ia;
}

class Matrix {
protected:
  double r[9] = { 1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0 };  //!< The array storing the coefficients of the matrix.
public:
  //! Empty.
  Matrix() {}
  explicit Matrix(const double&);
  explicit Matrix(const Matrix2&);
  explicit Matrix(const Vector&);
  explicit Matrix(const Vector&, const Vector&, const Vector&);
  explicit Matrix(const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);

  static const Matrix Null;
  static const Matrix Identity;

  // Functions to access Matrix class components  
  constexpr double& operator[] (int);
  constexpr double operator[] (int) const;

  constexpr double& operator() (int, int);
  constexpr double operator() (int, int) const;

  // Access to column vectors
  Vector C(int) const;
  Vector Row(int) const;

  // Unary operators
  Matrix operator+ () const { return *this; }
  Matrix operator- () const;

  // Assignment operators
  Matrix& operator+= (const Matrix&);
  Matrix& operator-= (const Matrix&);
  Matrix& operator*= (const Matrix&);
  Matrix& operator*= (double);
  Matrix& operator/= (double);

  // Binary operators
  friend Matrix operator+ (const Matrix&, const Matrix&);
  friend Matrix operator- (const Matrix&, const Matrix&);
  friend Matrix operator* (const Matrix&, const Matrix&);

  Vector operator*(const Vector&) const;

  friend Matrix operator* (const Matrix&, const double&);
  friend Matrix operator* (const double&, const Matrix&);

  friend Matrix operator/ (const Matrix&, const double&);
  friend Matrix operator/ (const double&, const Matrix&);

  static Matrix Rotation(const Vector&);
  static Matrix RotationMaya(const Vector&);
  static Matrix Rotation(const Vector&, const double&);
  static Matrix Rotation(const Vector&, const Vector&);
  static Matrix RotationCanonical(const Vector&);

  static Matrix RotationX(const double&);
  static Matrix RotationY(const double&);
  static Matrix RotationZ(const double&);

  static Matrix Symmetry(const Vector&);

  Vector GetRotationAngles() const;
  static Matrix Covariance(Vector*, int);

  static Matrix Lerp(const double&, const Matrix&, const Matrix&);

  static Matrix Rotation(RandomFast & = random239);

  friend std::ostream& operator<<(std::ostream&, const Matrix&);

  friend QString ToGLSL(const Matrix&);

  // Absolute value
  Matrix Abs() const;

  // Transpose
  Matrix T() const;

  // Adjoint
  Matrix Adjoint() const;

  // Determinant
  double Determinant() const;

  friend Matrix Inverse(const Matrix&);
  double Trace() const;

  // Algebra
  Cubic Characteristic() const;
  double SpectralNorm() const;
  double InfinityNorm() const;
  double FrobeniusNorm() const;

  // Singular values methods
  void SingularValueDecomposition(Matrix&, Vector&, Matrix&) const;

  void QDU(Matrix&, Vector&, Vector&) const;

  void ExtractAngleAxis(Vector&, double&) const;

  // Eigensolver, matrix must be symmetric
  void EigenSolveSymmetric(double[3], Vector[3]) const;
  Vector Eigen() const;

  void Tridiagonal(double[3], double[2]);
  int QLAlgorithm(double[3], double[3]);

  static Matrix Frenet(const Vector&, const Vector & = Vector::Z);

  static const Matrix RotationHalfPiZ;   //!< Half-pi rotation around vertical axis.
private:
  void Bidiagonalize(Matrix&, Matrix&);
  void GolubKahanStep(Matrix&, Matrix&);
private:
  static double epsilon; //!< Epsilon value used to check angles of a rotation matrix, and in the implementation of some algebra methods.
};

//! Direct access to the array of the matrix.
inline constexpr double& Matrix::operator[] (int i)
{
  return r[i];
}

//! Overloaded.
inline constexpr double Matrix::operator[] (int i) const
{
  return r[i];
}

/*!
\brief Get element (i,j) of the matrix.
\param i Row.
\param j Column.
*/
inline constexpr double& Matrix::operator() (int i, int j)
{
  return r[i + j + j + j];
}

//! Overloaded.
inline constexpr double Matrix::operator() (int i, int j) const
{
  return r[i + j + j + j];
}

/*!
\brief Creates a matrix with a set of double values.

Coefficients are given in column order.
*/
inline Matrix::Matrix(const double& a00, const double& a01, const double& a02, const double& a10, const double& a11, const double& a12, const double& a20, const double& a21, const double& a22)
{
  r[0] = a00;
  r[1] = a01;
  r[2] = a02;

  r[3] = a10;
  r[4] = a11;
  r[5] = a12;

  r[6] = a20;
  r[7] = a21;
  r[8] = a22;
}

/*!
\brief Creates a diagonal matrix.
\param a The diagonal value.
*/
inline Matrix::Matrix(const double& a)
{
  r[1] = r[2] = r[3] = r[5] = r[6] = r[7] = 0.0;
  r[0] = r[4] = r[8] = a;
}

/*!
\brief Creates a matrix form a lower dimension.
\param a The matrix.
*/
inline Matrix::Matrix(const Matrix2& a)
{
  r[0] = a[0];
  r[1] = a[1];
  r[3] = a[2];
  r[4] = a[3];
  r[2] = r[5] = r[6] = r[7] = 0.0;
  r[8] = 1.0;
}

/*!
\brief Create a diagonal matrix with diagonal terms set to the vector entries.
\param a %Vector of diagonal values.
*/
inline Matrix::Matrix(const Vector& a)
{
  r[1] = r[2] = r[3] = r[5] = r[6] = r[7] = 0.0;
  r[0] = a[0];
  r[4] = a[1];
  r[8] = a[2];
}

/*!
\brief Get the i-th column vector of the matrix.
\param i Index.
*/
inline Vector Matrix::C(int i) const
{
  return Vector(r[i * 3 + 0], r[i * 3 + 1], r[i * 3 + 2]);
}

/*!
\brief Get the i-th row vector of the matrix.
\param i Index.
*/
inline Vector Matrix::Row(int i) const
{
  return Vector(r[i], r[i + 3], r[i + 6]);
}

/*!
\brief Creates a matrix given column vectors.

\param a,b,c Column vectors.
*/
inline Matrix::Matrix(const Vector& a, const Vector& b, const Vector& c)
{
  r[0] = a[0];
  r[1] = a[1];
  r[2] = a[2];

  r[3] = b[0];
  r[4] = b[1];
  r[5] = b[2];

  r[6] = c[0];
  r[7] = c[1];
  r[8] = c[2];
}

/*!
\brief Transpose the matrix.

Note that if the matrix is a rotation matrix, then the transpose is equal to its inverse.
\sa T(const Matrix&)
*/
inline Matrix Matrix::T() const
{
  return Matrix(r[0], r[3], r[6], r[1], r[4], r[7], r[2], r[5], r[8]);
}

/*!
\brief Computes the determinant of the matrix.
*/
inline double Matrix::Determinant() const
{
  return r[0] * r[4] * r[8] + r[1] * r[5] * r[6] + r[2] * r[3] * r[7] - r[2] * r[4] * r[6] - r[1] * r[3] * r[8] - r[0] * r[5] * r[7];
}

//! Returns the opposite of a matrix -A.
inline Matrix Matrix::operator-() const
{
  return Matrix(-r[0], -r[1], -r[2], -r[3], -r[4], -r[5], -r[6], -r[7], -r[8]);
}

/*!
\brief Compute the trace (sum of diagonal terms) of a matrix.
*/
inline double Matrix::Trace() const
{
  return r[0] + r[4] + r[8];
}

/*!
\brief Right multiply by a vector.

Right multiplication by a vector is a member function and not a friend to avoid left/right
multiplication ambiguities since matrix multiplication is not commutative.

\param v %Vector.
\code
Vector v=Matrix::Rotation(Vector(0.0,Math::Pi/3.0,Math::Pi/4))*Vector(1.0,1.0,1.0);
\endcode
*/
inline Vector Matrix::operator*(const Vector& v) const
{
  return Vector(v[0] * r[0] + v[1] * r[3] + v[2] * r[6], v[0] * r[1] + v[1] * r[4] + v[2] * r[7], v[0] * r[2] + v[1] * r[5] + v[2] * r[8]);
}

//! Overloaded. 
inline Matrix operator-(const Matrix& u, const Matrix& v)
{
  return Matrix(u[0] - v[0], u[1] - v[1], u[2] - v[2], u[3] - v[3], u[4] - v[4], u[5] - v[5], u[6] - v[6], u[7] - v[7], u[8] - v[8]);
}

//! Overloaded. 
inline Matrix operator+(const Matrix& u, const Matrix& v)
{
  return Matrix(u[0] + v[0], u[1] + v[1], u[2] + v[2], u[3] + v[3], u[4] + v[4], u[5] + v[5], u[6] + v[6], u[7] + v[7], u[8] + v[8]);
}

//! Right multiply by a double. 
inline Matrix operator*(const Matrix& u, const double& a)
{
  return Matrix(a * u[0], a * u[1], a * u[2], a * u[3], a * u[4], a * u[5], a * u[6], a * u[7], a * u[8]);
}

//! Left multiply by a double. 
inline Matrix operator*(const double& a, const Matrix& u)
{
  return Matrix(a * u[0], a * u[1], a * u[2], a * u[3], a * u[4], a * u[5], a * u[6], a * u[7], a * u[8]);
}

//! Right divide by a double. 
inline Matrix operator/(const Matrix& u, const double& a)
{
  const double ia = 1.0 / a;
  return u * ia;
}

// Extended matrix
class Matrix4 {
protected:
  double r[16] = { 1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0 }; //!< Coefficients.
public:
  //! Empty
  Matrix4() {}
  explicit Matrix4(const double&);
  explicit Matrix4(const Vector&);
  explicit Matrix4(const Matrix&);
  explicit Matrix4(const Matrix&, const Vector&);
  explicit Matrix4(const Matrix&, const Vector&, const Vector&);
  explicit Matrix4(const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);

  static const Matrix4 Null; //!< Null matrix.
  static const Matrix4 Identity; //!< Identity matrix.

  // Functions to access class components
  constexpr double& operator[] (int i) { return r[i]; }
  constexpr const double& operator[] (int i) const { return r[i]; }

  constexpr double& operator() (int i, int j) { return r[i + (j << 2)]; }
  constexpr const double& operator() (int i, int j) const { return r[i + (j << 2)]; }

  // Unary operators
  Matrix4 operator+ () const { return *this; }
  Matrix4 operator- () const;

  // Assignment operators
  Matrix4& operator+= (const Matrix4&);
  Matrix4& operator-= (const Matrix4&);
  Matrix4& operator*= (const Matrix4&);
  Matrix4& operator*= (double);
  Matrix4& operator/= (double);

  // Binary operators
  friend Matrix4 operator+ (const Matrix4&, const Matrix4&);
  friend Matrix4 operator- (const Matrix4&, const Matrix4&);
  friend Matrix4 operator* (const Matrix4&, const Matrix4&);
  friend Matrix4 operator* (const Matrix4&, const double&);

  Vector operator* (const Vector&) const;

  friend QString ToGLSL(const Matrix4&);

  // Determinant
  double Determinant() const;

  double Trace() const;

  // Transpose
  Matrix4 T() const;
  friend Matrix4 Inverse(const Matrix4&);

  static Matrix4 LookAt(const Vector&, const Vector&, const double&);
  static Matrix4 Rotation(const Vector&, const double&);
  static Matrix4 Rotation(const Vector&);
  static Matrix4 Scale(const Vector&);
  static Matrix4 Translate(const Vector&);
  static Matrix4 Shear(const Vector&);

  // Float
  void Float(float[16]) const;

  // Get
  Matrix Sub() const;

  friend std::ostream& operator<<(std::ostream&, const Matrix4&);
public:
  static const Matrix4 Hermite; //!< Hermite matrix, usefull for bi-cubic Hermite interpolation.
protected:
  static const double epsilon; //!< Epsilon value used to check if the determinant of a matrix is null (used in Inverse(const Matrix4&)).
};

/*!
\brief Create a matrix given its coefficients.
*/
inline Matrix4::Matrix4(const double& a0, const double& a1, const double& a2, const double& a3, const double& a4, const double& a5, const double& a6, const double& a7, const double& a8, const double& a9, const double& a10, const double& a11, const double& a12, const double& a13, const double& a14, const double& a15)
{
  r[0] = a0;
  r[1] = a1;
  r[2] = a2;
  r[3] = a3;
  r[4] = a4;
  r[5] = a5;
  r[6] = a6;
  r[7] = a7;
  r[8] = a8;
  r[9] = a9;
  r[10] = a10;
  r[11] = a11;
  r[12] = a12;
  r[13] = a13;
  r[14] = a14;
  r[15] = a15;
}

/*!
\brief Compute the trace (sum of diagonal terms).
*/
inline double Matrix4::Trace() const
{
  return r[0] + r[5] + r[10] + r[15];
}

/*!
\brief Return the sub-matrix obtained by removing the last column and row.
*/
inline Matrix Matrix4::Sub() const
{
  return Matrix(Vector(r[0], r[1], r[2]), Vector(r[4], r[5], r[6]), Vector(r[8], r[9], r[10]));
}

/*!
\brief Multiplication by a real.
\param A matrix.
\param r Real.
*/
inline Matrix4 operator*(const Matrix4& A, const double& r)
{
  return Matrix4(A.r[0] * r, A.r[1] * r, A.r[2] * r, A.r[3] * r, A.r[4] * r, A.r[5] * r, A.r[6] * r, A.r[7] * r, A.r[8] * r, A.r[9] * r, A.r[10] * r, A.r[11] * r, A.r[12] * r, A.r[13] * r, A.r[14] * r, A.r[15] * r);
}
