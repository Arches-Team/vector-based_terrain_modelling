// Matrix

#include "libs/matrix.h"

/*!
\class Matrix4 matrix.h
\brief This class implements 4<SUP>2</SUP> matrix.

Components are stored in a single dimension array, starting from
element a<SUB>00</SUB> and sorting components by column.

The diagonal elements of a matrix A are A<SUB>00</SUB>=A[0], A<SUB>11</SUB>=A[5], A<SUB>22</SUB>=A[10], A<SUB>33</SUB>=A[15].
The terms corresponding to the translation vector are A[12], A[13], A[14].
\ingroup MathGroup
*/

const double Matrix4::epsilon = 1.0e-18;

const Matrix4 Matrix4::Null(0.0);

const Matrix4 Matrix4::Identity(1.0);

const Matrix4 Matrix4::Hermite(
  2.0, -3.0, 0.0, 1.0,
  -2.0, 3.0, 0.0, 0.0,
  1.0, -2.0, 1.0, 0.0,
  1.0, -1.0, 0.0, 0.0);

/*!
\brief Create an scaling matrix with the same diagonal terms.

The last term of the matrix is set to 1.0.
\param a Value of diagonal entries.
*/
Matrix4::Matrix4(const double& a) :r{ a,0.0,0.0,0.0,0.0,a,0.0,0.0,0.0,0.0,a,0.0,0.0,0.0,0.0,1.0 }
{
  // *this = Matrix4::Null;
  // r[0] = r[5] = r[10] = a;
  // r[15] = 1.0;
}

/*!
\brief Create a homogeneous diagonal matrix with diagonal terms
set to the vector entries.

Last diagonal entry is set to 1.0.
\param a Vector of diagonal entries.
*/
Matrix4::Matrix4(const Vector& a) :r{ a[0],0.0,0.0,0.0,0.0,a[1],0.0,0.0,0.0,0.0,a[2],0.0,0.0,0.0,0.0,1.0 }
{
  //  *this = Matrix4::Null;
  //  r[0] = a[0];
  //  r[5] = a[1];
  //  r[10] = a[2];
  //  r[15] = 1.0;
}

/*!
\brief Create an homogeneous matrix from a simple Matrix.

The translation and shear coefficients are set to 0.0.
\param a Matrix.
*/
Matrix4::Matrix4(const Matrix& a)
{
  // Rotation and scale
  r[0] = a[0];
  r[1] = a[1];
  r[2] = a[2];
  r[4] = a[3];
  r[5] = a[4];
  r[6] = a[5];
  r[8] = a[6];
  r[9] = a[7];
  r[10] = a[8];

  // Translation
  r[3] = r[7] = r[11] = 0.0;

  // Shear
  r[12] = r[13] = r[14] = 0.0;

  // Scale
  r[15] = 1.0;
}

/*!
\brief Create an homogeneous matrix from a simple Matrix and a translation vector

The shear coefficients are set to 0.0.
\param a Matrix.
\param t Translation vector.
*/
Matrix4::Matrix4(const Matrix& a, const Vector& t) :Matrix4(a, t, Vector::Null)
{
}

/*!
\brief Constructor from a Matrix, a translation and a shear vector.
\param a Matrix.
\param t Translation vector.
\param s Shear vector.
*/
Matrix4::Matrix4(const Matrix& a, const Vector& t, const Vector& s)
{
  // Rotation and scale
  r[0] = a[0];
  r[1] = a[1];
  r[2] = a[2];
  r[4] = a[3];
  r[5] = a[4];
  r[6] = a[5];
  r[8] = a[6];
  r[9] = a[7];
  r[10] = a[8];

  // Translation
  r[12] = t[0];
  r[13] = t[1];
  r[14] = t[2];

  // Shear
  r[3] = s[0];
  r[7] = s[1];
  r[11] = s[2];

  // Scale
  r[15] = 1.0;
}

/*!
\brief Convert to a generic 4&times;4 float matrix.
\param a Returned matrix.
*/
void Matrix4::Float(float a[16]) const
{
  for (int i = 0; i < 16; i++)
  {
    a[i] = float(r[i]);
  }
}

/*!
\brief Transpose a matrix.
*/
Matrix4 Matrix4::T() const
{
  return Matrix4(r[0], r[4], r[8], r[12], r[1], r[5], r[9], r[13], r[2], r[6], r[10], r[14], r[3], r[7], r[11], r[15]);
}

/*!
\brief Returns the opposite of a matrix -A.
*/
Matrix4 Matrix4::operator-() const
{
  Matrix4 n;
  for (int i = 0; i < 16; i++)
  {
    n.r[i] = -r[i];
  }
  return n;
}

/*!
\brief Multiplication.
\param u,v Argument matrix.
*/
Matrix4 operator*(const Matrix4& u, const Matrix4& v)
{
  Matrix4 a;
  for (int i = 0; i < 4; i++)
  {
    int k = i << 2;
    for (int j = 0; j < 4; j++)
    {
      a.r[k + j] = u.r[j] * v.r[k] + u.r[4 + j] * v.r[k + 1] + u.r[8 + j] * v.r[k + 2] + u.r[12 + j] * v.r[k + 3];
    }
  }
  return a;
}

/*!
\brief Destructive multiplication.
\param v %Matrix.
*/
Matrix4& Matrix4::operator*=(const Matrix4& v)
{
  Matrix4 r;

  r = (*this) * v;
  *this = r;
  return *this;
}

/*!
\brief Right multiply by a vector.
\param v %Vector.
*/
Vector Matrix4::operator*(const Vector& v) const
{
  double w = 1.0 / (v[0] * r[3] + v[1] * r[7] + v[2] * r[11] + r[15]);

  Vector res;
  for (int i = 0; i < 3; i++)
  {
    res[i] = (v[0] * r[i] + v[1] * r[4 + i] + v[2] * r[8 + i] + r[12 + i]) * w;
  }
  return res;
}

/*!
\brief Computes the inverse of a matrix A<SUP>-1</SUP>.

This function returns the null matrix if A cannot be inverted.
The threshold value involved in the singular matrix detection is
set to 10<SUP>-18</SUP>.

\param m Argument matrix.
*/
Matrix4 Inverse(const Matrix4& m)
{
  double d00 = m[5] * m[10] * m[15] + m[9] * m[14] * m[7] + m[13] * m[6] * m[11] - m[7] * m[10] * m[13] - m[11] * m[14] * m[5] - m[15] * m[6] * m[9];
  double d01 = m[1] * m[10] * m[15] + m[9] * m[14] * m[3] + m[13] * m[2] * m[11] - m[3] * m[10] * m[13] - m[11] * m[14] * m[1] - m[15] * m[2] * m[9];
  double d02 = m[1] * m[6] * m[15] + m[5] * m[14] * m[3] + m[13] * m[2] * m[7] - m[3] * m[6] * m[13] - m[7] * m[14] * m[1] - m[15] * m[2] * m[5];
  double d03 = m[1] * m[6] * m[11] + m[5] * m[10] * m[3] + m[9] * m[2] * m[7] - m[3] * m[6] * m[9] - m[7] * m[10] * m[1] - m[11] * m[2] * m[5];

  // Test if singular matrix
  double d = m[0] * d00 - m[4] * d01 + m[8] * d02 - m[12] * d03;

  if (fabs(d) < Matrix4::epsilon)
  {
    return Matrix4::Null;
  }

  double d10 = m[4] * m[10] * m[15] + m[8] * m[14] * m[7] + m[12] * m[6] * m[11] - m[7] * m[10] * m[12] - m[11] * m[14] * m[4] - m[15] * m[6] * m[8];
  double d11 = m[0] * m[10] * m[15] + m[8] * m[14] * m[3] + m[12] * m[2] * m[11] - m[3] * m[10] * m[12] - m[11] * m[14] * m[0] - m[15] * m[2] * m[8];
  double d12 = m[0] * m[6] * m[15] + m[4] * m[14] * m[3] + m[12] * m[2] * m[7] - m[3] * m[6] * m[12] - m[7] * m[14] * m[0] - m[15] * m[2] * m[4];
  double d13 = m[0] * m[6] * m[11] + m[4] * m[10] * m[3] + m[8] * m[2] * m[7] - m[3] * m[6] * m[8] - m[7] * m[10] * m[0] - m[11] * m[2] * m[4];

  double d20 = m[4] * m[9] * m[15] + m[8] * m[13] * m[7] + m[12] * m[5] * m[11] - m[7] * m[9] * m[12] - m[11] * m[13] * m[4] - m[15] * m[5] * m[8];
  double d21 = m[0] * m[9] * m[15] + m[8] * m[13] * m[3] + m[12] * m[1] * m[11] - m[3] * m[9] * m[12] - m[11] * m[13] * m[0] - m[15] * m[1] * m[8];
  double d22 = m[0] * m[5] * m[15] + m[4] * m[13] * m[3] + m[12] * m[1] * m[7] - m[3] * m[5] * m[12] - m[7] * m[13] * m[0] - m[15] * m[1] * m[4];
  double d23 = m[0] * m[5] * m[11] + m[4] * m[9] * m[3] + m[8] * m[1] * m[7] - m[3] * m[5] * m[8] - m[7] * m[9] * m[0] - m[11] * m[1] * m[4];

  double d30 = m[4] * m[9] * m[14] + m[8] * m[13] * m[6] + m[12] * m[5] * m[10] - m[6] * m[9] * m[12] - m[10] * m[13] * m[4] - m[14] * m[5] * m[8];
  double d31 = m[0] * m[9] * m[14] + m[8] * m[13] * m[2] + m[12] * m[1] * m[10] - m[2] * m[9] * m[12] - m[10] * m[13] * m[0] - m[14] * m[1] * m[8];
  double d32 = m[0] * m[5] * m[14] + m[4] * m[13] * m[2] + m[12] * m[1] * m[6] - m[2] * m[5] * m[12] - m[6] * m[13] * m[0] - m[14] * m[1] * m[4];
  double d33 = m[0] * m[5] * m[10] + m[4] * m[9] * m[2] + m[8] * m[1] * m[6] - m[2] * m[5] * m[8] - m[6] * m[9] * m[0] - m[10] * m[1] * m[4];

  // Inverse of determinant
  d = 1.0 / d;

  // Create inverse
  Matrix4 r;

  r[0] = d00 * d;  r[4] = -d10 * d; r[8] = d20 * d; r[12] = -d30 * d;
  r[1] = -d01 * d; r[5] = d11 * d; r[9] = -d21 * d; r[13] = d31 * d;
  r[2] = d02 * d; r[6] = -d12 * d; r[10] = d22 * d; r[14] = -d32 * d;
  r[3] = -d03 * d; r[7] = d13 * d; r[11] = -d23 * d; r[15] = d33 * d;
  return r;
  //   return Matrix4(d00, -d01, d02, -d03, -d10, d11, -d12, d13, d20, -d21, d22, -d23, -d30, d31, -d32, d33)*d;
}

/*!
\brief Destructive addition.
*/
Matrix4& Matrix4::operator+= (const Matrix4& u)
{
  for (int i = 0; i < 16; i++)
  {
    r[i] += u.r[i];
  }
  return *this;
}

/*!
\brief Destructive subtraction.
*/
Matrix4& Matrix4::operator-= (const Matrix4& u)
{
  for (int i = 0; i < 16; i++)
  {
    r[i] -= u.r[i];
  }
  return *this;
}

/*!
\brief Destructive scalar multiply.
*/
Matrix4& Matrix4::operator*= (double a)
{
  for (int i = 0; i < 16; i++)
  {
    r[i] *= a;
  }
  return *this;
}

/*!
\brief Destructive scalar divide.
*/
Matrix4& Matrix4::operator/= (double a)
{
  const double ia = 1.0 / a;
  for (int i = 0; i < 16; i++)
  {
    r[i] *= ia;
  }
  return *this;
}

/*!
\brief Overloaded.
\param u,v Argument matrixes.
*/
Matrix4 operator+ (const Matrix4& u, const Matrix4& v)
{
  return Matrix4(
    u[0] + v[0], u[1] + v[1], u[2] + v[2], u[3] + v[3],
    u[4] + v[4], u[5] + v[5], u[6] + v[6], u[7] + v[7],
    u[8] + v[8], u[9] + v[9], u[10] + v[10], u[11] + v[11],
    u[12] + v[12], u[13] + v[13], u[14] + v[14], u[15] + v[15]);
}

/*!
\brief Overloaded.
\param u,v Argument matrixes.
*/
Matrix4 operator- (const Matrix4& u, const Matrix4& v)
{
  return Matrix4(
    u[0] - v[0], u[1] - v[1], u[2] - v[2], u[3] - v[3],
    u[4] - v[4], u[5] - v[5], u[6] - v[6], u[7] - v[7],
    u[8] - v[8], u[9] - v[9], u[10] - v[10], u[11] - v[11],
    u[12] - v[12], u[13] - v[13], u[14] - v[14], u[15] - v[15]);
}

/*!
\brief Overloaded.
\param s Stream.
\param matrix The matrix.
*/
std::ostream& operator<<(std::ostream& s, const Matrix4& matrix)
{
  s << "Matrix4(";
  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      s << matrix(i, j);
      if (i * 4 + j != 15)
      {
        s << ',';
      }
    }
  }
  s << ')';
  return s;
}

/*!
\brief Return QString of the form 'mat4(m(0, 0), ..., m(3, 3))' that write the matrix in GLSL.
\author Hubert-Brierre Pierre

\param m %Matrix.
*/
QString ToGLSL(const Matrix4& m)
{
  QString mat_str = "mat4(" + QString::number(float(m[0]));
  for (int i = 1; i < 16; i++)
  {
    mat_str += ", " + QString::number(float(m[i]));
  }
  mat_str += ")";
  return mat_str;
}

/*!
\brief Create a viewing transformation matrix.
\param eye Eye location.
\param look Look at position.
\param twist Twisting angle from vertical.
*/
Matrix4 Matrix4::LookAt(const Vector& eye, const Vector& look, const double& twist)
{
  // Translate to eye
  Matrix4 torg = Matrix4::Translate(-eye);

  // Y axis rotation
  Matrix4 roty = Matrix4::Identity;

  Vector d = look - eye;
  double denom1 = sqrt((d[0] * d[0]) + (d[2] * d[2]));
  if (denom1 != 0.0)
  {
    double si = d[0] / denom1;
    double co = -d[2] / denom1;
    roty(0, 0) = co;
    roty(0, 2) = -si;
    roty(2, 0) = si;
    roty(2, 2) = co;
  }

  // X axis rotation
  Matrix4 rotx = Matrix4::Identity;
  double denom2 = sqrt((denom1 * denom1) + (d[1] * d[1]));
  if (denom2 != 0.0)
  {
    double si = -d[1] / denom2;
    double co = denom1 / denom2;
    rotx(1, 1) = co;
    rotx(1, 2) = si;
    rotx(2, 1) = -si;
    rotx(2, 2) = co;
  }

  // Z axis rotation
  Matrix4 rotz = Matrix4::Rotation(Vector(0.0, 0.0, -twist));

  // Concatenate matrices
  Matrix4 total = ((torg * roty) * rotx) * rotz;
  return total;
}

/*!
\brief Create a rotation matrix about an arbitrary axis.
\param v %Axis (should be normalized).
\param a Angle.
*/
Matrix4 Matrix4::Rotation(const Vector& v, const double& a)
{
  return Matrix4(Matrix::Rotation(v, a));
}

/*!
\brief Create a rotation matrix about the orthogonal axes.
\param u Vector of angles in radian.
*/
Matrix4 Matrix4::Rotation(const Vector& u)
{
  return Matrix4(Matrix::Rotation(u));
}

/*!
\brief Creates a scaling matrix.
\param u Scaling vector.
*/
Matrix4 Matrix4::Scale(const Vector& u)
{
  return Matrix4(u);
}

/*!
\brief Creates a translation matrix.
\param t Translation vector.
*/
Matrix4 Matrix4::Translate(const Vector& t)
{
  return Matrix4(Matrix::Identity, t);
}

/*!
\brief Creates a shear matrix.
\param s Shear vector.
*/
Matrix4 Matrix4::Shear(const Vector& s)
{
  return Matrix4(Matrix::Identity, Vector::Null, s);
}

/*!
\brief Compute the determinant of the matrix.

\sa Matrix::Determinant().
*/
double Matrix4::Determinant() const
{
  const Matrix4& M = *this;
  return M(0, 0) * Matrix(M(1, 1), M(1, 2), M(1, 3), M(2, 1), M(2, 2), M(2, 3), M(3, 1), M(3, 2), M(3, 3)).Determinant()
    - M(1, 0) * Matrix(M(0, 1), M(0, 2), M(0, 3), M(2, 1), M(2, 2), M(2, 3), M(3, 1), M(3, 2), M(3, 3)).Determinant()
    + M(2, 0) * Matrix(M(0, 1), M(0, 2), M(0, 3), M(1, 1), M(1, 2), M(1, 3), M(3, 1), M(3, 2), M(3, 3)).Determinant()
    - M(3, 0) * Matrix(M(0, 1), M(0, 2), M(0, 3), M(1, 1), M(1, 2), M(1, 3), M(2, 1), M(2, 2), M(2, 3)).Determinant();
}
