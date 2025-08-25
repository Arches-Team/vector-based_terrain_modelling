// Matrix 

#include "libs/matrix.h"
#include "libs/quadric.h"

/*!
\class Matrix2 matrix.h
\brief This class implements 2<SUP>2</SUP> matrix.

Operators have been overloaded to behave as expected. Thus it is possible to write:
\code
Matrix2 a = Matrix2(Vector2(1.0, 2.0));  // Scale
Matrix2 b = Matrix2::Rotation(Math::Pi / 4.0); // Rotation matrix
Vector2 u(2.0, -1.0); // Vector
Vector2 v = (a*b)*u;  // Matrix vector product
\endcode

\ingroup MathGroup
*/

double Matrix2::epsilon = 1.0e-8;

const Matrix2 Matrix2::RotationHalfPi(0.0, 1.0, -1.0, 0.0);
const Matrix2 Matrix2::RotationTwoThirdsPi(0.5, Math::Sqrt3/2.0, -Math::Sqrt3 / 2.0, 0.5);

/*!
\brief This static member defines the null matrix.
*/
const Matrix2 Matrix2::Null(0.0);

/*!
\brief This static member defines the identity matrix.
*/
const Matrix2 Matrix2::Identity(1.0);

/*!
\brief Computes the inverse of a matrix A<SUP>-1</SUP>.

This function returns the null matrix if A cannot be inversed.

The threshold value involved in the singular matrix detection
is set to 10<SUP>-18</SUP>.
*/
Matrix2 Matrix2::Inverse() const
{
  double e = Determinant();

  if (fabs(e) < 1.0e-18)
  {
    return Matrix2::Null;
  }

  return Matrix2(-r[3], -r[1], r[2], r[0]) / e;
}

/*!
\brief Multiplication.

\param u, v Input matrices.
*/
Matrix2 operator* (const Matrix2& u, const Matrix2& v)
{
  return Matrix2(u[0] * v[0] + u[2] * v[1], u[1] * v[0] + u[3] * v[1], u[0] * v[2] + u[2] * v[3], u[1] * v[2] + u[3] * v[3]);
}

/*!
\brief Destructive multiplication.
*/
Matrix2& Matrix2::operator*=(const Matrix2& M)
{
  Matrix2 r = (*this) * M;
  *this = r;
  return *this;
}

/*!
\brief Computes the inverse matrix and scales it by a double.
\param x %Real.
\param a %Matrix.
*/
Matrix2 operator/(const double& x, const Matrix2& a)
{
  return x * a.Inverse();
}

/*!
\brief Create a rotation matrix.

\param a Angle in radian.
*/
Matrix2 Matrix2::Rotation(const double& a)
{
  double c = cos(a);
  double s = sin(a);
  return Matrix2(c, s, -s, c);
}

/*!
\brief Create a random rotation matrix.

This is a simple function which is equivalent to the following inlined code:
\code
RandomFast r;
Matrix2 m=Matrix2::Rotation(2.0*Math::Pi*r.Uniform())
\endcode

\param ra %Random number generator.
*/
Matrix2 Matrix2::Rotation(RandomFast& ra)
{
  return Matrix2::Rotation(Math::TwoPi * ra.Uniform());
}

/*!
\brief Create a tensor matrix.

\param c,s Cosine and sine of the angle of the tensor.
*/
Matrix2 Matrix2::Tensor2(const double& c, const double& s)
{
  return Matrix2(c, s, s, -c);
}

/*!
\brief Create an anisotropic tensor matrix.

\param angle The angle of the tensor.
*/
Matrix2 Matrix2::Tensor2(const double& angle)
{
  double c = cos(2.0 * angle);
  double s = sin(2.0 * angle);
  return Matrix2(c, s, s, -c);
}

/*!
\brief Destructive addition operator.
\param u %Matrix.
*/
Matrix2& Matrix2::operator+=(const Matrix2& u)
{
  for (int i = 0; i < 4; i++)
  {
    r[i] += u.r[i];
  }
  return *this;
}

/*!
\brief Destructive subtraction operator.
\param u %Matrix.
*/
Matrix2& Matrix2::operator-=(const Matrix2& u)
{
  for (int i = 0; i < 4; i++)
  {
    r[i] -= u.r[i];
  }
  return *this;
}

/*!
\brief Destructive multiplication operator.
*/
Matrix2& Matrix2::operator*=(double a)
{
  for (int i = 0; i < 4; i++)
  {
    r[i] *= a;
  }
  return *this;
}
/*!
\brief Hadamart product, which is the element-wise product.

\paramu %Matrix.
*/
Matrix2 Matrix2::Scaled(const Matrix2& u) const
{
  return Matrix2(r[0] * u[0], r[1] * u[1], r[2] * u[2], r[3] * u[3]);
}

/*!
\brief Destructive division operator.
*/
Matrix2& Matrix2::operator/=(double a)
{
  const double ia = 1.0 / a;
  for (int i = 0; i < 4; i++)
  {
    r[i] *= ia;
  }
  return *this;
}

/*!
\brief Overloaded.
\param s Stream.
\param matrix The matrix.
*/
std::ostream& operator<<(std::ostream& s, const Matrix2& matrix)
{
  s << "Matrix2(" << matrix[0] << ',' << matrix[1] << ',' << matrix[2] << ',' << matrix[3] << ')';

  return s;
}

/*!
\brief Return QString of the form 'mat2(m(0, 0), ..., m(1, 1))' that write the matrix in GLSL.
\author Hubert-Brierre Pierre

\param m %Matrix.
*/
QString ToGLSL(const Matrix2& m)
{
  QString mat_str = "mat2(" + QString::number(float(m[0]));
  for (int i = 1; i < 4; i++)
    mat_str += ", " + QString::number(float(m[i]));
  mat_str += ")";
  return mat_str;
}


/*!
\brief Compute the characteristic polynomial.
*/
Quadric Matrix2::Characteristic() const
{
  return Quadric(1.0, -Trace(), Determinant());
}

/*!
\brief Compute the spectral norm of the matrix.
*/
double Matrix2::SpectralNorm() const
{
  const Matrix2& m = *this;

  Matrix2 A = m.T() * m;

  Quadric quadric(1.0, -A.Trace(), A.Determinant());

  double r[2];

#pragma warning(push)
#pragma warning(disable: 4189)  
  int n = quadric.Solve(r[0], r[1]);
#pragma warning(pop)
  double s = Math::Max(r[0], r[1]);

  return sqrt(s);
}

/*!
\brief Compute the infinity norm of the matrix.
*/
double Matrix2::InfinityNorm() const
{
  return Math::Max(fabs(r[0]), fabs(r[1]), fabs(r[2]), fabs(r[3]));
}

/*!
\brief Compute the Frobenius norm of the matrix.
*/
double Matrix2::FrobeniusNorm() const
{
  return sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2] + r[3] * r[3]);
}

/*!
\brief Compute the eigen values.
*/
Vector2 Matrix2::Eigen() const
{
  Vector2 a;
  Quadric q(1.0, -Trace(), Determinant());
  q.Solve(a[0], a[1]);

  return a;
}

/*!
\brief Compute the eigen values of the matrix.

Return eigen vectors as well.

\param a An array of two eigen values.
\param e The array containing the two eigen vectors.
*/
void Matrix2::EigenSolveSymmetric(double a[2], Vector2 e[2]) const
{
  Quadric q(1.0, -Trace(), Determinant());
  q.Solve(a[0], a[1]);

  if (r[1] != 0.0)
  {
    e[0] = Vector2(a[0] - r[3], r[1]);
    e[1] = Vector2(a[1] - r[3], r[1]);
  }
  else
  {
    /*
    e[0] = Vector2(r[2], a[0] - r[0]);
    e[1] = Vector2(r[2], a[1] - r[0]);
    */
    // Case where r[2] and r[1] are 0.0 (symmetric matrix) thus : identity
    e[0] = Vector2::X;
    e[1] = Vector2::Y;
  }
}
