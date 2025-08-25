// Matrix 

#include "libs/matrix.h"
#include "libs/random.h"
#include "libs/quaternion.h"

/*!
\class Matrix matrix.h
\brief This class implements 3<SUP>2</SUP> matrix.

Operators have been overloaded so as to behave as expected
with the Vector class.

\image html matrix3.png

Many functions have been coded inline out of efficiency.

The constructors of this class are used to create matrices
in the general case.
See the static member functions Matrix::Rotation() to create
different kinds of rotation matrices. Canonical symmetry matrices should
be created as follows:
\code
Matrix planar=Matrix(Vector(1.0,-1.0,1.0)); // Symmetry, plane y=0
Matrix origin=Matrix(Vector(-1.0,-1.0,-1.0)); // Symmetry around origin
\endcode
More general symmetries can be created using static member functions:
\code
Matrix planar=Matrix::Symmetry(Vector(2.0,-3.0,1.0)); // Symmetry, plane with normal vector as argument
\endcode

Components are stored in a single dimension array, sorting components
by column.

\ingroup MathGroup
*/

double Matrix::epsilon = 1.0e-8;

/*!
\brief Null matrix.
*/
const Matrix Matrix::Null(0.0);

/*!
\brief Identity matrix.
*/
const Matrix Matrix::Identity(1.0);

const Matrix Matrix::RotationHalfPiZ(0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0);

/*!
\brief Compute the inverse of a matrix A<SUP>-1</SUP>.

Recall that A<SUP>-1</SUP> can be defined as the transposed adjoint matrix divided by the determinant.

This function returns the null matrix if A cannot be inversed.

The threshold value involved in the singular matrix detection
is set to 10<SUP>-18</SUP>.

\param A %Matrix.
*/
Matrix Inverse(const Matrix& A)
{
  double e = A.Determinant();

  if (fabs(e) < 1.0e-18)
  {
    return Matrix::Null;
  }

  return A.Adjoint().T() / e;
}

/*!
\brief Multiplication.

\param u, v Input matrices.
*/
Matrix operator*(const Matrix& u, const Matrix& v)
{
  return Matrix(u[0] * v[0] + u[3] * v[1] + u[6] * v[2], u[1] * v[0] + u[4] * v[1] + u[7] * v[2], u[2] * v[0] + u[5] * v[1] + u[8] * v[2],
    u[0] * v[3] + u[3] * v[4] + u[6] * v[5], u[1] * v[3] + u[4] * v[4] + u[7] * v[5], u[2] * v[3] + u[5] * v[4] + u[8] * v[5],
    u[0] * v[6] + u[3] * v[7] + u[6] * v[8], u[1] * v[6] + u[4] * v[7] + u[7] * v[8], u[2] * v[6] + u[5] * v[7] + u[8] * v[8]);
}

/*!
\brief Destructive multiplication.
*/
Matrix& Matrix::operator*=(const Matrix& M)
{
  Matrix r = (*this) * M;
  *this = r;
  return *this;
}

/*!
\brief Computes the inverse matrix and scales it by a double.

This function calls Matrix::Inverse().
\param x %Real.
\param a %Matrix.
*/
Matrix operator/(const double& x, const Matrix& a)
{
  return x * Inverse(a);
}

/*!
\brief Create a rotation matrix given a
vector of angles that specifies the rotation
around each world coordinate axis.

Rotations are concatenated, starting by rotating
around the z axis, then y and eventually x.

\param v Vector of angles in radian.
*/
Matrix Matrix::Rotation(const Vector& v)
{
  return Matrix::RotationX(v[0]) * Matrix::RotationY(v[1]) * Matrix::RotationZ(v[2]);
}

/*!
\brief Create a rotation matrix given a
vector of angles that specifies the rotation
around each world coordinate axis.

Rotations are concatenated, starting by rotating
around the z axis, then y and eventually x.

\param v Vector of angles in radian.
*/
Matrix Matrix::RotationMaya(const Vector& v)
{
  return Matrix::RotationZ(v[2]) * Matrix::RotationY(v[1]) * Matrix::RotationX(v[0]);
}

/*!
\brief Create a rotation matrix around the x-axis.

This static member function is provided out of efficiency as
it is much faster than any other general Matrix::Rotation() member.

\param a Angle (in radian).
*/
Matrix Matrix::RotationX(const double& a)
{
  const double c = cos(a);
  const double s = sin(a);
  return Matrix(1.0, 0.0, 0.0, 0.0, c, s, 0.0, -s, c);
}

/*!
\brief Create a rotation matrix around the y-axis.

\sa Matrix::RotationX(const double&)

\param a Angle (in radian).
*/
Matrix Matrix::RotationY(const double& a)
{
  const double c = cos(a);
  const double s = sin(a);
  return Matrix(c, 0.0, -s, 0.0, 1.0, 0.0, s, 0.0, c);
}

/*!
\brief Create a rotation matrix around the z-axis.

\sa Matrix::RotationX(const double&)

\param a Angle (in radian).
*/
Matrix Matrix::RotationZ(const double& a)
{
  const double c = cos(a);
  const double s = sin(a);
  return Matrix(c, s, 0.0, -s, c, 0.0, 0.0, 0.0, 1.0);
}

/*!
\brief Create a planar symmetry matrix.

\param n Normal to the plane.
*/
Matrix Matrix::Symmetry(const Vector& n)
{
  Matrix r = Matrix::Rotation(Normalized(n), Vector::Z);
  return r * Matrix(Vector(1.0, 1.0, -1.0)) * r.T();
}

/*!
\brief Create a rotation matrix that rotates a normalized
vector into another one.

\param a, b Initial and final vector (should be normalized).
*/
Matrix Matrix::Rotation(const Vector& a, const Vector& b)
{
  Matrix matrix;

  Vector v = a / b;
  double e = a * b;

  // Almost identical vectors
  if (e > 1.0 - epsilon)
  {
    return Matrix::Identity;
  }
  // Almost opposite vectors
  else if (e < epsilon - 1.0)
  {
    Vector left(0.0, a[2], -a[1]);
    if (left * left < epsilon)
    {
      left[0] = -a[2]; left[1] = 0.0; left[2] = a[0];
    }

    Normalize(left);
    Vector up = left / a;
    // We have a coordinate system, i.e. a basis
    // M = (a, up, left), and we want to rotate to:    
    // N = (-a, up, -left). This is done with the matrix:
    // N*M^T where M^T is the transpose of M          
    double fxx = -a[0] * a[0]; double fyy = -a[1] * a[1]; double fzz = -a[2] * a[2];
    double fxy = -a[0] * a[1]; double fxz = -a[0] * a[2]; double fyz = -a[1] * a[2];

    double uxx = up[0] * up[0]; double uyy = up[1] * up[1]; double uzz = up[2] * up[2];
    double uxy = up[0] * up[1]; double uxz = up[0] * up[2]; double uyz = up[1] * up[2];

    double lxx = -left[0] * left[0]; double lyy = -left[1] * left[1]; double lzz = -left[2] * left[2];
    double lxy = -left[0] * left[1]; double lxz = -left[0] * left[2]; double lyz = -left[1] * left[2];

    // Symmetric matrix
    matrix(0, 0) = fxx + uxx + lxx; matrix(0, 1) = fxy + uxy + lxy; matrix(0, 2) = fxz + uxz + lxz;
    matrix(1, 0) = matrix(0, 1); matrix(1, 1) = fyy + uyy + lyy; matrix(1, 2) = fyz + uyz + lyz;
    matrix(2, 0) = matrix(0, 2); matrix(2, 1) = matrix(1, 2); matrix(2, 2) = fzz + uzz + lzz;
  }
  else
  {
    double h = (1.0 - e) / (v * v);
    double hvx = h * v[0];
    double hvz = h * v[2];
    double hvxy = hvx * v[1];
    double hvxz = hvx * v[2];
    double hvyz = hvz * v[1];
    matrix(0, 0) = e + hvx * v[0]; matrix(0, 1) = hvxy - v[2];    matrix(0, 2) = hvxz + v[1];
    matrix(1, 0) = hvxy + v[2];  matrix(1, 1) = e + h * v[1] * v[1]; matrix(1, 2) = hvyz - v[0];
    matrix(2, 0) = hvxz - v[1];  matrix(2, 1) = hvyz + v[0];    matrix(2, 2) = e + hvz * v[2];
  }
  return matrix;
}

/*!
\brief Create a rotation matrix that rotates (0,0,1) to a normalized vector.

This is the same as:
\code
Vector v;
Matrix r=Matrix::Rotation(Vector::Z,v);
\endcode

\param v Final vector (should be normalized).
*/
Matrix Matrix::RotationCanonical(const Vector& v)
{
  return Matrix::Rotation(Vector::Z, v);
}

/*!
\brief Create a rotation matrix about an arbitrary axis.

\code
Matrix R=Matrix::Rotation(Normalized(Vector(1,2,-1)),Math::DegreeToRadian(30.0));
\endcode

\param u Rotation axis, which should be unit.
\param a %Angle, in radian.
*/
Matrix Matrix::Rotation(const Vector& u, const double& a)
{
  double c = cos(a);
  double s = sin(a);
  double ux = u[0];
  double uy = u[1];
  double uz = u[2];

  double ic = 1.0 - c;
  return Matrix(c + ux * ux * ic,
    uy * ux * ic + uz * s,
    uz * ux * ic - uy * s,
    ux * uy * ic - uz * s,
    c + uy * uy * ic,
    uz * uy * ic + ux * s,
    ux * uz * ic + uy * s,
    uy * uz * ic - ux * s,
    c + uz * uz * ic);
}

/*!
\brief Destructive addition operator.
*/
Matrix& Matrix::operator+=(const Matrix& u)
{
  for (int i = 0; i < 9; i++)
  {
    r[i] += u.r[i];
  }
  return *this;
}

/*!
\brief Destructive subtraction operator.
*/
Matrix& Matrix::operator-=(const Matrix& u)
{
  for (int i = 0; i < 9; i++)
  {
    r[i] -= u.r[i];
  }
  return *this;
}

/*!
\brief Destructive multiplication operator.
*/
Matrix& Matrix::operator*=(double a)
{
  for (int i = 0; i < 9; i++)
  {
    r[i] *= a;
  }
  return *this;
}

/*!
\brief Destructive division operator.
*/
Matrix& Matrix::operator/=(double a)
{
  const double ia = 1.0 / a;
  for (int i = 0; i < 9; i++)
  {
    r[i] *= ia;
  }
  return *this;
}

/*!
\brief Compute the absolute value of the matrix.
*/
Matrix Matrix::Abs() const
{
  return Matrix(fabs(r[0]), fabs(r[1]), fabs(r[2]), fabs(r[3]), fabs(r[4]), fabs(r[5]), fabs(r[6]), fabs(r[7]), fabs(r[8]));
}

/*!
\brief Compute the adjoint, i.e, the comatrix of the matrix.

\sa Inverse()
*/
Matrix Matrix::Adjoint() const
{
  return Matrix(
    r[4] * r[8] - r[7] * r[5],
    -(r[3] * r[8] - r[6] * r[5]),
    r[3] * r[7] - r[6] * r[4],

    -(r[1] * r[8] - r[7] * r[2]),
    r[0] * r[8] - r[6] * r[2],
    -(r[0] * r[7] - r[6] * r[1]),

    r[1] * r[5] - r[4] * r[2],
    -(r[0] * r[5] - r[3] * r[2]),
    r[0] * r[4] - r[3] * r[1]);
}

/*!
\brief Overloaded.
\param s Stream.
\param matrix The matrix.
*/
std::ostream& operator<<(std::ostream& s, const Matrix& matrix)
{
  s << "Matrix(";
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      s << matrix(i, j);
      if (i * 3 + j != 8)
      {
        s << ',';
      }
    }
  }
  s << ')';
  return s;
}

/*!
\brief Compute the covariance matrix of a point cloud.

This matrix is useful for computing the principal directions of a point of clouds.
The principal directions are defined as the eigen vectors of the matrix.

\sa EigenSolveSymmetric
\param p Array of vertices.
\param n Number of vertices.
*/
Matrix Matrix::Covariance(Vector* p, int n)
{
  Vector c(0.0);
  // Compute center of cloud
  for (int i = 0; i < n; i++)
  {
    c += p[i];
  }
  c /= n;

  Matrix C(0.0);
  for (int i = 0; i < n; i++)
  {
    Vector e = p[i] - c;
    C[0] += e[0] * e[0];
    C[1] += e[0] * e[1];
    C[2] += e[0] * e[2];

    C[4] += e[1] * e[1];
    C[5] += e[1] * e[2];

    C[8] += e[2] * e[2];
  }
  // Some coefficients are the same
  C[3] = C[1];
  C[6] = C[2];
  C[7] = C[5];
  C /= n;
  return C;
}

/*!
\brief Compute the angles of the rotation matrix.
*/
Vector Matrix::GetRotationAngles() const
{
  // Angles
  double x, y, z;

  const Matrix& m = *this;

  if (m(0, 2) < +1.0)
  {
    if (m(0, 2) > -1.0)
    {
      y = asin(m(0, 2));
      x = atan2(-m(1, 2), m(2, 2));
      z = atan2(-m(0, 1), m(0, 0));
    }
    else // r02 = -1
    {
      // Not a unique solution: thetaZ - thetaX = atan2(r10,r11)
      y = -Math::HalfPi;
      x = -atan2(m(1, 0), m(1, 1));
      z = 0.0;
    }
  }
  else // r02 = +1
  {
    // Not a unique solution: thetaZ + thetaX = atan2(r10,r11)
    y = +Math::HalfPi;
    x = atan2(m(1, 0), m(1, 1));
    z = 0.0;
  }

  return Vector(x, y, z);
}

/*!
\brief Spherical interpolation of two rotation matrices.

This is an expensive function that first converts rotations to quaternions,
then interpolates the quaternions, before finaly reconstructing the corresponding
matrix.

\sa Quaternion Lerp(const double&, const Quaternion&, const Quaternion&);
\param t Interpolation parameter.
\param a,b Rotation matrix.
*/
Matrix Matrix::Lerp(const double& t, const Matrix& a, const Matrix& b)
{
  Quaternion qa = Quaternion(a);
  Quaternion qb = Quaternion(b);
  Quaternion q = Quaternion::Lerp(t, qa, qb);
  Matrix m = q.RotationMatrix();
  return m;
}

/*!
\brief Compute the Frenet frame.

Given tangent <B>t</B> and up vector <B>u</B>, computes the matrix (<B>t</B>,<B>t</B>^<B>u</B>,<B>u</B>)).
This function stores the tangent, the binormal and the normal vectors in the matrix.
The normalization of the vectors is performed internally.

\param t Tangent, which need not be normalized.
\param u Up vector.
*/
Matrix Matrix::Frenet(const Vector& t, const Vector& u)
{
  return Matrix(Normalized(t), Normalized(t / u), Normalized(u));
}

/*!
\brief Generate uniformly distributed random rotation matrix.

See Fast Random Rotation Matrices, J. Arvo. <I>Graphics Gems III</I>, 1992.

\param ra Fast random number generator.
*/
Matrix Matrix::Rotation(RandomFast& ra)
{
  double theta = ra.Uniform() * Math::TwoPi; // Rotation about the pole (Z).
  double phi = ra.Uniform() * Math::TwoPi; // For direction of pole deflection.
  double z = ra.Uniform() * 2.0;                 // For magnitude of pole deflection.

  // Compute a vector V used for distributing points over the sphere  via the reflection I - V t(V). 
  // This formulation of V will guarantee that if x[1] and x[2] are uniformly distributed, the reflected points will be uniform on the sphere.  
  // Note that V  has length sqrt(2) to eliminate the 2 in the Householder matrix. 
  double r = sqrt(z);
  double vx = sin(phi) * r;
  double vy = cos(phi) * r;
  double vz = sqrt(2.0 - z);

  // Compute the row vector S = t(V) * R, where R is a simple rotation by theta about the z-axis.  
  // No need to compute Sz since it's just Vz.                                                    
  double st = sin(theta);
  double ct = cos(theta);
  double sx = vx * ct - vy * st;
  double sy = vx * st + vy * ct;

  // Construct the rotation matrix : ( V t(V) - I ) R, 
  // which is equivalent to V S - R (last element == Vz * Vz - 1)
  return Matrix(
    vx * sx - ct, vx * sy - st, vx * vz,
    vy * sx + st, vy * sy - ct, vy * vz,
    vz * sx, vz * sy, 1.0 - z);
}

/*!
\brief Return QString of the form 'mat3(m(0, 0), ..., m(3, 3))' that write the matrix in GLSL.
\author Hubert-Brierre Pierre

\param m %Matrix.
*/
QString ToGLSL(const Matrix& m)
{
  QString mat_str = "mat3(" + QString::number(float(m[0]));
  for (int i = 1; i < 9; i++)
    mat_str += ", " + QString::number(float(m[i]));
  mat_str += ")";
  return mat_str;
}

