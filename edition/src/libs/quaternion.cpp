// Quaternion

#include "libs/quaternion.h"

const double Quaternion::Epsilon = 1e-03;

const Quaternion Quaternion::Null(0, 0, 0, 0);
const Quaternion Quaternion::Identity(1, 0, 0, 0);

/*!
\class Quaternion quaternion.h
\brief Quaternions.

\ingroup MathGroup

*/

/*!
\brief Create a quaternion from an axis vector and a rotation angle.

Axis should be unit length.

\param axis Normalized axis.
\param a Angle.
*/
Quaternion::Quaternion(const Vector& axis, const double& a)
{
  double s = sin(0.5 * a);
  w = cos(0.5 * a);
  x = s * axis[0];
  y = s * axis[1];
  z = s * axis[2];
}

/*!
\brief Create a quaternion from two (normalized) vectors representing the rotation.

Simply coded as:
\code
Vector a,b; // Should be unit
Quaternion q(a/b,acos(a*b));
\endcode
\param a,b %Vectors (should be unit).
*/
Quaternion::Quaternion(const Vector& a, const Vector& b)
{
  Vector c = a / b;
  x = c[0];
  y = c[1];
  z = c[2];
  w = acos(Clamp(a * b, -1.0, 1.0));
}

/*!
\brief Create a quaternion from Euler rotation angles.

Using this function is simpler than calling the constructor from a rotation matrix.

\param a,b,c Rotation angles around x, y and z axes (should be in radian).
*/
Quaternion Quaternion::FromAngles(const double& a, const double& b, const double& c)
{
  Quaternion qa(cos(a / 2.0), sin(a / 2.0), 0.0, 0.0);
  Quaternion qb(cos(b / 2.0), 0.0, sin(b / 2.0), 0.0);
  Quaternion qc(cos(c / 2.0), 0.0, 0.0, sin(c / 2.0));

  qa = qa * (1.0 / Norm(qa));
  qb = qb * (1.0 / Norm(qb));
  qc = qc * (1.0 / Norm(qc));

  return qc * qb * qa;
}

/*!
\brief Create a quaternion from a rotation matrix.

The algorithm can be found in Ken Shoemake's article in
1987 Siggraph course notes Quaternion Calculus and Fast Animation.
\param R Rotation matrix.
*/
Quaternion::Quaternion(const Matrix& R)
{
  double trace = R.Trace();

  if (trace > 0.0)
  {
    // |w|>1/2, may as well choose w>1/2
    double root = sqrt(trace + 1.0);  // 2w
    w = 0.5 * root;
    root = 0.5 / root;  // 1/(4w)
    x = (R[7] - R[5]) * root;
    y = (R[2] - R[6]) * root;
    z = (R[3] - R[1]) * root;
  }
  else
  {
    // |w|<=1/2
    static const int next[3] = { 1, 2, 0 };
    int i = 0;
    if (R[4] > R[0])
    {
      i = 1;
    }
    if (R[8] > R[3 * i + i])
    {
      i = 2;
    }
    int j = next[i];
    int k = next[j];

    double root = sqrt(R[3 * i + i] - R[3 * j + j] - R[3 * k + k] + 1.0);
    double* quat[3] = { &x, &y, &z };
    *quat[i] = 0.5 * root;
    root = 0.5 / root;
    w = (R[3 * k + j] - R[3 * j + k]) * root;
    *quat[j] = (R[3 * j + i] + R[3 * i + j]) * root;
    *quat[k] = (R[3 * k + i] + R[3 * i + k]) * root;
  }
}

/*!
\brief Transforms a quaternion into a rotation matrix.
*/
Matrix Quaternion::RotationMatrix() const
{
  double tx = 2.0 * x;
  double ty = 2.0 * y;
  double tz = 2.0 * z;
  double twx = tx * w;
  double twy = ty * w;
  double twz = tz * w;
  double txx = tx * x;
  double txy = ty * x;
  double txz = tz * x;
  double tyy = ty * y;
  double tyz = tz * y;
  double tzz = tz * z;

  return Matrix(
    1.0 - (tyy + tzz), txy - twz, txz + twy,
    txy + twz, 1.0 - (txx + tzz), tyz - twx,
    txz - twy, tyz + twx, 1.0 - (txx + tyy));
}

/*!
\brief Transform a quaternion into an rotation defined by its axis
and the angle of rotation angle.

The quaternion representing the rotation is : q=cos(a/2)+sin(a/2)*(x*i+y*j+z*k).

\param angle The angle.
\param axis Vector axis.
*/
void Quaternion::AngleAxis(double& angle, Vector& axis) const
{
  double n = x * x + y * y + z * z;
  if (n > 0.0)
  {
    angle = 2.0 * acos(w);
    double invlen = 1.0 / sqrt(n);
    axis = Vector(x, y, z) * invlen;
  }
  else
  {
    // Angle is 0 (mod 2*pi), so any axis will do
    angle = 0.0;
    axis = Vector::X;
  }
}

/*!
\brief Overloaded.
*/
Quaternion Quaternion::operator+ (const Quaternion& q) const
{
  return Quaternion(w + q.w, x + q.x, y + q.y, z + q.z);
}

/*!
\brief Overloaded.
*/
Quaternion Quaternion::operator- (const Quaternion& q) const
{
  return Quaternion(w - q.w, x - q.x, y - q.y, z - q.z);
}

/*!
\brief Multiplication of two quaternion.

Note that multiplication is not generally commutative.
*/
Quaternion Quaternion::operator* (const Quaternion& q) const
{
  return Quaternion(w * q.w - x * q.x - y * q.y - z * q.z,
    w * q.x + x * q.w + y * q.z - z * q.y,
    w * q.y + y * q.w + z * q.x - x * q.z,
    w * q.z + z * q.w + x * q.y - y * q.x);
}

/*!
\brief Overloaded.
*/
Quaternion Quaternion::operator- () const
{
  return Quaternion(-w, -x, -y, -z);
}

/*!
\brief Squared length of a quaternion.
\param q %Quaternion.
*/
double Norm(const Quaternion& q)
{
  return q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z;
}

/*!
\brief Compute the inverse of a quaternion.

This should be applied
to non-zero quaternions only. If quaternion happens to be null, then
the function return also null wich is an invalid result that flags
the error.
*/
Quaternion Quaternion::Inverse() const
{
  double n = w * w + x * x + y * y + z * z;
  if (n > 0.0)
  {
    n = 1.0 / n;
    return Quaternion(w * n, -x * n, -y * n, -z * n);
  }
  else
  {
    // Return an invalid result to flag the error
    return Null;
  }
}

/*!
\brief Compute the exponential of a quaternion.

If q=a*(x,y,z) where (x,y,z) is unit length, then
exp(q)=cos(a)+sin(a)*(x,y,z). If sin(a) is near zero,
use exp(a)=cos(a)+a*(x,y,z).
*/
Quaternion Quaternion::Exp() const
{
  double angle = sqrt(x * x + y * y + z * z);
  double sn = sin(angle);

  if (fabs(sn) >= Epsilon)
  {
    double c = sn / angle;
    return Quaternion(cos(angle), c * x, c * y, c * z);
  }
  else
  {
    return Quaternion(cos(angle), x, y, z);
  }
}

/*!
\brief Computes the logarithm of a quaternion.

If quaternion is of unit length, then q=cos(a)+sin(a).v, then log(q)=a.v.
*/
Quaternion Quaternion::Log() const
{
  if (fabs(w) < 1.0)
  {
    double angle = acos(w);
    double sn = sin(angle);
    if (fabs(sn) >= Epsilon)
    {
      double c = angle / sn;
      return Quaternion(0.0, c * x, c * y, c * z);
    }
  }
  return Quaternion(0.0, x, y, z);
}

/*!
\brief Compute the product between a quaternion and a vector.

This member calls Quaternion::RotationMatrix().
As such, the quaternion representation of a rotation
matrix requires less space than the matrix and more time to compute
the rotated vector. This is a typical space-time tradeoff.
*/
Vector Quaternion::operator* (const Vector& p) const
{
  return RotationMatrix() * p;
}

/*!
\brief Perform quaternion spherical interpolation, also known as slerping.

The dot product of argument quaternions should be positive.
\param p, q Argument quaternions.
\param t Interpolating value.
*/
Quaternion Quaternion::Lerp(const double& t, const Quaternion& p, const Quaternion& q)
{
  // Dot product
  double cs = p.w * q.w + p.x * q.x + p.y * q.y + p.z * q.z;

  double sn = sqrt(fabs(1.0 - cs * cs));
  if (fabs(sn) < Epsilon)
    return p;

  double angle = atan2(sn, cs);
  double invSn = 1.0 / sn;
  double c0 = sin((1.0 - t) * angle) * invSn;
  double c1 = sin(t * angle) * invSn;

  return c0 * p + c1 * q;
}

/*!
\brief Overloaded.
\param s Stream.
\param quaternion The quaternion.
*/
std::ostream& operator<<(std::ostream& s, const Quaternion& quaternion)
{
  s << "Quaternion (w: " << quaternion.w << ", x: " << quaternion.x << ", y: " << quaternion.y << ", z: " << quaternion.z << ")";
  return s;
}