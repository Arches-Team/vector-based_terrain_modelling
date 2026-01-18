// Maya  

#pragma once

#include "libs/evector.h"

class VectorFloat
{
public:
  float x, y, z; //!< Components.
public:
  explicit VectorFloat();
  explicit VectorFloat(const Vector&);
  explicit VectorFloat(float, float, float);
};

//! Create a vector.
inline VectorFloat::VectorFloat()
{
  x = y = z = 0.0f;
}

//! Create a vector.
inline VectorFloat::VectorFloat(const Vector& v) :x(v[0]), y(v[1]), z(v[2])
{
}

//! Create a vector.
inline VectorFloat::VectorFloat(float x, float y, float z) : x(x), y(y), z(z)
{
}

class Vector2Float
{
public:
  float x, y;
public:
  explicit Vector2Float();
  explicit Vector2Float(const Vector2&);
};

//! Create a vector.
inline Vector2Float::Vector2Float()
{
  x = y = 0.0f;
}

//! Create a vector.
inline Vector2Float::Vector2Float(const Vector2& v) :x(v[0]), y(v[1])
{
}


class Matrix4Float
{
public:
  float r[16] = { 1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,1.0,0.0,0.0,1.0 }; //!< Coefficients, set to identity matrix.
public:
  explicit Matrix4Float();
  explicit Matrix4Float(float);
  explicit Matrix4Float(float[16]);

  float& operator()(int, int);
  float operator()(int, int) const;

  friend Matrix4Float operator* (const Matrix4Float&, const Matrix4Float&);
  friend bool operator!= (const Matrix4Float&, const Matrix4Float&);

  static Matrix4Float LookAt(const Vector& eye, const Vector& at, const Vector& up);
  static Matrix4Float Perspective(float fovy, float zNear, float zFar, float width, float height);
  static Matrix4Float Orthographic(float left, float right, float bottom, float top, float zNear, float zFar);
  static Matrix4Float inverse(const Matrix4Float&);

public:
  static Matrix4Float Identity;
};

//! Create a matrix.
inline Matrix4Float::Matrix4Float()
{
}

//! Create a diagonal matrix.
inline Matrix4Float::Matrix4Float(float x)
{
	for (int i = 0; i < 16; i++)
		r[i] = 0.0f;
	r[0] = r[5] = r[10] = r[15] = x;
}

//! Create a matrix from coeficients.
inline Matrix4Float::Matrix4Float(float x[16])
{
	for (int i = 0; i < 16; i++)
		r[i] = x[i];
}

//! Accessor.
inline float& Matrix4Float::operator()(int i, int j)
{
	return r[i * 4 + j];
}

inline float Matrix4Float::operator()(int i, int j) const
{
	return r[i * 4 + j];
}

//! Matrix multiplication.
inline Matrix4Float operator*(const Matrix4Float& u, const Matrix4Float& v)
{
	Matrix4Float a;
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

//! Matrix comparison
inline bool operator!=(const Matrix4Float& u, const Matrix4Float& v)
{
	for (int i = 0; i < 16; i++)
		if (u.r[i] != v.r[i])
			return true;
	return false;

}
