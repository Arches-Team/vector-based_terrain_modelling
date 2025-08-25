// Noise

#pragma once


#include "libs/evector.h"

class RandomFast;

class InfiniteArray
{
protected:
  unsigned int o = 0; //!< Random offset.
protected:
  InfiniteArray();
  void Square(const Vector2&, int&, int&) const;
  void Square(const Vector2&, int&, int&, double&, double&) const;
  void Square(const Vector&, int&, int&, int&) const;
  void Square(const Vector&, int&, int&, int&, double&, double&, double&) const;
public:
  Vector2 CellVertex(int, int) const;
  Vector CellVertex(int, int, int) const;
protected:
  static RandomFast random; //!< Fast random number generator
};

/*!
\brief Compute the cubic cell coordinates.
\param p Point.
\param x,y Returned integer coordinates.

\sa InfiniteArray::Square(const Vector2& , int& , int& , double& , double& ) const
*/
inline void InfiniteArray::Square(const Vector2& p, int& x, int& y) const
{
  double int_x = floor(p[0]);
  double int_y = floor(p[1]);

  // Cell
  x = int(int_x);
  y = int(int_y);
}

/*!
\brief Compute the local coordinates inside a cubic cell.
\param p Point.
\param x,y Returned integer coordinates.
\param u,v Returned local coordinates.
*/
inline void InfiniteArray::Square(const Vector2& p, int& x, int& y, double& u, double& v) const
{
  double int_x = floor(p[0]);
  double int_y = floor(p[1]);

  u = p[0] - int_x;
  v = p[1] - int_y;

  // Cell
  x = int(int_x);
  y = int(int_y);
}

/*!
\brief Compute the cubic cell coordinates.
\param p Point.
\param x,y,z Returned integer coordinates.
*/
inline void InfiniteArray::Square(const Vector& p, int& x, int& y, int& z) const
{
  double int_x = floor(p[0]);
  double int_y = floor(p[1]);
  double int_z = floor(p[2]);

  // Cell
  x = int(int_x);
  y = int(int_y);
  z = int(int_z);
}

/*!
\brief Compute the local coordinates inside a cubic cell.
\param p Point.
\param x,y,z Returned integer coordinates.
\param u,v,w Returned local coordinates.
*/
inline void InfiniteArray::Square(const Vector& p, int& x, int& y, int& z, double& u, double& v, double& w) const
{
  double int_x = floor(p[0]);
  double int_y = floor(p[1]);
  double int_z = floor(p[2]);

  u = p[0] - int_x;
  v = p[1] - int_y;
  w = p[2] - int_z;

  // Cell
  x = int(int_x);
  y = int(int_y);
  z = int(int_z);
}

// Cellular noise
class CellularNoise2 :protected InfiniteArray
{
protected:
public:
  CellularNoise2();
  double Value(const Vector2&) const;
};

class CellularNoise :protected InfiniteArray
{
protected:
public:
  CellularNoise();
  double Value(const Vector&) const;
};


// Cellular noise
class VoroNoise :protected InfiniteArray
{
protected:
public:
  VoroNoise();
  double Value(const Vector2&) const;
};

