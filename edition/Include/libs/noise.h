// Noise

#pragma once

#include "libs/evector.h"
#include "libs/scalarfield.h"

class Noise :public AnalyticScalarField
{
protected:
  bool quintic = false; //!< Flag defining interpolation type: either cubic (false) or quintic (true).
public:
  Noise(bool = false);
  //! Empty
  ~Noise() {}

  double Value(const Vector&) const;
  double K() const;
  
  Vector Gradient(const Vector&) const;

  double AtGradient(const Vector&, Vector&) const;


  // Vector noise
  Vector AtVector(const Vector&) const;
protected:
  static double rtable[267]; //!< Random double table.
  static const short hashTable[4096]; //!< Hash table.
protected:
  // Hash code
  short Hash1d(int, int) const;
  short Hash2d(int, int) const;
  short Hash3d(int, int, int) const;

  double IncrSum(int, double, double, double, double) const;

  double IncrSumP(double*, double, double, double, double) const;

  double GradientVertex(double*, double, double, double) const;
  double GradientVertex(double*, double, double) const;
  double ValueVertex(double*) const;
  Vector GradientVertexGrad(double*) const;
};

class Noise2 :public AnalyticScalarField2
{
protected:
  bool quintic = false; //!< Flag defining interpolation type: either cubic (false) or quintic (true).
public:
  Noise2(bool = false);
  //! Empty
  ~Noise2() {}

  double Value(const Vector2&) const;
  double K() const;

protected:
  static double rtable[267]; //!< Random double table.
  static const short hashTable[4096]; //!< Hash table.
protected:
  // Hash code
  short Hash1d(int, int) const;
  short Hash2d(int, int) const;

  double IncrSum(int, double, double, double, double) const;
  double IncrSumP(double*, double, double, double, double) const;

  double GradientVertex(double*, double, double) const;
  double ValueVertex(double*) const;
};


// Simplex noise
class SimplexNoise :public AnalyticScalarField {
protected:
public:
  //! Empty.
  SimplexNoise() {}
  //! Empty.
  ~SimplexNoise() {}

  double Value(const Vector&) const;
  double K() const;

protected:
  static const Vector gradient[12]; //!< Array of gradients for 3D noise.
  static const int perm[512];    //!< Permutation table, 256 entries duplicated once to avoid modulo computations.
protected:
  static const double F3, G3; //!< Unskew factors for 3D case.
};


// Simplex noise
class SimplexNoise2 :public AnalyticScalarField2 {
protected:
public:
  //! Empty.
  SimplexNoise2() {}
  //! Empty.
  ~SimplexNoise2() {}

  double Value(const Vector2&) const;
  double K() const;

protected:
  static const Vector2 grad2[8];  //!< Array of gradients.
  static const int perm[512];    //!< Permutation table, 256 entries duplicated once to avoid modulo computations.
protected:
  static const double F2, G2; //!< Unskew factors for planar case.
};


// Simplex noise
class SimplexNoise4 {
protected:
public:
  //! Empty.
  SimplexNoise4() {}
  //! Empty.
  ~SimplexNoise4() {}

  double Value(const Vector&, const double&) const;
  double K() const;

protected:
  double dot(const int*, const double&, const double&, const double&, const double&) const;

protected:
  static const int grad4[32][4]; //!< Array of gradients for 4D noise.
  static const int perm[512];    //!< Permutation table, 256 entries duplicated once to avoid modulo computations.
  static const int simplex[64][4]; //!< Simplex data for 4D noise.
protected:
  static const double F4, G4; //!< Unskew factors for 3D case.
};


class ExponentialSimplexNoise :protected SimplexNoise2
{
protected:
public:
  ExponentialSimplexNoise() {}
  ~ExponentialSimplexNoise() {}
  double Value(const Vector2&) const;

  using AnalyticScalarField2::Sample;

protected:
  double Exponential(const double&) const;
protected:
  static const Vector2 g[8];
};
