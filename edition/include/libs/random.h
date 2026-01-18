// Core class

#pragma once

// Random number generator
class Random {
protected:
  int vec[48];
  int feed, tap, borrow;
  int V[571];
  int Y;
public:
  Random(int = 239);
  ~Random();

  void Seed(int);
  int Integer(int);

  // Uniform distributions
  double Uniform();
  double Uniform(const double&);
  double Uniform(const double&, const double&);

  //double Normal(const double&, const double&);
  double Exponential(const double&);

  int geometric(double);
  double NormalBoxMuller(const double&, const double&);
  double Angle(unsigned int);

public:
  static Random R239;
protected:
  int Int31();
  int SimpleInt31();
};

// Fast random number generator
class RandomFast
{
private:
  unsigned int x; //!< Variable for computing random values
public:
  RandomFast(unsigned int = 239);
  RandomFast(int);
  //! Empty
  ~RandomFast() {}

  void Seed(unsigned int);
  void Seed(int);

  unsigned int Integer();
  unsigned int Integer(int);
  double Uniform();
  double Uniform(const double&);
  double Uniform(const double&, const double&);

  unsigned int Poisson(const double&);
  int Pick(const double*, int);
};

// Xor and Shift random number generator
class RandomXorShift
{
private:
  unsigned long x, y, z; //!< Internal variables for computing random values.
public:
  RandomXorShift();
  unsigned long Integer();
};

class RandomMul
{
private:
  unsigned int x;//!< Variable for computing random values
public:
  RandomMul(unsigned int = 239);
  unsigned int Integer();
};
