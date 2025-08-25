
// Core class

#include <ostream>
#include <cmath>
#include "libs/random.h"

/*!
\class RandomFast random.h
\brief A fast linear congruential random number generator.

It uses the linear congruential generator
x<SUB>n+1</SUB> = 3039177861 x<SUB>n</SUB> (mod 232), introduced by Borosh and
Niederreiter. This generator is fast and, despite its simplicity, optimal with
respect to statistical independence of successive pseudo-random numbers.

\ingroup MathGroup

\sa RandomMul, RandomXorShift

*/

/*!
\brief Create a fast random number generator.
\param s Seed.
*/
RandomFast::RandomFast(unsigned int s)
{
  Seed(s);
}

/*!
\brief Create a fast random number generator.
\param s Seed.
*/
RandomFast::RandomFast(int s)
{
  Seed(s);
}

/*!
\brief Changer the seed of the random number generator.
\param s Seed.
*/
void RandomFast::Seed(unsigned int s)
{
  x = s;
}

/*!
\brief Changer the seed of the random number generator.
\param s Seed.
*/
void RandomFast::Seed(int s)
{
  x = (unsigned int)(s);
}

/*!
\brief Compute a random unsigned integer.
*/
unsigned int RandomFast::Integer()
{
  x *= 3039177861;
  return x;
}

/*!
\brief Compute a random unsigned integer within interval [0,a].
\param a Integer.
*/
unsigned int RandomFast::Integer(int a)
{
  return Integer() % a;
}

/*!
\brief Compute uniform distribution in [0, 1[.
*/
double RandomFast::Uniform()
{
  return double(Integer()) * (1.0 / double(4294967295U)); // Should be UINT_MAX
}

/*!
\brief Compute uniform distribution in [0, a[.
\param a Amplitude.
*/
double RandomFast::Uniform(const double& a)
{
  return a * Uniform();
}

/*!
\brief Compute uniform distribution in [a, b[.
\param a, b Amplitude interval.
*/
double RandomFast::Uniform(const double& a, const double& b)
{
  return a + (b - a) * Uniform();
}

/*!
\brief Poisson random sequence.
\param mu Mean value.
*/
unsigned int RandomFast::Poisson(const double& mu)
{
  double g = std::exp(-mu);
  unsigned int em = 0;
  double t = Uniform();
  while (t > g)
  {
    ++em;
    t *= Uniform();
  }
  return em;
}

/*!
\brief Randomly pick one item given an array of probabilities.
\param a %Array.
\param n Size.
*/
int RandomFast::Pick(const double* a, int n)
{
  double s = 0.0;
  for (int i = 0; i < n; i++)
  {
    s += a[i];
  }
  double r = Uniform(s);

  s = 0.0;
  for (int i = 0; i < n; i++)
  {
    s += a[i];
    if (r < s)
    {
      return i;
    }
  }
  return n - 1;
}
