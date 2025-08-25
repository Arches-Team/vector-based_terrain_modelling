// Shuffling algorithms

#include "libs/shuffle.h"
#include "libs/mathematics.h"

/*!
\class Shuffle shuffle.h
\brief %Shuffling algorithms for generation permutation tables.

\ingroup MathGroup
*/

RandomFast Shuffle::r;

/*!
\brief Compute a permutation table.

This is an implementation of a modified version of the Fisherâ€“Yates algorithm.

\param n Size.
*/
QVector<int> Shuffle::Table(int n)
{
  QVector<int> t(n);
  for (int i = 0; i < n; i++)
  {
    int j = r.Integer() % (i + 1); // Random integer such that 0 <= j <= i
    if (j != i)
    {
      t[i] = t[j];
    }
    t[j] = i;
  }
  return t;
}
