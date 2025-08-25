// Polynomials

#include "libs/polynomial.h"

/*!
\brief Create the Sturm sequences of a polynomial.

Recursively compute the successive derivatives of the polynomial, normalize the
leading coefficient and the polynomial accordingly and build the
sequence.

\param s Sequence, <I>i.e.</I> an array of polynomials.
*/
int Polynomial::Sturm(Polynomial* s) const
{
  Polynomial* sp;

  // Initialize sequence
  s[0] = *this;

  s[0].n = n;
  s[1].n = n - 1;

  // Compute the derivative and normalize the leading coefficient
  double f = fabs(s[0][n] * n);
  double* fp = &(s[1][0]);
  double* fc = &(s[0][1]);

  for (int i = 1; i <= n; i++)
  {
    *fp++ = *fc++ * i / f;
  }

  // Create the rest of the sequence
  for (sp = &(s[2]); modp(sp - 2, sp - 1, sp); sp++)
  {
    // Reverse the sign and normalize     
    f = -fabs((*sp)[sp->n]);

    for (fp = &((*sp)[sp->n]); fp >= &((*sp)[0]); fp--)
    {
      *fp /= f;
    }
  }
  (*sp)[0] = -(*sp)[0];

  return int(sp - s);
}

/*!
\brief Given a Sturm sequence, compute the visible roots.
*/
int Polynomial::VisibleRoots(int np, Polynomial* sseq, int& atzer, int& atpos)
{
  int atposinf = 0;
  int atzero = 0;
  Polynomial* s;

  // changes at positve infinity 

  double lf = sseq[0][sseq[0].n];

  for (s = sseq + 1; s <= sseq + np; s++)
  {
    double f = (*s)[s->n];

    if (lf == 0.0 || lf * f < 0)
    {
      atposinf++;
    }

    lf = f;
  }

  // Changes at zero 

  lf = sseq[0][0];

  for (s = sseq + 1; s <= sseq + np; s++)
  {
    double f = (*s)[0];

    if (lf == 0.0 || lf * f < 0.0)
    {
      atzero++;
    }

    lf = f;
  }

  atzer = atzero;
  atpos = atposinf;

  return (atzero - atposinf);
}

/*!
\brief Given a Sturm sequence, compute the number of sign changes for a real number.
\param np Order of the polynomial.
\param sseq Sturm sequence.
\param a Value.
*/
int Polynomial::Changes(int np, Polynomial* sseq, const double& a)
{
  // Changes
  int c = 0;

  double lf = sseq[0](a);

  for (Polynomial* s = sseq + 1; s <= sseq + np; s++)
  {
    double f = (*s)(a);

    if ((lf == 0.0) || (lf * f < 0.0))
    {
      c++;
    }
    lf = f;
  }

  return c;
}

/*!
\brief Seach the all the roots of a polynomial performing bissection over a given interval.

This function first tries the standard bissection method if there is only one root,
and performs more complex search otherwize.

\param a,b Interval.
\param np Order of the polynomial.
\param atmin, atmax ?
\param roots Array of roots.
*/
int Polynomial::Bissection(int np, double a, double b, int atmin, int atmax, double* roots)
{
  double mid = (a + b) / 2.0;

  if ((atmin - atmax) == 1)
  {
    // Try regula-falsa to find the root

    if (Bissection(a, b, *roots, epsilon))
    {
      return 1;
    }
    else
    {
      // That failed, so now find it by bisection 

      for (int i = 0; i < iterations; i++)
      {
        mid = (a + b) / 2.0;

        int atmid = Changes(np, this, mid);

        if (fabs(mid) > epsilon)
        {
          if (fabs((b - a) / mid) < epsilon)
          {
            roots[0] = mid;

            return 1;
          }
        }
        else
        {
          if (fabs(b - a) < epsilon)
          {
            roots[0] = mid;

            return 1;
          }
        }

        if ((atmin - atmid) == 0)
        {
          a = mid;
        }
        else
        {
          b = mid;
        }
      }

      // Bisection took too long-just return what we got

      roots[0] = mid;

      return 1;
    }
  }

  // There is more than one root in the interval : bisect to find the next intervals

  for (int i = 0; i < iterations; i++)
  {
    mid = (a + b) / 2.0;
    int atmid = Changes(np, this, mid);

    int n1 = atmin - atmid;
    int n2 = atmid - atmax;

    if ((n1 != 0) && (n2 != 0))
    {
      n1 = Bissection(np, a, mid, atmin, atmid, roots);
      n2 = Bissection(np, mid, b, atmid, atmax, &roots[n1]);

      return n1 + n2;
    }

    if (n1 == 0)
    {
      a = mid;
    }
    else
    {
      b = mid;
    }
  }

  // Took too long to bisect-just return what we got
  roots[0] = mid;

  return 1;
}

/*!
\brief Search the roots of a polynomial equation.

This function calls other root solving functions,
invoking Quartic::Quartic(), Cubic::Solve() or
Quadratic::Solve() as appropriate if the degree is lower
of equal than four (which means that closed form expressions of the
roots exist).

Otherwise, rely on a numeric root finder to isolate roots and converge
using the Sturm sequences.

\param roots Array for storing the roots.
\return The number of roots.
*/
int Polynomial::Solve(double* roots) const
{
  if (n < 5)
  {
    if (n == 4)
      return Quartic(c[4], c[3], c[2], c[1], c[0]).Solve(roots);
    else if (n == 3)
      return Cubic(c[3], c[2], c[1], c[0]).Solve(roots);
    else
      return Quadric(c[2], c[1], c[0]).Solve(roots);
  }
  if (n == 5)
  {
    if (c[5] == 0.0)
    {
      return Quartic(c[4], c[3], c[2], c[1], c[0]).Solve(roots);
    }
  }

  return SturmSolve(roots);
}

/*!
\brief Search the roots of a polynomial equation using the Sturm sequences.
\param roots Array for storing the roots.
\return The number of roots.
*/
int Polynomial::SturmSolve(double* roots) const
{
  Polynomial sseq[16];

  int atmin, atmax;

  // Build the sequence   
  int np = Sturm(&sseq[0]);

  int nroots = VisibleRoots(np, sseq, atmin, atmax);
  // Get the total number of visible roots  
  if (nroots == 0)
  {
    return 0;
  }

  // Bracket the roots
  const double min_value = 0.0;
  const double max_value = 1.0e7;

  atmin = Changes(np, sseq, min_value);
  atmax = Changes(np, sseq, max_value);

  nroots = atmin - atmax;

  if (nroots == 0)
  {
    return 0;
  }

  // perform the bisection   
  int temp = sseq->Bissection(np, min_value, max_value, atmin, atmax, roots);
  return temp;
}

/*!
\brief Search the roots of a polynomial equation over a given interval.
\param roots The array of roots.
\param a, b The search interval.
\return The number of roots.
*/
int Polynomial::Solve(double* roots, const double& a, const double& b) const
{
  int r = Solve(roots);
  int j = 0;

  for (int i = 0; i < r; i++)
  {
    if ((roots[i] > a) && (roots[i] < b))
    {
      roots[j] = roots[i];
      j++;
    }
  }

  return j;
}
