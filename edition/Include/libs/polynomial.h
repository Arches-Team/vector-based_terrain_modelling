// Polynomials

#pragma once

#include "libs/nonic.h"

// Polynomials
class Polynomial {
protected:
  int n = 0;   //!< The degree of the polynomial.
  static constexpr int MaxDegree = 11; //!< Maximum degree.
  double c[MaxDegree + 1] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; //!< %Array of coefficients.
public:
  Polynomial();
  //! Empty.
  ~Polynomial() {}

  explicit Polynomial(const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);
  explicit Polynomial(const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);
  explicit Polynomial(const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);
  explicit Polynomial(const double&, const double&, const double&, const double&, const double&, const double&, const double&);
  explicit Polynomial(const double&, const double&, const double&, const double&, const double&, const double&);
  explicit Polynomial(const double&, const double&, const double&, const double&, const double&);
  explicit Polynomial(const double&, const double&, const double&, const double&);
  explicit Polynomial(const double&, const double&, const double&);
  explicit Polynomial(const double&, const double&);
  explicit Polynomial(const double&);
  Polynomial(const Linear&);
  Polynomial(const Quadric&);
  Polynomial(const Cubic&);
  Polynomial(const Quartic&);
  Polynomial(const Quintic&);
  Polynomial(const Sextic&);
  Polynomial(const Septic&);
  Polynomial(const Octic&);
  Polynomial(const Nonic&);

  constexpr double& operator[] (int);
  constexpr double operator[] (int) const;

  //! Overloaded.
  Polynomial operator+ () const { return *this; }
  Polynomial operator- () const;

  Polynomial& operator+= (const Polynomial&);
  Polynomial& operator-= (const Polynomial&);
  Polynomial& operator*= (const Polynomial&);
  Polynomial& operator*= (const double&);
  Polynomial& operator/= (const double&);

  // Binary operators
  friend Polynomial operator+ (const Polynomial&, const Polynomial&);
  friend Polynomial operator- (const Polynomial&, const Polynomial&);

  friend Polynomial operator* (const Polynomial&, const Polynomial&);

  friend Polynomial operator* (const Polynomial&, const double&);
  friend Polynomial operator* (const double&, const Polynomial&);
  friend Polynomial operator/ (const Polynomial&, const double&);

  Polynomial Prime() const;
  Polynomial Compose(const Polynomial&) const;

  Polynomial Reversed() const;

  // Evaluates polynomial
  double operator()(const double&) const;

  friend std::ostream& operator<<(std::ostream&, const Polynomial&);

  // Solve
  int Solve(double*) const;
  int SturmSolve(double*) const;
  int Solve(double*, const double&, const double&) const;

  int Degree() const;
  void Check();
  void Range(double&, double&, const double& = 0.0, const double& = 1.0) const;

  static Polynomial Compose(const Quartic&, const Quadric&);
  int Bissection(double, double, double&, const double&) const;

  static const Polynomial Null;
private:
  int Bissection(int, double, double, int, int, double*);
  int Sturm(Polynomial*) const;
  static int VisibleRoots(int, Polynomial*, int&, int&);
  static int Changes(int, Polynomial*, const double&);
  static int modp(const Polynomial*, const Polynomial*, Polynomial*);
private:
  static const double tiny; //!< Very small \htmlonly\epsilon;\endhtmlonly, indicating coefficient limit.
  static const double epsilon; //!< Small \htmlonly\epsilon;\endhtmlonly, indicating tolerance around 0.0.
  static const int iterations; //!< Maximum number of iterations.
};

//! Creates a degree 9 polynomial.
inline Polynomial::Polynomial(const double& a9, const double& a8, const double& a7, const double& a6, const double& a5, const double& a4, const double& a3, const double& a2, const double& a1, const double& a0) { c[0] = a0; c[1] = a1; c[2] = a2; c[3] = a3; c[4] = a4; c[5] = a5; c[6] = a6; c[7] = a7; c[8] = a8; c[9] = a9; n = 9; }

//! Creates a degree 8 polynomial.
inline Polynomial::Polynomial(const double& a8, const double& a7, const double& a6, const double& a5, const double& a4, const double& a3, const double& a2, const double& a1, const double& a0) { c[0] = a0; c[1] = a1; c[2] = a2; c[3] = a3; c[4] = a4; c[5] = a5; c[6] = a6; c[7] = a7; c[8] = a8; n = 8; }

//! Creates a septic.
inline Polynomial::Polynomial(const double& a7, const double& a6, const double& a5, const double& a4, const double& a3, const double& a2, const double& a1, const double& a0) { c[0] = a0; c[1] = a1; c[2] = a2; c[3] = a3; c[4] = a4; c[5] = a5; c[6] = a6; c[7] = a7; n = 7; }

//! Creates a sextic.
inline Polynomial::Polynomial(const double& a6, const double& a5, const double& a4, const double& a3, const double& a2, const double& a1, const double& a0) { c[0] = a0; c[1] = a1; c[2] = a2; c[3] = a3; c[4] = a4; c[5] = a5; c[6] = a6; n = 6; }

//! Creates a quintic.
inline Polynomial::Polynomial(const double& a5, const double& a4, const double& a3, const double& a2, const double& a1, const double& a0) { c[0] = a0; c[1] = a1; c[2] = a2; c[3] = a3; c[4] = a4; c[5] = a5; n = 5; }

//! Creates a quartic.
inline Polynomial::Polynomial(const double& a4, const double& a3, const double& a2, const double& a1, const double& a0) { c[0] = a0; c[1] = a1; c[2] = a2; c[3] = a3; c[4] = a4; n = 4; }

//! Creates a cubic.
inline Polynomial::Polynomial(const double& a, const double& b, const double& a1, const double& d) { c[0] = d; c[1] = a1; c[2] = b; c[3] = a; n = 3; }

//! Creates a quadric.
inline Polynomial::Polynomial(const double& a, const double& b, const double& a0) { c[0] = a0; c[1] = b; c[2] = a; n = 2; }

//! Creates a linear polynomial.
inline Polynomial::Polynomial(const double& a, const double& b) { c[0] = b; c[1] = a; n = 1; }

//! Creates a constant polynomial.
inline Polynomial::Polynomial(const double& a0) { c[0] = a0; n = 0; }

/*!
\brief Creates a polynomial given a linear.
\param p %Linear polynomial.
*/
inline Polynomial::Polynomial(const Linear& p) :n(1)
{
  c[0] = p[0];
  c[1] = p[1];
}

/*!
\brief Creates a polynomial given a quadric.
\param p %Quadric polynomial.
*/
inline Polynomial::Polynomial(const Quadric& p) :n(2)
{
  c[0] = p[0];
  c[1] = p[1];
  c[2] = p[2];
}

/*!
\brief Creates a polynomial given a cubic.

\param p %Cubic polynomial.
*/
inline Polynomial::Polynomial(const Cubic& p) :n(3)
{
  c[0] = p[0];
  c[1] = p[1];
  c[2] = p[2];
  c[3] = p[3];
}

/*!
\brief Creates a polynomial given a quartic.
\param p %Quartic polynomial.
*/
inline Polynomial::Polynomial(const Quartic& p) : n(4)
{
  c[0] = p[0];
  c[1] = p[1];
  c[2] = p[2];
  c[3] = p[3];
  c[4] = p[4];
}

/*!
\brief Creates a polynomial given a quintic.
\param p %Quintic polynomial.
*/
inline Polynomial::Polynomial(const Quintic& p) :n(5)
{
  c[0] = p[0];
  c[1] = p[1];
  c[2] = p[2];
  c[3] = p[3];
  c[4] = p[4];
  c[5] = p[5];
}

/*!
\brief Creates a polynomial given a sextic.
\param p %Sextic polynomial.
*/
inline Polynomial::Polynomial(const Sextic& p) :n(6)
{
  c[0] = p[0];
  c[1] = p[1];
  c[2] = p[2];
  c[3] = p[3];
  c[4] = p[4];
  c[5] = p[5];
  c[6] = p[6];
}

/*!
\brief Creates a polynomial given a septic.
\param p %Septic polynomial.
*/
inline Polynomial::Polynomial(const Septic& p) :n(7)
{
  c[0] = p[0];
  c[1] = p[1];
  c[2] = p[2];
  c[3] = p[3];
  c[4] = p[4];
  c[5] = p[5];
  c[6] = p[6];
  c[7] = p[7];
}

/*!
\brief Creates a polynomial given an octic.
\param p %Octic polynomial.
*/
inline Polynomial::Polynomial(const Octic& p) :n(8)
{
  c[0] = p[0];
  c[1] = p[1];
  c[2] = p[2];
  c[3] = p[3];
  c[4] = p[4];
  c[5] = p[5];
  c[6] = p[6];
  c[7] = p[7];
  c[8] = p[8];
}

/*!
\brief Creates a polynomial given a nonic.
\param p %Nonic polynomial.
*/
inline Polynomial::Polynomial(const Nonic& p) :n(9)
{
  c[0] = p[0];
  c[1] = p[1];
  c[2] = p[2];
  c[3] = p[3];
  c[4] = p[4];
  c[5] = p[5];
  c[6] = p[6];
  c[7] = p[7];
  c[8] = p[8];
  c[9] = p[9];
}

//! Returns the degree of the polynomial.
inline int Polynomial::Degree() const
{
  return n;
}

//! Access to the cefficients of the polynomial.
inline constexpr double& Polynomial::operator[] (int i)
{
  return c[i];
}

//! Overloaded.
inline constexpr double Polynomial::operator[] (int i) const
{
  return c[i];
}

