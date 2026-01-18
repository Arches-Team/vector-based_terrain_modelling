// Cpu  

#pragma once


#define __Double__
//#define __Float__
//#define __Real__

// Use C++ real
#ifdef __Double__

#endif

#include <QtCore/QString>
#include <QtCore/QElapsedTimer>
#include <QtGui/QPainter>

// Use C++ float
#ifdef __Float__

#define double float

//! Compute the inverse square root.
inline float isqrt(float x)
{
  float y = 0.5f * x;
  int i = *(int*)&x;
  i = 0x5f3759df - (i >> 1);
  x = *(float*)&i;
  x *= (1.5f - y * x * x);
  x *= (1.5f - y * x * x);
  return x;
}

#endif

// Use C++ real class
#ifdef __Real__

#include "libs/mathematics.h"

#include <ostream>
//using namespace std;

class Real {
protected:
  double x;
public:
  Real() {}
  Real(const double&);
  ~Real() {}

  // Unary operators
  Real operator+ () const;
  Real operator- () const;

  // Assignment operators
  Real& operator+= (const Real&);
  Real& operator-= (const Real&);
  Real& operator*= (const Real&);
  Real& operator/= (const Real&);

  // Binary operators
  friend bool operator> (const Real&, const Real&);
  friend bool operator< (const Real&, const Real&);

  friend bool operator>= (const Real&, const Real&);
  friend bool operator<= (const Real&, const Real&);

  // Binary operators
  friend Real operator+ (const Real&, const Real&);
  friend Real operator- (const Real&, const Real&);
  friend Real operator* (const Real&, const Real&);
  friend Real operator/ (const Real&, const Real&);

  // Boolean functions
  friend bool operator==(const Real&, const Real&);
  friend bool operator!=(const Real&, const Real&);

  friend Real sqrt(const Real&);
  friend Real sin(const Real&);
  friend Real cos(const Real&);
  friend Real fabs(const Real&);
  friend Real asin(const Real&);
  friend Real acos(const Real&);

  friend Real exp(const Real&);
  friend Real log(const Real&);
  friend Real floor(const Real&);
  friend Real ceil(const Real&);

  friend std::ostream& operator<<(std::ostream&, const Real&);
  friend Real atan2(const Real&, const Real&);
  friend Real pow(const Real&, const Real&);

public:
  static void View();
protected:
  static unsigned long long n[10]; //!< Counters for arimthmetic operations
};

inline Real::Real(const double& x) :x(x)
{
}

// Unary operators

inline Real Real::operator+ () const
{
  return *this;
}

inline Real Real::operator- () const
{
  n[0]++;
  return Real(-x);
}

//! Destructive addition.
inline Real& Real::operator+= (const Real& u)
{
  n[0]++;
  x += u.x;
  return *this;
}

//! Destructive subtraction.
inline Real& Real::operator-= (const Real& u)
{
  n[0]++;
  x -= u.x;
  return *this;
}

//! Destructive multiplication.
inline Real& Real::operator*= (const Real& a)
{
  n[1]++;
  x *= a.x;
  return *this;
}

//! Destructive division by a scalar.
inline Real& Real::operator/= (const Real& a)
{
  n[2]++;
  x /= a.x;
  return *this;
}

//! Overloaded
inline bool operator> (const Real& u, const Real& v)
{
  Real::n[3]++;
  return (u.x > v.x);
}

//! Overloaded
inline bool operator< (const Real& u, const Real& v)
{
  Real::n[3]++;
  return (u.x < v.x);
}

//! Overloaded
inline bool operator>= (const Real& u, const Real& v)
{
  Real::n[3]++;
  return (u.x >= v.x);
}

//! Overloaded
inline bool operator<= (const Real& u, const Real& v)
{
  Real::n[3]++;
  return (u.x <= v.x);
}

//! Overloaded
inline Real operator+ (const Real& u, const Real& v)
{
  Real::n[0]++;
  return Real(u.x + v.x);
}

//! Overloaded
inline Real operator- (const Real& u, const Real& v)
{
  Real::n[0]++;
  return Real(u.x - v.x);
}

//! Overloaded
inline Real operator* (const Real& u, const Real& v)
{
  Real::n[1]++;
  return Real(u.x * v.x);
}

//! Overloaded
inline Real operator/ (const Real& u, const Real& v)
{
  Real::n[2]++;
  return Real(u.x / v.x);
}

//! Overloaded
inline bool operator== (const Real& u, const Real& v)
{
  Real::n[3]++;
  return (u.x == v.x);
}

//! Overloaded
inline bool operator!= (const Real& u, const Real& v)
{
  Real::n[3]++;
  return (u.x != v.x);
}

//! Square root
inline Real sqrt(const Real& r)
{
  Real::n[4]++;
  return sqrt(r.x);
}

//! Sine
inline Real sin(const Real& r)
{
  Real::n[5]++;
  return sin(r.x);
}

//! Cosine
inline Real cos(const Real& r)
{
  Real::n[5]++;
  return cos(r.x);
}

inline Real fabs(const Real& r)
{
  Real::n[7]++;
  return Real(fabs(r.x));
}

inline Real atan2(const Real& r, const Real& s)
{
  Real::n[5]++;
  return atan2(r.x, s.x);
}

inline Real pow(const Real& r, const Real& s)
{
  Real::n[6]++;
  return pow(r.x, s.x);
}

inline Real acos(const Real& r)
{
  Real::n[5]++;
  return acos(r.x);
}

inline Real asin(const Real& r)
{
  Real::n[5]++;
  return asin(r.x);
}

inline Real exp(const Real& r)
{
  Real::n[6]++;
  return exp(r.x);
}

inline Real log(const Real& r)
{
  Real::n[6]++;
  return log(r.x);
}

inline Real floor(const Real& r)
{
  Real::n[8]++;
  return floor(r.x);
}

inline Real ceil(const Real& r)
{
  Real::n[8]++;
  return ceil(r.x);
}

//! Overloaded output-stream operator.
inline std::ostream& operator<<(std::ostream& s, const Real& u)
{
  s << u.x;
  return s;
}

#endif

#include "libs/box.h"

class System
{
protected:
  static bool avx; //!< Advanced Vector Extensions flag.
public:
#ifdef _MSC_VER
  static const char* WinGetEnv(const char*);
#endif
  static QString GetEnv(const QString&);
  static QString GetArchesLib();
  static QString GetHeightFieldDir();
  static QString GetDesktop();
  static QString GetResource(const QString&, const QString & = QString(""));
  static QString Memory();

  static QString Elapsed(const QElapsedTimer&);
  static void ShowElapsed(const QElapsedTimer&);
  static QString LongInteger(int);
  static QString VeryLongInteger(long long);
  static QString DateTime();

  static void SavePdf(QGraphicsScene&, const QString&);
  static void SavePdf(QGraphicsScene&, const QString&, const Box2&);
  static void FlipHorizontal(QGraphicsScene&, const Box2&);
  static void FlipVertical(QGraphicsScene&, const Box2&);

  static QImage Rasterize(QGraphicsScene&, const Box2&, int);
  static void Compose(QImage&, const QImage&, QPainter::CompositionMode = QPainter::CompositionMode_SourceOver);
  static QImage Compose(const QImage&, const QImage&, QPainter::CompositionMode = QPainter::CompositionMode_SourceOver);

  static void SaveImage(const QImage&, const QString & = QString("../image"));

  static int CodeLines(const QString & = GetArchesLib(), bool = true, bool = true);

  static void SetAvx(bool);
  static bool Avx();
protected:
  static int NumberOfLines(const QString&);
  static int RecurseDirectory(const QString&, bool, bool);
};

class Timer {
protected:
  QElapsedTimer qtimer; //<! Qt
  int t; //!< Elapsed time.
public:
  Timer() { qtimer.restart(); t = 0; }
  int Elapsed() const { return t + qtimer.elapsed(); }
  void Restart() { qtimer.restart(); }
  void Start() { t = 0; qtimer.restart(); }
  void Stop() { t += qtimer.elapsed(); }
};

#include <iostream>
// using namespace std;

class Byte {
public:
  static void Set(unsigned long&, int);
};

inline void Byte::Set(unsigned long& u, int b)
{
  u = (1L << b) | u;
}
