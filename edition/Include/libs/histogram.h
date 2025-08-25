#pragma once

#include <iostream>
#include <QVector>
#include <QImage>

#include "libs/ia.h"


class Histogram {
protected:
  QVector<int> v;   //!< Counts.
  double a = 0.0, b = 1.0;  //!< Interval.
  double e = 1.0; // Internal factor for computing bins.
public:
  explicit Histogram(int, double = 0.0, double = 1.0);
  Histogram(int, double, double, const QVector<double>&);
  Histogram(int, const QVector<double>&);

  //! Empty.
  ~Histogram() {}
  //Histogram Cumulated() const;

  void Insert(const double&);

  double Select(const double&) const;

  int GetSize() const;
  Ia Range() const;

  QVector<int> GetValues() const;

  void Draw(const QString & = "Histogram") const;
  void ExportPDF(const QString & = "histogram.pdf") const;
  void ExportPNG(const QString & = "histogram.png") const;
  QImage CreateImage() const;

  friend std::ostream& operator<<(std::ostream&, const Histogram&);

  int MaxCount() const;
protected:
  int Index(const double&) const;
};

/*!
\brief Compute the entry of a given value.
\param x Real.
*/
inline int Histogram::Index(const double& x) const
{
  int k = Clamp(int((x - a) * e), 0, v.size() - 1);
  return k;
}

/*!
\brief Return the size of the histogram.
*/
inline int Histogram::GetSize() const
{
  return v.size();
}
/*!
\brief Return the range of values of the histogram.
*/
inline Ia Histogram::Range() const
{
  return Ia(a, b);
};

/*!
\brief Return the set containing the values.
*/
inline QVector<int> Histogram::GetValues() const
{
  return v;
}
