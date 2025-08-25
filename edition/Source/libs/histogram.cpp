//#include <QtCharts/QtCharts>
//#include <QtCharts/QLineSeries>

#include "libs/histogram.h"
#include "libs/mathematics.h"
#include "libs/chart.h"

/*!
\brief Create a histogram.
\param a,b Range.
\param size Number of bins.
*/
Histogram::Histogram(int size, double a, double b) :a(a), b(b), v(size, 0), e(size / (b - a))
{
}

/*!
\brief Create a histogram from a set.
\param values Set of values.
\param a,b Range.
\param size Number of bins.
*/
Histogram::Histogram(int size, double a, double b, const QVector<double>& values) :Histogram(size, a, b)
{
  for (const double& x : values)
  {
    v[Index(x)]++;
  }
}
/*!
\brief Add a value to the histogram.
\param x Real.
*/
void Histogram::Insert(const double& x)
{
  v[Index(x)]++;
}

/*!
\brief Return the size of the histogram.
\param size Number of bins.
\param values Set of values.
*/
Histogram::Histogram(int size, const QVector<double>& values) :Histogram(size, 0.0, 0.0)
{
  if (size == 0) return;

  a = values.at(0);
  b = values.at(0);

  for (const double& x : values)
  {
    Math::SetMinMax(x, a, b);
  }

  e = size / (b - a);

  for (const double& x : values)
  {
    v[Index(x)]++;
  }
}


/*!
\brief Draw.
*/
void Histogram::Draw(const QString& title) const
{

  DrawChart graph = DrawChart(title);
  graph.AddBars(v, a, b);
  graph.Display();

}

/*!
\brief Draw.
*/
void Histogram::ExportPDF(const QString& filename) const
{
  DrawChart graph = DrawChart("Slope Histogram");
  graph.AddBars(v, a, b);
  //graph.AddLine(histo);
  //graph.AddScatter(histo);
  graph.ExportPDF(filename);
}

/*!
\brief Draw.
*/
void Histogram::ExportPNG(const QString& filename) const
{
  DrawChart graph = DrawChart("Slope Histogram");
  //graph.AddBars(histo, min, max);
  graph.AddLine(v);
  //graph.AddScatter(histo);
  graph.ExportPNG(filename);
}


/*!
\brief Transform into a cumulated histogram
*/
/*
Histogram Histogram::Cumulated() const
{
  Histogram c(histo.size(), a, b);

  c.histo[0] = histo.at(0);
  for (int i = 1; i < histo.size(); i++)
  {
    c.histo[i] = c.histo[i - 1] + histo.at(i);
  }

  return c;
}
*/

/*!
\brief Compute the value in the histogram corresponding to a given percentage of the cumulated values.
\param t Percentage, in [0.0,1.0].
*/
double Histogram::Select(const double& t) const
{
  if (t <= 0.0)
    return a;

  double vt = Math::Lerp(a, b, t);
  int k = v.size() - 1;

  double c = 0.0;
  for (int i = 0; i < v.size(); i++)
  {
    c += v.at(i);
    if (c > vt)
    {
      k = i;
      break;
    }
  }
  return v.at(k); // Should use linear interpolation instead
}

/*!
\brief Overloaded.
\param s Stream.
\param h The histogram.
*/
std::ostream& operator<<(std::ostream& s, const Histogram& h)
{
  s << "Histogram(" << h.v.size() << ", " << h.a << ", " << h.b;

  for (int i : h.v)
  {
    s << ", " << std::endl << i;
  }
  s << ")" << std::endl;
  return s;
}

int Histogram::MaxCount() const
{
  int m = 0;
  for (int c : v)
  {
    if (c > m)
    {
      m = c;
    }
  }
  return m;
}

/*!
\brief Create a (small) image showing the histogram distribution.
*/
QImage Histogram::CreateImage() const
{
  constexpr const int width = 256;
  constexpr const int height = 256;
  QImage image(width, height, QImage::Format_ARGB32);

  double M = double(MaxCount());

  // White
  image.fill(QColor(255, 255, 255));

  // Horizontal separator
  for (int i = 0; i < width; i++)
  {
    image.setPixelColor(i, height / 2, QColor(245, 245, 245));
  }

  for (int i = 0; i < width; i++)
  {
    int k = int(double(v.size() - 1) * Math::Unit(i, width));

    int y = int((height - 1) * double(v[k]) / M);

    for (int j = height - 1 - y; j < height; j++)
    {
      image.setPixelColor(i, j, QColor(181, 175, 169));
    }

    image.setPixelColor(i, height - 1 - y, QColor(255, 175, 169));

  }

  return image;
}