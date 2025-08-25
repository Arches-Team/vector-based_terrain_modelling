// Fields
#include <immintrin.h>
#include <math.h>

#include <QtCore/QFile>
#include <QtCore/QTextStream>
#include <QtGui/QImageReader>

#include "libs/scalarfield.h"
#include "libs/vectorfield.h"
#include "libs/cubic.h"
#include "libs/cpu.h"
#include "libs/histogram.h"
#include "libs/quintic.h"
#include "libs/segment.h"
#include "libs/polygon.h"
#include "libs/mesh.h"


/*!
\class ScalarField2 scalarfield.h
\brief A base two-dimensional field of real values.

\ingroup StructureGroup
*/

/*!
\brief Create the field structure.
\param a Array representing the grid domain.
\param v Constant value of field.
*/
ScalarField2::ScalarField2(const Array2& a, const double& v) :Array2(a)
{
  field.fill(v, nx * ny);
}

/*!
\brief Create the field structure.
\param box The box.
\param x,y Size of the array.
\param v Constant value of field.
*/
ScalarField2::ScalarField2(const Box2& box, int x, int y, const double& v) :ScalarField2(Array2(box, x, y), v)
{
}

/*!
\brief Create the field structure.
\param box The box.
\param x,y Size of the array.
\param v Array of scalar values.
*/
ScalarField2::ScalarField2(const Box2& box, int x, int y, const QVector<double>& v) :Array2(box, x, y),field(v)
{
}

/*!
\brief Create a field structure from an image.
\param image The image.
\param v Range of values to scale the image, values will vary between 0.0 and v.
\param d Distance between samples.
\param grayscale Read image as grayscale if set to true, otherwize use color for a better accuracy.
*/
ScalarField2::ScalarField2(const QImage& image, const double& v, const double& d, bool grayscale) :ScalarField2(Box2(image.size()).Scaled(d), image, 0.0, v, grayscale)
{
}

/*!
\brief Create a field structure from an image.
\param box The box.
\param image The image.
\param a,b Range of values to scale the image.
\param grayscale Read image as grayscale if set to true, otherwize use color for a better accuracy.
*/
ScalarField2::ScalarField2(const Box2& box, const QImage& image, const double& a, const double& b, bool grayscale) :Array2(box, image.width(), image.height())
{
  // Set size
  field.resize(nx * ny);

  // Write Heightmap
  for (int i = 0; i < image.width(); i++)
  {
    for (int j = 0; j < image.height(); j++)
    {
      double t = 0.0;

      // Grayscale
      if (grayscale)
      {
        // Grayscale 16 bits
        if (image.format() == QImage::Format_Grayscale16 || image.format() == QImage::Format_RGBA64)
        {
          QColor thecolor = image.pixelColor(i, j);
          t = thecolor.blueF();
        }
        // Grayscale 8 bits
        else
        {
          QRgb color = image.pixel(i, j);
          t = double(qGray(color)) / 255.0;
        }
      }
      // Color 
      else
      {
        QRgb color = image.pixel(i, j);
        // Maximum value is 256^3-1
        int u = qRed(color) << 16 | qGreen(color) << 8 | qBlue(color);

        t = double(u) / (16777216.0 - 1.0);
      }
      field[VertexIndex(i, j)] = Math::Lerp(a, b, t);
    }
  }
}

/*!
\brief Create a field structure from an image.

This is a convenience constructor which is the same as:
\code
ScalarField2 s(Box2(0.5), QImage(name),0.0,1.0,true);
\endcode

\param name Name of the image.
\param a,b Range of values to scale the image.
*/
ScalarField2::ScalarField2(const QString& name, double a, double b) :ScalarField2(Box2(0.5), QImage(name), a, b, true)
{
}

/*!
\brief Loads a scalar field from an image.

The function uses a specific image reader with unlimited allocation limit.

\param box The box.
\param path Name of the image.
\param a,b Range of values to scale the image.
*/

ScalarField2 ScalarField2::Load(Box2 box, const QString& path, double a, double b)
{
  QImageReader reader;

  int l = QImageReader::allocationLimit();

  QImageReader::setAllocationLimit(0);
  reader.setFileName(path);
  QImage image = reader.read();
  QImageReader::setAllocationLimit(l);

  return ScalarField2(box, image, a, b);
}

/*!
\brief Import FAW Float 32 bits format that comes from World Machine.

Suppose that the domain is a square.
*/
ScalarField2 ScalarField2::LoadFromR32(const Box2& box, const QString& filename)
{
  QVector<double> values;
  int nx, ny;
  QFile file(filename);
  file.open(QIODevice::ReadOnly);
  QDataStream in(&file);
  in.setByteOrder(QDataStream::LittleEndian);
  in.setFloatingPointPrecision(QDataStream::SinglePrecision);
  float v;
  while (!file.atEnd()) {
    in >> v;
    values << v;
  }
  file.close();

  nx = (int)sqrt(values.size());
  ny = nx;
  ScalarField2 sf(box, nx, ny, values);
  return sf;
}

/*!
\brief Empty.
*/
ScalarField2::~ScalarField2()
{
  //std::cout << "~ScalarField2 "<<this << std::endl;
}

/*!
\brief Return the array representing the grid domain.
*/
Array2 ScalarField2::GetArray() const
{
  return Array2(Box2(a, b), nx, ny);
}

/*!
\brief Set all values in the field that are lower than epsilon to true zero.
\param e Epsilon coefficient.
*/
void ScalarField2::CutEpsilon(const double& e)
{
  // Process all elements
  for (int i = 0; i < nx * ny; i++)
  {
    if (field.at(i) < e)
    {
      field[i] = 0.0;
    }
  }
}

/*!
\brief Get the integral of the scalar field.

Compute the sum of all the elements and make the product by the area of small cell element.

\sa Sum()
*/
double ScalarField2::Integral() const
{
  double s = Sum();

  // Scale by the size of a cell element
  s *= CellArea();

  return s;
}

/*!
\brief Compute the average value of the elements in the scalar field.

\avx

\sa Sum()
*/
double ScalarField2::Average() const
{
  const int size = nx * ny;
  double s = Sum();
  s /= size;
  return s;
}


/*!
\brief Compute the sum of the elements in the scalar field.

\avx

\sa Integral(), Average()
*/
double ScalarField2::Sum() const
{
  const int size = nx * ny;

  double s = 0.0;

#ifdef _MSC_VER
  if (System::Avx() && field.size() > 8)
  {
    const double* p = field.data();

    __m256d avxa = _mm256_load_pd(p);
    __m256d avxb = _mm256_load_pd(p + 4);

    const int offset = size % 8;
    for (int i = 1; i < size / 8; i++)
    {
      __m256d avxi = _mm256_load_pd(p + i * 8);
      avxa = _mm256_add_pd(avxa, avxi);
      __m256d avxj = _mm256_load_pd(p + i * 8 + 4);
      avxb = _mm256_add_pd(avxb, avxj);
    }
    avxa = _mm256_add_pd(avxa, avxb);
    double e[4];
    _mm256_store_pd(e, avxa);
    s = e[0] + e[1] + e[2] + e[3];

    // Could try this instead of the previous four lines
    // avxa = _mm256_hadd_pd(avxa, avxa);
    // s= ((double*)&avxa)[0] + ((double*)&avxa)[2];

    for (int i = size - offset; i < size; i++)
    {
      s += field.at(i);
    }
  }
  else
#endif
  {
    // Process all elements

    for (int i = 0; i < size; i++)
    {
      s += field.at(i);
    }
  }
  return s;
}

/*!
\brief Perform a symmetry.
\param x, y Axes of symmetry.
*/
void ScalarField2::Symmetry(bool x, bool y)
{
  if ((x == false) && (y == false)) return;

  else if ((x == false) && (y == true))
  {
    for (int i = 0; i < nx; i++)
    {
      for (int j = 0; j < ny / 2; j++)
      {
        Math::Swap(field[VertexIndex(i, j)], field[VertexIndex(i, ny - 1 - j)]);
      }
    }
  }
  else if ((x == true) && (y == false))
  {
    for (int i = 0; i < nx / 2; i++)
    {
      for (int j = 0; j < ny; j++)
      {
        Math::Swap(field[VertexIndex(i, j)], field[VertexIndex(nx - 1 - i, j)]);
      }
    }
  }
  else
  {
    Symmetry(true, false);
    Symmetry(false, true);
  }
}

/*!
\brief Counter-clockwise rotation, keep the center of the domain unchanged.
*/
void ScalarField2::Rotate()
{
  // Square
  if (nx == ny)
  {
    for (int x = 0; x < nx / 2; x++)
    {
      // Consider elements in group of 4 in current square
      for (int y = x; y < ny - x - 1; y++)
      {
        // Store current cell in
        // temp variable
        double temp = at(x, y);

        // Move values from right to top
        (*this)(x, y) = (*this)(y, nx - 1 - x);

        // Move values from bottom to right
        (*this)(y, nx - 1 - x) = (*this)(nx - 1 - x, nx - 1 - y);

        // Move values from left to bottom
        (*this)(nx - 1 - x, nx - 1 - y) = (*this)(nx - 1 - y, x);

        // Assign temp to left
        (*this)(nx - 1 - y, x) = temp;
      }
    }
  }
  else
  {
    Vector2 c = Box2::Center();
    Vector2 d = Box2::Diagonal();
    Math::Swap(d[0], d[1]);

    // Change box
    a = c - d;
    b = c + d;

    // Create scalar field
    ScalarField2 rotated(Box2(a, b), ny, nx);

    // Transform data
    for (int x = 0; x < nx; x++)
    {
      for (int y = 0; y < ny; y++)
      {
        rotated(y, x) = (*this)(x, y);
      }
    }
    *this = rotated;
  }
}

/*!
\brief Compute the norm of the scalar field.

\avx

\sa Average()
*/
double ScalarField2::ScaledNorm() const
{
  const int size = nx * ny;

#ifdef _MSC_VER
  if (System::Avx() && field.size() > 8)
  {
    double s = 0.0;

    const double* p = field.data();

    __m256d avxa = _mm256_load_pd(p);
    __m256d avxb = _mm256_load_pd(p + 4);

    const int offset = size % 8;
    for (int i = 1; i < size / 8; i++)
    {
      __m256d avxi = _mm256_load_pd(p + i * 8);
      avxi = _mm256_mul_pd(avxi, avxi);
      avxa = _mm256_add_pd(avxa, avxi);
      __m256d avxj = _mm256_load_pd(p + i * 8 + 4);
      avxj = _mm256_mul_pd(avxj, avxj);
      avxb = _mm256_add_pd(avxb, avxj);
    }
    avxa = _mm256_add_pd(avxa, avxb);
    double e[4];
    _mm256_store_pd(e, avxa);
    s = e[0] + e[1] + e[2] + e[3];

    for (int i = size - offset; i < size; i++)
    {
      s += Math::Sqr(field.at(i));
    }

    // Scale by the size of a cell element
    s = sqrt(s) / size;

    return s;
  }
  else
#endif
  {
    double s = 0.0;
    // Process all elements

    for (int i = 0; i < size; i++)
    {
      s += Math::Sqr(field.at(i));
    }
    // Scale by the size of a cell element
    s = sqrt(s) / size;

    return s;
  }
}
/*!
\brief Get the range of the field.

\avx

\param a,b Returned minimum and maximum.
*/
void ScalarField2::GetRange(double& a, double& b) const
{
  const int size = nx * ny;
  // Escape
  if (size == 0)
  {
    a = b = 0.0;
    return;
  }

#ifdef _MSC_VER
  if (System::Avx() && field.size() > 8)
  {
    const double* p = field.data();

    __m256d avxa = _mm256_load_pd(p);
    __m256d avxb = avxa;

    const int offset = size % 4;
    for (int i = 1; i < size / 4; i++)
    {
      __m256d avxi = _mm256_load_pd(p + i * 4);
      avxa = _mm256_min_pd(avxa, avxi);
      avxb = _mm256_max_pd(avxb, avxi);
    }
    double e[4];
    _mm256_store_pd(e, avxa);
    a = Math::Min(e[0], e[1], e[2], e[3]);

    _mm256_store_pd(e, avxb);
    b = Math::Max(e[0], e[1], e[2], e[3]);

    for (int i = size - offset; i < size; i++)
    {
      double x = field.at(i);
      if (x < a)
      {
        a = x;
      }
      else if (x > b)
      {
        b = x;
      }
    }
  }
  else
#endif
  {
    a = field.at(0);
    b = a;

    for (int i = 1; i < field.size(); i++)
    {
      double x = field.at(i);
      if (x < a)
      {
        a = x;
      }
      else if (x > b)
      {
        b = x;
      }
    }
  }
}


/*!
\brief Compute the histogram of the scalar field.
\param n Discretization.
*/
Histogram ScalarField2::GetHistogram(int n) const
{
  return Histogram(n, field);
}

/*!
\brief Compute the histogram of the scalar field.
\param a,b Interval.
\param n Discretization.
*/
Histogram ScalarField2::GetHistogram(int n, double a, double b) const
{
  return Histogram(n, a, b, field);
}

/*!
\brief Compute a histogram.

Compute the range of values taken by the scalar field, and compute the histogram.
\param n Discretization.
*/
QVector<std::pair<double, int>> ScalarField2::OldHistogram(int n) const
{
  QVector<std::pair<double, int>> histogram;
  double min, max;
  GetRange(min, max);
  double range = max - min;

  for (int i = 0; i < n; i++)
  {
    histogram.append(std::pair<double, int>((i + 1.0) * range / n, 0));
  }

  for (double val : field)
  {
    for (int i = 0; i < n; i++)
    {
      if (val <= histogram[i].first)
      {
        histogram[i].second++;
        break;
      }
    }
  }

  return histogram;
}

/*!
\brief Compute a cumulative histogram.

\sa Histogram
\param n Discretization.
*/
QVector<std::pair<double, int>> ScalarField2::CumulativeHistogram(int n) const
{
  QVector<std::pair<double, int>> histogram = OldHistogram(n);
  for (int i = 1; i < n; i++)
  {
    histogram[i].second += histogram[i - 1].second;
  }

  return histogram;
}

/*!
\brief Compute a normalized cumulative histogram.

\sa Histogram
\param n Discretization.
*/
QVector<std::pair<double, double>> ScalarField2::CumulativeNormedHistogram(int n) const
{
  QVector<std::pair<double, int>> histogram = CumulativeHistogram(n);
  QVector<std::pair<double, double>> histogramNormed;
  for (int i = 0; i < n; i++)
  {
    histogramNormed.append(std::pair<double, double>(histogram[i].first, (double)histogram[i].second / (double)field.size()));
  }
  return histogramNormed;
}

/*!
\brief Create an image from the field.
\param grayscale Export as grayscale if set to true, and to color otherwise (provides a better precision).
*/
QImage ScalarField2::CreateImage(bool grayscale) const
{
  double a, b;
  GetRange(a, b);

  // Check if the scalarfield is constant, in that case modify the range to avoid division by 0 errors
  if (a == b)
  {
    b = a + 1.0;
  }
  return CreateImage(a, b, grayscale);
}

/*!
\brief Create an image from the field using a palette.
\param palette The palette.
\param symmetric If set to true, adapt the range interval so that shading should be performed with median alignment.
\param transparent If set to true, the lowest values of the scalar field will be transparent.
*/
QImage ScalarField2::CreateImage(const AnalyticPalette& palette, bool symmetric, bool transparent) const
{
  double a, b;
  GetRange(a, b);
  // Check if the scalarfield is constant, in that case modify the range to avoid division by 0 errors
  if (a == b)
  {
    b = a + 1.0;
  }
  if (symmetric)
  {
    double m = Math::Max(fabs(a), fabs(b));
    a = -m;
    b = m;
  }

  return CreateImage(a, b, palette, transparent);
}

/*!
\brief Create an image from the field using a palette.
\param palette The palette.
*/
QImage ScalarField2::CreateImage(const Palette& palette) const
{
  double a, b;
  GetRange(a, b);
  // Check if the scalarfield is constant, in that case modify the range to avoid division by 0 errors
  if (a == b)
  {
    b = a + 1.0;
  }

  // Define image
  QImage image(nx, ny, QImage::Format_ARGB32);

  // Write Heightmap
  for (int i = 0; i < image.width(); i++)
  {
    for (int j = 0; j < image.height(); j++)
    {
      double x = field.at(VertexIndex(i, j));
      double y = Linear::Step(x, a, b);

      QColor color = palette.GetColor(y).GetQt();

      image.setPixel(i, j, color.rgb());
    }
  }

  return image;
}

/*!
\brief Create an image from the field using a palette.
\param palette The palette.
*/
QImage ScalarField2::CreateImage(const GenericPalette& palette) const
{
  double a, b;
  GetRange(a, b);
  return CreateImage(a, b, palette);
}


/*!
\brief Create an image from the field using a palette.
\param a,b Range of elevation that will be mapped to image scale.
\param palette The palette.
\param transparent If set to true, the lowest values of the scalar field will be transparent.
*/
QImage ScalarField2::CreateImage(const double& a, const double& b, const GenericPalette& palette, bool transparent) const
{
  // Define image
  QImage image(nx, ny, QImage::Format_ARGB32);

  // Write Heightmap
  for (int i = 0; i < image.width(); i++)
  {
    for (int j = 0; j < image.height(); j++)
    {
      double x = field.at(VertexIndex(i, j));
      double y = Linear::Step(x, a, b);

      QColor color = palette.GetColor(y).GetQt();

      // Transparency
      if (transparent)
      {
        if (y < 0.01)
        {
          color.setAlpha(0);
        }
      }
      image.setPixel(i, j, color.rgba());
    }
  }
  return image;
}


/*!
\brief Create an image from the field using a palette.
\param a,b Range of elevation that will be mapped to image scale.
\param palette The palette.
\param transparent If set to true, the lowest values of the scalar field will be transparent.
*/
QImage ScalarField2::CreateImage(const double& a, const double& b, const AnalyticPalette& palette, bool transparent) const
{
  // Define image
  QImage image(nx, ny, QImage::Format_ARGB32);

  // Write Heightmap
  for (int i = 0; i < image.width(); i++)
  {
    for (int j = 0; j < image.height(); j++)
    {
      double x = field.at(VertexIndex(i, j));
      double y = Linear::Step(x, a, b);

      QColor color = palette.GetColor(y).GetQt();

      // Transparency
      if (transparent)
      {
        if (y < 0.001)
        {
          color.setAlpha(0);
        }
      }
      image.setPixel(i, j, color.rgba());
    }
  }
  return image;
}

/*!
\brief Create an image from the field.
\param a,b Range of elevation that will be mapped to image scale.
\param grayscale Export as grayscale if set to true, color otherwise.
*/
QImage ScalarField2::CreateImage(const double& a, const double& b, bool grayscale) const
{
  // Write Heightmap
  if (grayscale)
  {
    QImage image(nx, ny, QImage::Format_Grayscale16);
    for (int j = 0; j < image.height(); j++)
    {
      quint16* dst = reinterpret_cast<quint16*>(image.bits() + j * image.bytesPerLine());
      for (int i = 0; i < image.width(); i++)
      {
        double x = field.at(VertexIndex(i, j));
        double y = Linear::Step(x, a, b);

        unsigned short pixelval = ushort(y * 65535.0);

        dst[i] = pixelval;
      }
    }
    return image;
  }
  else
  {
    // Define image
    QImage image(nx, ny, QImage::Format_ARGB32);

    // Write Heightmap
    for (int i = 0; i < image.width(); i++)
    {
      for (int j = 0; j < image.height(); j++)
      {
        double x = field.at(VertexIndex(i, j));
        double y = Linear::Step(x, a, b);

        QColor color;

        int c = int(y * (256.0 * 256.0 * 256.0 - 1.0));
        int cr = (c >> 16) & 255;
        int cv = (c >> 8) & 255;
        int cb = c & 255;

        color = QColor(cr, cv, cb);

        image.setPixel(i, j, color.rgb());
      }
    }
    return image;
  }

}

/*!
\brief Export the scalarfield as a PGM file in 16 bits.
\param url File name.
*/
bool ScalarField2::ExportPGM(const QString& url) const
{
  double a, b;
  GetRange(a, b);
  return ExportPGM(url, a, b);
}

/*!
\brief Export the scalarfield as a PGM file in 16 bits.
\param url File name.
\param a,b Interval range, values in the field will be clamped to [a,b] and rescaled 16 bits range.
*/
bool ScalarField2::ExportPGM(const QString& url, const double& a, const double& b) const
{
  QFile data(url);
  if (!data.open(QFile::WriteOnly))
  {
    return false;
  }

  QTextStream out(&data);
  out << "P2" << "\n";
  out << "# Created by LibCore" << "\n";
  out << nx << " " << ny << "\n";
  out << "65535" << "\n";

  int cpt = 0;

  for (int j = 0; j < ny; j++)
  {
    for (int i = 0; i < nx; i++)
    {
      double x = field.at(VertexIndex(i, j));
      double y = Linear::Step(x, a, b);
      int value = int(y * 65535.0);
      out << value << " ";
      cpt++;
      if (cpt % 10 == 0)
        out << "\n";
    }
  }

  out.flush();
  data.close();

  return true;
}

/*!
\brief Get the field value with world coordinate system.

This function computes a bi-linear interpolation of the values.

\sa Math::Bilinear
\param p Point position (should be strictly inside the box domain).
*/
double ScalarField2::Value(const Vector2& p) const
{
  double u, v;
  int i, j;
  CellInteger(p, i, j, u, v);

  // Test position
  if (!InsideCellIndex(i, j))
  {
    return 0.0;
  }

  return Math::Bilinear(at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1), u, v);
}

/*!
\brief Get the field value with world coordinate system.

This function computes a triangular interpolation of the values.

\param p Point (should be strictly within bounding box of the domain).
*/
double ScalarField2::Triangular(const Vector2& p) const
{
  double u, v;
  int i, j;
  CellInteger(p, i, j, u, v);

  if (!InsideCellIndex(i, j))
  {
    return 0.0;
  }

  if (u + v < 1.0)
  {
    return (1 - u - v) * at(i, j) + u * at(i + 1, j) + v * at(i, j + 1);
  }
  else
  {
    return (u + v - 1) * at(i + 1, j + 1) + (1 - v) * at(i + 1, j) + (1 - u) * at(i, j + 1);
  }
}

/*!
\brief Get the field value with world coordinate system.

This function computes a bicubic interpolation of the values, with null derivatives.

\sa Math::BiCubic

\param p Point (should be strictly within bounding box of the domain).
*/
double ScalarField2::BiCubic(const Vector2& p) const
{
  double u, v;
  int i, j;
  CellInteger(p, i, j, u, v);

  if (!InsideCellIndex(i, j))
  {
    return 0.0;
  }

  return Math::BiCubic(u, v, at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
}

/*!
\brief Compute the bicubic interpolation.
\author Lois Paulin.
\param p Point.
*/
double ScalarField2::BiCubicValue(const Vector2& p) const
{
  double u, v;
  int i, j;
  CellInteger(p, i, j, u, v);

  //if (!InsideCellIndex(i, j))
  //  return -1.0;

  double a00 = InsideVertexIndex(i, j) ? at(i, j) : 0.0;
  double a01 = InsideVertexIndex(i, j + 1) ? at(i, j + 1) : a00;
  double a10 = InsideVertexIndex(i + 1, j) ? at(i + 1, j) : a00;
  double a11 = InsideVertexIndex(i + 1, j + 1) ? at(i + 1, j + 1) : a00;

  double am10 = InsideVertexIndex(i - 1, j) ? at(i - 1, j) : a00;
  double am11 = InsideVertexIndex(i - 1, j + 1) ? at(i - 1, j + 1) : a01;

  double a0m1 = InsideVertexIndex(i, j - 1) ? at(i, j - 1) : a00;
  double a1m1 = InsideVertexIndex(i + 1, j - 1) ? at(i + 1, j - 1) : a10;

  double a02 = InsideVertexIndex(i, j + 2) ? at(i, j + 2) : a01;
  double a12 = InsideVertexIndex(i + 1, j + 2) ? at(i + 1, j + 2) : a11;

  double a20 = InsideVertexIndex(i + 2, j) ? at(i + 2, j) : a10;
  double a21 = InsideVertexIndex(i + 2, j + 1) ? at(i + 2, j + 1) : a11;

  double am1m1 = InsideVertexIndex(i - 1, j - 1) ? at(i - 1, j - 1) : am10;
  double am12 = InsideVertexIndex(i - 1, j + 2) ? at(i - 1, j + 2) : am11;
  double a2m1 = InsideVertexIndex(i + 2, j - 1) ? at(i + 2, j - 1) : a20;
  double a22 = InsideVertexIndex(i + 2, j + 2) ? at(i + 2, j + 2) : a21;

  am10 = Cubic::Interpolation(v, am1m1, am10, am11, am12);
  a00 = Cubic::Interpolation(v, a0m1, a00, a01, a02);
  a10 = Cubic::Interpolation(v, a1m1, a10, a11, a12);
  a20 = Cubic::Interpolation(v, a2m1, a20, a21, a22);

  return Cubic::Interpolation(u, am10, a00, a10, a20);
}

/*!
\brief Subdivide the scalar field.

Simply doubles the resolution and interpolate the values from the existing ones.
*/
void ScalarField2::Subdivide()
{
  // Extended scalar field
  ScalarField2 extended(Box2(a, b), 2 * nx - 1, 2 * ny - 1);

  // Copy old values at right location 
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      extended(2 * i, 2 * j) = at(i, j);
    }
  }

  // Interpolate : square
  for (int i = 0; i < nx - 1; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      // Term (i,j) maps to vertical interpolation 
      extended(2 * i + 1, 2 * j) = 0.5 * (at(i, j) + at(i + 1, j));
    }
  }
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny - 1; j++)
    {
      // Term (i,j) maps to horizontal interpolation 
      extended(2 * i, 2 * j + 1) = 0.5 * (at(i, j) + at(i, j + 1));
    }
  }

  // Interpolate : diamond
  for (int i = 1; i < extended.nx; i += 2)
  {
    for (int j = 1; j < extended.ny; j += 2)
    {
      extended(i, j) = 0.5 * (extended.at(i - 1, j) + extended.at(i + 1, j));
    }
  }

  // Update field
  *this = extended;
}

/*!
\brief Sets the entire field with a constant value.
\param s Scalar.
*/
void ScalarField2::Fill(const double& s)
{
  field.fill(s, nx * ny);
}

/*!
\brief Change the resolution of the scalar field.

Note that because the box should be strictly inside the original domain, this
function is not the same as:
\code
ScalarField field(Box2(Vector2(-2.0),Vector(3.0),5,5); // Original field
ScalarField s=Sample(field.GetBox(),12,12);
\endcode
\param x,y Sampling size.
\param bicubic Bicubic flag, set to true to use bicubic interpolation.
*/
ScalarField2 ScalarField2::SetResolution(int x, int y, bool bicubic) const
{
  // Sampled scalar field
  ScalarField2 sampled(Box2(a, b), x, y);

  // Corners
  sampled(0, 0) = at(0, 0);
  sampled(0, y - 1) = at(0, ny - 1);
  sampled(x - 1, 0) = at(nx - 1, 0);
  sampled(x - 1, y - 1) = at(nx - 1, ny - 1);

  if (bicubic == false)
  {
    // Borders (use linear interpolation)
    for (int i = 1; i < x - 1; i++)
    {
      double tx = (nx - 1) * (i / double(x - 1));
      int x0 = int(floor(tx));
      int x1 = int(ceil(tx));

      sampled(i, 0) = Math::Lerp(at(x0, 0), at(x1, 0), tx - x0);
      sampled(i, y - 1) = Math::Lerp(at(x0, ny - 1), at(x1, ny - 1), tx - x0);
    }

    for (int j = 1; j < y - 1; j++)
    {
      double ty = (ny - 1) * (j / double(y - 1));
      int y0 = int(floor(ty));
      int y1 = int(ceil(ty));

      sampled(0, j) = Math::Lerp(at(0, y0), at(0, y1), ty - y0);
      sampled(x - 1, j) = Math::Lerp(at(nx - 1, y0), at(nx - 1, y1), ty - y0);
    }

    // Interior
    for (int i = 1; i < x - 1; i++)
    {
      for (int j = 1; j < y - 1; j++)
      {
        sampled(i, j) = Value(sampled.ArrayVertex(i, j));
      }
    }
  }
  else
  {
    // Edges
    for (int i = 1; i < x - 1; i++)
    {
      double tx = (nx - 1) * (i / double(x - 1));
      int x0 = int(floor(tx));
      int x1 = int(ceil(tx));

      sampled(i, 0) = Math::Lerp(at(x0, 0), at(x1, 0), tx - x0);
      sampled(i, y - 1) = Math::Lerp(at(x0, ny - 1), at(x1, ny - 1), tx - x0);
    }

    for (int j = 1; j < y - 1; j++)
    {
      double ty = (ny - 1) * (j / double(y - 1));
      int y0 = int(floor(ty));
      int y1 = int(ceil(ty));

      sampled(0, j) = Math::Lerp(at(0, y0), at(0, y1), ty - y0);
      sampled(x - 1, j) = Math::Lerp(at(nx - 1, y0), at(nx - 1, y1), ty - y0);
    }

    // Interior
    for (int i = 1; i < x - 1; i++)
    {
      for (int j = 1; j < y - 1; j++)
      {
        sampled(i, j) = BiCubicValue(sampled.ArrayVertex(i, j));
      }
    }
  }
  return sampled;
}

/*!
\brief Resample a rectangular region in the scalar field.
\param box Rectangular region, should be strictly inside the domain.
\param x,y Sampling size.
\param bicubic Bicubic flag, set to true to use bicubic interpolation.
*/
ScalarField2 ScalarField2::Sample(const Box2& box, int x, int y, bool bicubic) const
{
  // Sampled scalar field
  ScalarField2 sampled(box, x, y);

  // Sample the domain
  for (int i = 0; i < x; i++)
  {
    for (int j = 0; j < y; j++)
    {
      if (bicubic == false)
      {
        sampled(i, j) = Value(sampled.ArrayVertex(i, j));
      }
      else
      {
        sampled(i, j) = BiCubicValue(sampled.ArrayVertex(i, j));
      }
    }
  }
  return sampled;
}


/*!
\brief Sample along a thick segment.
\param segment The segment.
\param radius Distance in the orthogonal direction to the segment.

Discretization is adapted from the resolution of the scalar field.
*/
ScalarField2 ScalarField2::Sample(const Segment2& segment, const double& radius) const
{
  const double Epsilon = 1.0e-4;
  double length = segment.Length() + Epsilon;
  int nx = int(length / celldiagonal[0]);

  int ny = int(nx * radius / length);
  return Sample(segment, radius, nx + 1, ny + 1);
}

/*!
\brief Sample along a thick segment.
\param segment The segment.
\param radius Distance in the orthogonal direction to the segment.
\param nx,ny Discretization.
*/
ScalarField2 ScalarField2::Sample(const Segment2& segment, const double& radius, int nx, int ny) const
{
  ScalarField2 s(Box2(segment.Length(), radius), nx, ny);

  Vector2 y = radius * segment.GetAxis().Orthogonal();

  for (int i = 0; i < nx; i++)
  {
    Vector2 p = segment.VertexAt(Math::Unit(i, nx));
    for (int j = 0; j < ny; i++)
    {
      Vector2 q = p + y * (2.0 * Math::Unit(j, ny) - 1.0);
      if (Inside(q))
      {
        s(i, j) = Value(q);
      }
    }
  }
  return s;
}

/*!
\brief Crops a rectangular region in the scalar field.

Note that the corners of the region of the scalar field are (0,0) and (nx-1,ny-1).

\param a,b Rectangular region.
*/
ScalarField2 ScalarField2::Crop(const QPoint& a, const QPoint& b) const
{
  // Vertices of the domain
  Vector2 va = ArrayVertex(a.x(), a.y());
  Vector2 vb = ArrayVertex(b.x(), b.y());

  // Size
  int x = b.x() - a.x() + 1;
  int y = b.y() - a.y() + 1;

  // Sampled scalar field
  ScalarField2 sampled(Box2(va, vb), x, y);

  // Copy
  for (int i = 0; i < x; i++)
  {
    for (int j = 0; j < y; j++)
    {
      sampled(i, j) = at(a.x() + i, a.y() + j);
    }
  }
  return sampled;
}

/*!
\brief Lower the resolution of the scalarfield.
\param f Lowering factor.
*/
ScalarField2 ScalarField2::DownSample(int f) const
{
  ScalarField2 res(Box2(a, b), nx / f, ny / f);

  for (int i = 0; i < nx / f; i++)
  {
    for (int j = 0; j < ny / f; j++)
    {
      double sij = 0.0;
      for (int ti = i * f; ti < i * f + f; ti++)
      {
        for (int tj = j * f; tj < j * f + f; tj++)
        {
          sij += at(ti, tj);
        }
      }
      res(i, j) = sij / (f * f);
    }
  }

  return res;
}

/*!
\brief Upsampling of the scalarfield.
\param f Upsampling factor.
\param bicubic Bicubic flag, set to true to use bicubic interpolation.
\sa ScalarField2::SetResolution
*/
ScalarField2 ScalarField2::UpSample(int f, bool bicubic) const
{
  return SetResolution(nx * f, ny * f, bicubic);
}

/*!
\brief Return the field gradient at a given array vertex.
\param i,j Integer coordinates of the array vertex.
\sa at(int,int)
*/
double ScalarField2::Value(int i, int j) const
{
  return field.at(VertexIndex(i, j));
}

/*!
\brief Compute the gradient at a given array vertex.

\param i,j Integer coordinates of the array vertex.
*/
Vector2 ScalarField2::Gradient(int i, int j) const
{
  Vector2 n;

  // Gradient along x axis
  if (i == 0)
  {
    n[0] = (at(i + 1, j) - at(i, j)) * inversecelldiagonal[0];
  }
  else if (i == nx - 1)
  {
    n[0] = (at(i, j) - at(i - 1, j)) * inversecelldiagonal[0];
  }
  else
  {
    n[0] = (at(i + 1, j) - at(i - 1, j)) * 0.5 * inversecelldiagonal[0];
  }

  // Gradient along y axis
  if (j == 0)
  {
    n[1] = (at(i, j + 1) - at(i, j)) * inversecelldiagonal[1];
  }
  else if (j == ny - 1)
  {
    n[1] = (at(i, j) - at(i, j - 1)) * inversecelldiagonal[1];
  }
  else
  {
    n[1] = (at(i, j + 1) - at(i, j - 1)) * 0.5 * inversecelldiagonal[1];
  }

  return n;
}

/*!
\brief Compute the smooth gradient at a given array vertex using a second order approximation.

\param i,j Integer coordinates of the array vertex.
*/
Vector2 ScalarField2::GradientSmooth(int i, int j) const
{
  Vector2 g;

  if (InsideVertexIndex(i, j, 2))
  {
    g[0] = -at(i - 2, j - 1) - 2.0 * at(i - 1, j - 1) + 2.0 * at(i + 1, j - 1) + at(i + 2, j - 1)
      - 2.0 * at(i - 2, j) - 4.0 * at(i - 1, j) + 4.0 * at(i + 1, j) + 2.0 * at(i + 2, j)
      - at(i - 2, j + 1) - 2.0 * at(i - 1, j + 1) + 2.0 * at(i + 1, j + 1) + at(i + 2, j + 1);

    g[1] = -at(i - 1, j - 2) - 2.0 * at(i - 1, j - 1) + 2.0 * at(i - 1, j + 1) + at(i - 1, j + 2)
      - 2.0 * at(i, j - 2) - 4.0 * at(i, j - 1) + 4.0 * at(i, j + 1) + 2.0 * at(i, j + 2)
      - at(i + 1, j - 2) - 2.0 * at(i + 1, j - 1) + 2.0 * at(i + 1, j + 1) + at(i + 1, j + 2);

    g *= inversecelldiagonal[0] / 32.0;
  }
  else
  {
    g = Gradient(i, j);
  }
  return g;
}

/*!
\brief Compute the local neighborhood in the one-ring aroung a point.
\param i, j Point.
*/
Matrix ScalarField2::Local(int i, int j) const
{
  // Borders
  int ia = (i == 0) ? 0 : i - 1;
  int ib = (i == nx - 1) ? nx - 1 : i + 1;

  int ja = (j == 0) ? 0 : j - 1;
  int jb = (j == ny - 1) ? ny - 1 : j + 1;

  // Local matrix
  Matrix a;
  a(0, 0) = at(ia, ja);
  a(1, 0) = at(i, ja);
  a(2, 0) = at(ib, ja);

  a(1, 0) = at(ia, j);
  a(1, 1) = at(i, j);
  a(1, 2) = at(ib, j);

  a(2, 0) = at(ia, jb);
  a(2, 1) = at(i, jb);
  a(2, 2) = at(ib, jb);

  return a;
}

/*!
\brief Compute the hessian at a given array vertex.

\param i,j Integer coordinates of the array vertex.
*/
Matrix2 ScalarField2::Hessian(int i, int j) const
{
  Vector2 d = celldiagonal.Scaled(celldiagonal);

  // Compute local neighborhood
  Matrix a = Local(i, j);

  return Matrix2(
    (a(2, 1) + a(0, 1) - 2.0 * a(1, 1)) / d[0],
    (a(1, 2) + a(1, 0) - 2.0 * a(1, 1)) / d[0],
    ((a(2, 2) - a(0, 2)) - (a(2, 0) - a(0, 0))) / (4.0 * d[0]));
}

/*!
\brief Compute the gradient field.
*/
VectorField2 ScalarField2::Gradient() const
{
  VectorField2 v(GetBox(), nx, ny);
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      v(i, j) = Gradient(i, j);
    }
  }
  return v;
}

/*!
\brief Compute the smooth gradient field.
*/
VectorField2 ScalarField2::GradientSmooth() const
{
  VectorField2 v(GetBox(), nx, ny);
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      v(i, j) = GradientSmooth(i, j);
    }
  }
  return v;
}

/*!
\brief Compute the %Laplacian at a given sample.

Returns a discrete approximation of Laplaceâ€™s differential operator applied to the field.
The implementation is equivalent to the 4*del2(U) Matlab function, except on the edges of the field.

\sa Laplacian()
\author Mathieu Gaillard
\param i,j Integer coordinates of the sample.
*/
double ScalarField2::Laplacian(int i, int j) const
{
  double laplacian = 0.0;

  // Divergence along x axis
  if (i == 0)
  {
    laplacian += (at(i, j) - 2.0 * at(i + 1, j) + at(i + 2, j)) / (celldiagonal[0] * celldiagonal[0]);
  }
  else if (i == nx - 1)
  {
    laplacian += (at(i, j) - 2.0 * at(i - 1, j) + at(i - 2, j)) / (celldiagonal[0] * celldiagonal[0]);
  }
  else
  {
    laplacian += (at(i + 1, j) - 2.0 * at(i, j) + at(i - 1, j)) / (celldiagonal[0] * celldiagonal[0]);
  }

  // Divergence along y axis
  if (j == 0)
  {
    laplacian += (at(i, j) - 2.0 * at(i, j + 1) + at(i, j + 2)) / (celldiagonal[1] * celldiagonal[1]);
  }
  else if (j == ny - 1)
  {
    laplacian += (at(i, j) - 2.0 * at(i, j - 1) + at(i, j - 2)) / (celldiagonal[1] * celldiagonal[1]);
  }
  else
  {
    laplacian += (at(i, j + 1) - 2.0 * at(i, j) + at(i, j - 1)) / (celldiagonal[1] * celldiagonal[1]);
  }

  return laplacian;
}

/*!
\brief Compute the %Laplacian field.

Returns a discrete approximation of Laplaceâ€™s differential operator applied to the field.
The implementation is equivalent to the 4*del2(U) Matlab function, except on the edges of the field.

\sa Laplacian(int, int)
\author Mathieu Gaillard
*/
ScalarField2 ScalarField2::Laplacian() const
{
  ScalarField2 laplacian_field(GetBox(), nx, ny);

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      laplacian_field(i, j) = Laplacian(i, j);
    }
  }

  return laplacian_field;
}

/*!
\brief Compute the gradient norm scalar field.
*/
ScalarField2 ScalarField2::GradientNorm() const
{
  // Scalar field of the same size
  ScalarField2 n(GetBox(), nx, ny);

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      n(i, j) = Norm(Gradient(i, j));
    }
  }
  return n;
}

/*!
\brief Compute the gradient norm scalar field.
*/
ScalarField2 ScalarField2::GradientSmoothNorm() const
{
  // Scalar field of the same size
  ScalarField2 n(GetBox(), nx, ny);

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      n(i, j) = Norm(GradientSmooth(i, j));
    }
  }
  return n;
}

/*!
\brief Compute the logarithm of the scalar field.
*/
ScalarField2 ScalarField2::Ln() const
{
  // Scalar field of the same size
  ScalarField2 r(GetBox(), nx, ny);

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      r(i, j) = log(at(i, j));
    }
  }
  return r;
}

/*!
\brief Compute the square root of the scalar field.

\avx
*/
ScalarField2 ScalarField2::Sqrt() const
{
  // Scalar field of the same size
  ScalarField2 r(GetBox(), nx, ny);

  const int size = nx * ny;

#ifdef _MSC_VER
  if (System::Avx())
  {
    double* sp = r.field.data();
    const double* p = field.data();

    const int offset = size % 4;
    for (int i = 0; i < size / 4; i++)
    {
      __m256d avx = _mm256_load_pd(p + i * 4);
      __m256d sqrtavx = _mm256_sqrt_pd(avx);
      _mm256_store_pd(sp + i * 4, sqrtavx);
    }

    for (int i = size - offset; i < size; i++)
    {
      r[i] = sqrt(at(i));
    }
  }
  else
#endif
  {

    for (int i = 0; i < size; i++)
    {
      r[i] = sqrt(at(i));
    }
  }

  return r;
}

/*!
\brief Compute the quartic root of the scalar field.

Direct quartic root computation is more efficicent that calling ScalarField2::Sqrt() twice.

\avx
*/
ScalarField2 ScalarField2::Qurt() const
{
  // Scalar field of the same size
  ScalarField2 r(GetBox(), nx, ny);

  const int size = nx * ny;

  if (System::Avx())
  {
    double* sp = r.field.data();
    const double* p = field.data();

    const int offset = size % 4;
    for (int i = 0; i < size / 4; i++)
    {
      __m256d avx = _mm256_load_pd(p + i * 4);
      __m256d sqrtavx = _mm256_sqrt_pd(avx);
      sqrtavx = _mm256_sqrt_pd(sqrtavx);
      _mm256_store_pd(sp + i * 4, sqrtavx);
    }

    for (int i = size - offset; i < size; i++)
    {
      r[i] = Math::Sqrt4(at(i));
    }
  }
  else
  {
    for (int i = 0; i < size; i++)
    {
      r[i] = Math::Sqrt4(at(i));
    }
  }

  return r;
}

/*!
\brief Compute the quartic root of the scalar field.

Direct quartic root computation is more efficicent that calling ScalarField2::Sqrt() twice.

\avx
*/
ScalarField2 ScalarField2::Cbrt() const
{
  // Scalar field of the same size
  ScalarField2 r(GetBox(), nx, ny);

  const int size = nx * ny;

  // not supported by GCC
#ifdef _MSC_VER
  if (System::Avx())
  {
    double* sp = r.field.data();
    const double* p = field.data();

    const int offset = size % 4;
    for (int i = 0; i < size / 4; i++)
    {
      __m256d avx = _mm256_load_pd(p + i * 4);
      __m256d cbrtavx = _mm256_cbrt_pd(avx);
      _mm256_store_pd(sp + i * 4, cbrtavx);
    }

    for (int i = size - offset; i < size; i++)
    {
      r[i] = cbrt(at(i));
    }
  }
  else
#endif
  {
    for (int i = 0; i < size; i++)
    {
      r[i] = cbrt(at(i));
    }
  }

  return r;
}

/*!
\brief Compute the square root of the scalar field.

This function changes the values of the  field.

\sa ScalarField2::Sqrt() const
\avx
*/
void ScalarField2::Sqrted()
{
  const int size = nx * ny;

#ifdef _MSC_VER
  if (System::Avx())
  {
    double* p = field.data();

    const int offset = size % 4;
    for (int i = 0; i < size / 4; i++)
    {
      __m256d avx = _mm256_load_pd(p + i * 4);
      __m256d sqrtavx = _mm256_sqrt_pd(avx);
      _mm256_store_pd(p + i * 4, sqrtavx);
    }

    for (int i = size - offset; i < size; i++)
    {
      field[i] = sqrt(at(i));
    }
  }
  else
#endif
  {
    for (int i = 0; i < size; i++)
    {
      field[i] = sqrt(at(i));
    }
  }
}

/*!
\brief Compute the Lipschitz constant of the elevation function.

Note that this is equivalent to the following code, with the
difference that this function does not store the norm of the gradient in memory.
\code
ScalarField2 n=f.GradientNorm();
double k=0.0;
for (int i=0;i<f.VertexSize();i++)
{
k=Math::Max(k,f.at(i));
}
\endcode
*/
double ScalarField2::K() const
{
  double k = 0.0;

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      k = Math::Max(k, Norm(Gradient(i, j)));
    }
  }
  return k;
}

/*!
\brief Normalize the values of a scalar field to unit interval.

Compute range using GetRange() and apply the affine transformation to map values to [0,1].

\sa Unitize()
\avx
*/
void ScalarField2::Normalize()
{
  const int size = nx * ny;
  const int offset = size % 4;

  // Escape
  if (size == 0)
    return;


  double a, b;
  GetRange(a, b);
  if (a == b)
  {
    field.fill(1.0);
    return;
  }
#ifdef _MSC_VER
  if (System::Avx())
  {
    double* p = field.data();

    __m256d avxiba = _mm256_set1_pd(1.0 / (b - a));
    __m256d avxa = _mm256_set1_pd(a);

    for (int i = 0; i < size / 4; i++)
    {
      __m256d avx = _mm256_load_pd(p + i * 4);
      // Could use _mm256_fmadd_ps
      avx = _mm256_sub_pd(avx, avxa);
      avx = _mm256_mul_pd(avx, avxiba);
      _mm256_store_pd(p + i * 4, avx);
    }

    for (int i = size - offset; i < size; i++)
    {
      field[i] = (field[i] - a) / (b - a);
    }
  }
  else
#endif
  {
    for (int i = 0; i < field.size(); i++)
    {
      field[i] = (field[i] - a) / (b - a);
    }
  }
}

/*!
\brief Linearly scale the values of a scalar field to unit interval.

This function differs from ScalarField2::Normalize() which maps values to [0,1]. Instead,
we divide scalard field values by the maximum of the absolute value of the scalar field,
which maps values to [-1,1] and preserves 0.

\sa Normalize()

\avx
*/
void ScalarField2::Unitize()
{
  double a, b;
  GetRange(a, b);
  if (a == b)
  {
    field.fill(1.0);
    return;
  }

  const double ab = 1.0 / Math::Max(fabs(a), fabs(b));

#ifdef _MSC_VER
  if (System::Avx())
  {
    double* p = field.data();

    const int size = nx * ny;
    const int offset = size % 4;

    __m256d avxiba = _mm256_set1_pd(ab);

    for (int i = 0; i < size / 4; i++)
    {
      __m256d avx = _mm256_load_pd(p + i * 4);
      // Could use _mm256_fmadd_ps
      avx = _mm256_mul_pd(avx, avxiba);
      _mm256_store_pd(p + i * 4, avx);
    }

    for (int i = size - offset; i < size; i++)
    {
      field[i] = field[i] * ab;
    }
  }
  else
#endif
  {
    for (int i = 0; i < field.size(); i++)
    {
      field[i] = field[i] * ab;
    }
  }
}

/*!
\brief Clamp the values of a scalar field.
\param a,b Interval.
*/
void ScalarField2::Clamp(const double& a, const double& b)
{
  for (int i = 0; i < field.size(); i++)
  {
    field[i] = Math::Clamp(field.at(i), a, b);
  }
}


/*!
\brief Compute the arctangent of the values of a scalar field, may be usefull for converting slope to angles.

\avx
*/
void ScalarField2::Atan()
{
  // not supported by gcc
#ifdef _MSC_VER
  if (System::Avx())
  {
    double* p = field.data();

    const int size = nx * ny;
    const int offset = size % 4;

    for (int i = 0; i < size / 4; i++)
    {
      __m256d avx = _mm256_load_pd(p + i * 4);
      // Could use _mm256_fmadd_ps
      avx = _mm256_atan_pd(avx);
      _mm256_store_pd(p + i * 4, avx);
    }

    for (int i = size - offset; i < size; i++)
    {
      field[i] = atan(field[i]);
    }
  }
  else
#endif
  {

    for (int i = 0; i < field.size(); i++)
    {
      field[i] = atan(field.at(i));
    }
  }
}

/*!
\brief Add a constant to the values the scalar field.
\param c Constant.
*/
ScalarField2& ScalarField2::operator+=(const double& c)
{
  for (int i = 0; i < field.size(); i++)
  {
    field[i] += c;
  }
  return *this;
}

/*!
\brief Scale the values of a scalar field.

\param s Scaling factor.
*/
ScalarField2& ScalarField2::operator*=(const double& s)
{
  for (int i = 0; i < field.size(); i++)
  {
    field[i] *= s;
  }
  return *this;
}

/*!
\brief Subtraction.
\param s Scalar field.

\avx
*/
ScalarField2& ScalarField2::operator-=(const ScalarField2& s)
{
#ifdef _MSC_VER
  if (System::Avx())
  {
    const double* sp = s.field.data();
    double* p = field.data();

    const int size = nx * ny;
    const int offset = size % 4;
    for (int i = 0; i < size / 4; i++)
    {
      __m256d avx = _mm256_load_pd(p + i * 4);
      __m256d savx = _mm256_load_pd(sp + i * 4);
      __m256d subavx = _mm256_sub_pd(avx, savx);
      _mm256_store_pd(p + i * 4, subavx);
    }

    for (int i = size - offset; i < size; i++)
    {
      field[i] -= s.field.at(i);
    }
  }
  else
#endif
  {
    int n = field.size();
    double* fieldp = field.data();
    const double* sfieldp = s.field.data();

    for (int i = 0; i < n; i++)
    {
      fieldp[i] -= sfieldp[i];
    }
  }
  return *this;
}

/*!
\brief Addition.
\param s Scalar field.
\avx
*/
ScalarField2& ScalarField2::operator+=(const ScalarField2& s)
{
#ifdef _MSC_VER
  if (System::Avx())
  {
    const double* sp = s.field.data();
    double* p = field.data();

    const int size = nx * ny;
    const int offset = size % 4;
    for (int i = 0; i < size / 4; i++)
    {
      __m256d avx = _mm256_load_pd(p + i * 4);
      __m256d savx = _mm256_load_pd(sp + i * 4);
      __m256d addavx = _mm256_add_pd(avx, savx);
      _mm256_store_pd(p + i * 4, addavx);
    }

    for (int i = size - offset; i < size; i++)
    {
      field[i] += s.field.at(i);
    }
  }
  else
#endif
  {
    int n = field.size();
    double* fieldp = field.data();
    const double* sfieldp = s.field.data();

    for (int i = 0; i < n; i++)
    {
      fieldp[i] += sfieldp[i];
    }
  }
  return *this;
}

/*!
\brief Multiplication.
\param s Scalar field.
*/
ScalarField2& ScalarField2::operator*=(const ScalarField2& s)
{
  for (int i = 0; i < field.size(); i++)
  {
    field[i] *= s.at(i);
  }
  return *this;
}

/*!
\brief Power.
\param s Real.
*/
void ScalarField2::Pow(const double& s)
{
  for (int i = 0; i < field.size(); i++)
  {
    field[i] = pow(field[i], s);
  }
}

/*!
\brief Inversion (take the opposite of each field value)
*/
void ScalarField2::Invert()
{
  for (int i = 0; i < field.size(); i++)
  {
    field[i] = -field[i];
  }
}

/*!
\brief Negate (invert 0 and 1 values of a binary field)
*/
void ScalarField2::Negate()
{
  for (int i = 0; i < field.size(); i++)
  {
    if (field[i] > 0.5)
      field[i] = 0.0;
    else
      field[i] = 1.0;
  }
}



/*!
\brief Perform a linear step over the values of a scalar field.

The function Linear::Step() is applied to all the values of the field.
\sa Linear::Step()
\sa ScalarField2::Normalize()

\param a,b Interval.
*/
void ScalarField2::Step(const double& a, const double& b)
{
  // The step could use _mm256_fmadd_pd with field.at(i), 1.0/(b-a) and -a/(b-a)
  // or _mm256_fadd_pd with -a and then _mm256_fmull_pd with 1.0/(b-a)
  // Finally, clamp to [0,1] by using _mm256_max_pd and _mm256_min_pd

  for (int i = 0; i < field.size(); i++)
  {
    field[i] = Linear::Step(field.at(i), a, b);
  }
}

/*!
\brief Compute the standard deviation.
\param m Mean.
*/
double ScalarField2::StandardDeviation(const double& m) const
{
  double s = 0.0;
#ifdef _MSC_VER
  if (System::Avx())
  {
    const double* p = field.data();

    const int size = nx * ny;
    const int offset = size % 4;

    __m256d avxm = _mm256_set1_pd(m);
    __m256d avxs = _mm256_setzero_pd();

    for (int i = 0; i < size / 4; i++)
    {
      __m256d avx = _mm256_load_pd(p + i * 4);
      avx = _mm256_sub_pd(avx, avxm);
      avx = _mm256_mul_pd(avx, avx);
      avxs = _mm256_add_pd(avxs, avx);
    }

    double e[4];
    _mm256_store_pd(e, avxs);
    s = e[0] + e[1] + e[2] + e[3];

    for (int i = size - offset; i < size; i++)
    {
      double x = field.at(i) - m;
      s += x * x;
    }
  }
  else
#endif
  {
    for (int i = 0; i < field.size(); i++)
    {
      double x = field.at(i) - m;
      s += x * x;
    }
  }

  return sqrt(s / field.size());
}

/*!
\brief Compute the standard deviation.

\avx
*/
double ScalarField2::StandardDeviation() const
{
  // Mean
  double m = Average();
  return StandardDeviation(m);
}

/*!
\brief Linear interpolation between two scalar fields.

\param a,b Scalar fields, should have the same resolution.
\param t Interpolant.
\avx
*/
void ScalarField2::Lerp(const ScalarField2& a, const ScalarField2& b, const double& t)
{
  // Coefficients
  const double u = 1.0 - t;
  const double v = t;

  // Size
  const int size = field.size();

  // Pointers to data
  double* p = field.data();
  const double* ap = a.field.data();
  const double* bp = b.field.data();

#ifdef _MSC_VER
  if (System::Avx())
  {
    __m256d avxu = _mm256_set1_pd(u);
    __m256d avxv = _mm256_set1_pd(v);

    const int offset = size % 4;
    for (int i = 0; i < size / 4; i++)
    {
      __m256d avxa = _mm256_load_pd(ap + i * 4);
      __m256d avxb = _mm256_load_pd(bp + i * 4);
      avxa = _mm256_mul_pd(avxa, avxu);
      avxb = _mm256_mul_pd(avxb, avxv);
      avxb = _mm256_add_pd(avxa, avxb);
      _mm256_store_pd(p + i * 4, avxb);
    }

    for (int i = size - offset; i < size; i++)
    {
      p[i] = u * ap[i] + v * bp[i];
    }
  }
  else
#endif
  {
    for (int i = 0; i < size; i++)
    {
      p[i] = u * ap[i] + v * bp[i];
    }
  }
}

/*!
\brief Interpolation between two scalar fields with an alpha scalar field.

Compute (1-t) a + t b.

\param a,b Scalar fields, should have the same resolution.
\param t Alpha.
\avx
*/
void ScalarField2::Lerp(const ScalarField2& a, const ScalarField2& b, const ScalarField2& t)
{
  // Size
  const int size = field.size();

  // Pointers to data
  double* p = field.data();
  const double* ap = a.field.data();
  const double* bp = b.field.data();
  const double* tp = t.field.data();

#ifdef _MSC_VER
  if (System::Avx())
  {
    // Constant
    __m256d avxu = _mm256_set1_pd(1.0);

    const int offset = size % 4;
    for (int i = 0; i < size / 4; i++)
    {
      __m256d avxa = _mm256_load_pd(ap + i * 4);
      __m256d avxb = _mm256_load_pd(bp + i * 4);

      // Coefficient t
      __m256d avxt = _mm256_load_pd(tp + i * 4);

      // Coefficient 1-t
      __m256d avxv = _mm256_sub_pd(avxu, avxt);

      avxa = _mm256_mul_pd(avxa, avxu);
      avxb = _mm256_mul_pd(avxb, avxv);
      avxb = _mm256_add_pd(avxa, avxb);
      _mm256_store_pd(p + i * 4, avxb);
    }

    for (int i = size - offset; i < size; i++)
    {
      p[i] = (1.0 - tp[i]) * ap[i] + tp[i] * bp[i];
    }
  }
  else
#endif
  {
    for (int i = 0; i < size; i++)
    {
      p[i] = (1.0 - tp[i]) * ap[i] + tp[i] * bp[i];
    }
  }
}

/*!
\brief Scales the scalar field to a given range interval.

\param a,b Interval.
*/
void ScalarField2::SetRange(const double& a, const double& b)
{
  double x, y;
  GetRange(x, y);
  if (x == y)
  {
    field.fill(a, field.size());
  }
  else
  {
    const double c = (b - a) / (y - x);
    for (int i = 0; i < field.size(); i++)
    {
      field[i] = a + c * (field.at(i) - x);
    }
  }
}

/*!
\brief Transform the field into a binary (0, 1 values) field using a threshold.

\param t Threshold.
*/
void ScalarField2::Binarize(const double& t)
{
  for (int i = 0; i < field.size(); i++)
  {
    if (field[i] >= t)
      field[i] = 1.0;
    else
      field[i] = 0.0;
  }
}

/*!
\brief Perform a smooth step over the values of a scalar field.

Either the function Cubic::Step() or Quintic::Step() is applied to all the values of the field.
\sa Cubic::Step(), Quintic::Step()
\param a,b Interval.
\param s Smoothing Boolean value, apply a quintic smoothing kernel if set to true.
*/
void ScalarField2::SmoothStep(const double& a, const double& b, bool s)
{
  for (int i = 0; i < field.size(); i++)
  {
    // Cubic smoothing
    if (s == false)
    {
      field[i] = Cubic::SmoothStep(field.at(i), a, b);

    }
    // Quintic smoothing
    else
    {
      field[i] = Quintic::SmoothStep(field.at(i), a, b);
    }
  }
}

/*!
\brief Compute the absolute value of a scalar field.
*/
void ScalarField2::Abs()
{
  for (int i = 0; i < field.size(); i++)
  {
    field[i] = fabs(field[i]);
  }
}

/*!
\brief Translate the domain of the scalar field.

\param t Translation vector.
*/
void ScalarField2::Translate(const Vector2& t)
{
  Box2::Translate(t);
}

/*!
\brief Scale the domain of the scalar field.

\param s Scaling factor.
*/
void ScalarField2::Scale(const Vector2& s)
{
  Array2::Scale(s);
}

/*!
\brief Scale the scalar field.

This function also scales the values, it is different from ScalarField2::Scale(const Vector2&);

\param s Scaling factor.
*/
void ScalarField2::Scale(const double& s)
{
  Array2::Scale(s);
  for (int i = 0; i < nx * ny; i++)
  {
    field[i] *= s;
  }
}

/*!
\brief Add material with gaussian distribution.
\param center Center of the distribution.
\param radius Radius
\param height Maximum height of the distribution.
*/
void ScalarField2::Gaussian(const Vector2& center, const double& radius, const double& height)
{
  // Compute modification Area
  const QRect area = VertexIntegerArea(Circle2(center, radius));

  // Compute thickness
  for (int y = area.y(); y <= area.y() + area.height(); y++)
  {
    for (int x = area.x(); x <= area.x() + area.width(); x++)
    {
      // Distance between central point and current point
      double u = SquaredNorm(center - ArrayVertex(x, y));

      if (u < radius * radius)
      {
        field[VertexIndex(x, y)] += height * Cubic::Smooth(u, radius * radius);
      }
    }
  }
}

/*!
\brief Add a value to the scalar field with diffusion.

The value is distributed among the neighboring vertices.
\param p Point.
\param x Value.
*/
void ScalarField2::Add(const Vector2& p, const double& x)
{
  double u, v;
  int i, j;
  CellInteger(p, i, j, u, v);

  if (!InsideCellIndex(i, j))
    return;

  field[VertexIndex(i, j)] += x * (1.0 - u) * (1.0 - v);
  field[VertexIndex(i + 1, j)] += x * u * (1.0 - v);
  field[VertexIndex(i + 1, j + 1)] += x * u * v;
  field[VertexIndex(i, j + 1)] += x * (1.0 - u) * v;
}

/*!
\brief Copy the values of the edges to a new field.

\param s New scalar field.
*/
void ScalarField2::CopyEdges(ScalarField2& s) const
{
  for (int i = 0; i < nx; i++)
  {
    s(i, 0) = at(i, 0);
  }
  for (int i = 0; i < nx; i++)
  {
    s(i, ny - 1) = at(i, ny - 1);
  }
  for (int j = 1; j < ny - 1; j++)
  {
    s(0, j) = at(0, j);
  }
  for (int j = 1; j < ny - 1; j++)
  {
    s(nx - 1, j) = at(nx - 1, j);
  }
}

/*!
\brief Compute the median filtered scalar field.

Edges are not processed.
*/
ScalarField2 ScalarField2::MedianFilter() const
{
  ScalarField2 filtered(GetArray());

  CopyEdges(filtered);
  //   Move window through all elements of the image
  for (int x = 1; x < nx - 1; x++)
  {
    for (int y = 1; y < ny - 1; y++)
    {
      //   Pick up window elements
      int k = 0;
      double window[9];
      for (int j = x - 1; j < x + 2; j++)
      {
        for (int i = y - 1; i < y + 2; i++)
        {
          window[k++] = at(i, j);
        }
      }
      //   Order elements (only half of them)
      for (int j = 0; j < 5; ++j)
      {
        //   Find position of minimum element
        int min = j;
        for (int l = j + 1; l < 9; ++l)
        {
          if (window[l] < window[min])
            min = l;
        }

        Math::Swap(window[j], window[min]);
      }
      // Get result - the middle element
      filtered(x, y) = window[4];
    }
  }
  return filtered;
}

/*!
\brief Flatten the scalar field.
\param center Center.
\param radius Radius.
\param scaling Scaling of the effect.
*/
void ScalarField2::Flatten(const Vector2& center, const double& radius, const double& scaling)
{
  // Compute modification Area
  const QRect area = VertexIntegerArea(Circle2(center, radius));

  // Elevation
  double pz = Value(center);

  for (int y = area.y(); y <= area.y() + area.height(); y++)
  {
    for (int x = area.x(); x <= area.x() + area.width(); x++)
    {
      // Distance between central point and current point
      double u = SquaredNorm(center - ArrayVertex(x, y));

      if (u < radius * radius)
      {
        double a = Cubic::Smooth(u, radius * radius);

        // We modify the field proportionnaly to the height difference
        field[VertexIndex(x, y)] += (pz - field[VertexIndex(x, y)]) * a * scaling;
      }
    }
  }
}

/*!
\brief Multiply the scalar field by a scaling factor locally (not equivalent to Flatten())
\param center Center.
\param radius Radius.
\param scaling Scaling of the effect.
*/
void ScalarField2::Multiply(const Vector2& center, const double& radius, const double& scaling)
{
  // Compute modification Area
  const QRect area = VertexIntegerArea(Circle2(center, radius));

  for (int y = area.y(); y <= area.y() + area.height(); y++)
  {
    for (int x = area.x(); x <= area.x() + area.width(); x++)
    {
      // Distance between central point and current point
      double u = SquaredNorm(center - ArrayVertex(x, y));

      if (u < radius * radius)
      {
        double a = Cubic::Smooth(u, radius * radius);

        // We modify the field proportionnaly to the height difference
        field[VertexIndex(x, y)] *= a * scaling + (1 - a); // interpolates scaling and 1
      }
    }
  }
}

/*!
\brief Level the scalar field.

\param center Center.
\param radius Radius.
\param z Elevation.
*/
void ScalarField2::Level(const Vector2& center, const double& radius, const double& z)
{
  // Compute modification Area
  const QRect area = VertexIntegerArea(Circle2(center, radius));

  // Compute thickness
  for (int y = area.y(); y <= area.y() + area.height(); y++)
  {
    for (int x = area.x(); x <= area.x() + area.width(); x++)
    {
      // Distance between central point and current point
      double u = SquaredNorm(center - ArrayVertex(x, y));

      if (u < radius * radius)
      {
        field[VertexIndex(x, y)] = Math::Lerp(field[VertexIndex(x, y)], z, Cubic::Smooth(u, radius * radius));
      }
    }
  }
}

/*!
\brief Create the set of scalar points.
*/
QVector<ScalarPoint2> ScalarField2::GetScalarPoints() const
{
  // Create array
  QVector<ScalarPoint2> e(nx * ny);

  int k = 0;
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      e[k++] = ScalarPoint2(QPoint(i, j), at(i, j));
    }
  }
  return e;
}

/*!
\brief Compute and return an array of 4 byte floats.
*/
FloatArray ScalarField2::GetAsFloats() const
{
  FloatArray a(nx * ny);
  for (int i = 0; i < nx * ny; i++)
  {
    a[i] = at(i);
  }
  return a;
}

/*!
\brief Compute the iso-contour of the scalar field.
\param T Threshold.
\param closed Set to true if polylines should be closed at the border of the domain.
*/
SegmentSet2 ScalarField2::LineSegments(const double& T, bool closed) const
{
  // Vertices
  QVector<Vector2> v;
  v.reserve(nx * ny / 4);
  int nv = 0;

  // Edges
  QVector<int> eax(nx * ny, -1);
  QVector<int> eay(nx * ny, -1);

  // Compute straddling edges 
  for (int i = 0; i < nx - 1; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      // Different signs
      if (!Math::SameSign(at(i, j) - T, at(i + 1, j) - T))
      {
        v.append(Dichotomy(ArrayVertex(i, j), ArrayVertex(i + 1, j), at(i, j) - T, at(i + 1, j) - T, celldiagonal[0], T, epsilon));
        eax[VertexIndex(i, j)] = nv;
        nv++;
      }
    }
  }

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny - 1; j++)
    {
      // Different signs
      if (!Math::SameSign(at(i, j) - T, at(i, j + 1) - T))
      {
        v.append(Dichotomy(ArrayVertex(i, j), ArrayVertex(i, j + 1), at(i, j) - T, at(i, j + 1) - T, celldiagonal[1], T, epsilon));
        eay[VertexIndex(i, j)] = nv;
        nv++;
      }
    }
  }

  // Array for edge vertices
  QVector<int> edges;

  // There may be no more than nv segments
  edges.reserve(2 * nv);

  // Create edges for interior cells
  for (int i = 0; i < nx - 1; i++)
  {
    for (int j = 0; j < ny - 1; j++)
    {
      int squareindex = 0;
      if (at(i, j) - T > 0.0) squareindex |= 1;
      if (at(i + 1, j) - T > 0.0) squareindex |= 2;
      if (at(i, j + 1) - T > 0.0) squareindex |= 4;
      if (at(i + 1, j + 1) - T > 0.0) squareindex |= 8;

      switch (squareindex)
      {
      case 0:
        break;
      case 1:
        edges.append(eay[VertexIndex(i, j)]);
        edges.append(eax[VertexIndex(i, j)]);
        break;
      case 2:
        edges.append(eay[VertexIndex(i + 1, j)]);
        edges.append(eax[VertexIndex(i, j)]);
        break;
      case 3:
        edges.append(eay[VertexIndex(i + 1, j)]);
        edges.append(eay[VertexIndex(i, j)]);
        break;
      case 4:
        edges.append(eay[VertexIndex(i, j)]);
        edges.append(eax[VertexIndex(i, j + 1)]);
        break;
      case 5:
        edges.append(eax[VertexIndex(i, j)]);
        edges.append(eax[VertexIndex(i, j + 1)]);
        break;
      case 6:
        edges.append(eay[VertexIndex(i, j)]);
        edges.append(eax[VertexIndex(i, j + 1)]);
        edges.append(eay[VertexIndex(i + 1, j)]);
        edges.append(eax[VertexIndex(i, j)]);
        break;
      case 7:
        edges.append(eay[VertexIndex(i + 1, j)]);
        edges.append(eax[VertexIndex(i, j + 1)]);
        break;
      case 8:
        edges.append(eax[VertexIndex(i, j + 1)]);
        edges.append(eay[VertexIndex(i + 1, j)]);
        break;
      case 9:
        edges.append(eax[VertexIndex(i, j)]);
        edges.append(eay[VertexIndex(i, j)]);
        edges.append(eax[VertexIndex(i, j + 1)]);
        edges.append(eay[VertexIndex(i + 1, j)]);
        break;
      case 10:
        edges.append(eax[VertexIndex(i, j + 1)]);
        edges.append(eax[VertexIndex(i, j)]);
        break;
      case 11:
        edges.append(eax[VertexIndex(i, j + 1)]);
        edges.append(eay[VertexIndex(i, j)]);
        break;
      case 12:
        edges.append(eay[VertexIndex(i, j)]);
        edges.append(eay[VertexIndex(i + 1, j)]);
        break;
      case 13:
        edges.append(eax[VertexIndex(i, j)]);
        edges.append(eay[VertexIndex(i + 1, j)]);
        break;
      case 14:
        edges.append(eay[VertexIndex(i, j)]);
        edges.append(eax[VertexIndex(i, j)]);
        break;
      case 15:
        break;
      };
    }
  }

  // Escape if no boundary
  if (closed == false)
    return SegmentSet2(v, edges);

  // Border vertexes
  const int nb = VertexBorderSize();

  QVector<int> vxy(nb);

  for (int k = 0; k < nb; k++)
  {
    if (at(VertexBorder(k)) - T > 0.0)
    {
      v.append(ArrayVertex(VertexBorder(k).x(), VertexBorder(k).y()));
      vxy[k] = nv;
      nv++;
    }
  }

  // Create edges for boundary

  // Horizontal 
  for (int i = 0; i < nx - 1; i++)
  {
    // Straddling edge
    if (eax[VertexIndex(i, 0)] >= 0)
    {
      if (at(i, 0) - T > 0.0)
      {
        edges.append(vxy[VertexBorderIndex(i, 0)]);
        edges.append(eax[VertexIndex(i, 0)]);
      }
      else
      {
        edges.append(eax[VertexIndex(i, 0)]);
        edges.append(vxy[VertexBorderIndex(i + 1, 0)]);
      }
    }
    else
    {
      if ((at(i, 0) - T > 0.0) && (at(i + 1, 0) - T > 0.0))
      {
        edges.append(vxy[VertexBorderIndex(i, 0)]);
        edges.append(vxy[VertexBorderIndex(i + 1, 0)]);
      }
    }
  }

  for (int i = 0; i < nx - 1; i++)
  {
    // Straddling edge
    if (eax[VertexIndex(i, ny - 1)] >= 0)
    {
      if (at(i, ny - 1) - T > 0.0)
      {
        edges.append(vxy[VertexBorderIndex(i, ny - 1)]);
        edges.append(eax[VertexIndex(i, ny - 1)]);
      }
      else
      {
        edges.append(eax[VertexIndex(i, ny - 1)]);
        edges.append(vxy[VertexBorderIndex(i + 1, ny - 1)]);
      }
    }
    else
    {
      if ((at(i, ny - 1) - T > 0.0) && (at(i + 1, ny - 1) - T > 0.0))
      {
        edges.append(vxy[VertexBorderIndex(i, ny - 1)]);
        edges.append(vxy[VertexBorderIndex(i + 1, ny - 1)]);
      }
    }
  }

  // Vertical 
  for (int j = 0; j < ny - 1; j++)
  {
    // Straddling edge
    if (eay[VertexIndex(0, j)] >= 0)
    {
      if (at(0, j) - T > 0.0)
      {
        edges.append(vxy[VertexBorderIndex(0, j)]);
        edges.append(eay[VertexIndex(0, j)]);
      }
      else
      {
        edges.append(eay[VertexIndex(0, j)]);
        edges.append(vxy[VertexBorderIndex(0, j + 1)]);
      }
    }
    else
    {
      if ((at(0, j) - T > 0.0) && (at(0, j + 1) - T > 0.0))
      {
        edges.append(vxy[VertexBorderIndex(0, j)]);
        edges.append(vxy[VertexBorderIndex(0, j + 1)]);
      }
    }
  }

  for (int j = 0; j < ny - 1; j++)
  {
    // Straddling
    if (eay[VertexIndex(nx - 1, j)] >= 0)
    {
      if (at(nx - 1, j) - T > 0.0)
      {
        edges.append(vxy[VertexBorderIndex(nx - 1, j)]);
        edges.append(eay[VertexIndex(nx - 1, j)]);
      }
      else
      {
        edges.append(eay[VertexIndex(nx - 1, j)]);
        edges.append(vxy[VertexBorderIndex(nx - 1, j + 1)]);
      }
    }
    else
    {
      if ((at(nx - 1, j) - T > 0.0) && (at(nx - 1, j + 1) - T > 0.0))
      {
        edges.append(vxy[VertexBorderIndex(nx - 1, j)]);
        edges.append(vxy[VertexBorderIndex(nx - 1, j + 1)]);
      }
    }
  }
  return SegmentSet2(v, edges);

}

/*!
\brief Compute the intersection between a segment and the scalar field.

\param a,b End vertices of the segment straddling the surface.
\param va,vb Field function value at those end vertices.
\param length Distance between vertices.
\param epsilon Precision.
\param T Threshold value.
\return Point.
*/
Vector2 ScalarField2::Dichotomy(Vector2 a, Vector2 b, double va, double vb, double length, const double& T, const double& epsilon) const
{
  int ia = va > 0 ? 1 : -1;

  // Get an accurate first guess
  Vector2 c = Vector2::Solve(a, b, va, vb);

  while (length > epsilon)
  {
    double vc = Value(c) - T;
    int ic = vc > 0.0 ? 1 : -1;
    if (ia + ic == 0)
    {
      b = c;
    }
    else
    {
      ia = ic;
      a = c;
    }
    length *= 0.5;
    c = 0.5 * (a + b);
  }
  return c;
}

/*!
\brief Compute the memory size of a scalarfield.
*/
int ScalarField2::Memory() const
{
  return int(sizeof(ScalarField2)) + int(sizeof(double)) * field.size();
}

/*!
\brief Applies a function to the scalar field.

\param function The function.

\code
ScalarField2 s;
s.SelfOp([&](double x) { return x-pow(x, 0.7); });
\endcode
*/
void ScalarField2::SelfOp(const std::function<double(double)>& function)
{
  const int size = nx * ny;
  for (int k = 0; k < size; k++)
  {
    field[k] = function(field[k]);
  }
}

/*!
\brief Applies a function to the scalar field.

\param function The function.
\return A scalarfield with computed values.

\sa ScalarField2::SelfOp(const std::function<double(double)>&)
*/
ScalarField2 ScalarField2::Op(const std::function<double(double)>& function) const
{
  ScalarField2 ret(GetBox(), nx, ny);

  const int size = nx * ny;
  for (int k = 0; k < size; k++)
  {
    ret[k] = function(field[k]);
  }

  return ret;
}

/*!
\brief Applies a function to the scalar field.

\param function The function.
\return A scalarfield with computed values.
*/
ScalarField2 ScalarField2::Op(const std::function<double(const ScalarField2&, int, int)>& function) const
{
  ScalarField2 ret(GetBox(), nx, ny);

  for (int j = 0; j < nx; j++)
  {
    for (int i = 0; i < ny; i++)
    {
      ret(i, j) = function(*this, i, j);
    }
  }

  return ret;
}

/*!
\brief Compute PSNR and RMSE between two fields.

PSNR is computed assuming the data has a normalized [0,1] range.

\param s Second scalar field.
\param psnr Returned PSNR.
\param rmse Returned RMSE.
*/
void ScalarField2::SignalError(const ScalarField2& s, double& psnr, double& rmse) const
{
  double mse;
  SignalError(s, psnr, mse, rmse);
}

/*!
\brief Compute PSNR, MSE and RMSE between two fields.

\param s Second scalar field.
\param psnr Returned PSNR.
\param mse  Returned MSE.
\param rmse Returned RMSE.
*/
void ScalarField2::SignalError(const ScalarField2& s, double& psnr, double& mse, double& rmse) const
{
  double x = 0.0;
  for (int k = 0; k < field.size(); k++)
  {
    x += Math::Sqr(field.at(k) - s.field.at(k));
  }
  mse = x / field.size();
  rmse = sqrt(mse);
  psnr = 20.0 * log10(1.0 / rmse);
}

static double* dt1d(double* f, int n)
{
  double* d = new double[n];
  int* v = new int[n];
  double* z = new double[n + 1];
  int k = 0;

  v[0] = 0;
  z[0] = -Math::Infinity;
  z[1] = Math::Infinity;
  for (int q = 1; q <= n - 1; q++)
  {
    double s = ((f[q] + q * q) - (f[v[k]] + v[k] * v[k])) / double(2 * q - 2 * v[k]);
    while (s <= z[k])
    {
      k--;
      s = ((f[q] + q * q) - (f[v[k]] + v[k] * v[k])) / double(2 * q - 2 * v[k]);
    }
    k++;
    v[k] = q;
    z[k] = s;
    z[k + 1] = Math::Infinity;
  }

  k = 0;
  for (int q = 0; q <= n - 1; q++)
  {
    while (z[k + 1] < q)
    {
      k++;
    }
    d[q] = (q - v[k]) * (q - v[k]) + f[v[k]];
  }

  delete[] v;
  delete[] z;
  return d;
}

/*
\brief Compute the squared distance field from a scalar field.

Points with values below the threshold value are considered as inside, others as outside.
\param threshold Threshold value defining inside and outside.
*/
ScalarField2 ScalarField2::DistanceTransform(const double& threshold) const
{
  double* f = new double[max(nx, ny)];
  ScalarField2 dist(GetArray());

  for (int x = 0; x < nx; x++)
  {
    for (int y = 0; y < ny; y++)
    {
      dist(x, y) = at(x, y) < threshold ? at(x, y) : Math::Large;
    }
  }

  // Transform along columns
  for (int x = 0; x < nx; x++)
  {
    for (int y = 0; y < ny; y++)
    {
      f[y] = dist(x, y);
    }
    double* d = dt1d(f, ny);
    for (int y = 0; y < ny; y++)
    {
      dist(x, y) = d[y];
    }
    delete[] d;
  }

  // Transform along rows
  for (int y = 0; y < ny; y++)
  {
    for (int x = 0; x < nx; x++)
    {
      f[x] = dist(x, y);
    }
    double* d = dt1d(f, nx);
    for (int x = 0; x < nx; x++)
    {
      dist(x, y) = d[x];
    }
    delete[] d;
  }

  return dist;
}

/*!
\brief Compute the max filtered scalar field.
\param w Window size.
*/
ScalarField2 ScalarField2::MaxFilter(int w) const
{
  ScalarField2 filtered(GetArray());

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      double v = at(i, j);
      for (int di = -w; di <= w; di++)
      {
        for (int dj = -w; dj <= w; dj++)
        {
          if (InsideVertexIndex(i + di, j + dj))
          {
            v = Math::Max(v, at(i + di, j + dj));
          }
        }
      }
      filtered(i, j) = v;
    }
  }

  return filtered;
}

/*!
\brief Compute the min filtered scalar field.
\param w Window size.
*/
ScalarField2 ScalarField2::MinFilter(int w) const
{
  ScalarField2 filtered(GetArray());

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      double v = at(i, j);
      for (int di = -w; di <= w; di++)
      {
        for (int dj = -w; dj <= w; dj++)
        {
          if (InsideVertexIndex(i + di, j + dj))
          {
            v = Math::Min(v, at(i + di, j + dj));
          }
        }
      }
      filtered(i, j) = v;
    }
  }

  return filtered;
}

/*!
\brief Perform a Gaussian blur.
\param r Radius of blur.

\sa ScalarField2::GaussianBlur(int)
*/
ScalarField2 ScalarField2::GaussianBlur(const double& r) const
{
  int n = r / Array2::celldiagonal[0];
  return GaussianBlur(n);
}

/*!
\brief Perform a Gaussian blur.
\param r Radius (in cells) of blur.
*/
ScalarField2 ScalarField2::GaussianBlur(int r) const
{
  // Fixed sigma to half kernel radius
  double sigma = 0.5 * double(r);
  double sigma2 = sigma * sigma;
  double norm = 1.0 / (sqrt(Math::TwoPi) * sigma);

  // Kernel
  std::vector<double> kernel(r + 1);
  for (int i = 0; i <= r; i++)
  {
    kernel[i] = norm * exp((-0.5 * double(i) * double(i)) / sigma2);
  }

  QVector<double> resx(nx * ny);

  // Filter x direction
  for (int y = 0; y < ny; y++)
  {
    for (int x = 0; x < nx; x++)
    {
      double v = kernel[0] * at(x, y);
      double w = kernel[0];
      for (int k = 1; k <= r; k++)
      {
        if (InsideVertexIndex(x + k, y))
        {
          v += kernel[k] * at(x + k, y);
          w += kernel[k];
        }
        if (InsideVertexIndex(x - k, y))
        {
          v += kernel[k] * at(x - k, y);
          w += kernel[k];
        }
      }
      resx[VertexIndex(x, y)] = v / w;
    }
  }

  ScalarField2 res(Box2(a, b), nx, ny);

  // Filter y direction
  for (int x = 0; x < nx; x++)
  {
    for (int y = 0; y < ny; y++)
    {
      double v = kernel[0] * resx[VertexIndex(x, y)];
      double w = kernel[0];
      for (int k = 1; k <= r; k++)
      {
        if (InsideVertexIndex(x, y + k))
        {
          v += kernel[k] * resx[VertexIndex(x, y + k)];
          w += kernel[k];
        }
        if (InsideVertexIndex(x, y - k))
        {
          v += kernel[k] * resx[VertexIndex(x, y - k)];
          w += kernel[k];
        }
      }
      res(x, y) = v / w;
    }
  }

  return res;
}


/*!
\brief Polygonize the scalar field into a surface mesh using the marching square algorithm.

\return The mesh.
*/
Mesh2 ScalarField2::Polygonize() const
{
  QVector<Vector2> vertex;
  QVector<int> triangle;
  vertex.reserve(nx * ny / 16 + 4);
  triangle.reserve(nx * ny / 16 + 4);

  int* eax = new int[nx - 1];
  int* ebx = new int[nx - 1];
  int* ecy = new int[nx];

  int* evax = new int[nx];
  int* evbx = new int[nx];

  int nv = 0;

  // Compute grid vertexes along lower horizontal edge
  for (int i = 0; i < nx; i++)
  {
    if (at(i, 0) > 0.0)
    {
      vertex.append(ArrayVertex(i, 0));
      evax[i] = nv;
      nv++;
    }
  }

  // Compute vertex along lower horizontal edge
  for (int i = 0; i < nx - 1; i++)
  {
    if (!((at(i, 0) < 0.0) == !(at(i + 1, 0) >= 0.0)))
    {
      vertex.append(Vector2::Solve(ArrayVertex(i, 0), ArrayVertex(i + 1, 0), at(i, 0), at(i + 1, 0)));
      eax[i] = nv;
      nv++;
    }
  }

  for (int j = 0; j < ny - 1; j++)
  {
    // Compute grid vertexes along upper horizontal edge
    for (int i = 0; i < nx; i++)
    {
      if (at(i, j + 1) > 0.0)
      {
        vertex.append(ArrayVertex(i, j + 1));
        evbx[i] = nv;
        nv++;
      }
    }
    // Compute vertex along upper horizontal edge
    for (int i = 0; i < nx - 1; i++)
    {
      if (!((at(i, j + 1) < 0.0) == !(at(i + 1, j + 1) >= 0.0)))
      {
        vertex.append(Vector2::Solve(ArrayVertex(i, j + 1), ArrayVertex(i + 1, j + 1), at(i, j + 1), at(i + 1, j + 1)));
        ebx[i] = nv;
        nv++;
      }
    }
    // Compute vertex along vertical edges
    for (int i = 0; i < nx; i++)
    {
      if (!((at(i, j) < 0.0) == !(at(i, j + 1) >= 0.0)))
      {
        vertex.append(Vector2::Solve(ArrayVertex(i, j), ArrayVertex(i, j + 1), at(i, j), at(i, j + 1)));
        ecy[i] = nv;
        nv++;
      }
    }

    // Configuration
    for (int i = 0; i < nx - 1; i++)
    {
      int squareindex = 0;
      if (at(i, j) > 0.0)
      {
        squareindex |= 1;
      }
      if (at(i + 1, j) > 0.0)
      {
        squareindex |= 2;
      }
      if (at(i, j + 1) > 0.0) {
        squareindex |= 4;
      }
      if (at(i + 1, j + 1) > 0.0)
      {
        squareindex |= 8;
      }

      // Create mesh
      switch (squareindex)
      {
      case 0:
        break;
      case 1:

        break;
      case 2:
        break;
      case 3:
        break;
      case 4:
        break;
      case 5:
        break;
      case 6:
        break;
      case 7:
        break;
      case 8:
        break;
      case 9:
        break;
      case 10:
        break;
      case 11:
        break;
      case 12:
        break;
      case 13:
        break;
      case 14:
        triangle.append(eax[i]);
        triangle.append(evax[i + 1]);
        triangle.append(evbx[i + 1]);
        triangle.append(evax[i + 1]);
        triangle.append(evbx[i + 1]);
        triangle.append(ecy[i]);
        triangle.append(evax[i + 1]);
        triangle.append(evbx[i + 1]);
        triangle.append(evbx[i]);
        break;
      case 15:
        triangle.append(evax[i]);
        triangle.append(evax[i + 1]);
        triangle.append(evbx[i]);
        triangle.append(evax[i + 1]);
        triangle.append(evbx[i + 1]);
        triangle.append(evbx[i]);
        break;
      };

    }
    Math::Swap(eax, ebx);
    Math::Swap(evax, evbx);
  }

  delete[]eax;
  delete[]ebx;
  delete[]ecy;

  return Mesh2(vertex, triangle);
}

/*!
\brief Compute the summed area table of a scalar field.
*/
ScalarField2 ScalarField2::SummedAreaTable() const
{
  ScalarField2 sat(*this);

  for (int i = 1; i < nx; i++)
  {
    sat(i, 0) += sat(i - 1, 0);
  }

  for (int j = 1; j < ny; j++)
  {
    sat(0, j) += sat(0, j - 1);
  }

  for (int i = 1; i < nx; i++)
  {
    for (int j = 1; j < ny; j++)
    {
      sat(i, j) += sat(i - 1, j) + sat(i, j - 1) - sat(i - 1, j - 1);
    }
  }
  return sat;
}

/*!
\brief Overloaded.
\param s Stream.
\param scalar The scalar field.
*/
std::ostream& operator<<(std::ostream& s, const ScalarField2& scalar)
{
  s << Array2(scalar) << std::endl;

  for (int j = scalar.ny - 1; j > -1; j--)
  {
    for (int i = 0; i < scalar.nx; i++)
    {
      s << scalar.at(i, j) << ' ';
    }
    s << std::endl;
  }
  return s;
}

/*!
\brief Adaptive bilateral filtering.
\param spatial_sigma Width (in number of cells) of spatial gaussian kernel.
\param range_sigma Intensity gaussian kernel.
\author Hugo Schott
*/
void ScalarField2::AdaptiveBilateralFiltering(int spatial_sigma, double range_sigma)
{
  ScalarField2 new_map = ScalarField2(*this, 0.0);

  // Window size
  int w = 3 * spatial_sigma;

  // Inverse squared values
  double is = 1.0 / (2.0 * spatial_sigma * spatial_sigma);
  double ir = 1.0 / (2.0 * range_sigma * range_sigma);

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      QPoint q = QPoint(i, j);

      double norm = 0.0;
      double value = 0.0;

      for (int x = -w; x < w; x++)
      {
        for (int y = -w; y < w; y++)
        {
          QPoint p = QPoint(q.x() + x, q.y() + y);
          if (!InsideVertexIndex(p)) continue;

          // Spatial
          double spatial_distance = x * x + y * y;
          double spatial_coeff = exp(-spatial_distance * is);

          // Range
          double range_distance = at(q) - at(p);
          double range_coeff = exp(-range_distance * range_distance * ir);

          double coeff = spatial_coeff * range_coeff;
          norm += coeff;
          value += coeff * at(p);
        }
      }

      if (norm == 0.) norm = 1.0;
      norm = 1. / norm;

      new_map(q) = norm * value;
    }
  }

  *this = new_map;
}

/*!
\brief Adaptive bilateral filtering.
\param w Width of spatial gaussian kernel.
\param s Intensity gaussian kernel.
\author Hugo Schott
*/
void ScalarField2::AdaptiveBilateralFiltering(double w, double s)
{
  AdaptiveBilateralFiltering(int(w * inversecelldiagonal[0]), s);
}