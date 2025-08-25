// Fields

#include "libs/arrayinteger.h"
#include "libs/cubic.h"

/*!
\class Array2I scalarfield.h
\brief A base two-dimensional array of integers.

\ingroup StructureGroup
*/

/*!
\brief Create the field structure.
\param a The array.
\param v Default value of field.
*/
Array2I::Array2I(const Array2& a, int v) :Array2(a)
{
  field.fill(v, nx * ny);
}

/*!
\brief Create the field structure.
\param box The box.
\param x,y Size of the array.
\param v Default value of field.
*/
Array2I::Array2I(const Box2& box, int x, int y, int v) :Array2I(Array2(box, x, y), v)
{
}

/*!
\brief Create the field structure.
\param box The box.
\param x,y Size of the array.
\param v Array of scalar values.
*/
Array2I::Array2I(const Box2& box, int x, int y, const QVector<int>& v) :Array2(box, x, y)
{
  field = v;
}

/*!
\brief Create a field structure from an image.
\param box The box.
\param image The image.
\param grayscale Read image as grayscale if set to true, otherwize use color for a better accuracy.
*/
Array2I::Array2I(const Box2& box, const QImage& image, bool grayscale) :Array2(box, image.width(), image.height())
{
  // Set size
  field.resize(nx * ny);

  // Write Heightmap
  for (int i = 0; i < image.width(); i++)
  {
    for (int j = 0; j < image.height(); j++)
    {
      double t;
      QRgb color = image.pixel(i, j);

      // Greyscale
      if (grayscale)
      {
        t = qGray(color);
      }
      // Color 
      else
      {
        t = double(qRed(color) << 16 | qGreen(color) << 8 | qBlue(color));
      }
      field[VertexIndex(i, j)] = t;
    }
  }
}

/*!
\brief Create a field structure from a scalar field.
\param s The scalar field.
\param a,b Integer range, the scalar field will be normalized and scale to interval [a,b].
*/
Array2I::Array2I(const ScalarField2& s, int a, int b):Array2(s)
{
  double u, v;
  s.GetRange(u, v);
  int size = nx * ny;
  field.resize(size);

  if (u == v)
  {
    field.fill(a, size);
    return;
  }

  for (int i = 0; i < size; i++)
  {
    field[i] = a + int((b - a) * Linear::Step(s.at(i), u, v));
  }
}

/*!
\brief Empty.
*/
Array2I::~Array2I()
{
}

/*!
\brief Get the range of the field.
\param a,b Returned minimum and maximum.
*/
void Array2I::GetRange(int& a, int& b) const
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

/*!
\brief Sets the entire field with a constant value.
\param s Scalar.
*/
void Array2I::Fill(int s)
{
  field.fill(s, nx * ny);
}

/*!
\brief Return the field gradient at a given array vertex.
\param i,j Integer coordinates of the vertex.
\sa at(int,int)
*/
int Array2I::Value(int i, int j) const
{
  return field.at(VertexIndex(i, j));
}

/*!
\brief Translate the domain of the scalar field.

\param t Translation vector.
*/
void Array2I::Translate(const Vector2& t)
{
  Box2::Translate(t);
}

/*!
\brief Scale the domain of the scalar field.

\param s Scaling factor.
*/
void Array2I::Scale(const Vector2& s)
{
  Array2::Scale(s);
}

/*!
\brief Create an image from the field.
\param grayscale Export as grayscale if set to true, and to color otherwise (provides a better precision).
*/
QImage Array2I::CreateImage(bool grayscale) const
{
  int a, b;
  GetRange(a, b);

  return CreateImage(a, b, grayscale);
}

QImage Array2I::CreateImage(const LookupPalette& palette) const
{
  // Define image
  QImage image(nx, ny, QImage::Format_ARGB32);

  // Write Heightmap
  for (int i = 0; i < image.width(); i++)
  {
    for (int j = 0; j < image.height(); j++)
    {
      image.setPixel(i, j, palette.GetColor(at(i, j)).GetQt().rgb());
    }
  }
  return image;
}

/*!
\brief Create an image from the field.
\param a,b Range of elevation that will be mapped to image scale.
\param grayscale Export as grayscale if set to true, color otherwise.
*/
QImage Array2I::CreateImage(int a, int b, bool grayscale) const
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
      // Greyscale: dynamic range limited
      if (grayscale)
      {
        int c = int(y * 255.0);
        color = QColor(c, c, c);
      }
      else
      {
        int c = int(y * (256.0 * 256.0 * 256.0 - 1.0));
        int red = (c >> 16) & 255;
        int green = (c >> 8) & 255;
        int blue = c & 255;

        color = QColor(red, green, blue);
      }
      image.setPixel(i, j, color.rgb());
    }
  }
  return image;
}

/*
\brief Cropping.
\param a,b Points defining the cropping region.
*/
Array2I Array2I::Crop(const QPoint& a, const QPoint& b) const
{
  // Vertices of the domain
  Vector2 va = ArrayVertex(a);
  Vector2 vb = ArrayVertex(b);

  // Size
  int x = b.x() - a.x() + 1;
  int y = b.y() - a.y() + 1;

  // Sampled integer array
  Array2I sampled(Box2(va, vb), x, y);

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
\brief Overloaded.
\param s Stream.
\param a The array.
*/
std::ostream& operator<<(std::ostream& s, const Array2I& a)
{
  s << Array2(a) << std::endl;

  for (int j = a.ny - 1; j > -1; j--)
  {
    for (int i = 0; i < a.nx; i++)
    {
      s << a.at(i, j) << ' ';
    }
    s << std::endl;
  }
  return s;
}
