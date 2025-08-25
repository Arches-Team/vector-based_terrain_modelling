// Vector fields
#include <QtGui/QPainter>
#include <QtWidgets/QGraphicsScene>
#include <immintrin.h>

#include "libs/vectorfield.h"
#include "libs/triangle.h"
#include "libs/cubic.h"
#include "libs/curveset.h"
#include "libs/curve.h"
#include "libs/cpu.h"
#include "libs/sampling.h"
#include "libs/segment.h"
#include "libs/curvepoint.h"


/*!
\class VectorField2 vectorfield.h
\brief A base two-dimensional field of Vector2 values.

\ingroup StructureGroup
*/

/*!
\brief Create the field structure.
\param box The box.
\param x,y Size of the array.
\param v Default value of field.
*/
VectorField2::VectorField2(const Box2& box, int x, int y, const Vector2& v) :Array2(box, x, y)
{
  field.fill(v, nx * ny);
}

/*!
\brief Create the field structure from an analytic vector field.
\param box The box.
\param x,y Size of the array.
\param a The analytic vector field.
*/
VectorField2::VectorField2(const Box2& box, int x, int y, const AnalyticVectorField2& a) :Array2(box, x, y)
{
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      field[VertexIndex(i, j)] = a.Value(ArrayVertex(i, j));
    }
  }
}

/*!
\brief Compute the orthogonal vector field.

Note that this is the same as rotating the field, but more efficient:
\code
VectorField2 v;
v.Rotate(Matrix2(Math::Pi/2)).
\endcode
*/
VectorField2 VectorField2::Orthogonal() const
{
  VectorField2 r(Box2(a, b), nx, ny);
  for (int i = 0; i < nx * ny; i++)
  {
    r[i] = field[i].Orthogonal();
  }
  return r;
}

/*!
\brief Compute the rotated vector field.
\param r %Rotation matrix.
*/
void VectorField2::Rotate(const Matrix2& r)
{
  for (int i = 0; i < nx * ny; i++)
  {
    field[i] = r * field[i];
  }
}

/*!
\brief Linear interpolation between two vector fields.

\param a,b Scalar fields, should have the same resolution.
\param t Interpolant.
\avx
*/
void VectorField2::Lerp(const VectorField2& a, const VectorField2& b, const double& t)
{
  // Size
  const int size = field.size();

  // Coefficients
  const double u = 1.0 - t;
  const double v = t;

  // Pointers to data
  Vector2* p = field.data();
  const Vector2* ap = a.field.data();
  const Vector2* bp = b.field.data();

  for (int i = 0; i < size; i++)
  {
    p[i] = u * ap[i] + v * bp[i];
  }
}

/*!
\brief Multiply the vector field by a scaling factor locally.

At the center, the scaling will be applied and at the boundary it will be equal to 1.
\param center Center.
\param radius Radius.
\param scaling Scaling of the effect.
*/
void VectorField2::Multiply(const Vector2& center, const double& radius, const double& scaling)
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
\brief Adds to the vector field a radial component centered at center and goind towards the center (with scaling being positive) and towards the exterior of the center (with scaling being negative)
\param center Center.
\param radius Radius.
\param scaling Scaling of the effect.
*/
void VectorField2::Radial(const Vector2& center, const double& radius, const double& scaling)
{
  // Compute modification Area
  const QRect area = VertexIntegerArea(Circle2(center, radius));

  for (int y = area.y(); y <= area.y() + area.height(); y++)
  {
    for (int x = area.x(); x <= area.x() + area.width(); x++)
    {
      // Distance between central point and current point
      Vector2 direction = center - ArrayVertex(x, y);
      double u = Norm(direction);
      direction /= u;

      if (u < radius)
      {
        double a = Cubic::Smooth(abs(0.5 * radius - u), 0.5 * radius);

        // We modify the field proportionnaly to the height difference
        field[VertexIndex(x, y)] += scaling * a * direction;
      }
    }
  }
}

/*!
\brief Clone part of the vector field from source to destination positions, with a given radius. The destination is blended with a decreasing weight towards the disk boundary where it is null.
\param center Destination position.
\param source Source position.
\param strokebegin Stroke starting position.
\param radius Scaling of the effect.
*/
void VectorField2::Clone(const Vector2& center, const Vector2& source, const Vector2& strokebegin, const double& radius)
{
  // Compute modification Area
  const QRect area = VertexIntegerArea(Circle2(center, radius));
  QPoint sourcexy = VertexInteger(source);
  int x0 = sourcexy.x();
  int y0 = sourcexy.y();
  QPoint destxy = VertexInteger(center);

  QPoint beginxy = VertexInteger(strokebegin);
  int xb = beginxy.x();
  int yb = beginxy.y();

  for (int y = area.y(); y <= area.y() + area.height(); y++)
  {
    for (int x = area.x(); x <= area.x() + area.width(); x++)
    {
      // Distance between central point and current point
      double u = SquaredNorm(center - ArrayVertex(x, y));

      if (u < radius * radius)
      {
        double a = Cubic::Smooth(u, radius * radius);
        a = sqrt(a);
        // We modify the field proportionnaly to the height difference
        field[VertexIndex(x, y)] = (1. - a) * field[VertexIndex(x, y)] + a * field[VertexIndex(x0 + x - xb, y0 + y - yb)]; // blend
      }
    }
  }
}

/*!
\brief Clone part of the vector field from source to destination positions.

The destination is blended with a decreasing weight towards the disk boundary where it is null.
\param center Target position.
\param lassoedZone Selected area.
\param borderBlendTreshold Border threshold.
\param gradientSource Source.
*/
void VectorField2::CloneZone(const Vector2& center, const QPolygon& lassoedZone, int borderBlendTreshold, const VectorField2* gradientSource)
{
  const QRect area = lassoedZone.boundingRect();
  const QPoint sourceCorner = area.topLeft();
  const QPoint destCorner(center[0] - area.width() / 2, center[1] - area.height() / 2);

  ScalarField2 distField(Array2(Box2(Vector2(0, 0), Vector2(area.width(), area.height())), area.width(), area.height()), 0.0);
  for (int y = 0; y <= area.height(); y++)
  {
    for (int x = 0; x <= area.width(); x++)
    {
      if (lassoedZone.containsPoint(sourceCorner + QPoint(x, y), Qt::OddEvenFill))
      { //the grid contains 0 inside the polygon, and 0 outsides 
        distField(x, y) = 1.0;
      }
    }
  }
  //computes distanceField
  distField = distField.DistanceTransform(0.5);   //the bounds of the polygone is where it goes from 0 to 1, so 0.5 is the limit
  distField.Sqrted();

  distField.SetRange(0.0, 1.0);   //set all the distances to border between 0 and 1
  double blendThreshold = (double)borderBlendTreshold / 100.0; //if the point is less than distToBlend to the border, then it is blended with the original gradient, else it is just cloned as is

  for (int y = 0; y <= area.height(); y++)   //relative position of the point in the bounding rect of the zone to clone
  {
    for (int x = 0; x <= area.width(); x++)
    {
      if (lassoedZone.containsPoint(sourceCorner + QPoint(x, y), Qt::OddEvenFill))
      {  //absolute position of the point to clone
        QPoint sourcePoint;  //the points coordonates in the field object

        if (gradientSource != nullptr) {
          sourcePoint = gradientSource->VertexInteger(Vector2(double(sourceCorner.x() + x), double(sourceCorner.y() + y)));
        }
        else
        {
          sourcePoint = VertexInteger(Vector2(sourceCorner.x() + (double)x, sourceCorner.y() + (double)y));
        }
        QPoint destPoint = VertexInteger(Vector2(destCorner.x() + (double)x, destCorner.y() + (double)y));

        if (InsideVertexIndex(destPoint) && ((gradientSource != nullptr) ? gradientSource->InsideVertexIndex(sourcePoint) : InsideVertexIndex(sourcePoint)))
        {
          double blendCoeff = ((distField(x, y) <= blendThreshold) ? (distField(x, y) / blendThreshold) : 1.0);
          //if the point is close enough to the border, we'll blend its gradient with the new gradient's value, else we just take the new gradients value
          blendCoeff = Cubic::Smooth(blendCoeff * blendCoeff);
          blendCoeff = sqrt(blendCoeff);
          if (gradientSource != nullptr)
          {
            field[VertexIndex(destPoint.x(), destPoint.y())] = (1. - blendCoeff) * field[VertexIndex(destPoint.x(), destPoint.y())] + (blendCoeff)*gradientSource->field[gradientSource->VertexIndex(sourcePoint.x(), sourcePoint.y())]; //blend
          }
          else
          {
            field[VertexIndex(destPoint.x(), destPoint.y())] = (1. - blendCoeff) * field[VertexIndex(destPoint.x(), destPoint.y())] + (blendCoeff)*field[VertexIndex(sourcePoint.x(), sourcePoint.y())]; //blend
          }
        }
      }
    }
  }
}

/*!
\brief clones the gradient from the sourceLine along the destLine, cloning all points that are less than halfWidth away from the curves
*/
void VectorField2::CloneBrokenLine(const QVector<Vector2>& sourceLine, const QVector<Vector2>& destLine, int halfWidth, const VectorField2* gradientSource)
{
  QuadricCurve2Set sourceCurve(sourceLine);
  QuadricCurve2Set destCurve(destLine);
  double ratio = sourceCurve.GetLength() / destCurve.GetLength(); //the 2 curve probably don't have the same length
  //QRect destRect = destCurve.GetBox().GetQtRect().toRect();
  for (int ab = 0; ab < destCurve.GetLength() * 10; ab++)
  {
    //for each point of the destCurve
    double a = double(ab) / 10.0;
    for (int d = -halfWidth; d < halfWidth; d++)
    {      //for every point on the normal to the curve in a, who is less than halfwidth away from the curve
      double u;   //abscisse curviligne of the sub-curve of the set, where the point of the curveSet with the curviligne abscisse "a" is 

      int i = destCurve.U(a, u);//the index of the curve of the Set where the projection is
      Vector2 onCurvePoint = destCurve(i)(u);
      Vector2 tangente = destCurve(i).Tangent(u);
      Vector2 normale = Normalized(tangente.Orthogonal());
      Vector2 point = onCurvePoint + normale * d;   //point to get the clone gradient value from

      int si; //i but for the source curve
      double su;  //u but for the source curve
      si = sourceCurve.U(a * ratio, su);
      Vector2 onSourceCurvePoint = sourceCurve(si)(su);
      Vector2 sourceTangente = sourceCurve(si).Tangent(su);
      Vector2 sourceNormale = Normalized(sourceTangente.Orthogonal());
      Vector2 pointToClone = onSourceCurvePoint + sourceNormale * d;   //point to get the clone gradient value from

      //convert points into vertex for field
      QPoint sourcePoint = (gradientSource != nullptr) ? gradientSource->VertexInteger(pointToClone) : VertexInteger(pointToClone);
      QPoint destPoint = VertexInteger(point);

      if (InsideVertexIndex(destPoint) && ((gradientSource != nullptr) ? gradientSource->InsideVertexIndex(sourcePoint) : InsideVertexIndex(sourcePoint)))
      {
        //we need to rotate the gradient we copy by the difference of angles of normale and sourceNormale to match the new orientation
        double angle = sourceTangente.Angle(tangente);
        Matrix2 rotation = Matrix2::Rotation(angle);

        double blendCoeff = (double)(abs(d)) / halfWidth; //blend old and new gradient
        blendCoeff = Cubic::SmoothStep(blendCoeff, 0.4, 1.);
        /*blendCoeff = Cubic::Smooth(blendCoeff * blendCoeff);
        blendCoeff = sqrt(blendCoeff);*/
        if (gradientSource != nullptr)
        {
          field[VertexIndex(destPoint.x(), destPoint.y())] = blendCoeff * field[VertexIndex(destPoint.x(), destPoint.y())] + (1.0 - blendCoeff) * rotation * gradientSource->field[gradientSource->VertexIndex(sourcePoint.x(), sourcePoint.y())];
        }
        else
        {
          field[VertexIndex(destPoint.x(), destPoint.y())] = blendCoeff * field[VertexIndex(destPoint.x(), destPoint.y())] + (1.0 - blendCoeff) * rotation * field[VertexIndex(sourcePoint.x(), sourcePoint.y())];
        }
      }
    }
  }
}

/*!
\brief Draw a vector field.
\param scene Graphics scene.
\param pen The pen.
\param brush The brush.
*/
void VectorField2::Draw(QGraphicsScene& scene, const QPen& pen, const QBrush& brush) const
{
  Vector2 va, vb;
  GetRange(va, vb);
  double s = Math::Max(Norm(va), Norm(vb));

  // Field
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      Vector2 v = at(i, j);
      if (v != Vector2::Null)
      {
        Vector2 p = Array2::ArrayVertex(i, j);
        Vector2 q = (v / s).Scaled(celldiagonal);
        Segment2(p, p + q).DrawArrow(scene, 0.1 * NormInfinity(celldiagonal), pen, brush);
      }
    }
  }
}


/*!
\brief Create an image from the field, showing trajectory with particle and streamlines.

\author Axel Paris.

\param dt Time step for tracing particle trajectories.
*/
QImage VectorField2::CreateImageParticles(const double& dt) const
{
  Box2 box = Box2(a, b);

  QVector<QuadricCurve2> curves;

  DiscTile poisson = DiscTile(box.Size()[0], box.Size()[0] / 150.0, box.Width() * 2.0);

  for (int i = 0; i < poisson.Size(); i++)
  {
    QVector<Vector2> pts;
    Vector2 p = poisson.Vertex(i);
    int j = 0;
    while (box.Inside(p) && j < 10)
    {
      pts.push_back(p);
      Vector2 w = Value(p);
      if (w == Vector2::Null)
        break;
      p = p + w * dt;
      j++;
    }
    if (pts.size() <= 1)
      continue;
    curves.append(QuadricCurve2::Bezier(pts));
  }

  QGraphicsScene s;
  s.addPixmap(QPixmap::fromImage(GetNorm().CreateImage(AnalyticPalette(3))));
  for (int i = 0; i < curves.size(); i++)
  {
    curves[i].Draw(s, QPen(QColor(255, 255, 255, 120)));
  }
  return System::Rasterize(s, box, box.Width());
}

/*!
\brief Create an image from the field.
*/
QImage VectorField2::CreateImage() const
{
  QGraphicsScene gs(0, 0, nx + 1, ny + 1);
  QPen pen;
  pen.setColor(QColor(255, 0, 0));
  pen.setWidthF(0.1);

  // Boundingbox
  Array2::GetBox().Draw(gs, pen);

  // Field
  Draw(gs, pen);

  QImage im(2048, 2048, QImage::Format_ARGB32);

  QPainter p(&im);
  gs.render(&p, QRectF(0, 0, 2048, 2048), gs.itemsBoundingRect());

  return im;
}


/*!
\brief Create an image from the field.
*/
QImage VectorField2::CreateStreamLines() const
{
  // Compute norm
  ScalarField2 n = GetNorm();
  double sa, sb;

  // Get range
  n.GetRange(sa, sb);
  if (sa == sb) { sb += 1.0; }

  // Crate base image with norm
  QImage image = n.CreateImage(sa, sb, AnalyticPalette(5));

  double t = Norm(celldiagonal);

  for (int i = 0; i < nx * ny / 16; i++)
  {
    Vector2 p = RandomInside();

    QPoint q = VertexInteger(p);
    Vector2 speed = Value(p) / sb;
    Vector2 pj = p;
    for (int j = 0; j < 50; j++)
    {
      pj += speed * t;
      QPoint qj = VertexInteger(pj);
      if (qj == q)
      {
        continue;
      }
      if (!InsideVertexIndex(qj)) { break; }
      Color background = Color(image.pixel(q));
      Color c = Color::Lerp(Norm(speed) * Math::Unit(j, 50), background, background + Color(0.2, 0.5, 0.5));
      image.setPixelColor(q, c.GetQt());
      q = qj;
    }
  }
  return image;
}

/*!
\brief Return the field vector at a given array vertex.
\param i,j Integer coordinates of the vertex.
\sa at(int,int)
*/
Vector2 VectorField2::Value(int i, int j) const
{
  return field.at(VertexIndex(i, j));
}

/*!
\brief Compute the divergence field.
\sa Divergence(int, int)
\author Mathieu Gaillard
*/
ScalarField2 VectorField2::Divergence() const
{
  ScalarField2 d(GetBox(), nx, ny);

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      d(i, j) = Divergence(i, j);
    }
  }

  return d;
}

/*!
\brief Compute the divergence at a given sample.
\sa Divergence()
\author Mathieu Gaillard
\param i,j Integer coordinates of the sample.
*/
double VectorField2::Divergence(int i, int j) const
{
  double div = 0.0;

  // Divergence along x axis
  if (i == 0)
  {
    div += (at(i + 1, j)[0] - at(i, j)[0]) * inversecelldiagonal[0];
  }
  else if (i == nx - 1)
  {
    div += (at(i, j)[0] - at(i - 1, j)[0]) * inversecelldiagonal[0];
  }
  else
  {
    div += (at(i + 1, j)[0] - at(i - 1, j)[0]) * 0.5 * inversecelldiagonal[0];
  }

  // Divergence along y axis
  if (j == 0)
  {
    div += (at(i, j + 1)[1] - at(i, j)[1]) * inversecelldiagonal[1];
  }
  else if (j == ny - 1)
  {
    div += (at(i, j)[1] - at(i, j - 1)[1]) * inversecelldiagonal[1];
  }
  else
  {
    div += (at(i, j + 1)[1] - at(i, j - 1)[1]) * 0.5 * inversecelldiagonal[1];
  }

  return div;
}

/*!
\brief Add two vector fields.
\param v Argument vector field.
\avx
*/
void VectorField2::Add(const VectorField2& v)
{
#ifdef _MSC_VER
  if (System::Avx())
  {
    const double* sp = (double*)(v.field.data());
    double* p = (double*)(field.data());

    const int size = 2 * nx * ny;
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
      p[i] += sp[i];
    }
  }
  else
#endif
  {
    int n = 2 * field.size();
    double* fieldp = (double*)(field.data());
    const double* sfieldp = (double*)(v.field.data());

    for (int i = 0; i < n; i++)
    {
      fieldp[i] += sfieldp[i];
    }
  }
}

/*!
\brief Subtract two vector fields.
\param v Argument vector field.
\avx
*/
void VectorField2::Sub(const VectorField2& v)
{
#ifdef _MSC_VER
  if (System::Avx())
  {
    const double* sp = (double*)(v.field.data());
    double* p = (double*)(field.data());

    const int size = 2 * nx * ny;
    const int offset = size % 4;
    for (int i = 0; i < size / 4; i++)
    {
      __m256d avx = _mm256_load_pd(p + i * 4);
      __m256d savx = _mm256_load_pd(sp + i * 4);
      __m256d addavx = _mm256_sub_pd(avx, savx);
      _mm256_store_pd(p + i * 4, addavx);
    }

    for (int i = size - offset; i < size; i++)
    {
      p[i] -= sp[i];
    }
  }
  else
#endif
  {
    int n = 2 * field.size();
    double* fieldp = (double*)(field.data());
    const double* sfieldp = (double*)(v.field.data());

    for (int i = 0; i < n; i++)
    {
      fieldp[i] -= sfieldp[i];
    }
  }
}

/*!
\brief Subtract two vector fields.
\param s Scalar.
\avx
*/
void VectorField2::Mul(const double& s)
{
#ifdef _MSC_VER
  if (System::Avx())
  {
    double* p = (double*)(field.data());

    const int size = 2 * nx * ny;
    const int offset = size % 4;

    const double ps[4] = { s,s,s,s };
    __m256d savx = _mm256_load_pd(ps);
    for (int i = 0; i < size / 4; i++)
    {
      __m256d avx = _mm256_load_pd(p + i * 4);
      __m256d addavx = _mm256_mul_pd(avx, savx);
      _mm256_store_pd(p + i * 4, addavx);
    }

    for (int i = size - offset; i < size; i++)
    {
      p[i] *= s;
    }
  }
  else
#endif
  {
    const int n = 2 * field.size();
    double* fieldp = (double*)(field.data());

    for (int i = 0; i < n; i++)
    {
      fieldp[i] *= s;
    }
  }
}

/*!
\brief Compute the curl of the vector field.
\sa Curl(int, int)
*/
ScalarField2 VectorField2::Curl() const
{
  ScalarField2 c(GetBox(), nx, ny);

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      c(i, j) = Curl(i, j);
    }
  }

  return c;
}

/*!
\brief Compute the curl at a given sample.
\sa Curl()
\param i,j Integer coordinates of the sample.
*/
double VectorField2::Curl(int i, int j) const
{
  double c = 0.0;

  // Divergence along x axis
  if (i == 0)
  {
    c += (at(i + 1, j)[1] - at(i, j)[1]) / celldiagonal[0];
  }
  else if (i == nx - 1)
  {
    c += (at(i, j)[1] - at(i - 1, j)[1]) / celldiagonal[0];
  }
  else
  {
    c += (at(i + 1, j)[1] - at(i - 1, j)[1]) / (2.0 * celldiagonal[0]);
  }

  // Divergence along y axis
  if (j == 0)
  {
    c -= (at(i, j + 1)[0] - at(i, j)[0]) / celldiagonal[1];
  }
  else if (j == ny - 1)
  {
    c -= (at(i, j)[0] - at(i, j - 1)[0]) / celldiagonal[1];
  }
  else
  {
    c -= (at(i, j + 1)[0] - at(i, j - 1)[0]) / (2.0 * celldiagonal[1]);
  }

  return c;
}

/*!
\brief Normalize the vector field, i.e. set the range of norm of the vectors [0,1].
*/
void VectorField2::Normalize()
{
  // Compute norm of the vectors in the field
  ScalarField2 n = GetNorm();

  double a, b;
  n.GetRange(a, b);
  if (a == b)
  {
    if (a == 0.0)
    {
    }
    else
    {
      a = 1.0 / a;
      for (int i = 0; i < nx * ny; i++)
      {
        field[i] *= a;
      }
    }
  }
  else
  {
    for (int i = 0; i < nx * ny; i++)
    {
      field[i] *= Linear::Affine(n.at(i), a, b);
    }
  }
}

/*!
\brief Set all vectors to unit, i.e. set the norm of the vectors to 1.

\sa VectorField2::Normalize()
*/
void VectorField2::Unit()
{
  for (int i = 0; i < nx * ny; i++)
  {
    field[i] /= Norm(field[i]);
  }
}

/*!
\brief Compute the range of the vector field.
\param a,b Returned minimum and maximum vector values.
*/
void VectorField2::GetRange(Vector2& a, Vector2& b) const
{
  a = b = field.at(0);
  for (int i = 1; i < nx * ny; i++)
  {
    Vector2::SetMinMax(field.at(i), a, b);
  }
}

/*!
\brief Compute the range of norm the vector field.
\param a,b Returned minimum and maximum norms.
*/
void VectorField2::GetNormRange(double& a, double& b) const
{
  a = b = Norm(field.at(0));
  for (int i = 1; i < nx * ny; i++)
  {
    double x = Norm(field.at(i));
    Math::SetMinMax(x, a, b);
  }
}

/*!
\brief Compute the average norm.
*/
double VectorField2::AverageNorm() const
{
  double an = 0;
  const int n = nx * ny;
  for (int i = 0; i < n; i++)
  {
    an += Norm(field.at(i));
  }
  return an / double(n);
}

/*!
\brief Compute the norm of the vector field.
*/
ScalarField2 VectorField2::GetNorm() const
{
  ScalarField2 n(Array2::GetBox(), nx, ny);
  for (int i = 0; i < nx * ny; i++)
  {
    n[i] = Norm(field.at(i));
  }
  return n;
}

/*!
\brief Compute the angle of the vector field.
*/
ScalarField2 VectorField2::GetAngle() const
{
  ScalarField2 sf(Box2(a, b), nx, ny);

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      Vector2 v = at(i, j);
      if (v == Vector2::Null)
      {
        sf(i, j) = 0.0;
      }
      else
      {
        sf(i, j) = v.Angle();
      }
    }
  }

  return sf;
}

/*!
\brief Create and return the scalar field with the x or y components.
\param i Index.
*/
ScalarField2 VectorField2::GetComponent(int i) const
{
  ScalarField2 xy(Box2(a, b), nx, ny);
  for (int j = 0; j < nx * ny; j++)
  {
    xy[j] = field[j][i];
  }
  return xy;
}

/*!
\brief Create the angle image.
*/
QImage VectorField2::CreateAngleImage() const
{
  QImage image(nx, ny, QImage::Format_ARGB32);

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      Vector2 v = at(i, j);
      if (v == Vector2::Null)
      {
        image.setPixel(i, j, QColor(0, 0, 0).rgb());
      }
      else
      {
        v = Normalized(v);
        double a = Math::Clamp(0.5 * (1.0 + v.Angle() / Math::Pi), 0.0, 1.0);
        image.setPixel(i, j, (Color(0.5) + 0.5 * Color::Wheel(a)).GetQt().rgb());
      }
    }
  }

  return image;
}

/*!
\brief Create the angle image.
*/
QImage VectorField2::CreateBlueRedImage() const
{
  QImage image(nx, ny, QImage::Format_ARGB32);

  double a, b;
  GetNormRange(a, b);

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      Vector2 v(0.0);
      if (a == b)
      {
        if (a == 0.0)
        {
        }
        else
        {
          v = at(i, j) / a;
        }
      }
      else
      {
        v = at(i, j) * (Norm(at(i, j)) - a) / (b - a);
      }

      double red = 0.5 * (1.0 + v[0]);
      double green = 0.5 * (1.0 + v[1]);
      image.setPixel(i, j, (Color(0.25) + 0.75 * Color(red, green, 1.0)).GetQt().rgb());

    }
  }

  return image;
}

/*!
\brief Create the gradient image (similar to CreateBlueRedImage() with different color settings and 16 bits depth)
*/
QImage VectorField2::CreateGradientImage(const double& max) const
{
  QImage image(nx, ny, QImage::Format_RGBA64);

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      Vector2 v = at(i, j) / max;

      double red = 0.5 * (1.0 + v[0]);
      double green = 0.5 * (1.0 + v[1]);
      image.setPixel(i, j, (Color(red, green, 1.0)).GetQt().rgb());
    }
  }

  return image;
}

/*!
\brief Get the field value at a given point.

This function relies on a bi-linear interpolation of the vectors.

\sa Math::Bilinear
\param p Point.
*/
Vector2 VectorField2::Value(const Vector2& p) const
{
  double u, v;
  int i, j;
  CellInteger(p, i, j, u, v);

  // Test position
  if (!InsideCellIndex(i, j))
    return Vector2::Null;

  return Vector2::Bilinear(at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1), u, v);
}

/*!
\brief Destructive sum of two vector fields.
\param f %Vector field.
*/
VectorField2& VectorField2::operator += (const VectorField2& f)
{
  for (int i = 0; i < field.size(); i++)
  {
    field[i] += f.at(i);
  }
  return *this;
}

/*!
\brief Destructive difference of two vector fields.
\param f %Vector field.
*/
VectorField2& VectorField2::operator -= (const VectorField2& f)
{
  for (int i = 0; i < field.size(); i++)
  {
    field[i] -= f.at(i);
  }
  return *this;
}

/*!
\brief Set the whole vector field to a constant vector.
\param v %Vector.
*/
void VectorField2::Set(const Vector2& v)
{
  field.fill(v, nx * ny);
}

/*!
\brief Compute the new position of a point in the vector field using forward Euler integration.

Direct forward Euler integration simply returns <B>p</B> + t v(<B>p</B>), error is bounded by O(t<SUP>2</SUP>).

Improved Eulerâ€™s method use a prediction-correction scheme, error is bounded by O(t<SUP>3</SUP>).
\param p Point.
\param t Time step.
\param s Scheme: false is simple Euler method, true for prediction-correction.
*/
Vector2 VectorField2::Euler(const Vector2& p, const double& t, bool s) const
{
  if (s == false)
  {
    return p + t * Value(p);
  }
  else
  {
    Vector2 vp = Value(p);
    Vector2 q = p + 0.5 * t * vp;
    return p + 0.5 * t * (vp + Value(q));
  }
}

/*!
\brief Compute the new position of a point in the vector field using Runge Kutta 4 integration.

Error is bounded by O(t<SUP>4</SUP>).

\param p Point.
\param t Time step.
*/
Vector2 VectorField2::RungeKutta(const Vector2& p, const double& t) const
{
  Vector2 k1 = t * Value(p);
  Vector2 v1 = p + 0.5 * k1;

  Vector2 k2 = t * Value(v1);
  Vector2 v2 = p + 0.5 * k2;

  Vector2 k3 = t * Value(v2);
  Vector2 v3 = p + k3;

  Vector2 k4 = t * Value(v3);

  return p + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
}


/*!
\brief Compute trajectory of a point using Euler integration.

\param p Point.
\param t Time step.
\param n Number of steps.
\param s Scheme: false is simple Euler method, true for prediction-correction.
*/
PointCurve2 VectorField2::EulerSteps(const Vector2& p, const double& t, int n, bool s) const
{
  PointCurve2 c;

  Vector2 q = p;

  c.Append(q);
  for (int i = 1; i < n; i++)
  {
    q = Euler(q, t, s);
    if (!Inside(q))
    {
      break;
    }
    c.Append(q);
  }
  return c;
}

/*!
\brief Compute line integral convolution image.

This is the implementation of Cabral, Brian and Leedom, Leith Casey. Imaging Vector Fields Using Line Integral Convolution. Proceedings of SIGGRAPH, 263â€“270, 1993.

Stalling, Detlev and Hege, Hans-Christian. Fast and resolution independent line integral convolution. Proceedings of SIGGRAPH, 249â€“256, 1995.

The vector field should be converted to unit direction vector using VectorField2::Unit.

\sa VectorField2::Unit().

\param n Number of steps.
*/
ScalarField2 VectorField2::LineIntegral(int n) const
{
  ScalarField2 line(*this);

  static Random random;

  // Kernel and integral
  QVector<double> kernel(n);
  double integral = 0;
  for (int k = 0; k < n; k++)
  {
    double a = 1.0 - Cubic::Smooth(Math::Unit(k, n));
    kernel[k] = a;
    integral += a;
  }
  for (int k = 0; k < n; k++)
  {
    kernel[k] /= integral;
  }

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      // White noise
      double x = random.Uniform();

      // Forward
      PointCurve2 curve = EulerSteps(ArrayVertex(i, j), 0.5, n, false);
      for (int k = 0; k < curve.Size(); k++)
      {
        Vector2 pk = curve[k];
        double a = kernel[k];
        line.Add(pk, x * a);
      }
      // Backward tracking
      curve = EulerSteps(ArrayVertex(i, j), -0.5, n, false);

      for (int k = 0; k < curve.Size(); k++)
      {
        Vector2 pk = curve[k];
        double a = kernel[k];
        line.Add(pk, x * a);
      }
    }
  }

  return line;
}

/*!
\brief Perform a Gaussian blur.
\param r Radius of blur.

\sa ScalarField2::GaussianBlur(int)
*/
VectorField2 VectorField2::GaussianBlur(const double& r) const
{
  int n = r / Array2::celldiagonal[0];
  return GaussianBlur(n);
}

/*!
\brief Perform a Gaussian blur.
\param r Radius (in cells) of blur.
*/
VectorField2 VectorField2::GaussianBlur(int r) const
{
  // Fixed sigma to half kernel radius
  double sigma = 0.5 * double(r);
  double sigma2 = sigma * sigma;
  double norm = 1.0 / (sqrt(2.0 * Math::Pi) * sigma);

  // Kernel
  std::vector<double> kernel(r + 1);
  for (int i = 0; i <= r; i++)
  {
    kernel[i] = norm * exp((-0.5 * double(i) * double(i)) / sigma2);
  }

  QVector<Vector2> resx(nx * ny);

  // Filter x direction
  for (int y = 0; y < ny; y++)
  {
    for (int x = 0; x < nx; x++)
    {
      Vector2 v = kernel[0] * at(x, y);
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

  VectorField2 res(Box2(a, b), nx, ny);

  // Filter y direction
  for (int x = 0; x < nx; x++)
  {
    for (int y = 0; y < ny; y++)
    {
      Vector2 v = kernel[0] * resx[VertexIndex(x, y)];
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