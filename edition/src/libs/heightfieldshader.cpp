// Heightfield

#include "libs/gputerrainsimu.h"
#include "libs/heightfieldshader.h"
#include "libs/vectorfield.h"
#include "libs/cubic.h"
#include "libs/cpu.h"

void Normalize(Vector&);

/*!
\class HeightFieldShader heightfieldshader.h
\brief Relief shading class.

The texture image representing the shaded terrain for Hugo's 2024 paper:
 \code
 return LaplacianShading(Palette::HugoShading, 1.6, 0.6);
\endcode

*/

/*!
\brief Create a shader.
\param field Reference to a heightfield.
*/
HeightFieldShader::HeightFieldShader(const HeightField& field) :HeightField(field)
{
}



/*!
\brief Create a texture image representing the shaded terrain.
*/
QImage HeightFieldShader::ReliefLighting() const
{
  double ra, rb;

  // Get elevation range
  GetRange(ra, rb);

  // Modify range by a 5%
  ra -= (rb - ra) * 0.05;
  rb += (rb - ra) * 0.05;

  // Shading
  QImage shading(nx, ny, QImage::Format_ARGB32);

  ScalarField2 light = Light(Normalized(Vector(-1.0, 0.5, 2.5)));

  // Height field
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      // Normalized height
      double z = Linear::Step(at(i, j), ra, rb);

      // Height shading: color is a linear interpolation of height colors
      Color cz = Color::Lerp(z, Color(0.75, 0.725, 0.70), Color(0.95, 0.925, 0.90));

      double s = light(i, j);
      s *= s;
      s *= s;

      Color c = cz * (1.0 + s) / 2.0;

      shading.setPixel(i, j, c.GetQt().rgb());
    }
  }
  return shading;
}

/*!
\brief Diffuse lighting with enhanced curvature shaded relief.

Create a texture representing the shaded terrain.
*/
QImage HeightFieldShader::Relief(bool enhanced) const
{
  double ra, rb;

  // Get elevation range
  GetRange(ra, rb);

  // Modify range by a 5%
  ra -= (rb - ra) * 0.05;
  rb += (rb - ra) * 0.05;

  // Shading
  QImage shading = CreateEmptyImage();

  // Light
  ScalarField2 light = Light(Normalized(Vector(-1.0, 0.5, 2.5)));

  // Laplacian for details
  ScalarField2 laplacian = Laplacian();

  laplacian.Step(-0.20, 0.20);

  // Height field
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      // Normalized height
      double z = Linear::Step(at(i, j), ra, rb);

      // Height shading: color is a linear interpolation of height colors
      Color cz = Color::Lerp(z, Color(0.75, 0.725, 0.70), Color(0.95, 0.925, 0.90));

      double s = light(i, j);
      s *= s;
      s *= s;

      double la = 0.0;

      if (enhanced == true)
      {
        la = 1.0 - 2.0 * laplacian(i, j);

        la *= fabs(la);
      }
      Color c = cz * ((1.0 + s) / 2.0 + 0.15 * la);

      c = c.Clamp();
      shading.setPixel(i, j, c.GetQt().rgb());
    }
  }
  return shading;
}




// usefull function for the next method
double ImageShadeRadianceScalingWarp(double impfunc, double beta) {
  const double alpha = 0.1;
  double expbeta = exp(-beta);
  return (alpha * expbeta + impfunc * (1.0 - alpha - alpha * expbeta)) / (alpha + impfunc * (expbeta - alpha - alpha * expbeta));
}

// **** WEIGHT FUNCTIONS ****
double ImageShadeRadianceScalingsilhouetteWeight(double s) {
  const double ts = 0.07;
  const double t2 = 0.9 + ts;
  const double t1 = t2 - 0.01;

  return Cubic::SmoothStep(t1, t2, max(1.0 - s, 0.9));
}

double ImageShadeRadianceScalingTanh(double c, double en) {
  double cmax = en * 15.0;
  const double tanhmax = 3.11622;

  double x = ((c * cmax * 1.0) / tanhmax);
  double e = exp(-2.0 * x);
  double t = Clamp((1.0 - e) / (1.0 + e), -1.0, 1.0);

  return t;
}

double ImageShadeRadianceScalingCurvature(double w, Vector h, double e) {
  double c = ImageShadeRadianceScalingTanh(-(h[0] + h[1]) / 2.0, e);
  return  c * max(w - 0.5, 0.0); //invert ? -c * max(w - 0.5, 0.0) :
}

/*!
\brief Create a scalarfield representing the shaded terrain with an irradiance scaling effect.
\param pos is the position of the lights typical range is [0.1-1.0]
       dir1 = [pos,pos,0.5]
       dir2 = [-2*pos,-2*pos,0.5]
\param factor factor of radiance scaling typical values [0.5,4.0]
Returns a normalized scalarfield that can be used to enhance high curvature crests
*/
ScalarField2 HeightFieldShader::RadianceScaling(const double& pos, const double& factor) const
{
  // result
  ScalarField2 shading(Box2(a, b), nx, ny);

  HeightField hf = *this;
  hf.Unity();

  Vector lightdir1 = Normalized(Vector(pos, pos, 0.5));
  Vector lightdir2 = Normalized(Vector(-2.0 * pos, -2.0 * pos, 0.5));

  VectorField2 g = hf.Gradient();

  // Height field
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      // "hessian" computation
      Vector2 gxp = i < nx - 1 ? g(i + 1, j) : Vector2(0.0);
      Vector2 gxm = i > 0 ? g(i - 1, j) : Vector2(0.0);
      double xx = gxp[0] - gxm[0];
      double xy = gxp[1] - gxm[1];

      Vector2 gyp = j < ny - 1 ? g(i, j + 1) : Vector2(0.0);
      Vector2 gym = j > 0 ? g(i, j - 1) : Vector2(0.0);
      double yx = gyp[0] - gym[0];
      double yy = gyp[1] - gym[1];

      Vector h(xx, yy, (xy + yx) / 2.0);
      double w = 0.98;// ImageShadeRadianceScalingWeight(); // this function should calculate the average of depth value differences with respect to the central point and on which we apply the silhouette function

      double c = ImageShadeRadianceScalingCurvature(w, h, factor);

      Vector n = hf.Normal(i, j);

      // shading with light1
      double cosineTerm1 = Math::Clamp(n * lightdir1, 0.0, 1.0);
      double cosineTerm2 = Math::Clamp(n * lightdir2, 0.0, 1.0);


      // final radiance
      double radiance = Math::Clamp(0.5 * cosineTerm1 * ImageShadeRadianceScalingWarp(cosineTerm1, c) + 0.5 * cosineTerm2 * ImageShadeRadianceScalingWarp(cosineTerm2, c), 0., 1.);
      //double radiance = Math::Clamp(abs(c), 0., 1.);// 0.6*cosineTerm1 + 0.6*cosineTerm2;

      shading(i, j) = radiance;
    }
  }
  return shading;
}


/*!
\brief Create a QImage representing the shaded terrain with an irradiance scaling effect with two colors
*/
QImage HeightFieldShader::ImageShadeRadianceScaling(const Color& color1, const Color& color2) const
{
  ScalarField2 rs = RadianceScaling(0.3, 1.0);

  // Shading
  QImage shading(nx, ny, QImage::Format_ARGB32);
  // Height field
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      Color mix = rs(i, j) * color1 + (1.0 - rs(i, j)) * color2;
      shading.setPixelColor(i, j, mix.GetQt());
    }
  }
  return shading;
}


/*!
\brief Create a texture image representing the shaded terrain.
*/
QImage HeightFieldShader::Mitsuba() const
{
  // Shading
  QImage shading(nx, ny, QImage::Format_ARGB32);

  // Height field
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      Vector n = Normal(i, j);
      double s = n * Normalized(Vector(-2.0, 1.0, 4.0));
      s = 0.5 * (1.0 + s);
      s *= s;
      s *= s;

      // Normal shading: color is a combination of cool and cold colors according to the orientation
      Color cn = Color::Lerp(s, Color(50, 65, 85), Color::White); //52, 78, 118

      Color c = 0.15 * Color(0.975, 0.975, 0.975) + 0.85 * cn;
      shading.setPixel(i, j, c.GetQt().rgb());
    }
  }
  return shading;
}

/*!
* \brief Create a texture image representing the shaded terrain with a Laplacian shading effect.
* \param palette Palette used to color the terrain.
* \param gamma Gamma correction.
* \param lapamp Laplacian enhancement factor
* \return QImage representing the shaded terrain.
* */
QImage HeightFieldShader::LaplacianShading(const Palette& palette, double gamma, double lapamp) const
{
  QElapsedTimer timer;
  timer.restart();// change resolution


  HeightField heightfield = *this;

  if (heightfield.GetSizeX() < 2048) // we don't want to downscale big maps !
  {
    heightfield.Smooth(); // optional
    heightfield = heightfield.SetResolution(2048, 2048, true);
  }
  else if (heightfield.GetSizeX() > 4096) // if it's too big, downscale
  {
    heightfield = heightfield.SetResolution(4096, 4096, true);
  }
  std::cout << "resize: "; System::ShowElapsed(timer);

  // breach (optional)
  //heightfield.CompleteBreach();
  std::cout << heightfield.GetSizeX() << std::endl;
  // Define texture from palette
  QImage texture = heightfield.ScalarField2::CreateImage(palette);

  std::cout << "palette: "; System::ShowElapsed(timer);

  // fractlap analysis
  GPUHeightFieldAnalysis analysis;
  ScalarField2 fractlap = analysis.FractionalLaplacian(heightfield, 20, 0.5, 10);

  std::cout << "fraclapl: "; System::ShowElapsed(timer);
  double a, b;
  fractlap.GetRange(a, b);
  fractlap.Scale(1.0 / max(b, abs(a))); // normalize and keep it centered
  std::cout << "centered: ";  System::ShowElapsed(timer);

  // final postprocess, adjust gamma, saturation, and brightness using fractlap
#pragma omp parallel for
  for (int i = 0; i < texture.width(); i++)
    for (int j = 0; j < texture.height(); j++)
    {
      double lap = fractlap(i, j);
      double e = 0.7;
      if (lap > 0.)
        lap = pow(lap, e);
      else
        lap = -pow(-lap, e);

      QColor qcolor = texture.pixelColor(i, j);
      Color color(qcolor);

      // apply gamma correction
      color = Color(pow(color[0], gamma), pow(color[1], gamma), pow(color[2], gamma));

      qcolor = QColor(255. * color[0], 255. * color[1], 255. * color[2]);
      float h, s, v;
      qcolor.getHsvF(&h, &s, &v);

      s = Math::Clamp(s - 0.05f, 0.0, 1.0); // decrease saturation
      v = Math::Clamp(v + lap * lapamp, 0., 1.); // use laplacian to darken and lighten depending wether it's concave or convex // Gundy: used to be 0.6, was simply too much

      qcolor.setHsvF(h, s, v);

      texture.setPixelColor(i, j, qcolor);
    }
  std::cout << "adjust: "; System::ShowElapsed(timer);
  return texture;
}


QImage HeightFieldShader::BrownShading() const
{
  double h_min, h_max;
  GetRange(h_min, h_max);
  return BrownShading(h_min, h_max);
}


QImage HeightFieldShader::BrownShading(double h_min, double h_max) const
{
  // see heightfield_raytrace.glsl in LibMaya for comments on the texture computation
  QImage texture = CreateEmptyImage();
  ScalarField2 f_lap = Laplacian();// = GPUHeightFieldAnalysis().FractionalLaplacian(*this, 4.0, 0.25, 3);
  f_lap.Normalize();

  // parameters
  Color white(255, 255, 255);
  Color black(0, 0, 0);
  Color brown(125, 70, 45);
  double alpha = 2.0;
  double aa = 0.0;
  double bb = 1.0;
  Vector light_dir(0.25, 0.25, 2.5);
  light_dir = Normalized(light_dir);

  for (int i = 0; i < texture.width(); i++)
  {
    for (int j = 0; j < texture.height(); j++)
    {

      // elevation color
      double h = at(i, j);
      double linearStep = (h - h_min) / (h_max - h_min);
      linearStep = std::min(std::max(linearStep, 0.), 1.);
      double t = alpha * linearStep;
      t = std::min(std::max(t, aa), bb);
      Color color = Color::Lerp(t, white, brown);

      // slope color
      Vector normal = Normal(i, j);
      t = 1. - normal[2];
      t = pow(t, 2.);
      color = Color::Lerp(t, color, black);

      // laplacian color
      t = 0.4 * f_lap.at(i, j);
      color = Color::Lerp(t, color, white);

      // lighting
      double s = normal * light_dir;
      s = 0.5 * (1.0 + s);
      s *= s;
      color = color * s;

      texture.setPixelColor(i, j, color.GetQt());
    }
  }

  return texture;
}


/*!
\brief Create a texture image representing the iso contours of the terrain.

Contours are created from the origin elevation.
\param color Color of the minor iso-lines.
\param bold Color of the major iso-lines.
\param height Height range between two major iso-lines.
\param n Number of minor iso-lines between two major iso-lines.
*/
QImage HeightFieldShader::IsoLines(const QColor& color, const QColor& bold, const double& height, int n) const
{
  QImage image(nx, ny, QImage::Format::Format_ARGB32);

  // Set transparent image
  image.fill(QColor(0, 0, 0, 0));
  double dh = height / n;

  // Scan 
  for (int i = 0; i < image.width(); i++)
  {
    for (int j = 0; j < image.height() - 1; j++)
    {
      double a = at(i, j);
      double b = at(i, j + 1);

      int na = int(a / dh);
      int nb = int(b / dh);
      if (na != nb)
      {
        QColor c = color;
        if ((nb % n == 0) || (na % n == 0))
        {
          c = bold;
        }
        // Alpha

        //double z = Math::Max(na,nb) * dh;
        //double wa = fabs(z - b)/(b-a);
        //double wb = fabs(z - a) / (b - a);

        //QColor cwa = c;
        //cwa.setAlphaF(wa);
        //QColor cwb = c;
        //cwb.setAlphaF(wb);

        //image.setPixel(i, j, cwa.rgba());
        //image.setPixel(i, j + 1, cwb.rgba());

        image.setPixel(i, j, c.rgba());
        image.setPixel(i, j + 1, c.rgba());

      }
    }
  }

  for (int i = 0; i < image.width() - 1; i++)
  {
    for (int j = 0; j < image.height(); j++)
    {
      double a = at(i, j);
      double b = at(i + 1, j);

      int na = int(a / dh);
      int nb = int(b / dh);
      if (na != nb)
      {
        if ((nb % n == 0) || (na % n == 0))
        {
          image.setPixel(i, j, bold.rgb());
          image.setPixel(i + 1, j, bold.rgb());
        }
        else
        {
          image.setPixel(i, j, color.rgb());
          image.setPixel(i + 1, j, color.rgb());
        }
      }
    }
  }
  return image;
}
/*

QImage image(nx, ny, QImage::Format::Format_ARGB32);

// Set transparent image
image.fill(QColor(0, 0, 0, 0));
double dh = height / n;

// Scan
for (int i = 0; i < image.width(); i++)
{
  for (int j = 0; j < image.height() - 1; j++)
  {
    double zij = at(i, j);

    QPoint pij(i, j);

    for (int k = 0; k < 8; k++)
            {
      QPoint q = Next(pij, k);
      // Outside domain
      if (!InsideVertexIndex(q)) continue;
      double z = at(q);

      int na = int(zij / dh);
      int nb = int(z / dh);
      if (na == nb)
        continue;

        QColor c = color;
      if ((nb % n == 0) || (na % n == 0))
      {
        c = bold;
      }
      double z0 = Math::Max(na, nb) * dh;

      // Alpha
      double t=Linear::Solve(0.0, length[k], zij - z0, z - z0);

    }


  }
}
*/

/*!
\brief Create a texture image representing the iso contours of the terrain with an anti-aliasing algorithm in altitude space (resulting in larger lines in flat areas)

Contours are created from the origin elevation.
\param color Color of the minor iso-lines.
\param bold Color of the major iso-lines.
\param height Height range between two major iso-lines.
\param n Number of minor iso-lines between two major iso-lines.
*/
QImage HeightFieldShader::IsoLinesAA(const QColor& color, const QColor& bold, const double& height, int n) const
{
  QImage image(nx, ny, QImage::Format::Format_ARGB32);

  // Set transparent image
  QColor bg(255, 255, 255, 255);
  image.fill(bg);

  // Scan 
  for (int i = 0; i < image.width(); i++)
  {
    for (int j = 0; j < image.height() - 1; j++)
    {
      double a = at(i, j);
      double a2 = a / height * 10.; // normalize to a 10 interval to fit the following function
      double a1 = a / height * n * 10.; // normalize to a 10 interval to fit the following function

      double alpha1 = Math::Clamp(5 * abs(fmod(a1, 10.) - 5.) - 23.5, 0., 1.);
      double alpha2 = Math::Clamp(5 * abs(fmod(a2, 10.) - 5.) - 23.5, 0., 1.);
      QColor c = color;
      double alpha = alpha1;
      if (alpha2 > alpha1) {
        c = bold;
        alpha = alpha2;
      }
      QColor mix(alpha * c.red() + (1.0 - alpha) * bg.red(), alpha * c.green() + (1.0 - alpha) * bg.green(), alpha * c.blue() + (1.0 - alpha) * bg.blue(), alpha * c.alpha() + (1.0 - alpha) * bg.alpha());
      image.setPixel(i, j, mix.rgb());
    }
  }
  return image;
}



/*!
\brief Diffuse lighting with enhanced curvature shaded relief.

Create a texture representing the shaded terrain.
*/
QImage HeightFieldShader::ShadedRelief() const
{
  double ra, rb;

  // Get elevation range
  GetRange(ra, rb);

  // Modify range by a 5%
  ra -= (rb - ra) * 0.05;
  rb += (rb - ra) * 0.05;

  // Shading
  QImage shading = CreateEmptyImage();

  // Light
  ScalarField2 light = Light(Normalized(Vector(-1.0, 0.5, 2.5)));

  ScalarField2 fractional = GPUHeightFieldAnalysis().FractionalLaplacian(*this, 2.0, 0.75, 50);
  fractional.Unitize();

  double u, v;
  fractional.GetRange(u, v);
  ScalarField2 accessibility = GPUHeightFieldAnalysis().Accessibility(*this, 25.0 * Norm(CellDiagonal()), 64);

  // Laplacian for details
  ScalarField2 laplacian = Laplacian();
  laplacian.Step(-0.20, 0.20);



  // Height field
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      // Normalized height
      double z = Linear::Step(at(i, j), ra, rb);

      // Height shading: color is a linear interpolation of height colors
      Color cz = Color::Lerp(z, Color(0.75, 0.725, 0.70), Color(0.95, 0.925, 0.90));

      double s = light(i, j);
      s *= s;
      s *= s;

      Color c = cz * (1.0 + s) / 2.0;
      double la = 1.0 - 2.0 * laplacian(i, j);

      la *= fabs(la);


      c *= 0.70 + 0.30 * ((1.0 + fractional(i, j)) / 2.0) + 0.05 * accessibility(i, j) + 0.05 * la;

      shading.setPixel(i, j, c.GetQt().rgb());
    }
  }
  return shading;
}


/*!
\brief Cartographic-like shading with green/yellow colors, using fractional laplacian to emphasize ridges.

Rivers are detected after a breaching and shaded in blue.

*/
QImage HeightFieldShader::CartographicGreen() const
{
  // this has to be adjusted depending on the terrain size and dynamic
  ScalarField2 sf = *this;
  sf.Normalize();
  sf.Scale(3000. / 65535.0);

  // rivers are computed at 1024 resolution, to avoid having too thin lines
  ScalarField2 rivers;
  HeightField hf1024(sf);
  hf1024.CompleteBreach();

  ScalarField2 flowmap = hf1024.StreamAreaWeighted(2.0);
  flowmap.Normalize();
  rivers = flowmap.SetResolution(2048, 2048);
  sf = sf.SetResolution(2048, 2048);

  HeightField heightfield(sf);
  heightfield.CompleteBreach();
  ScalarField2 sfSlope = heightfield.AverageSlope();
  sfSlope.Normalize();

  // automatically compute 80% of the culumated histogram to keep the data usable
  QVector<QPair<double, double>> histo;
  histo = sfSlope.CumulativeNormedHistogram(5000);
  int i = 0;
  while (histo[i].second < 0.8)
    i++;
  double bound = histo[i].first;
  sfSlope.Clamp(0.0, bound); // StreamAreaWeighted [1; 30k]
  sfSlope.Normalize();

  ScalarField2 sfArea = heightfield.StreamAreaWeighted(2.0);
  flowmap = sfArea;
  flowmap.Normalize();

  // automatically compute 50% of the culumated histogram to keep the data usable
  histo = flowmap.CumulativeNormedHistogram(10000);
  i = 0;
  while (histo[i].second < 0.5)
    i++;
  bound = histo[i].first;
  flowmap.Clamp(0.0, bound);
  flowmap.Normalize();
  flowmap.SelfOp([&](double x) { return pow(x, 0.7); });

  QImage image(System::GetArchesLib() + QString("/Gradients/cartography.png"));

  QImage texture = flowmap.CreateImage(Palette(image));
  GPUHeightFieldAnalysis analysis;

  // parameters could be adjusted
  ScalarField2 fractlap = analysis.FractionalLaplacian(heightfield, 10, 0.5, 10);
  fractlap.Unitize();
  ScalarField2 slopeField = Slope();

  for (i = 0; i < texture.width(); i++)
    for (int j = 0; j < texture.height(); j++)
    {
      double lap = fractlap(i, j);
      double e = 0.7;
      if (lap > 0.)
        lap = pow(lap, e);
      else
        lap = -pow(-lap, e);

      QColor qcolor = texture.pixelColor(i, j);
      Color color(qcolor);

      double slope = sfSlope(i, j);
      double grass = pow(Math::Clamp(slope * 2.0 - 0.3, 0., 1.), 0.3);

      // insert grass with green color only on convex parts
      if (lap > -0.2)
        color = grass * color + (1.0 - grass) * Color(149, 169, 122); // insert grass on gentle slopes

      double riv = pow(rivers(i, j), 0.74);
      color = riv * Color(15, 0, 196) + (1.0 - riv) * color; // insert rivers
      color = color.Clamp();

      // apply gamma correction
      double gamma = 1.0;
      color = Color(pow(color[0], gamma), pow(color[1], gamma), pow(color[2], gamma));
      qcolor = QColor(255. * color[0], 255. * color[1], 255. * color[2]);
      float h, s, v;
      qcolor.getHsvF(&h, &s, &v);

      // rise saturation
      s = Math::Clamp(s + 0.25);

      // use laplacian to darken and lighten depending wether it's concave or convex
      v = Math::Clamp(v + lap * 0.8);

      qcolor.setHsvF(h, s, v);
      texture.setPixelColor(i, j, qcolor);
    }
  return texture;
}

/*!
\brief Compute an image representing a river of varying with by using the steepest drainage area.

\param t Stream area threshold: river will exist only for drainage area above this value.
\param size Largest rive size (in pixels).
*/
QImage HeightFieldShader::RiverWidth(const double& t, int size) const
{
  // River map
  ScalarField2 alpha(Box2(a, b), nx, ny, 0.0);
  ScalarField2 tint(Box2(a, b), nx, ny, 0.0);

  ScalarField2 sa = StreamAreaSteepest();
  sa.SelfOp([&](double x) { return x > t ? pow(x - t, 0.34) : 0.0; });

  sa.Normalize();

  const double r = Norm(celldiagonal);

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      double a = sa(i, j);

      // Skip if empty
      if (a == 0.0) continue;

      // Width
      double w = 1.0 + size * a;

      // Number of cells
      int n = int(w) + 1;
      QRect area(i - n, j - n, 2 * n + 1, 2 * n + 1);

      // Limit to domain
      QRect mask(0, 0, nx - 1, ny - 1);
      area = area.intersected(mask);

      Vector2 center = ArrayVertex(i, j);
      double re = w * r;
      double ri = (w - 1.0) * r;

      for (int y = area.y(); y <= area.y() + area.height(); y++)
      {
        for (int x = area.x(); x <= area.x() + area.width(); x++)
        {
          // Distance between central point and current point
          double u = Norm(center - ArrayVertex(x, y));

          double f = 1.0 - Cubic::SmoothStep(u, ri, re);

          if (f > alpha[VertexIndex(x, y)])
          {
            alpha[VertexIndex(x, y)] = f;
            tint[VertexIndex(x, y)] = a;
          }
        }
      }
    }
  }

  QImage image(nx, ny, QImage::Format_ARGB32);

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      // Interpolate opaque colors
      Color c = Color::Lerp(tint(i, j), Color(0.1, 0.8, 0.9, 1.0), Color(0.1, 0.4, 0.5, 1.0));
      // And set alpha
      c[3] = alpha(i, j);
      // Shift to white if transparent
      if (alpha(i, j) == 0.0) { c = Color(1.0, 1.0, 1.0, 0.0); }
      image.setPixelColor(i, j, c.GetQt());
    }
  }

  return image;
}



/*!
\brief Cartographic-like shading with green/yellow colors.
*/
QImage HeightFieldShader::CartographicGreenYellow() const
{
  // Define palette
  Palette anchorPalette({ Color(0.54, 0.54, 0.40), Color(0.59, 0.60, 0.44), Color(0.67, 0.64, 0.47), Color(0.76, 0.69, 0.54), Color(0.79, 0.74, 0.59), Color(0.84, 0.79, 0.61), Color(0.95, 0.86, 0.70), Color(0.97, 0.97, 0.97) });

  QImage texture = CreateImage(Palette(anchorPalette));

  // parameters could be adjusted
  GPUHeightFieldAnalysis analysis;
  ScalarField2 fractlap = analysis.FractionalLaplacian(*this, 10, 0.5, 10);

  fractlap.Unitize();

  for (int i = 0; i < texture.width(); i++)
  {
    for (int j = 0; j < texture.height(); j++)
    {
      double lap = fractlap(i, j);

      QColor qcolor = texture.pixelColor(i, j);
      Color color(qcolor);

      qcolor = QColor(255. * color[0], 255. * color[1], 255. * color[2]);
      float h, s, v;
      qcolor.getHsvF(&h, &s, &v);

      /*
      Color cc(qcolor);
            cc=cc.Darken(lap * 0.15);
            cc.Saturate(-0.1);
            qcolor = cc.GetQt();
      */

      s = Math::Clamp(s - 0.1); // decrease saturation
      // Use laplacian to darken and lighten depending whether it's concave or convex
      v = Math::Clamp(0.95 * v + lap * 0.15);

      qcolor.setHsvF(h, s, v);

      texture.setPixelColor(i, j, qcolor);
    }

  }
  return texture;
}
