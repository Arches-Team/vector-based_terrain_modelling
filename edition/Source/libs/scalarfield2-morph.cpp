// Fields

#include "libs/scalarfield.h"
#include "libs/cubic.h"

/*!
\brief Convolution of the field using a kernel.

\param t Array of scalar values.
\param s Size (odd number) of the square kernel.
*/
ScalarField2 ScalarField2::Convolution(double t[], int s)
{
  int radius = (s - 1) / 2;
  ScalarField2 result(Box2(a, b), nx, ny);

  // room for optimization here
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      double v = 0.0;
      for (int ii = -radius; ii <= radius; ii++)
      {
        for (int jj = -radius; jj <= radius; jj++)
        {
          if (i + ii >= 0 && i + ii < nx && j + jj >= 0 && j + jj < ny) // f...ing borders
            v += at(i + ii, j + jj) * t[ii + radius + s * (jj + radius)];
        }
      }
      result(i, j) = v;
    }
  }
  return result;
}

/*!
\brief Compute the skeleton using successive erosion/dilatations.

\param n Maximum number of iterations, which is related to the maximum thickness of the shape.
*/
ScalarField2 ScalarField2::MorphSkeleton(int n)
{
  ScalarField2 previous(*this);
  ScalarField2 current(*this);
  ScalarField2 eroded;
  ScalarField2 skeleton(Box2(a, b), nx, ny);

  for (int i = 0; i < n; i++)
  {
    QImage im3 = current.CreateImage(true);
    im3.save("init.png");
    current.MorphErode();
    QImage im = current.CreateImage(true);
    im.save("erode.png");
    eroded = current;
    current.MorphDilate();
    QImage im2 = current.CreateImage(true);
    im2.save("redilate.png");
    previous -= current;

    QImage im4 = previous.CreateImage(true);
    im4.save("diff.png");
    skeleton += previous;
    previous = eroded;
    current = eroded;
  }
  return skeleton;
}

void RotateKernel(double k[]) {
  double r[9];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      r[i + j * 3] = k[j + (2 - i) * 3];
    }
  }
  for (int i = 0; i < 9; i++)
    k[i] = r[i];
}

/*!
\brief Compute the skeleton using thinning.

\param n Maximum number of iterations, which is related to the maximum thickness of the shape.
*/
ScalarField2 ScalarField2::MorphSkeletonConnected(int n)
{
  double thin1[9] = { -1, -1, -1, 0, 1, 0, 1, 1, 1 };
  double thin2[9] = { 0, -1, -1, 1, 1, -1, 0, 1, 0 };

  ScalarField2 skeleton(*this);

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      skeleton.MorphThin(thin1, 3);
      skeleton.MorphThin(thin2, 3);
      RotateKernel(thin1);
      RotateKernel(thin2);
    }
  }
  return skeleton;
}

/*!
\brief Dilate the binary field using a D4 structuring element.
*/
void ScalarField2::MorphDilate()
{
  double square3[9] = { 0,1,0,1,1,1,0,1,0 };
  ScalarField2 conv = Convolution(square3, 3);
  conv.Binarize(0.9);
  field = conv.field; // is this correct?
}

/*!
\brief Erode the binary field using a D4 structuring element.
*/
void ScalarField2::MorphErode()
{
  double square3[9] = { 0,1,0,1,1,1,0,1,0 };
  ScalarField2 conv = Convolution(square3, 3);
  conv.Binarize(4.9);
  field = conv.field; // is this correct?
}

/*!
\brief Remove the ends of lines by authorizing at least 3 pixels.
*/
void ScalarField2::MorphRemoveEnds()
{
  double square3[9] = { 1,1,1,1,1,1,1,1,1 };
  ScalarField2 conv = Convolution(square3, 3);
  conv.Binarize(2.9);
  this->operator*=(conv);
  Binarize(0.5);
}

/*!
\brief Perform a hit and miss operation on the field.
\param k is the kernel, containing 1 and -1 for hit and miss resp. and 0 for "don't care"
\param s Size of the kernel (should be odd)
*/
void ScalarField2::MorphHitAndMiss(double k[], int s)
{
  int taille = s * s;
  double* hit = new double[taille];
  double* miss = new double[taille];
  double nhit = 0.0;
  double nmiss = 0.0;
  for (int i = 0; i < taille; i++)
  { // build the hit and miss kernels
    if (k[i] > 0.5) {
      hit[i] = 1.0;
      miss[i] = 0.0;
      nhit += 1.0;
    }
    else if (k[i] < -0.5)
    {
      hit[i] = 0.0;
      miss[i] = 1.0;
      nmiss += 1.0;
    }
    else
    {
      hit[i] = 0.0;
      miss[i] = 0.0;
    }
  }
  ScalarField2 hitf = Convolution(hit, s);

  ScalarField2 neg(*this);
  neg.Negate();
  ScalarField2 missf = neg.Convolution(miss, s);

  hitf.Binarize(nhit - 0.1);
  missf.Binarize(nmiss - 0.1);
  hitf *= missf; // both hit and miss must have a 1.0 value
  hitf.Binarize(0.5);
  field = hitf.field;
  delete[]hit;
  delete[]miss;
}

/*!
\brief Thinning by using a kernel.
\param k Kernel, containing 1 and -1 for hit and miss respectively, and 0 for <I>don't care</I>.
\param s Size of the kernel (should be odd).
*/
void ScalarField2::MorphThin(double k[], int s)
{
  ScalarField2 ham(*this);
  ham.MorphHitAndMiss(k, s);
  ham.Negate();
  (*this) *= ham;
}
