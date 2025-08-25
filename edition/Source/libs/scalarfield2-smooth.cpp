
// Fields

#include "libs/scalarfield.h"
#include "libs/cubic.h"

/*!
\brief Smooth the scalar field within a given region.

The smoothing will be applied using the cubic falloff function Cubic::Smooth(const double&, const double&)
\param center Center.
\param r Radius.
*/
void ScalarField2::Smooth(const Vector2& center, double r)
{
  // Vertex
  QPoint q = VertexInteger(center);
  const int x = q.x();
  const int y = q.y();

  // Calcul de la zone de modification (area)
  const QRect area = VertexIntegerArea(Box2(center, r));

  const int vx = area.width() + 1;
  const int vy = area.height() + 1;

  QVector<double> h;
  h.resize(vx * vy);

  // Compute smoothing
  for (int i = 0; i < vx; i++)
  {
    for (int j = 0; j < vy; j++)
    {
      h[i + j * vx] = SmoothPoint(area.x() + i, area.y() + j);
    }
  }

  // Squared radius
  r *= r;

  // Update
  for (int i = area.x(); i <= area.x() + area.width(); i++)
  {
    for (int j = area.y(); j <= area.y() + area.height(); j++)
    {
      // Distance between central point and current point
      double u = SquaredNorm(center - ArrayVertex(x, y));

      if (u < r)
      {
        // Scaling set to half the cubic value
        const double a = 0.5 * Cubic::Smooth(u, r);
        field[VertexIndex(i, j)] = (1.0 - a) * at(i, j) + a * h[i - area.x() + (j - area.y()) * vx];
      }
    }
  }
}

/*!
\brief Set the border of the scalar field to a constant value.
\param c Constant.
*/
void ScalarField2::SetBorder(const double& c)
{
  ScalarField2& s = *this;

  for (int i = 0; i < nx; i++)
  {
    s(i, 0) = c;
    s(i, ny - 1) = c;
  }

  for (int i = 0; i < ny; i++)
  {
    s(0, i) = c;
    s(nx - 1, i) = c;
  }
}

/*
\brief Compute a smoothed value at a given point in the scalar field.

The function uses a 3<SUP>2</SUP> approximation of the Gaussian kernel.
\param x, y Integer coordinates of the point.
*/
double ScalarField2::SmoothPoint(int x, int y) const
{
  int k = VertexIndex(x, y);

  if (x == 0)
  {
    // Corner
    if (y == 0)
    {
      return (4.0 * at(k) + 2.0 * at(k + 1) + 2.0 * at(k + ny)) / 8.0;
    }
    // Corner
    else if (y == ny - 1)
    {
      return (4.0 * at(k) + 2.0 * at(k + 1) + 2.0 * at(k - ny)) / 8.0;
    }
    else
    {
      return (2.0 * at(k - ny) + 4.0 * at(k) + 2.0 * at(k + ny) + at(k + 1 - ny) + 2.0 * at(k + 1) + at(k + 1 + ny)) / 12.0;
    }
  }
  else if (x == nx - 1)
  {
    // Corner
    if (y == 0)
    {
      return (4.0 * at(k) + 2.0 * at(k - 1) + 2.0 * at(k + ny)) / 8.0;
    }
    // Corner
    else if (y == ny - 1)
    {
      return (4.0 * at(k) + 2.0 * at(k - 1) + 2.0 * at(k - ny)) / 8.0;
    }
    else
    {
      return (2.0 * at(k - ny) + 4.0 * at(k) + 2.0 * at(k + ny) + at(k - 1 - ny) + 2.0 * at(k - 1) + at(k - 1 + ny)) / 12.0;
    }
  }
  else
  {
    if (y == 0)
    {
      return (4.0 * at(k) + 2.0 * at(k - 1) + 2.0 * at(k + 1) + 2.0 * at(k + nx) + at(k - 1 + nx) + at(k + 1 + nx)) / 12.0;

    }
    else if (y == ny - 1)
    {
      return (4.0 * at(k) + 2.0 * at(k - 1) + 2.0 * at(k + 1) + 2.0 * at(k - nx) + at(k - 1 - nx) + at(k + 1 - nx)) / 12.0;
    }
    // Center
    else
    {
      return (4.0 * at(k) + 2.0 * at(k - 1) + 2.0 * at(k + 1) + 2.0 * at(k - ny) + 2.0 * at(k + ny) + at(k - 1 - ny) + at(k + 1 - ny) + at(k - 1 + ny) + at(k + 1 + ny)) / 16.0;
    }
  }
}

/*!
\brief Applies several smoothing steps to the scalar field.

\sa Smooth()
\param n Number of smoothing steps.
*/
void ScalarField2::Smooth(int n)
{
  for (int i = 0; i < n; i++)
  {
    Smooth();
  }
}

/*!
\brief Smooth the scalar field using a discrete gaussian kernel.

The function uses a 3<SUP>2</SUP> approximation of the Gaussian kernel.
*/
void ScalarField2::Smooth()
{
  QVector<double> smoothed;
  smoothed.resize(nx * ny);

  int k;

  // Smooth center
  for (int i = 1; i < nx - 1; i++)
  {
    for (int j = 1; j < ny - 1; j++)
    {
      k = VertexIndex(i, j);
      smoothed[k] = (4.0 * at(k) + 2.0 * at(k - 1) + 2.0 * at(k + 1) + 2.0 * at(k - nx) + 2.0 * at(k + nx) + at(k - 1 - nx) + at(k + 1 - nx) + at(k - 1 + nx) + at(k + 1 + nx)) / 16.0;
    }
  }

  // Smooth edges
  for (int i = 1; i < nx - 1; i++)
  {
    k = VertexIndex(i, 0);
    smoothed[k] = (4.0 * at(k) + 2.0 * at(k - 1) + 2.0 * at(k + 1) + 2.0 * at(k + nx) + at(k - 1 + nx) + at(k + 1 + nx)) / 12.0;

    k = VertexIndex(i, ny - 1);
    smoothed[k] = (4.0 * at(k) + 2.0 * at(k - 1) + 2.0 * at(k + 1) + 2.0 * at(k - nx) + at(k - 1 - nx) + at(k + 1 - nx)) / 12.0;
  }

  for (int j = 1; j < ny - 1; j++)
  {
    k = VertexIndex(0, j);
    smoothed[k] = (2.0 * at(k - nx) + 4.0 * at(k) + 2.0 * at(k + nx) + at(k + 1 - nx) + 2.0 * at(k + 1) + at(k + 1 + nx)) / 12.0;

    k = VertexIndex(nx - 1, j);
    smoothed[k] = (2.0 * at(k - nx) + 4.0 * at(k) + 2.0 * at(k + nx) + at(k - 1 - nx) + 2.0 * at(k - 1) + at(k - 1 + nx)) / 12.0;
  }

  // Corners
  k = VertexIndex(0, 0);
  smoothed[k] = (4.0 * at(k) + 2.0 * at(k + 1) + 2.0 * at(k + nx) + 1.0 * at(k + nx + 1)) / 9.0;

  k = VertexIndex(nx - 1, 0);
  smoothed[k] = (4.0 * at(k) + 2.0 * at(k - 1) + 2.0 * at(k + nx) + 1.0 * at(k + nx - 1)) / 9.0;

  k = VertexIndex(0, ny - 1);
  smoothed[k] = (4.0 * at(k) + 2.0 * at(k + 1) + 2.0 * at(k - nx) + 1.0 * at(k - nx + 1)) / 9.0;

  k = VertexIndex(nx - 1, ny - 1);
  smoothed[k] = (4.0 * at(k) + 2.0 * at(k - 1) + 2.0 * at(k - nx) + 1.0 * at(k - nx - 1)) / 9.0;

  // Center
  field = smoothed;
}

/*!
\brief Small smoothing of the scalar field using a discrete gaussian kernel.

The function uses a 3<SUP>2</SUP> approximation of the Gaussian kernel.
*/
void ScalarField2::SmoothSmall()
{
  QVector<double> smoothed;
  smoothed.resize(nx * ny);

  int k;

  // Smooth center
  for (int i = 1; i < nx - 1; i++)
  {
    for (int j = 1; j < ny - 1; j++)
    {
      k = VertexIndex(i, j);

      smoothed[k] = (186.0 * at(k) + 3.0 * at(k - 1) + 3.0 * at(k + 1) + 3.0 * at(k - nx) + 3.0 * at(k + nx) + 2.0 * at(k - 1 - nx) + 2.0 * at(k + 1 - nx) + 2.0 * at(k - 1 + nx) + 2.0 * at(k + 1 + nx)) / 206.0;
    }
  }

  // Smooth edges
  for (int i = 1; i < nx - 1; i++)
  {
    k = VertexIndex(i, 0);
    smoothed[k] = (186.0 * at(k) + 3.0 * at(k - 1) + 3.0 * at(k + 1) + 3.0 * at(k + nx) + 2.0 * at(k - 1 + nx) + 2.0 * at(k + 1 + nx)) / 199.0;

    k = VertexIndex(i, ny - 1);
    smoothed[k] = (186.0 * at(k) + 3.0 * at(k - 1) + 3.0 * at(k + 1) + 3.0 * at(k - nx) + 2.0 * at(k - 1 - nx) + 2.0 * at(k + 1 - nx)) / 199.0;
  }

  for (int j = 1; j < ny - 1; j++)
  {
    k = VertexIndex(0, j);
    smoothed[k] = (3.0 * at(k - nx) + 186.0 * at(k) + 3.0 * at(k + nx) + 2.0 * at(k + 1 - nx) + 3.0 * at(k + 1) + 2.0 * at(k + 1 + nx)) / 199.0;

    k = VertexIndex(nx - 1, j);
    smoothed[k] = (3.0 * at(k - nx) + 186.0 * at(k) + 3.0 * at(k + nx) + 2.0 * at(k - 1 - nx) + 3.0 * at(k - 1) + 2.0 * at(k - 1 + nx)) / 199.0;
  }

  // Corners
  k = VertexIndex(0, 0);
  smoothed[k] = (186.0 * at(k) + 3.0 * at(k + 1) + 3.0 * at(k + nx) + 2.0 * at(k + nx + 1)) / 194.0;

  k = VertexIndex(nx - 1, 0);
  smoothed[k] = (186.0 * at(k) + 3.0 * at(k - 1) + 3.0 * at(k + nx) + 2.0 * at(k + nx - 1)) / 194.0;

  k = VertexIndex(0, ny - 1);
  smoothed[k] = (186.0 * at(k) + 3.0 * at(k + 1) + 3.0 * at(k - nx) + 2.0 * at(k - nx + 1)) / 194.0;

  k = VertexIndex(nx - 1, ny - 1);
  smoothed[k] = (186.0 * at(k) + 3.0 * at(k - 1) + 3.0 * at(k - nx) + 2.0 * at(k - nx - 1)) / 194.0;

  // Center
  field = smoothed;
}

/*!
\brief Filter with median filter.
*/
void ScalarField2::MedianFilter()
{
  QVector<double> filtered;
  filtered.resize(nx * ny);

  // Smooth center
  for (int i = 1; i < nx - 1; i++)
  {
    for (int j = 1; j < ny - 1; j++)
    {
      int k = VertexIndex(i, j);

      std::vector<double> list;
      list.push_back(at(k));
      list.push_back(at(k - 1));
      list.push_back(at(k + 1));
      list.push_back(at(k - nx));
      list.push_back(at(k + nx));
      list.push_back(at(k - 1 - nx));
      list.push_back(at(k + 1 - nx));
      list.push_back(at(k - 1 + nx));
      list.push_back(at(k + 1 + nx));
      sort(list.begin(), list.end());

      filtered[k] = list[4];
    }
  }

  // Smooth edges
  for (int i = 1; i < nx - 1; i++)
  {
    int k = VertexIndex(i, 0);

    std::vector<double> listA;
    listA.push_back(at(k));
    listA.push_back(at(k - 1));
    listA.push_back(at(k + 1));
    listA.push_back(at(k + nx));
    listA.push_back(at(k - 1 + nx));
    listA.push_back(at(k + 1 + nx));
    sort(listA.begin(), listA.end());
    filtered[k] = listA[2];

    k = VertexIndex(i, ny - 1);

    std::vector<double> listB;
    listB.push_back(at(k));
    listB.push_back(at(k - 1));
    listB.push_back(at(k + 1));
    listB.push_back(at(k - nx));
    listB.push_back(at(k - 1 - nx));
    listB.push_back(at(k + 1 - nx));
    sort(listB.begin(), listB.end());
    filtered[k] = listB[2];
  }

  for (int j = 1; j < ny - 1; j++)
  {
    int k = VertexIndex(0, j);

    std::vector<double> listA;
    listA.push_back(at(k));
    listA.push_back(at(k + 1));
    listA.push_back(at(k - nx));
    listA.push_back(at(k + nx));
    listA.push_back(at(k + 1 - nx));
    listA.push_back(at(k + 1 + nx));
    sort(listA.begin(), listA.end());
    filtered[k] = listA[2];

    k = VertexIndex(nx - 1, j);

    std::vector<double> listB;
    listB.push_back(at(k));
    listB.push_back(at(k - 1));
    listB.push_back(at(k - nx));
    listB.push_back(at(k + nx));
    listB.push_back(at(k - 1 - nx));
    listB.push_back(at(k - 1 + nx));
    sort(listB.begin(), listB.end());
    filtered[k] = listB[2];
  }

  // Corners
  int k = VertexIndex(0, 0);
  std::vector<double> list;
  list.push_back(at(k));
  list.push_back(at(k + 1));
  list.push_back(at(k + nx));
  list.push_back(at(k + 1 + nx));
  sort(list.begin(), list.end());
  filtered[k] = list[1];

  k = VertexIndex(nx - 1, 0);
  list.clear();
  list.push_back(at(k));
  list.push_back(at(k - 1));
  list.push_back(at(k + nx));
  list.push_back(at(k - 1 + nx));
  sort(list.begin(), list.end());
  filtered[k] = list[1];

  k = VertexIndex(0, ny - 1);
  list.clear();
  list.push_back(at(k));
  list.push_back(at(k + 1));
  list.push_back(at(k - nx));
  list.push_back(at(k + 1 - nx));
  sort(list.begin(), list.end());
  filtered[k] = list[1];

  k = VertexIndex(nx - 1, ny - 1);
  list.clear();
  list.push_back(at(k));
  list.push_back(at(k - 1));
  list.push_back(at(k - nx));
  list.push_back(at(k - 1 - nx));
  sort(list.begin(), list.end());
  filtered[k] = list[1];

  // Center
  field = filtered;
}

/*!
\brief Blurs the scalar field.

The function uses a 3<SUP>2</SUP> blur kernel.
*/
void ScalarField2::Blur()
{
  QVector<double> smoothed;
  smoothed.resize(nx * ny);

  int k;

  // Smooth center
  for (int i = 1; i < nx - 1; i++)
  {
    for (int j = 1; j < ny - 1; j++)
    {
      k = VertexIndex(i, j);

      smoothed[k] = (at(k) + at(k - 1) + at(k + 1) + at(k - nx) + at(k + nx) + at(k - 1 - nx) + at(k + 1 - nx) + at(k - 1 + nx) + at(k + 1 + nx)) / 9.0;
    }
  }

  // Smooth edges
  for (int i = 1; i < nx - 1; i++)
  {
    k = VertexIndex(i, 0);
    smoothed[k] = (at(k) + at(k - 1) + at(k + 1) + at(k + nx) + at(k - 1 + nx) + at(k + 1 + nx)) / 6.0;

    k = VertexIndex(i, ny - 1);
    smoothed[k] = (at(k) + at(k - 1) + at(k + 1) + at(k - nx) + at(k - 1 - nx) + at(k + 1 - nx)) / 6.0;
  }

  for (int j = 1; j < ny - 1; j++)
  {
    k = VertexIndex(0, j);
    smoothed[k] = (at(k - nx) + at(k) + at(k + nx) + at(k + 1 - nx) + at(k + 1) + at(k + 1 + nx)) / 6.0;

    k = VertexIndex(nx - 1, j);
    smoothed[k] = (at(k - nx) + at(k) + at(k + nx) + at(k - 1 - nx) + at(k - 1) + at(k - 1 + nx)) / 6.0;
  }

  // Corners
  k = VertexIndex(0, 0);
  smoothed[k] = (at(k) + at(k + 1) + at(k + nx)) / 3.0;

  k = VertexIndex(nx - 1, 0);
  smoothed[k] = (at(k) + at(k - 1) + at(k + nx)) / 3.0;

  k = VertexIndex(0, ny - 1);
  smoothed[k] = (at(k) + at(k + 1) + at(k - nx)) / 3.0;

  k = VertexIndex(nx - 1, ny - 1);
  smoothed[k] = (at(k) + at(k - 1) + at(k - nx)) / 3.0;

  // Center
  field = smoothed;
}

/*!
\brief Applies several blurring steps to the scalar field.

\sa Blur()
\param n Number of blurring steps.
*/
void ScalarField2::Blur(int n)
{
  for (int i = 0; i < n; i++)
  {
    Blur();
  }
}
void ScalarField2::Debug() const
{
  std::cout << field.size() << "   / " << nx << ", " << ny << std::endl;
}