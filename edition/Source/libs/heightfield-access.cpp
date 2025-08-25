// Heightfield

#include "libs/heightfield.h"

#include "libs/hemisphere.h"

/*!
\brief Compute the Lipschitz constant of the signed distance field defined as f(x,y,z)=z-h(x,y).
*/
double HeightField::K() const
{
  // Lipschitz
  double lipschitz = ScalarField2::K();

  return sqrt(1.0 + lipschitz * lipschitz);
}

/*!
\brief Compute the accessibility.

\param r Radius.
\param n Number of rays.
*/
ScalarField2 HeightField::Accessibility(const double& r, int n) const
{
  double epsilon = 0.05;

  // Lipschitz
  double lipschitz = K();

  Box box = GetBox();

  ScalarField2 sf(Box2(a, b), nx, ny, 1.0);

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      Vector p = Vertex(i, j) + Vector(0.0, 0.0, epsilon);
      Vector normal = Normal(i, j);

      int hit = 0;
      for (int k = 0; k < n; k++)
      {
        Vector direction = HemiSphere::RandomDirection(normal);
        Ray ray(p, direction);
        double t;
        Vector q;
        if (Intersect(ray, t, q, box, lipschitz, r, epsilon / 2.0))
        {
          hit++;
        }
      }
      sf(i, j) = 1.0 - double(hit) / double(n);
    }
  }

  return sf;
}

/*!
\brief Compute the direct lighting.

\param u Light direction.
\param cosine Boolean, set to true to use scalar product between normal and light, false to use (1+c)/2.
*/
ScalarField2 HeightField::Light(const Vector& u, bool cosine) const
{
  ScalarField2 sf(Box2(a, b), nx, ny, 1.0);

  if (cosine == true)
  {
    for (int i = 0; i < nx; i++)
    {
      for (int j = 0; j < ny; j++)
      {
        double light = Normal(i, j) * u;

        light = Math::Max(0.0, light);


        sf(i, j) = light;
      }
    }
  }
  else
  {
    for (int i = 0; i < nx; i++)
    {
      for (int j = 0; j < ny; j++)
      {
        double light = Normal(i, j) * u;

        light = 0.5 * (1.0 + light);

        sf(i, j) = light;
      }
    }

  }

  return sf;
}

/*!
\brief Create a data structure for DX computation.
\param hf Heightfield.
*/
TerrainDX::TerrainDX(const HeightField& hf) :Array2(hf.GetArray()), terrain(hf)
{
  dirx = Array2I(hf.GetArray());
  accumulaton = ScalarField2(hf.GetArray());

  // Compute flow direction.

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      if (!(i == 0 || j == 0 || i == nx - 1 || j == ny - 1))
      {
        // Recherche de la pente minimum
        double cHeight = terrain.at(i, j);
        double min = 0.0;
        for (int dir = 1; dir < 256; dir <<= 1)
        {

          QPoint q = CodeToDir(dir);

          if (!(q.x() == 0 && q.y() == 0))
          {
            double dz = terrain.at(i + q.x(), j + q.y()) - cHeight;
            double s = dz / sqrt(double(q.x() * q.x() + q.y() * q.y()));

            if (s < min)
            {
              min = s;
              dirx(i, j) = dir;
            }
          }
        }
      }
    }
  }

  // Compute area accumulation
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      int it = i;
      int jt = j;

      int cpt = 0;
      // tant que la direction n'est pas nulle
      while (dirx(it, jt) != 0)
      {
        accumulaton(it, jt)++;

        int c = dirx(it, jt);
        QPoint q = CodeToDir(c);

        it += q.x();
        jt += q.y();

        cpt++;
      }
    }
  }
}




/*!
\brief Compute flow.

Uses a power law to compute flow from drainage area.
*/
double TerrainDX::getMeanFlow(int i, int j) const
{
  double area = accumulaton.Cell(0).Area();
  return 0.42 * pow(double(accumulaton.at(i, j)) * area, 0.69);
}

/*!
\brief Get the coded list direction of river input.
\return List of river direction.
*/
QVector<int> TerrainDX::getRiverInput(int i, int j) const
{
  QVector<int> list;

  for (int dir = 1; dir < 256; dir <<= 1)
  {
    QPoint q = CodeToDir(dir);
    if (isRiver(i + q.x(), j + q.y()))
    {
      int c;
      c = dirx.at(i + q.x(), j + q.y());
      QPoint qq = CodeToDir(c);

      if (q.x() + qq.x() == 0 && q.y() + qq.y() == 0)
        list.append(c);
    }
  }
  return list;
}

/*!
\brief Get the coded list direction of flow input.
\return List of flow direction.
*/
bool TerrainDX::isSource(int i, int j) const
{
  QVector<int> list;
  if (isRiver(i, j))
  {
    for (int dir = 1; dir < 256; dir <<= 1)
    {
      QPoint q = CodeToDir(dir);

      if (isRiver(i + q.x(), j + q.y()))
      {
        int c = dirx.at(i + q.x(), j + q.y());
        QPoint qq = CodeToDir(c);

        if (q.x() + qq.x() == 0 && q.y() + qq.y() == 0)
          list.append(c);
      }
    }
    return list.size() == 0;
  }
  else return false;
}

/*!
\brief Get the coded list direction of flow input.
\return List of flow direction.
*/
bool TerrainDX::isRiverBorder(int i, int j) const
{
  QVector<int> list;
  if (isRiver(i, j))
  {
    int c = getRiverOutput(i, j);

    QPoint q = CodeToDir(c);

    if ((i + q.x() == 0) || (j + q.y() == 0) || (i + q.x() == nx - 1) || (j + q.y() == ny - 1))
      return true;

    return false;
  }
  else return false;
}

/*!
\brief Get the coded list direction of flow input.
\return List of flow direction.
*/
QVector<int> TerrainDX::getFlowInput(int i, int j) const
{
  QVector<int> list;

  for (int dir = 1; dir < 256; dir <<= 1)
  {
    int c;
    QPoint q = CodeToDir(dir);

    c = dirx.at(i + q.x(), j + q.y());
    QPoint qq = CodeToDir(c);
    if (q.x() + qq.x() == 0 && q.y() + qq.y() == 0)
      list.append(c);
  }
  return list;
}

/*!
\brief Get the coded direction of the river flow output.
\return Flow direction.
*/
int TerrainDX::getRiverOutput(int i, int j) const
{
  return dirx.at(i, j);
}



