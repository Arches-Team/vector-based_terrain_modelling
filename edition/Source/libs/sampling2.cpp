#include <QtWidgets/QGraphicsTextItem>
#include <QtWidgets/QGraphicsScene>

#include "libs/sampling.h"
#include "libs/cubic.h"
#include "libs/arrayinteger.h"
#include "libs/shuffle.h"
#include "libs/scalarfield.h"

/*!
\class DiscTile sampling.h
\brief Poisson tiling.

Can be used for importance sampling
from Template Poisson disk tiles, Ares Lagae, Philip DutrÃ©.
See Ares Lagae, Philip DutrÃ©. A procedural object distribution function. <I>ACM Transactions on Graphics</I>. <B>24</B>(4), 1442â€“1461.

\ingroup PlanarGroup
*/

/*!
\brief Initialize a poisson disc tiling.
*/
DiscTile::DiscTile()
{
}

/*!
\brief Create a Poisson sampling in a square domain.
\param s Size of the side of the square.
\param r Radius of the samples.
\param v Set of points.
*/
DiscTile::DiscTile(const double& s, const double& r, const QVector<Vector2>& v) :s(s), r(r), e(r* r * 4.02), p(v)
{
}

/*!
\brief Create a poisson sampling in a square domain.
\param s Size of the side of the square.
\param r Radius of the samples.
\param n Maximum number of samples in the dart throwing process. Note that the number of samples will be different.
\param ra %Random number generator.
*/
DiscTile::DiscTile(const double& s, const double& r, int n, Random& ra) :s(s), r(r), e(r* r * 4.02)
{
  DiscTile::Generate(n, ra);
}

/*!
\brief Create a poisson sampling in a square domain.
\param s Size of the side of the square.
\param r Radius of the samples.
\param n Maximum number of samples in the dart throwing process. Note that the number of samples will be different.
\param p Set of samples that will be in the tile.
\param ra %Random number generator.
*/
DiscTile::DiscTile(const double& s, const double& r, int n, const QVector<Vector2>& p, Random& ra) :s(s), r(r), e(r* r * 4.02), p(p)
{
  DiscTile::Generate(n, ra);
}


/*!
\brief Scramble the points inside the generated tile.

Use a modified version of the Fisher-Yates algorithm.
\sa Shuffle::Table
*/
void DiscTile::Scramble()
{
  QVector<int> shuffle = Shuffle::Table(p.size());

  // New array
  QVector<Vector2> q(p.size());

  for (int i = 0; i < p.size(); i++)
  {
    q[i] = p[shuffle.at(i)];
  }
  p = q;
}

/*!
\brief Perform a relaxation step over the Poisson-Disc distribution.
\param c Relaxation coefficient, it the fraction of the radius of the Poisson distribution.
\param n Number of steps.
*/
void DiscTile::Relaxation(int n, const double& c)
{
  for (int i = 0; i < n; i++)
  {
    Relaxation(c);
  }
}

/*!
\brief Perform a relaxation step over the Poisson-Disc distribution.
\param c Relaxation coefficient, it is the fraction of the Poisson radius.
*/
void DiscTile::Relaxation(const double& c)
{
  // Size of the square
  const double t = 2.0 * r / sqrt(2.0);
  int mod = int(s / t) + 1;

  // Grid for speeding up queries
  QVector<int> grid(mod * mod);
  grid.fill(-1, mod * mod);

  for (int k = 0; k < p.size(); k++)
  {
    int xk = (p.at(k)[0] * mod / s);
    int yk = (p.at(k)[1] * mod / s);
    grid[xk * mod + yk] = k;
  }

  // Displacements
  QVector<Vector2> d(p.size());

  // Compute displacements
  for (int k = 0; k < p.size(); k++)
  {
    Vector2 eps = Vector2::Null;
    int xk = (p.at(k)[0] * mod / s);
    int yk = (p.at(k)[1] * mod / s);
    for (int i = -3; i < 4; i++)
    {
      for (int j = -3; j < 4; j++)
      {
        int ox = 0;
        int oy = 0;

        if (i + xk < 0) { ox = -1; }
        if (i + xk >= mod) { ox = 1; }
        if (j + yk < 0) { oy = -1; }
        if (j + yk >= mod) { oy = 1; }

        int gi = (i + xk + mod) % mod;
        int gj = (j + yk + mod) % mod;
        int ng = grid.at(gi * mod + gj);
        if ((ng != -1) && (ng != k))
        {
          Vector2 dis = p.at(ng) + Vector2(ox * s, oy * s) - p.at(k);
          double dd = dis * dis;

          // Move by a small amount of radius, smoothed according to distance
          eps += -r * c * Normalized(dis) * Cubic::SmoothCompact(dd, 16.0 * r * r);
        }
      }
    }
    // Move
    d[k] = eps;
  }

  // Apply displacements
  for (int k = 0; k < p.size(); k++)
  {
    p[k] += d[k];

    // Keep point within tile
    p[k] = Vector2::Mod(p[k], Vector2(s, s));
  }
}

/*!
\brief Check if a disc intersect existing discs inside the tile, without toric constraints.
\param q Point.
*/
bool DiscTile::Check(const Vector2& q) const
{
  bool c = false;
  for (int j = 0; j < p.size(); j++)
  {
    if (SquaredNorm(p.at(j) - q) < e)
    {
      c = true;
      break;
    }
  }
  return c;
}

/*!
\brief Check if a disc intersect existing discs inside the tile, with or without toric constraints.

This complexity is O(n), thus the function is expensive should it be called many times.

\param q Point.
\param t Toric constraint, set to true as default.
*/
bool DiscTile::Intersect(const Vector2& q, bool t) const
{
  // Simple test with discs inside the tile if toric flag is false
  if (t == false)
  {
    return Check(q);
  }

  bool c = Check(q);

  // Occurs in toric space if the point is away from the boundary
  if (fabs(GetBox().Signed(q)) > 2.0 * r)
  {
    return c;
  }

  // Toric space, moreover the point is close to the boundary

  if (fabs(q[0] - 0.0) < e)
  {
    c |= Check(q + Vector2(s, 0.0));
    if (c == true) return c;

    // Bottom left
    if (fabs(q[1] - 0.0) < e)
    {
      c |= Check(q + Vector2(s, s));
      c |= Check(q + Vector2(0.0, s));
    }
    // Top left
    else if (fabs(q[1] - s) < e)
    {
      c |= Check(q + Vector2(s, -s));
      c |= Check(q + Vector2(0.0, -s));
    }
    // Left 
    else
    {
    }
  }
  else if (fabs(q[0] - s) < e)
  {
    c |= Check(q + Vector2(-s, 0.0));
    if (c == true) return c;

    // Bottom right
    if (fabs(q[1] - 0.0) < e)
    {
      c |= Check(q + Vector2(-s, s));
      c |= Check(q + Vector2(0.0, s));
    }
    // Top right
    else if (fabs(q[1] - s) < e)
    {
      c |= Check(q + Vector2(-s, -s));
      c |= Check(q + Vector2(0.0, -s));
    }
    // Right 
    else
    {
    }
  }
  else
  {
    // Bottom 
    if (fabs(q[1] - 0.0) < e)
    {
      c |= Check(q + Vector2(s, 0.0));
    }
    // Top
    else
    {
      c |= Check(q + Vector2(-s, 0.0));
    }
  }
  return c;

}


/*!
\brief Generate the poisson sampling of the domain.
\param n Number of darts thrown on the domain. Note that the number of generated samples will be less than n.
\param random %Random number generator.
*/
void DiscTile::Generate(int n, Random& random)
{
  // Size of the square
  double t = 2.0 * r / sqrt(2.0);
  int mod = int(s / t) + 1;

  QVector<int> grid(mod * mod, -1);

  // Locate existing samples into the grid
  for (int k = 0; k < p.size(); k++)
  {
    double x = p.at(k)[0];
    double y = p.at(k)[1];

    int xk = (x * mod / s);
    int yk = (y * mod / s);

    grid[xk * mod + yk] = k;
  }
  static const QPoint set[25] = {
    QPoint(-2,-2),QPoint(-1,-2), QPoint(0,-2), QPoint(1,-2), QPoint(2,-2),
    QPoint(-2,-1),QPoint(-1,-1), QPoint(0,-1), QPoint(1,-1), QPoint(2,-1),
    QPoint(-2, 0),QPoint(-1, 0), QPoint(0, 0), QPoint(1, 0), QPoint(2, 0),
    QPoint(-2,+1),QPoint(-1,+1), QPoint(0,+1), QPoint(1,+1), QPoint(2,+1),
    QPoint(-2,+2),QPoint(-1,+2), QPoint(0,+2), QPoint(1,+2), QPoint(2,+2)
  };
  // Try to sample with n discs
  for (int k = 0; k < n; k++)
  {
    double x = random.Uniform(s);
    double y = random.Uniform(s);

    int xk = (x * mod / s);
    int yk = (y * mod / s);

    // Detect collisions with existing discs
    bool c = false;

    for (int l = 0; l < 25; l++)
    {
      int i = set[l].x();
      int j = set[l].y();
      int ox = 0;
      int oy = 0;

      if (i + xk < 0) { ox = -1; }
      if (i + xk >= mod) { ox = 1; }
      if (j + yk < 0) { oy = -1; }
      if (j + yk >= mod) { oy = 1; }

      int gi = (i + xk + 2 * mod) % mod;
      int gj = (j + yk + 2 * mod) % mod;
      int ng = grid.at(gi * mod + gj);
      if (ng != -1)
      {
        if (SquaredNorm(p.at(ng) + Vector2(ox * s, oy * s) - Vector2(x, y)) < e)
        {
          c = true;
          break;
        }
      }
      if (c == true) break;

    }

    if (c == false)
    {
      grid[xk * mod + yk] = p.size();
      p.append(Vector2(x, y));
    }
  }
}

/*!
\brief Draw the distribution.
\param scene The graphics scene.
*/
void DiscTile::Draw(QGraphicsScene& scene) const
{
  QFont serifFont("Calibri", 10, QFont::Light);

  QPen pen(QColor(200, 150, 150));
  QBrush brush(QColor(240, 240, 240));
  pen.setWidthF(s * 0.002);

  double e = s * 0.25;

  scene.setSceneRect(-e, -e, s + 2.0 * e, s + 2.0 * e);
  scene.addRect(0.0, 0.0, s, s, pen);

  pen.setWidthF(s * 0.002);

  Box2 box = GetBox();
  for (int k = 0; k < p.size(); k++)
  {
    double d = fabs(box.Signed(p.at(k)));
    if (d < r)
    {
      pen.setColor(QColor(150, 100, 130));
      brush.setColor(QColor(230, 180, 210));
    }
    else if (d < 2.0 * r)
    {
      pen.setColor(QColor(130, 100, 130));
      brush.setColor(QColor(210, 180, 210));
    }
    else
    {
      pen.setColor(QColor(100, 100, 130));
      brush.setColor(QColor(180, 180, 210));
    }

    scene.addEllipse(p.at(k)[0] - r, p.at(k)[1] - r, 2.0 * r, 2.0 * r, pen, brush);

    QString num = QString("%1").arg(k);
    QGraphicsTextItem* item = scene.addText(num, serifFont);
    item->setScale(0.025 * 2.0 * r);
    item->setPos(p.at(k)[0] - 0.5 * r, p.at(k)[1] - 0.5 * r);
  }
}

/*!
\brief Sample a domain according to an input importance function.
\param b The box.
\param f Importance function.
*/
QVector<Vector2> DiscTile::Sample(const Box2& b, double(*f)(const Vector2&)) const
{
  QVector<Vector2> samples;
  Box2 tile = GetBox();

  // Get range of tile across the domain
  const QRect rect = b.TileRange(tile);

  for (int j = rect.y(); j < rect.y() + rect.height(); j++)
  {
    for (int i = rect.x(); i < rect.x() + rect.width(); i++)
    {
      Box2 tiled = tile.Tile(i, j);
      for (int k = 0; k < p.size(); k++)
      {
        Vector2 pk = p.at(k) + tiled[0];
        if (b.Inside(pk))
        {
          // Call to the importance function
          double importance = (*f)(pk);

          // Not enough important : discarded
          if (!(importance > Math::Unit(k, p.size())))
            continue;

          samples.append(pk);
        }
      }
    }
  }
  return samples;
}

/*!
\brief Sample a domain according to an input importance function.
\param s Scalar field defining the importance function.
*/
QVector<Vector2> DiscTile::Sample(const ScalarField2& s) const
{
  QVector<Vector2> samples;
  Box2 tile = GetBox();

  Box2 box = s.Array2::GetBox();

  // Get range of tile across the domain
  const QRect rect = box.TileRange(tile);

  for (int j = rect.y(); j < rect.y() + rect.height(); j++)
  {
    for (int i = rect.x(); i < rect.x() + rect.width(); i++)
    {
      Box2 tiled = tile.Tile(i, j);
      for (int k = 0; k < p.size(); k++)
      {
        Vector2 pk = p[k] + tiled[0];
        if (box.Inside(pk))
        {
          // Call to the importance function
          double importance = s.Value(pk);

          // Not enough important : discarded : use strict > for comparison so that 0 importance does not generate sample
          if (!(importance > Math::Unit(k, p.size())))
            continue;

          samples.append(pk);
        }
      }
    }
  }
  return samples;
}

/*!
\brief Sample a domain according to an input importance function.
\param box The box.
*/
QVector<Vector2> DiscTile::Sample(const Box2& box) const
{
  QVector<Vector2> samples;
  Box2 tile = GetBox();

  // Get range of tile across the domain
  const QRect rect = box.TileRange(tile);

  for (int j = rect.y(); j < rect.y() + rect.height(); j++)
  {
    for (int i = rect.x(); i < rect.x() + rect.width(); i++)
    {
      Box2 tiled = tile.Tile(i, j);
      for (int k = 0; k < p.size(); k++)
      {
        Vector2 pk = p[k] + tiled[0];
        if (box.Inside(pk))
        {
          samples.append(pk);
        }
      }
    }
  }
  return samples;
}

/*
\brief Sample a domain according to an input importance function.
\param circle The circle.
*/
QVector<Vector2> DiscTile::Sample(const Circle2& circle) const
{
  Box2 box = circle.GetBox();

  QVector<Vector2> samples;
  Box2 tile = GetBox();

  // Get range of tile across the domain
  const QRect rect = box.TileRange(tile);

  for (int j = rect.y(); j < rect.y() + rect.height(); j++)
  {
    for (int i = rect.x(); i < rect.x() + rect.width(); i++)
    {
      Box2 tiled = tile.Tile(i, j);
      for (int k = 0; k < p.size(); k++)
      {
        Vector2 pk = p.at(k) + tiled[0];
        if (circle.Inside(pk))
        {
          samples.append(pk);
        }
      }
    }
  }
  return samples;
}

/*!
\brief Performs a Poisson disc check on a set of point and a candidate position.

Returns true if the candidate intersects with the distribution, false otherwise.
\param p Candidate position.
\param s Set of point.
\param r Poisson disc radius.
*/
bool DiscTile::Check(const Vector2& p, const QVector<Vector2>& s, double r)
{
  // Square radius
  r *= r;
  for (int i = 0; i < s.size(); i++)
  {
    if (SquaredNorm(p - s.at(i)) < r)
    {
      return true;
    }
  }
  return false;
}

/*!
\brief Compute the packing density.

\sa DiscTile::Packing()
*/
double DiscTile::Packing() const
{
  return p.size() * Circle::Area(r) / s * s;
}

/*!
\brief Overloaded.
\param s Stream.
\param disctile The disc tile.
*/
std::ostream& operator<<(std::ostream& s, const DiscTile& disctile)
{
  s << "DiscTile(" << disctile.s << ", " << disctile.r << ") // n = " << disctile.p.size();
  for (const Vector2& v : disctile.p)
  {
    s << v << std::endl;
  }
  return s;
}
