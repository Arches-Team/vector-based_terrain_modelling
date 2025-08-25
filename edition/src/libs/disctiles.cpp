#include <QtWidgets/QGraphicsTextItem>
#include <QtWidgets/QGraphicsScene>

#include "libs/sampling.h"
#include "libs/cubic.h"
#include "libs/arrayinteger.h"
#include "libs/draw.h"
#include "libs/cpu.h"
#include "libs/vectorset.h"
#include "libs/octogon.h"

/*
\class DiscTiles sampling.h
\brief Poisson disc sampling with two tiles and aperiodic tiling.
*/

/*!
\brief Tile selection.
\param x,y Integer coordinates of the tile.
*/
int DiscTiles::Select(int x, int y) const
{
  int a = (x * 139 + y);
  // Specific hashing
  a = (a ^ 61) ^ (a >> 16);
  a = a + (a << 3);
  a = a ^ (a >> 4);
  a = a * 0x27d4eb2d;
  a = a ^ (a >> 15);
  return a % 2;
}

/*!
\brief Return the number of samples.
\param t Tile.
*/
int DiscTiles::Size(int t) const
{
  return tiles[t].Size();
}

/*!
\brief Return the box.
*/
Box2 DiscTiles::GetBox() const
{
  return Box2(Vector2::Null, Vector2(s));
}

/*!
\brief Getter on the i-th sample.
\param t Tile.
\param i Sample index.
*/
Vector2 DiscTiles::Vertex(int t, int i) const
{
  return tiles[t].Vertex(i);
}

/*!
\brief Compute the i-th sample for a displaced tile.
\param i Sample index, should be lower than the size of tile at (x,y).
\param x,y Integer coordinates of the tile.
*/
Vector2 DiscTiles::Vertex(int x, int y, int i) const
{
  int c = Select(x, y);

  return tiles[c].Vertex(x, y, i);

}

/*!
\brief Create a Poisson sampling in a square domain.
\param s Size of the side of the square.
\param r Radius of the samples.
\param n Maximum number of samples in the dart throwing process. Note that the number of samples will be different.
\param ra %Random number generator.
*/
DiscTiles::DiscTiles(const double& s, const double& r, int n, Random& ra) : s(s), r(r), e(r* r * 4.02)
{
  // Generate first tile
  tiles[0] = DiscTile(s, r, n, ra);

  // Select points on the boundary
  QVector<Vector2> p;
  Box2 box = tiles[0].GetBox();
  for (int i = 0; i < tiles[0].Size(); i++)
  {
    // Select only those points that are close to the boundary
    Vector2 pi = tiles[0].Vertex(i);
    if (fabs(box.Signed(pi)) < 3.0 * r)
    {
      p.append(pi);
    }
  }

  // Create second tile with boundary points as constraints
  tiles[1] = DiscTile(s, r, n, p, ra);

  // At this stage, we scramble the distribution as we started from the boundary discs
  tiles[1].Scramble();
}


/*!
\brief Sample a domain according to an input importance function.
\param s Scalar field defining the importance function.
*/
QVector<Vector2> DiscTiles::Sample(const ScalarField2& s) const
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
      int c = Select(i, j);

      const DiscTile& selected = tiles[c];
      for (int k = 0; k < selected.Size(); k++)
      {
        Vector2 pk = selected.Vertex(k) + tiled[0];
        if (box.Inside(pk))
        {
          // Call to the importance function
          double importance = s.Value(pk);

          // Not enough important : discarded : use strict > for comparison so that 0 importance does not generate sample
          if (!(importance > Math::Unit(k, selected.Size())))
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
QVector<Vector2> DiscTiles::Sample(const Box2& box) const
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
      int c = Select(i, j);

      const DiscTile& selected = tiles[c];
      for (int k = 0; k < selected.Size(); k++)
      {
        Vector2 pk = selected.Vertex(k) + tiled[0];
        if (box.Inside(pk))
        {
          samples.append(pk);
        }
      }
    }
  }
  return samples;
}


/*!
\brief Compute the average packing density.

\sa DiscTile::Packing()
*/
double DiscTiles::Packing() const
{
  return (tiles[0].Packing() + tiles[1].Packing()) / 2.0;
}


void DiscTiles::DrawDebug() const
{
  QGraphicsScene scene;
  tiles[0].Draw(scene);
  Draw::CreateVector(scene, QString("../tile1.pdf"));
  scene.clear();
  tiles[1].Draw(scene);
  Draw::CreateVector(scene, QString("../tile2.pdf"));

  scene.clear();

  QPen pen(QColor(200, 150, 150));
  QBrush brush(QColor(240, 240, 240));
  pen.setWidthF(s * 0.005);

  Box2 a = GetBox().Extended(2.0 * s);

  scene.setSceneRect(a[0][0], a[0][1], a.Diagonal()[0], a.Diagonal()[1]);
  a.Draw(scene, pen);

  QVector<Vector2> p = Sample(a);

  for (int i = 0; i < p.size(); i++)
  {

    Disc2(p.at(i), r).Draw(scene, pen);

  }

  Draw::CreateVector(scene, QString("../tile3.pdf"));

}


// -------------------------------------------------------------------------------

/*
\class DiscTileCorner sampling.h
\brief Poisson disc sampling with sixteen tiles and aperiodic tiling.

Tiles are created using the Corner Cube algorithm, see:

A. Lagae, P. DutrÃ©. An alternative for Wang tiles: colored edges versus colored corners.
<i>ACM Transactions on Graphics<i> <b>25</b> (4), 1442-1459.

A. Lagae, P. DutrÃ©. A procedural object distribution function. <i>ACM Transactions on Graphics<i> <b>24</b> (4), 2005.

*/


inline int SelectVertex(int x, int y)
{
  int a = (x * 139 + y);
  // Specific hashing
  a = (a ^ 61) ^ (a >> 16);
  a = a + (a << 3);
  a = a ^ (a >> 4);
  a = a * 0x27d4eb2d;
  a = a ^ (a >> 15);
  return a % 2;
}

/*!
\brief Tile selection.
\param x,y Integer coordinates of the tile.
*/
int DiscTileCorner::Select(int x, int y) const
{
  return SelectVertex(x, y) + (SelectVertex(x + 1, y) << 1) + (SelectVertex(x, y + 1) << 2) + (SelectVertex(x + 1, y + 1) << 3);
}

/*!
\brief Return the number of samples.
\param t Tile.
*/
int DiscTileCorner::Size(int t) const
{
  return tiles[t].Size();
}

/*!
\brief Return the box.
*/
Box2 DiscTileCorner::GetBox() const
{
  return Box2(Vector2::Null, Vector2(s));
}

/*!
\brief Getter on the i-th sample.
\param t Tile.
\param i Sample index.
*/
Vector2 DiscTileCorner::Vertex(int t, int i) const
{
  return tiles[t].Vertex(i);
}

/*!
\brief Compute the i-th sample for a displaced tile.
\param i Sample index, should be lower than the size of tile at (x,y).
\param x,y Integer coordinates of the tile.
*/
Vector2 DiscTileCorner::Vertex(int x, int y, int i) const
{
  int c = Select(x, y);

  return tiles[c].Vertex(x, y, i);
}


/*!
\brief Create a Poisson sampling in a square domain.
\param s Size of the side of the square.
\param r Radius of the samples.
\param n Maximum number of samples in the dart throwing process. Note that the number of samples will be different.
\param ra %Random number generator.
*/
DiscTileCorner::DiscTileCorner(const double& s, const double& r, int n, Random& ra) : s(s), r(r), e(r* r * 4.02)
{
  // Vertexes sub-tiles
  VectorSet2 vertex[2];

  const double rad = (std::sqrt(4 + 2 * std::sqrt(2)) / 2) * 2 * r;
  Octogon2 corner(rad);
  const double rad_insc = corner.InscribedRadius();

  // Need to do a poisson sampling over an entire tile, and then cut out the excess points, otherwise the density is too high
  Box2 cbb(s, s);
  vertex[0] = VectorSet2(cbb.Poisson(r, n, ra)).Cut(corner);
  vertex[1] = VectorSet2(cbb.Poisson(r, n, ra)).Cut(corner);

  // Edges sub-tiles
  VectorSet2 edge[8];

  // Need to do a poisson sampling over 2 tiles, and then cut out the excess points, otherwise the density is too high
  Box2 horizontal = Box2(s - rad_insc * 2, 2 * r).Translated(Vector2(s / 2.0, 0.0));
  int ne = 2 * n;
  Box2 hbb = Box2(2 * s, s).Translated(Vector2(s / 2.0, 0.0));
  edge[0] = VectorSet2(hbb.Poisson(r, ne, (vertex[0] + vertex[0].Translated(Vector2(s, 0.0))).Get(), false, ra)).Cut(horizontal); // H00 
  edge[1] = VectorSet2(hbb.Poisson(r, ne, (vertex[1] + vertex[0].Translated(Vector2(s, 0.0))).Get(), false, ra)).Cut(horizontal); // H10 
  edge[2] = VectorSet2(hbb.Poisson(r, ne, (vertex[0] + vertex[1].Translated(Vector2(s, 0.0))).Get(), false, ra)).Cut(horizontal); // H01 
  edge[3] = VectorSet2(hbb.Poisson(r, ne, (vertex[1] + vertex[1].Translated(Vector2(s, 0.0))).Get(), false, ra)).Cut(horizontal); // H11 

  Box2 vertical = Box2(2 * r, s - rad_insc * 2).Translated(Vector2(0.0, s / 2.0));
  Box2 vbb = Box2(s, 2 * s).Translated(Vector2(0.0, s / 2.0));
  edge[4] = VectorSet2(vbb.Poisson(r, ne, (vertex[0] + vertex[0].Translated(Vector2(0.0, s))).Get(), false, ra)).Cut(vertical); // V00
  edge[5] = VectorSet2(vbb.Poisson(r, ne, (vertex[1] + vertex[0].Translated(Vector2(0.0, s))).Get(), false, ra)).Cut(vertical); // V10
  edge[6] = VectorSet2(vbb.Poisson(r, ne, (vertex[0] + vertex[1].Translated(Vector2(0.0, s))).Get(), false, ra)).Cut(vertical); // V01
  edge[7] = VectorSet2(vbb.Poisson(r, ne, (vertex[1] + vertex[1].Translated(Vector2(0.0, s))).Get(), false, ra)).Cut(vertical); // V11

  // Tiles generation
  QVector<Vector2> middle[16];
  cbb = cbb.Translated(Vector2(s / 2.0)).Extended(-r);

  int nm = int(n * (s * s - corner.Area() - 2 * horizontal.Area()) / (s * s));
  for (int i = 0; i < 16; i++)
  {
    VectorSet2 border;

    // Corners
    border += vertex[i & 1];
    border += vertex[(i & 2) >> 1].Translated(Vector2(s, 0.0));
    border += vertex[(i & 4) >> 2].Translated(Vector2(0.0, s));
    border += vertex[(i & 8) >> 3].Translated(Vector2(s, s));

    // Horizontal edges
    border += edge[((i & 1) >> 0) + ((i & 2) >> 0)];
    border += edge[((i & 4) >> 2) + ((i & 8) >> 2)].Translated(Vector2(0.0, s));

    // Vertical edges
    border += edge[4 + ((i & 1) >> 0) + ((i & 4) >> 1)];
    border += edge[4 + ((i & 2) >> 1) + ((i & 8) >> 2)].Translated(Vector2(s, 0.0));

    // Generate the center, and cut out the border (can't increase border density)
    VectorSet2 center = VectorSet2(cbb.Poisson(r, nm, border.Get(), false, ra)).Cut(cbb);
    center = center + border;

    // Preserve only points within the tile
    center = center.Cut(GetBox());

    // Generate the tile
    tiles[i] = DiscTile(s, r, center.Get());
  }

  // At this stage, we scramble the distribution as we started from the boundary discs
  for (int i = 0; i < 16; i++)
  {
    tiles[i].Scramble();
  }
}


/*!
\brief Sample a domain according to an input importance function.
\param s Scalar field defining the importance function.
*/
QVector<Vector2> DiscTileCorner::Sample(const ScalarField2& s) const
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
      int c = Select(i, j);

      const DiscTile& selected = tiles[c];
      for (int k = 0; k < selected.Size(); k++)
      {
        Vector2 pk = selected.Vertex(k) + tiled[0];
        if (box.Inside(pk))
        {
          // Call to the importance function
          double importance = s.Value(pk);

          // Not enough important : discarded : use strict > for comparison so that 0 importance does not generate sample
          if (!(importance > Math::Unit(k, selected.Size())))
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
QVector<Vector2> DiscTileCorner::Sample(const Box2& box) const
{
  QVector<Vector2> samples;
  Box2 tile = GetBox();

  // Get range of tile across the domain
  const QRect rect = box.TileRange(tile);

  for (int j = rect.y(); j < rect.y() + rect.height(); j++)
  {
    for (int i = rect.x(); i < rect.x() + rect.width(); i++)
    {
      Box2 offset = tile.Tile(i, j);
      int c = Select(i, j);

      const DiscTile& selected = tiles[c];
      for (int k = 0; k < selected.Size(); k++)
      {
        Vector2 pk = selected.Vertex(k) + offset[0];
        if (box.Inside(pk))
        {
          samples.append(pk);
        }
      }
    }
  }
  return samples;
}

/*!
\brief Compute the average packing density.

\sa DiscTileCorner::Packing()
*/
double DiscTileCorner::Packing() const
{
  double packing = 0.0;

  for (int i = 0; i < 16; i++)
  {
    packing += tiles[i].Packing();
  }
  return packing / 16.0;
}

