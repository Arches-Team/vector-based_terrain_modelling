// Polygons

#include "libs/polygon.h"

// Used for drawing
#include <QtGui/QPolygonF>

/*!
\class Polygons2 polygon.h
\brief Set of disjoint polygons.
\ingroup PlanarGroup
*/

/*!
\brief Empty.
*/
Polygons2::Polygons2()
{
}

/*!
\brief Create a set of polygons.
\param p Set of polygons.
*/
Polygons2::Polygons2(const QVector<Polygon2>& p) :poly(p)
{
}

/*!
\brief Create a set of polygons by projecting polygons in space.
\param p Set of polygons.
*/
Polygons2::Polygons2(const QVector<Polygonal>& p)
{
  poly.reserve(p.size());
  for (int i = 0; i < p.size(); i++)
  {
    poly.append(Polygon2(p.at(i)));
  }
}

/*!
\brief Add a polygon to the set.
\param p %Polygon.
*/
void Polygons2::Add(const Polygon2& p)
{
  poly.append(p);
}

/*!
\brief Draw a polygon.
\param scene Graphics scene.
\param pen The pen.
\param brush The brush, should the triangle be filled.
*/
void Polygons2::Draw(QGraphicsScene& scene, const QPen& pen, const QBrush& brush) const
{
  for (const Polygon2& p:poly)
  {
    p.Draw(scene, pen, brush);
  }
}

/*!
\brief Test if a point is inside the union of a set of polygons.

Polygons should be non intersecting. Polygons inside others will be considered as holes.

\param p Point.
*/
bool Polygons2::Inside(const Vector2& p) const
{
  int n = 0;

  for (int j = 0; j < poly.size(); j++)
  {
    if (poly.at(j).Inside(p))
    {
      n++;
    }
  }

  if (n % 2 == 0)
    return false;
  else
    return true;
}

/*!
\brief Compute the squared distance between a point and a set of polygons.
\param p Point.
*/
double Polygons2::RC(const Vector2& p) const
{
  double r = poly.at(0).RC(p);
  for (int i = 1; i < poly.size(); i++)
  {
    r = Math::Min(r, poly.at(i).RC(p));
  }
  return r;
}