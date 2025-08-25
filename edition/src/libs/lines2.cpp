// Set of line segments

#include <QtWidgets/QGraphicsScene>

#include "libs/segment.h"
#include "libs/circle.h"

/*!
\class SegmentSet2 segment.h
\brief Set of line segments in the plane.

Segments are stored as set of vertexes, and an index array.

This structure is effective for generating the contours of a ScalarField2.

\sa ScalarField2::LineSegments

\ingroup PlanarGroup
*/

/*!
\brief Empty.
*/
SegmentSet2::SegmentSet2()
{
}

/*!
\brief Create a set of line segments from a list of vertices and a indexes.

Indices should be a multiple of two; the number of segments is derived
from the size of the array.

\sa Mesh2

\param v Set of vertices.
\param i Indexes that represent the segments.
*/
SegmentSet2::SegmentSet2(const QVector<Vector2>& v, const QVector<int>& i) :vertices(v), indices(i)
{
}

/*!
\brief Get the set of line segments.
*/
QVector<Segment2>  SegmentSet2::GetSegments() const
{
  const int n = vertices.size() / 2;
  QVector<Segment2> s(n);
  for (int i = 0; i < n; i++)
  {
    s[i] = GetSegment(i);
  }
  return s;
}


/*!
\brief Draw a set of line segments.
\param scene Graphics scene.
\param pen The pen.
*/
void SegmentSet2::Draw(QGraphicsScene& scene, const QPen& pen) const
{
  const int n = indices.size() / 2;

  QPen joining = pen;
  joining.setCapStyle(Qt::RoundCap); // Instead of FlatCap that yields discontinuities
  joining.setJoinStyle(Qt::RoundJoin);

  for (int i = 0; i < n; i++)
  {
    Segment2 s = GetSegment(i);
    s.Draw(scene, joining);
    // GetSegment(i).DrawArrow(scene, 20.0, pen);
     //Circle2(s.Vertex(0), 15.).Draw(scene);
     //Circle2(s.Vertex(1), 15.).Draw(scene);
  }
  /*
  for (int i = 0; i < vertices.size(); i++)
   {
     Circle2(vertices.at(i), 10.).Draw(scene,QPen(QColor(200,50,0)));
   }
 */
}

#include "libs/polygon.h"

#include "libs/curvepoint.h"
/*!
\brief Get the subset of closed polygons.

Detects all loops in the structure and create closed polygons, open polylines are not returned.
*/
Polygons2 SegmentSet2::GetPolygons() const
{
  // Set of closed polygons
  Polygons2 ps;

  // Set of open curves
  QVector <PointCurve2> pc;

  QVector<int> queue = indices;

  while (!queue.empty())
  {
    // Take two elements that form the start of a polygon or a curve
    int iib = queue.last();
    int iia = queue.last();
    queue.pop_back();
    queue.pop_back();

    // Create a tentative loop of indexes
    QVector<int> loop;
    loop.append(iia);
    loop.append(iib);

    bool finished = false;
    while (!finished)
    {
      bool attached = false;
      // Find if another segment can be attached
      for (int j = 0; j < queue.size(); j += 2)
      {
        int ja = queue.at(j);
        int jb = queue.at(j + 1);

        // Attach from this side
        if (ja == loop.at(0))
        {
          attached = true;
          if (jb == loop.at(loop.size() - 1))
          {
            // Both ends match: close polygon
            ps.Add(Polygon2(vertices, loop));
            finished = true;
          }
          else
          {
            // Only attach
            loop.prepend(jb);
          }
          // Remove candidates : they have been attached
          queue.remove(j, 2);
          break;
        }
        // Or at the end
        else if (jb == loop.at(0))
        {
          attached = true;
          if (ja == loop.at(loop.size() - 1))
          {
            // Both ends match: close polygon
            ps.Add(Polygon2(vertices, loop));
            finished = true;
          }
          else
          {
            // Only attach
            loop.prepend(ja);
          }
          // Remove candidates : they have been attached
          queue.remove(j, 2);
          break;
        }
      }
      // We traversed all the set and could not attach any segment: this means the loop cannot be closed
      if (attached == false)
      {
        pc.append(PointCurve2(vertices, loop));
        finished = true;
      }
    }
  }
  return ps;
}
