// Particles

#include "libs/particle.h"

/*!
\class ParticleSet2 particle.h
\brief Core particle set in the plane.

The particles have the same radius.
*/

/*!
\brief Create an empty set of particles.
\param r Radius.
*/
ParticleSet2::ParticleSet2(const double& r) :r(r)
{
}

/*!
\brief Create a set of particles.
\param r Radius.
\param p First particle.
*/
ParticleSet2::ParticleSet2(const Vector2& p, const double& r) : r(r)
{
  points.append(p);
}

/*!
\brief Create a set of particles.
\param r Radius.
\param s Set of particles.
*/
ParticleSet2::ParticleSet2(const QVector<Vector2>& s, const double& r) : points(s), r(r)
{
}

/*!
\brief Add a new particle to the set.
\param p Point.
*/
void ParticleSet2::Append(const Vector2& p)
{
  points.append(p);
}

/*!
\brief Draw the set of particles.
\param scene Graphics scene.
\param pen The pen.
\param brush The brush.
*/
void ParticleSet2::Draw(QGraphicsScene& scene, const QPen& pen, const QBrush& brush) const
{
  for (int i = 0; i < points.size(); i++)
  {
    At(i).Draw(scene, pen, brush);
  }
}

/*!
\brief Get the k-th particle.
\param k Index.
*/
Circle2 ParticleSet2::At(int k) const
{
  return Circle2(points.at(k), r);
}

/*!
\brief Check the intersection between the cluster and a circle.
\param circle The circle.
*/
bool ParticleSet2::Intersect(const Circle2& circle) const
{
  const Vector2 c = circle.Center();
  double h = Math::Sqr(2.0 * r);
  for (int i = 0; i < points.size(); i++)
  {
    if (SquaredNorm(points.at(i) - c) < h)
    {
      return true;
    }
  }
  return false;
}

/*!
\brief Check the intersection between the cluster and a box.
\param box The box.
*/
bool ParticleSet2::Intersect(const Box2& box) const
{
  double rr = Math::Sqr(r);
  for (int i = 0; i < points.size(); i++)
  {
    if (box.R(points.at(i)) < rr)
    {
      return true;
    }
  }
  return false;
}

/*!
\brief Computes the bounding circle of the centers of the particles.

The exact embedding circle can be computed as:
\code
ParticleSet2 particles;
Circle2 c=particles.GetCircle().Extended(particles.Radius());
\endcode
*/
Circle2 ParticleSet2::GetCircle() const
{
  return Circle2(points);
}

/*!
\brief Computes the bounding box of the centers of the particles.

\sa GetCircle()
*/
Box2 ParticleSet2::GetBox() const
{
  return Box2(points);
}
