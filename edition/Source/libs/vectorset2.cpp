// Vector sets

#include "libs/vectorset.h"
#include "libs/random.h"

/*!
\class VectorSet2 vectorset.h
\brief A simple set of vectors in the plane.

The class provides range-based for loops.

\sa VectorSet

\ingroup PlanarGroup
*/

/*!
\brief Create a set of vectors.
\param s The set.
*/
VectorSet2::VectorSet2(const QVector<Vector2>& s) :v(s)
{
}

/*!
\brief Convert a series of displacements into positions.

Only the first vector will be kept as origin.
*/
VectorSet2 VectorSet2::VectorToPoint() const
{
  if (v.size() == 0)
    return VectorSet2();

  QVector<Vector2> p(v.size());

  p[0] = v[0];
  for (int i = 1; i < v.size(); i++)
  {
    p[i] = p[i - 1] + v[i];
  }
  return VectorSet2(p);
}

/*!
\brief Reverse the order of the elements.
*/
void VectorSet2::Reverse()
{
  const int n = v.size();
  for (int i = 0; i < n / 2; i++)
  {
    Swap(v[i], v[n - 1 - i]);
  }
}

/*!
\brief Compute the barycenter.
*/
Vector2 VectorSet2::Barycenter() const
{
  Vector2 g(0.0);

  // Escape
  if (v.size() == 0)
  {
    return g;
  }

  for (int i = 0; i < v.size(); i++)
  {
    g += v.at(i);
  }
  return g / v.size();
}

/*!
\brief Compute the bounding box.

Simply calls Box::Box(const QVector&).
*/
Box2 VectorSet2::GetBox() const
{
  return Box2(v);
}

/*!
\brief Compute the convex hull of the set.
*/
Convex2 VectorSet2::GetHull() const
{
  return Convex2::Hull(v);
}

/*!
\brief Rotate all the points.
\param r Rotation matrix.
*/
void VectorSet2::Rotate(const Matrix2& r)
{
  for (int i = 0; i < v.size(); i++)
  {
    v[i] = r * v.at(i);
  }
}

/*!
\brief Rotate all the points.
\param a Rotation angle.
*/
void VectorSet2::Rotate(const double& a)
{
  Rotate(Matrix2::Rotation(a));
}

/*!
\brief Translate all the points.
\param t Translation vector.
*/
void VectorSet2::Translate(const Vector2& t)
{
  for (int i = 0; i < v.size(); i++)
  {
    v[i] += t;
  }
}

/*!
\brief Scale all the points.
\param s Scaling factor.
*/
void VectorSet2::Scale(const double& s)
{
  for (int i = 0; i < v.size(); i++)
  {
    v[i] *= s;
  }
}

/*!
\brief Scale all the points.
\param s Scaling vector.
*/
void VectorSet2::Scale(const Vector2& s)
{
  for (int i = 0; i < v.size(); i++)
  {
    v[i] *= s;
  }
}

/*!
\brief Compute the translated set.
\param t Translation vector.
*/
VectorSet2 VectorSet2::Translated(const Vector2& t) const
{
  VectorSet2 w = *this;
  w.Translate(t);

  return w;
}

/*!
\brief Compute the transformed set.
\param frame The frame.
*/
VectorSet2 VectorSet2::Transformed(const Frame2& frame) const
{
  VectorSet2 w = *this;
  for (Vector2& p : w.v)
  {
    p = frame.Transform(p);
  }

  return w;
}

/*!
\brief Compute the rotated set.
\param a Rotation matrix.
*/
VectorSet2 VectorSet2::Rotated(const double& a) const
{
  VectorSet2 w = *this;
  w.Rotate(a);

  return w;
}

/*!
\brief Compute the scaled set.
\param s Scaling vector.
*/
VectorSet2 VectorSet2::Scaled(const Vector2& s) const
{
  VectorSet2 w = *this;
  w.Scale(s);

  return w;
}

/*!
\brief Return the subset of points that lie inside the box.
\param box The box.
*/
VectorSet2 VectorSet2::Cut(const Box2& box) const
{
  VectorSet2 w;
  for (int i = 0; i < v.size(); i++)
  {
    if (box.Inside(v.at(i)))
    {
      w.Append(v.at(i));
    }
  }
  return w;
}

/*!
\brief Return the subset of points that lie inside the octogon.
\param octogon The octogon.
*/
VectorSet2 VectorSet2::Cut(const Octogon2& octogon) const
{
  VectorSet2 w;
  for (int i = 0; i < v.size(); i++)
  {
    if (octogon.Inside(v.at(i)))
    {
      w.Append(v.at(i));
    }
  }
  return w;
}

/*
\brief Add a point to the set.
\param p Point.
*/
void VectorSet2::Append(const Vector2& p)
{
  v.append(p);
}

/*
\brief Add a set of points to the set.
\param s Set of point.
*/
void VectorSet2::Append(const QVector<Vector2>& s)
{
  v.append(s);
}

/*
\brief Add a set of points to the set.
\param s Set of point.
*/
void VectorSet2::Append(const VectorSet2& s)
{
  v.append(s.v);
}

/*!
\brief Randomly move the points.

Move points using a random direction.

\param epsilon Moving distance.
*/
void VectorSet2::Vibrate(const double& epsilon)
{
  for (Vector2& c : v)
  {
    c += epsilon * Vector2::Polar(Math::TwoPi * Random::R239.Uniform());
  }
}

/*!
\brief Overloaded output-stream operator.
\param set Set.
\param s Stream.
*/
std::ostream& operator<<(std::ostream& s, const VectorSet2& set)
{
  s << "VectorSet2(" << std::endl;
  if (set.v.size() > 0)
  {
    s << std::endl;
  }
  for (int i = 0; i < set.v.size() - 1; i++)
  {
    s << ' ' << set.At(i) << ',' << std::endl;
  }
  if (set.v.size() > 0)
  {
    s << ' ' << set.At(set.v.size() - 1) << ')' << std::endl;
  }
  else
  {
    s << ')' << std::endl;
  }

  return s;
}

/*!
\brief Compute the index of the nearest point to the argument point.
\param p %Point.
*/
int VectorSet2::Nearest(const Vector2& p) const
{
  int nearest = 0;
  double r = SquaredNorm(p - v.at(0));
  for (int i = 1; i < v.size(); i++)
  {
    double ri = SquaredNorm(p - v.at(i));
    if (ri < r)
    {
      r = ri;
      nearest = i;
    }
  }
  return nearest;
}

/*!
\brief Compute the set of indexes of the k-nearest points to the argument point.
\param p %Point.
\param n Number of points.
*/
QVector<int> VectorSet2::NearestIndexes(const Vector2& p, int n) const
{
  // Create nearest set and initialize to first point
  QVector<int> nearest(n, 0);

  // Initialize distance array to the distance to first point
  // Distances will be stored in increasing order
  QVector<double> r(n, Math::Infinity);

  for (int i = 0; i < v.size(); i++)
  {
    double ri = SquaredNorm(p - v.at(i));

    int k;
    for (k = n - 1; k >= 0; k--)
    {
      if (ri > r.at(k)) break;
    }

    if (k == n - 1) continue;

    // Shift distances and nearest set to preserve sorting 
    for (int j = k + 1; j < n - 1; j++)
    {
      r[j + 1] = r[j];
      nearest[j + 1] = nearest[j];
    }
    // Update distances and set
    r[k + 1] = ri;
    nearest[k + 1] = i;

  }
  return nearest;
}

/*!
\brief Compute the set of indexes of the k-nearest points to the argument point in the set.
\param ip %Index of the point.
\param n Number of points.
\param d Maximum distance threshold.
*/
QVector<int> VectorSet2::NearestIndexes(int ip, int n, const double& d) const
{
  // Create nearest set and initialize to first point
  QVector<int> nearest;

  // Initialize distance array to the distance to first point
  // Distances will be stored in increasing order
  QVector<double> r;

  const Vector2 p = v.at(ip);

  for (int i = 0; i < v.size(); i++)
  {
    // Skip argument point
    if (i == ip) continue;

    double ri = SquaredNorm(p - v.at(i));

    // Skip if distance is above threshold
    if (ri > d) continue;

    int k;
    for (k = nearest.size() - 1; k >= 0; k--)
    {
      if (ri > r.at(k)) break;
    }

    // Here, k+1 is the index where the new index should inserted

    // Only add a new (void) element if we haven't reached the maximum
    if (nearest.size() < n)
    {
      r.append(Math::Infinity);
      nearest.append(-1);
    }

    // Shift distances and nearest set to preserve sorting 
    for (int j = k + 1; j < nearest.size() - 1; j++)
    {
      r[j + 1] = r[j];
      nearest[j + 1] = nearest[j];
    }
    // Update distances and set
    r[k + 1] = ri;
    nearest[k + 1] = i;

  }
  return nearest;
}

/*!
\brief Compute the k-nearest points to the argument point.
\param p %Point.
\param n Number of points.
*/
QVector<Vector2> VectorSet2::Nearest(const Vector2& p, int n) const
{
  // Create nearest set and initialize to first point
  QVector<Vector2> nearest(n);

  // Initialize distance array to the distance to first point
  // Distances will be stored in increasing order
  QVector<double> r(n, Math::Infinity);

  for (int i = 0; i < v.size(); i++)
  {
    double ri = SquaredNorm(p - v.at(i));

    int k;
    for (k = n - 1; k >= 0; k--)
    {
      if (ri > r.at(k)) break;
    }

    if (k == n - 1) continue;

    // Shift distances and nearest set to preserve sorting
    for (int j = k + 1; j < n - 1; j++)
    {
      r[j + 1] = r[j];
      nearest[j + 1] = nearest[j];
    }

    // Update distances and set
    r[k + 1] = ri;
    nearest[k + 1] = v.at(i);

  }
  return nearest;
}

/*!
\brief Append a vector set to another one.
\param b Added set.
*/
VectorSet2 VectorSet2::operator+(const VectorSet2& b) const
{
  QVector<Vector2> ab = v;
  ab.append(b.v);
  return VectorSet2(ab);
}

/*!
\brief Append a vector set.
\param b Added set.
*/
void VectorSet2::operator+=(const VectorSet2& b)
{
  v.append(b.v);
}
