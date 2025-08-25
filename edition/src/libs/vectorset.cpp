// Vector sets

#include "libs/vectorset.h"

#include <QtCore/QRegularExpression>

/*!
\class VectorSet vectorset.h
\brief A simple set of vectors.

This class implements several simple functions to handle point sets:
affine transforms, bounding box computation.

The class provides range-based for loops:
\code
  VectorSet v; // Set of vectors

  for (Vector& e : v) // For all elements
  {
    e *= 2.0 * e + Vector(3.0, 2.0, -1);
  }
\endcode


\ingroup KernelGroup
*/

/*!
\brief Create a set of vectors.
\param s The set.
*/
VectorSet::VectorSet(const QVector<Vector>& s) :v(s)
{
}

/*!
\brief Create a set of vectors.
\param s The set.
*/
VectorSet::VectorSet(const VectorSet2& s) 
{
  const int n = s.Size();

  v.resize(n);
  for (int i = 0; i < n; i++)
  {
    v[i] = s.At(i).ToVector(0.0);
  }
}

/*!
\brief Create a set of vectors.
\param s The set.
\param n %Size.
*/
VectorSet::VectorSet(Vector* s, unsigned int n)
{
  v.resize(n);
  for (unsigned int i = 0; i < n; i++)
  {
    v[i] = s[i];
  }
}

/*!
\brief Convert a series of displacements into positions.

Only the first vector will be kept as origin.
*/
VectorSet VectorSet::VectorToPoint() const
{
  if (v.size() == 0)
    return VectorSet();

  QVector<Vector> p(v.size());

  p[0] = v[0];
  for (int i = 1; i < v.size(); i++)
  {
    p[i] = p[i - 1] + v[i];
  }
  return VectorSet(p);
}

/*!
\brief Reverse the order of the elements.
*/
void VectorSet::Reverse()
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
Vector VectorSet::Barycenter() const
{
  Vector g(0.0);

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
Box VectorSet::GetBox() const
{
  return Box(v);
}

/*!
\brief Rotate all the points.
\param r Rotation matrix.
*/
void VectorSet::Rotate(const Matrix& r)
{
  for (int i = 0; i < v.size(); i++)
  {
    v[i] = r * v.at(i);
  }
}

/*!
\brief Translate all the points.
\param t Translation vector.
*/
void VectorSet::Translate(const Vector& t)
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
void VectorSet::Scale(const double& s)
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
void VectorSet::Scale(const Vector& s)
{
  for (int i = 0; i < v.size(); i++)
  {
    v[i] *= s;
  }
}

/*!
\brief Transform all the points.
\param f %Frame.
*/
void VectorSet::Transform(const Frame& f)
{
  for (int i = 0; i < v.size(); i++)
  {
    v[i] = f.Transform(v[i]);
  }
}

/*!
\brief Transform all the points.
\param f %Frame.
*/
void VectorSet::Transform(const FrameScaled& f)
{
  for (int i = 0; i < v.size(); i++)
  {
    v[i] = f.Transform(v[i]);
  }
}

/*!
\brief Add a point to the set.
\param p Point.
*/
void VectorSet::Append(const Vector& p)
{
  v.append(p);
}

/*!
\brief Remove redundant points that are within epsilon distance.
\param e Epsilon.
*/
void VectorSet::Clean(const double& e)
{
  for (int i = 0; i < v.size(); i++)
  {
    for (int j = i + 1; j < v.size(); j++)
    {
      if (Vector::Equal(v.at(i), v.at(j), e))
      {
        v.erase(v.begin() + j);
        j--;
      }
    }
  }
}

/*!
\brief Overloaded output-stream operator.
\param set Set.
\param s Stream.
*/
std::ostream& operator<<(std::ostream& s, const VectorSet& set)
{
  s << "VectorSet(" << std::endl;
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
int VectorSet::Nearest(const Vector& p) const
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
QVector<int> VectorSet::NearestIndexes(const Vector& p, int n) const
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
\brief Compute the k-nearest points to the argument point.
\param p %Point.
\param n Number of points.
*/
QVector<Vector> VectorSet::Nearest(const Vector& p, int n) const
{
  // Create nearest set and initialize to first point
  QVector<Vector> nearest(n);

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
\brief Compute the set of indexes of the k-nearest points to the argument point in the set.
\param ip %Index of the point.
\param n Number of points.
\param d Maximum distance threshold.
*/
QVector<int> VectorSet::NearestIndexes(int ip, int n, const double& d) const
{
  // Create nearest set and initialize to first point
  QVector<int> nearest;

  // Initialize distance array to the distance to first point
  // Distances will be stored in increasing order
  QVector<double> r;

  const Vector p = v.at(ip);

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

#include <QtCore/QFile>
#include <QtCore/QTextStream>

/*!
\brief Extract vertexes from an .obj file.
\param name File name.
*/
void VectorSet::Load(const QString& name)
{
  QFile data(name);

  if (!data.open(QFile::ReadOnly))
    return;
  QTextStream in(&data);

  // Set of regular expressions : Vertex, Normal, Triangle
  QRegularExpression rexv("v\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)");

  while (!in.atEnd())
  {
    QString line = in.readLine();
    QRegularExpressionMatch match = rexv.match(line);
    if (match.hasMatch())
    {
      Vector q = Vector(match.captured(1).toDouble(), match.captured(2).toDouble(), match.captured(3).toDouble());
      v.append(q);
    }
  }
  data.close();
}