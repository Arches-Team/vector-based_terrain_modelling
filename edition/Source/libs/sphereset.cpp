// Sphere set

#include "libs/sphereset.h"

/*!
\class SphereSet sphereset.h
\brief A simple set of spheres.

This class implements functions to compute the
distance between a point and a set of spheres, or the detect whether
two set of spheres intersect.

\ingroup KernelGroup

*/

/*!
\brief Creates a sphere set given one sphere.
\param s %Sphere.
*/
SphereSet::SphereSet(const Sphere& s) 
{
  spheres.append(s);
}

/*!
\brief Creates a set of spheres.
\param s Set.
*/
SphereSet::SphereSet(const QVector<Sphere>& s) :spheres(s)
{
}

/*!
\brief Creates a set of spheres and keep only the subset of spheres that intersect the second argument sphere.
\param ss Set.
\param s %Sphere.
*/
SphereSet::SphereSet(const SphereSet& ss, const Sphere& s)
{
  for (auto i = ss.spheres.begin(); i != ss.spheres.end(); i++)
  {
    if (i->Intersect(s))
    {
      spheres << *i;
    }
  }
}

/*!
\brief Merge with another set of spheres.
\param s SphereSet.
*/
void SphereSet::Append(const SphereSet& s)
{
  spheres << s.spheres;
}

/*!
\brief returns true if the set is empty
*/
bool SphereSet::IsEmpty() const
{
  return spheres.isEmpty();
}

/*!
\brief Compute the squared distance between a point and the set of spheres.
\param p The point.
*/
double SphereSet::R(const Vector& p) const
{
  if (spheres.size() == 0) return 0.0;

  double a = spheres.at(0).R(p);
  for (int i = 1; i < spheres.size(); i++)
  {
    double r = spheres.at(i).R(p);
    a = Math::Min(a, r);
  }
  return a;
}

/*!
\brief Compute the signed distance between a point and the set of spheres.
\param p The point.
*/
double SphereSet::Signed(const Vector& p) const
{
  if (spheres.size() == 0) return 0.0;

  double a = spheres.at(0).R(p);
  for (int i = 1; i < spheres.size(); i++)
  {
    double r = spheres.at(i).Signed(p);
    a = Math::Min(a, r);
  }
  return a;
}

/*!
\brief Compute the bounding box of the set of spheres.
*/
Box SphereSet::GetBox() const
{
  if (spheres.size() == 0) return Box::Null;

  Box box = spheres.at(0).GetBox();
  for (int i = 1; i < spheres.size(); i++)
  {
    box = Box(box, spheres.at(i).GetBox());
  }
  return box;
}

/*!
\brief Check if a point is inside or outside the set of sphere.
\param p The point.
*/
bool SphereSet::Inside(const Vector& p) const
{
  for (int i = 0; i < spheres.size(); i++)
  {
    if (spheres.at(i).Inside(p))
    {
      return true;
    }
  }
  return false;
}

/*!
\brief Computes the normal vector between a point and the sphere set.

Simply project point onto the spheres, and return the corresponding
Euclidean distance vector.
\param p Point.
*/
Vector SphereSet::Normal(const Vector& p) const
{
  // Escape
  if (spheres.size() == 0) return Vector::Null;

  Vector n = p - spheres.at(0).Center();
  double r = n * n;
  int k = 0;

  // Inside first sphere
  if (r < Math::Sqr(spheres.at(0).Radius()))
  {
    return Vector::Null;
  }

  for (int i = 1; i < spheres.size(); i++)
  {
    Vector ni = p - spheres.at(i).Center();
    double ri = ni * ni;

    // Inside sphere i
    if (ri < Math::Sqr(spheres.at(i).Radius()))
    {
      return Vector::Null;
    }
    else
    {
      if (ri < r)
      {
        r = ri;
        n = ni;
        k = i;
      }
    }
  }

  // We avoided the square root until the last step
  r = sqrt(r);
  n *= (1.0 - spheres.at(k).Radius() / r);

  return n;
}

/*!
\brief Rotates a sphere set.

\param r Rotation matrix.
*/
void SphereSet::Rotate(const Matrix& r)
{
  for (int i = 0; i < spheres.size(); i++)
  {
    spheres[i].Rotate(r);
  }
}

/*!
\brief Translates a sphere set.

\param t Translation vector.
*/
void SphereSet::Translate(const Vector& t)
{
  for (int i = 0; i < spheres.size(); i++)
  {
    spheres[i].Translate(t);
  }
}

/*!
\brief Uniformly scales a sphere-set.

\param s Scaling factor.
*/
void SphereSet::Scale(const double& s)
{
  for (int i = 0; i < spheres.size(); i++)
  {
    spheres[i].Scale(s);
  }
}

/*!
\brief Transforms a sphere set given a frame transformation.
\param t Transformation.
*/
void SphereSet::Transform(const Frame& t)
{
  for (int i = 0; i < spheres.size(); i++)
  {
    spheres[i] = spheres.at(i).Transformed(t);
  }
}

/*!
\brief Returns the number of Spheres.
*/
int SphereSet::Size(void) const
{
  return spheres.size();
}

/*!
\brief Returns the number of Spheres.
*/
Sphere SphereSet::GetSphere(int i) const
{
  return spheres[i];
}

/*!
\brief Removes spheres that are close to or embedded in the sphere set.
\param t Threshold value.
*/
void SphereSet::RemoveDuplicates(const double& t)
{
  QVector<Sphere> set;

  const int n = spheres.size();
  for (int i = 0; i < n; i++)
  {
    bool keep = true;
    for (int j = i + 1; j < n; j++)
    {
      if (Norm(spheres.at(i).Center() - spheres.at(j).Center()) - spheres[j].Radius() + spheres[i].Radius() < t)
      {
        keep = false;
        break;
      }
    }
    if (keep)
      set.append(spheres.at(i));
  }
  spheres = set;
}

/*!
\brief Check if two sphere-sets intersect.

\param set Set of spheres.
*/
bool SphereSet::Intersect(const SphereSet& set) const
{
  for (int i = 0; i < spheres.size(); i++)
  {
    for (int j = 0; j < set.spheres.size(); j++)
    {
      if (spheres.at(i).Intersect(set.spheres.at(j)))
      {
        return true;
      }
    }
  }

  return false;
}

/*!
\brief Compute the signed distance between two sphere-sets.

If two spheres intersect, the result will be negative. 
If one of the set is empty, returns 0.0.

\param set Set of spheres.
*/
double SphereSet::R(const SphereSet& set) const
{
  if (set.spheres.size() == 0 || spheres.size() == 0)
    return 0.0;

  // Minimum distance
  double e = set.spheres.at(0).R(spheres[0]);

  for (int i = 0; i < spheres.size(); i++)
  {
    for (int j = 0; j < set.spheres.size(); j++)
    {
      double d = set.spheres.at(j).R(spheres[i]);
      if (d < e)
      {
        e = d;
      }
    }
  }

  return e;
}