// Pyramid

#include "libs/pyramid.h"
#include "libs/triangle.h"

/*!
\class Pyramid pyramid.h

\brief An axis aligned pyramid.

\ingroup ExtendedKernelGroup
*/

// Array of vertex indexes
const int Pyramid::face[6][3] = {
  { 0, 2, 1 },{ 1, 2, 3 },{ 0, 1, 4 },{ 1, 3, 4 },{ 3, 2, 4 },{ 2, 0, 4 }
};

// Array of vertexes 
const Vector Pyramid::vertex[5] = {
  Vector(-1.0,-1.0,0.0),
  Vector(1.0,-1.0,0.0),
  Vector(-1.0, 1.0,0.0),
  Vector(1.0, 1.0,0.0),
  Vector(0.0,0.0,1.0)
};

/*!
\brief Create a pyramid.
\param s Half base-side size.
\param z Height.
*/
Pyramid::Pyramid(const double& s, const double& z) :Pyramid(Vector::Null, s, z)
{
}

/*!
\brief Create a pyramid.
\param c Center.
\param s Half base-side size.
\param z Height.
*/
Pyramid::Pyramid(const Vector& c, const double& s, const double& z) :center(c), a(s), height(z)
{
  double edgelength = sqrt(a * a + height * height);
  double longedgelength = sqrt(2.0 * a * a + height * height);
  nt = Vector(height, 0.0, a) / edgelength;
  neaa = Vector(-a, -a, height) / longedgelength;
}

/*!
\brief Create a pyramid.
\param c Center.
\param s Half base-side size, which will be height as well.
*/
Pyramid::Pyramid(const Vector& c, const double& s) :Pyramid(c, s, s)
{
}

/*!
\brief Compute the volume of the pyramid.
*/
double Pyramid::Volume() const
{
  return a * a * height / 6.0;
}

/*!
\brief Compute the area of the pyramid.

It is eight times the area of the small triangle, plus surface of the base square
*/
double Pyramid::Area() const
{
  return 8.0 * a * sqrt(a * a + height * height);
}

/*!
\brief Check if a point is inside the pyramid.
\param p Point.
*/
bool Pyramid::Inside(const Vector& p) const
{
  Vector q = p - center;

  // Base plane
  if (q[2] < 0.0) return false;

  // Take advantage of symmetry
  q = Vector(fabs(q[0]), fabs(q[1]), q[2]);

  Vector aaq = q - Vector(a, a, 0.0);

  // First oblique plane
  if (aaq * Vector(height, 0.0, a) > 0.0)
  {
    return false;
  }
  // Second oblique plane
  if (aaq * Vector(0.0, height, a) > 0.0)
  {
    return false;
  }

  return true;
}

/*!
\brief Compute the normal vector between a point and the pyramid.
\param p Point.
*/
Vector Pyramid::Normal(const Vector& p) const
{
  Vector q = p - center;

  // Take advantage of axial symmetry
  Vector symmetry = Vector(q[0] < 0.0 ? -1.0 : 1.0, q[1] < 0.0 ? -1.0 : 1.0, 1.0);
  q = Vector(fabs(q[0]), fabs(q[1]), q[2]);

  // Position against diagonal point
  Vector aaq = q - Vector(a, a, 0.0);

  Vector n = Vector::Null;

  // Below: take distance to square
  if (q[2] < 0.0)
  {
    if (aaq[0] < 0.0)
    {
      // Square
      if (aaq[1] < 0.0)
      {
        n = Vector(0.0, 0.0, q[2]);
      }
      // Edge
      else
      {
        n = Vector(0.0, aaq[1], q[2]);
      }
    }
    else
    {
      // Edge
      if (aaq[1] < 0.0)
      {
        n = Vector(aaq[0], 0.0, q[2]);
      }
      // Vertex
      else
      {
        n = Vector(aaq[0], aaq[1], q[2]);
      }
    }
    n = n.Scaled(symmetry);
    return n;
  }

  // Above: take distance to the pyramid
  bool swap = false;
  // Diagonal symmetry

  if (q[1] > q[0])
  {
    Math::Swap(q[0], q[1]);
    Math::Swap(aaq[0], aaq[1]);
    swap = true;
  }

  // Normal to face
  // Some of the coefficients could be saved to save cross product and divisions
  Vector ntaa(nt[2], 0.0, -nt[0]);
  Vector ntaaz = nt / neaa;

  // Inside
  double dt = aaq * nt;
  if (dt < 0.0)
  {
    // Inside, normal kept to null
    return n;
  }

  if (aaq * ntaa > 0.0)
  {
    if (aaq[1] < 0.0)
    {
      // Edge Ax
      n = aaq;
      n[1] = 0.0;
    }
    else
    {
      // Vertex A
      n = aaq;
    }
  }
  else
  {
    Vector qz = q - Vector(0.0, 0.0, height);

    if (qz * neaa > 0.0)
    {
      // Vertex T
      n = qz;
    }
    else
    {
      if ((aaq * ntaaz < 0.0))
      {
        // Edge
        n = aaq - (aaq * neaa) * neaa;
      }
      else
      {
        // Face
        n = dt * nt;
      }
    }
  }
  if ((aaq[1] > 0.0) && (aaq * neaa < 0.0))
  {
    n = aaq;
  }

  // Apply inverse diagonal symmetry if needed
  if (swap == true)
  {
    Math::Swap(n[0], n[1]);
  }

  n = n.Scaled(symmetry);

  return n;
}

/*!
\brief Compute the signed distance between a point and the pyramid.
\param p Point.
*/
double Pyramid::Signed(const Vector& p) const
{
  Vector q = p - center;

  // Take advantage of axial symmetry
  Vector symmetry = Vector(q[0] < 0.0 ? -1.0 : 1.0, q[1] < 0.0 ? -1.0 : 1.0, 1.0);
  q = Vector(fabs(q[0]), fabs(q[1]), q[2]);

  // Position against diagonal point
  Vector aaq = q - Vector(a, a, 0.0);

  Vector n = Vector::Null;

  // Below: take distance to square
  if (q[2] < 0.0)
  {
    if (aaq[0] < 0.0)
    {
      // Square
      if (aaq[1] < 0.0)
      {
        n = Vector(0.0, 0.0, q[2]);
      }
      // Edge
      else
      {
        n = Vector(0.0, aaq[1], q[2]);
      }
    }
    else
    {
      // Edge
      if (aaq[1] < 0.0)
      {
        n = Vector(aaq[0], 0.0, q[2]);
      }
      // Vertex
      else
      {
        n = Vector(aaq[0], aaq[1], q[2]);
      }
    }
    return Norm(n);
  }

  // Above: take distance to the pyramid
  bool swap = false;
  // Diagonal symmetry

  if (q[1] > q[0])
  {
    Math::Swap(q[0], q[1]);
    Math::Swap(aaq[0], aaq[1]);
    swap = true;
  }

  // Normal to face
  // Some of the coefficients could be saved to save cross product and divisions
  const Vector ntaa(nt[2], 0.0, -nt[0]);
  const Vector ntaaz = nt / neaa;

  // Inside
  const double dt = aaq * nt;
  if (dt < 0.0)
  {
    // Inside: compute the maximum negative distance
    return Math::Max(dt, -q[2]);
  }

  if (aaq * ntaa > 0.0)
  {
    if (aaq[1] < 0.0)
    {
      // Edge Ax
      n = aaq;
      n[1] = 0.0;
    }
    else
    {
      // Vertex A
      n = aaq;
    }
  }
  else
  {
    Vector qz = q - Vector(0.0, 0.0, height);

    if (qz * neaa > 0.0)
    {
      // Vertex T
      n = qz;
    }
    else
    {
      if ((aaq * ntaaz < 0.0))
      {
        // Edge
        n = aaq - (aaq * neaa) * neaa;
      }
      else
      {
        // Face
        n = dt * nt;
      }
    }
  }
  if ((aaq[1] > 0.0) && (aaq * neaa < 0.0))
  {
    n = aaq;
  }

  // Inverse diagonal symmetry not needed for signed distance
  return Norm(n);
}

/*!
\brief Compute the squared distance between a point and the pyramid.
\param p The point.
*/
double Pyramid::R(const Vector& p) const
{
  //  const Vector n = Normal(p);
  //  return n * n;
  Vector q = p - center;

  // Take advantage of axial symmetry
  Vector symmetry = Vector(q[0] < 0.0 ? -1.0 : 1.0, q[1] < 0.0 ? -1.0 : 1.0, 1.0);
  q = Vector(fabs(q[0]), fabs(q[1]), q[2]);

  // Position against diagonal point
  Vector aaq = q - Vector(a, a, 0.0);

  Vector n = Vector::Null;

  // Below: take distance to square
  if (q[2] < 0.0)
  {
    if (aaq[0] < 0.0)
    {
      // Square
      if (aaq[1] < 0.0)
      {
        return q[2] * q[2];
      }
      // Edge
      else
      {
        return aaq[1] * aaq[1] + q[2] * q[2];
      }
    }
    else
    {
      // Edge
      if (aaq[1] < 0.0)
      {
        return aaq[0] * aaq[0] + q[2] * q[2];
      }
      // Vertex
      else
      {
        n = Vector(aaq[0], aaq[1], q[2]);
        return n * n;
      }
    }
  }

  // Above: take distance to the pyramid
  bool swap = false;
  // Diagonal symmetry

  if (q[1] > q[0])
  {
    Math::Swap(q[0], q[1]);
    Math::Swap(aaq[0], aaq[1]);
    swap = true;
  }

  // Normal to face
  // Some of the coefficients could be saved to save cross product and divisions
  Vector ntaa(nt[2], 0.0, -nt[0]);
  Vector ntaaz = nt / neaa;

  // Inside
  double dt = aaq * nt;
  if (dt < 0.0)
  {
    // Inside, normal kept to null
    return 0.0;
  }

  if (aaq * ntaa > 0.0)
  {
    if (aaq[1] < 0.0)
    {
      // Edge Ax
      n = aaq;
      n[1] = 0.0;
    }
    else
    {
      // Vertex A
      n = aaq;
    }
  }
  else
  {
    Vector qz = q - Vector(0.0, 0.0, height);

    if (qz * neaa > 0.0)
    {
      // Vertex T
      n = qz;
    }
    else
    {
      if ((aaq * ntaaz < 0.0))
      {
        // Edge
        n = aaq - (aaq * neaa) * neaa;
      }
      else
      {
        // Face
        n = dt * nt;
      }
    }
  }
  if ((aaq[1] > 0.0) && (aaq * neaa < 0.0))
  {
    n = aaq;
  }

  // Inverse diagonal symmetry not needed for signed distance
  return n * n;
}

/*!
\brief %Normalized normal to the i-th face.
\param i Index.
*/
Vector Pyramid::Normal(int i) const
{
  return Triangle(Vertex(i, 0), Vertex(i, 1), Vertex(i, 2)).Normal();
}

/*!
\brief Overloaded.
\param s Stream.
\param pyramid The pyramid.
*/
std::ostream& operator<<(std::ostream& s, const Pyramid& pyramid)
{
  s << "Pyramid(" << pyramid.center << ',' << pyramid.a << ',' << pyramid.height << ")";
  return s;
}

