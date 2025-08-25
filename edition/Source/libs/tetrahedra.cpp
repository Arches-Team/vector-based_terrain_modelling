// Tetrahedra

#include "libs/tetrahedra.h"
#include "libs/plane.h"
#include "libs/cpu.h"

/*!
\class Tetrahedra tetrahedra.h
\brief A simple tetrahedron.

The class stores the vertices, the normals of the triangles.
It does not store the unit vectors of the edges.

\image html tetrahedron.png

This class also implements some useful functions, such as the
intersection with a ray. The member function GetSphere() used to
compute the bounding sphere of a tetrahedron simply calls the
right Sphere constructor.

\ingroup ExtendedKernelGroup
*/


/*!
\brief Creates a tetrahedron.

Note that vertices should be ordered so that triangle normals
should point outward.

\param a,b,c,d Vertices.
*/
Tetrahedra::Tetrahedra(const Vector& a, const Vector& b, const Vector& c, const Vector& d) :p{ a,b,c,d }
{
  // Edges: could be pre-processed in the constructor
  Vector e01 = Normalized(p[1] - p[0]);
  Vector e02 = Normalized(p[2] - p[0]);
  Vector e03 = Normalized(p[3] - p[0]);
  Vector e12 = Normalized(p[2] - p[1]);
  Vector e13 = Normalized(p[3] - p[1]);
  Vector e23 = Normalized(p[3] - p[2]);

  // Normal 012
  n[0] = Normalized(e01 / e02);

  // Normal 031
  n[1] = Normalized(e03 / e01);

  // Normal 023
  n[2] = Normalized(e02 / e03);

  // Normal 132
  n[3] = Normalized(e13 / e12);
}

// Could have been : 0214 3241 3264 3145 3465 3675

/*!
\brief Returns the k-th tetrahedron of the box.

The box is divided into six tetrahedra.

\param box The box.
\param k Index of the tetrahedron.
*/
Tetrahedra::Tetrahedra(const Box& box, int k) :Tetrahedra(box.Vertex(split[k][0]), box.Vertex(split[k][1]), box.Vertex(split[k][2]), box.Vertex(split[k][3]))
{
}

/*!
\brief Compute the sphere bounding the argument tetrahedron.

This is convenience function, simply calling the constructor of Sphere with four vertices.

\sa Sphere(const Vector&, const Vector&, const Vector&, const Vector&)
*/
Sphere Tetrahedra::GetSphere() const
{
  return Sphere(p[0], p[1], p[2], p[3]);
}

/*!
\brief Returns the i-th plane of the tetrahedron.
\param i Index.
*/
Plane Tetrahedra::GetPlane(int i) const
{
  return Plane(Normal(i), Vertex(i));
}

/*!
\brief Compute the intersection between tetrahedra and a ray.

Sorted intersection depths are returned if intersection occurs.
\param ray The ray.
\param ta, tb Intersection depths.
*/
bool Tetrahedra::Intersect(const Ray& ray, double& ta, double& tb) const
{
  // Check all triangles
  for (int i = 0; i < 4; i++)
  {
    int k = 0;
    // Check all planes
    double x = n[i] * ray.Direction();
    if ((x < epsilon) && (x > -epsilon))
      continue;

    // Intersection
    double t = (n[i] * (p[i] - ray.Origin())) / x;

    // Point
    Vector q = ray(t);

    int j;
    // Check against other three planes
    for (j = 0; j < 4; j++)
    {
      if (j != i)
      {
        if ((q - p[i]) * n[i] > 0.0) break;
      }
    }
    // Valid intersection detected
    if (j == 4)
    {
      if (k == 0)
      {
        ta = t;
      }
      else
      {
        tb = t;
      }
      k++;
    }

    // Escape if all intersections have been detected
    if (k == 2)
    {
      Math::Sort(ta, tb);
      return true;
    }
  }

  return false;
}

/*!
\brief Compute the intersection between tetrahedra and a ray.

Sorted intersection depths are returned if intersection occurs.
\param ray The ray.
\param ta, tb Returned intersection depths.
\param na, nb Returned normals at intersection points.
*/
bool Tetrahedra::Intersect(const Ray& ray, double& ta, double& tb, Vector& na, Vector& nb) const
{
  // Check all triangles
  for (int i = 0; i < 4; i++)
  {
    int k = 0;
    // Check all planes
    double x = n[i] * ray.Direction();
    if ((x < epsilon) && (x > -epsilon))
      continue;

    // Intersection
    double t = (n[i] * (p[i] - ray.Origin())) / x;

    // Point
    Vector q = ray(t);

    int j;
    // Check against other three planes
    for (j = 0; j < 4; j++)
    {
      if (j != i)
      {
        if ((q - p[i]) * n[i] > 0.0) break;
      }
    }
    // Valid intersection detected
    if (j == 4)
    {
      if (k == 0)
      {
        ta = t;
        na = n[i];
      }
      else
      {
        tb = t;
        nb = n[i];
      }
      k++;
    }

    // Escape if all intersections have been detected
    if (k == 2)
    {
      if (tb < ta)
      {
        t = tb;
        tb = ta;
        ta = t;
        Vector n = na;
        na = nb;
        nb = n;
      }
      return true;
    }
  }

  return false;
}

/*!
\brief Check if a point lies inside the tetrahedron.
\param p Point.
*/
bool Tetrahedra::Inside(const Vector& p) const
{
  // Check against planes
  for (int i = 0; i < 4; i++)
  {
    if ((p - Tetrahedra::p[i]) * n[i] > (epsilon * epsilon))
    {
      return false;
    }
  }
  return true;
}

/*!
\brief Compute the volume of a tetrahedron.

The volume may be defined as one-sixth the volume
of a parallelepiped whose three main edges emanate from a single vertex.
*/
double Tetrahedra::Volume() const
{
  return fabs(Matrix(p[1] - p[0], p[2] - p[0], p[3] - p[0]).Determinant()) / 6.0;
}

/*!
\brief Compute the surface area of a tetrahedron.

The surface is defined as the sum of the
surface areas of its four triangle faces. Other
techniques define the surface area as:
<BR>A=(u+v+w).(v-w+u).(w-u+v).(u-v+w)
<BR>Here u, v and w denote the three edges length.
*/
double Tetrahedra::Area() const
{
  double u = Norm(p[1] - p[0]);
  double v = Norm(p[2] - p[0]);
  double w = Norm(p[3] - p[0]);
  return (u + v + w) * (v - w + u) * (w - u + v) * (u - v + w);
}

/*!
\brief Compute the bounding box of a tetrahedron.
*/
Box Tetrahedra::GetBox() const
{
  Vector a = Vector::Min(p[0], Vector::Min(p[1], Vector::Min(p[2], p[3])));
  Vector b = Vector::Max(p[0], Vector::Max(p[1], Vector::Max(p[2], p[3])));

  return Box(a, b);
}

/*!
\brief Compute if the tetrahedra is inside a half plane.
\param plane The plane.
*/
bool Tetrahedra::Inside(const Plane& plane) const
{
  return plane.Inside(p[0]) && plane.Inside(p[1]) && plane.Inside(p[2]) && plane.Inside(p[3]);
}

/*!
\brief Overloaded.
\param s Stream.
\param t The tetrahedron.
*/
std::ostream& operator<<(std::ostream& s, const Tetrahedra& t)
{
  s << "Tetrahedra(" << t.p[0] << ',' << t.p[1] << ',' << t.p[2] << ',' << t.p[3] << ')' << std::endl;
  return s;
}

/*!
\brief Generate a random point inside the tetrahedron.

The method proposed here is a generalization of one of the techniques
used by Turk to generate a random point inside a triangle. Generate a
random point in a parallelogram and reflectit around the center of the
parallelogram.

After C. Rocchini and P. Cignoni, Generating random points in a tetrahedron,
<I>Journal of graphics tools</I>, <B>5</B>(4):9-12, 2000.

\param r %Random number generator.
*/
Vector Tetrahedra::RandomInside(Random& r) const
{
  double s = r.Uniform();
  double t = r.Uniform();
  double u = r.Uniform();
  // Cut'n fold the cube into a prism
  if (s + t > 1.0)
  {
    s = 1.0 - s;
    t = 1.0 - t;
  }
  // Cut'n fold the prism into a tetrahedron
  if (t + u > 1.0)
  {
    double tu = u;
    u = 1.0 - s - t;
    t = 1.0 - tu;
  }
  else if (s + t + u > 1.0)
  {
    double tu = u;
    u = s + t + u - 1.0;
    s = 1 - t - tu;
  }

  return (1.0 - s - t - u) * p[0] + s * p[1] + t * p[2] + u * p[3];
}

/*!
\brief Computes the vector distance between the tetrahedron and a point.
\param x Point.
*/
Vector Tetrahedra::Normal(const Vector& x) const
{
  // Edges: could be pre-processed in the constructor
  Vector e01 = Normalized(p[1] - p[0]);
  Vector e02 = Normalized(p[2] - p[0]);
  Vector e03 = Normalized(p[3] - p[0]);
  Vector e12 = Normalized(p[2] - p[1]);
  Vector e13 = Normalized(p[3] - p[1]);
  Vector e23 = Normalized(p[3] - p[2]);

  Vector ne012 = e01 / n[0];
  Vector ne120 = e12 / n[0];
  Vector ne201 = -(e02 / n[0]);

  // Face 031
  Vector ne013 = (n[1] / e01);
  Vector ne130 = (n[1] / e13);
  Vector ne301 = -(n[1] / e03);

  // Face 023
  Vector ne023 = -(n[2] / e02);
  Vector ne230 = -(n[2] / e23);
  Vector ne302 = (n[2] / e03);

  // Face 132
  Vector ne123 = (n[3] / e12);
  Vector ne231 = (n[3] / e23);
  Vector ne312 = -(n[3] / e13);

  unsigned long u = 0;

  //   0-3: above Face 
  //  4-15: beyond edge (for example, edge 01 for face 012) : 013 023 123 012 021 120 013 031 130 023 032 230 123 132 231
  // 16-27: edge are of face 

  Vector pa0 = x - p[0];
  Vector pa1 = x - p[1];
  Vector pa2 = x - p[2];
  Vector pa3 = x - p[3];

  // 4 Faces
  double pa0n0 = pa0 * n[0];
  double pa0n1 = pa0 * n[1];
  double pa0n2 = pa0 * n[2];
  double pa1n3 = pa1 * n[3];

  if (pa0n0 > 0.0) Byte::Set(u, 0);
  if (pa0n1 > 0.0) Byte::Set(u, 1);
  if (pa0n2 > 0.0) Byte::Set(u, 2);
  if (pa1n3 > 0.0) Byte::Set(u, 3);

  // 3 Edges for 4 faces
  double pa0e01 = pa0 * e01;
  double pa0e02 = pa0 * e02;
  double pa0e03 = pa0 * e03;

  double pa1e01 = pa1 * e01;
  double pa1e12 = pa1 * e12;
  double pa1e13 = pa1 * e13;

  double pa2e02 = pa2 * e02;
  double pa2e12 = pa2 * e12;
  double pa2e23 = pa2 * e23;

  double pa3e03 = pa3 * e03;
  double pa3e13 = pa3 * e13;
  double pa3e23 = pa3 * e23;

  // Inside
  if (u == 0) return Vector::Null;

  // Face [012]
  if (pa0 * ne012 > 0.0)
  {
    Byte::Set(u, 4);
    if (pa0e01 > 0.0 && pa1e01 <= 0.0) Byte::Set(u, 16);
  }
  if (pa2 * ne201 > 0.0)
  {
    Byte::Set(u, 5);
    if (pa0e02 > 0.0 && pa2e02 <= 0.0) Byte::Set(u, 17);
  }
  if (pa1 * ne120 > 0.0)
  {
    Byte::Set(u, 6);
    if (pa1e12 > 0.0 && pa2e12 <= 0.0) Byte::Set(u, 18);
  }

  // Face [013]
  if (pa0 * ne013 > 0.0)
  {
    Byte::Set(u, 7);
    if (pa0e01 > 0.0 && pa1e01 <= 0.0) Byte::Set(u, 19);
  }
  if (pa3 * ne301 > 0.0)
  {
    Byte::Set(u, 8);
    if (pa0e03 > 0.0 && pa3e03 <= 0.0) Byte::Set(u, 20);
  }
  if (pa1 * ne130 > 0.0)
  {
    Byte::Set(u, 9);
    if (pa1e13 > 0.0 && pa3e13 <= 0.0) Byte::Set(u, 21);
  }

  // face 023
  if (pa0 * ne023 > 0.0)
  {
    Byte::Set(u, 10);
    if (pa0e02 > 0.0 && pa2e02 <= 0.0) Byte::Set(u, 22);
  }
  if (pa3 * ne302 > 0.0)
  {
    Byte::Set(u, 11);
    if (pa0e03 > 0.0 && pa3e03 <= 0.0) Byte::Set(u, 23);
  }
  if (pa2 * ne230 > 0.0)
  {
    Byte::Set(u, 12);
    if (pa2e23 > 0.0 && pa3e23 <= 0.0) Byte::Set(u, 24);
  }

  // face 123
  if (pa1 * ne123 > 0.0)
  {
    Byte::Set(u, 13);
    if (pa1e12 > 0.0 && pa2e12 <= 0.0) Byte::Set(u, 25);
  }
  if (pa3 * ne312 > 0.0)
  {
    Byte::Set(u, 14);
    if (pa1e13 > 0.0 && pa3e13 <= 0.0) Byte::Set(u, 26);
  }
  if (pa2 * ne231 > 0.0)
  {
    Byte::Set(u, 15);
    if (pa2e23 > 0.0 && pa3e23 <= 0.0) Byte::Set(u, 27);
  }

  // Faces
  if ((u & 1) > 0 && ((~u) & 112) == 112) return pa0n0 * n[0]; // face 012 : b0, !b4, !b5, !b6
  if ((u & 2) > 0 && ((~u) & 896) == 896) return pa0n1 * n[1]; // face 013 : b1, !b7, !b8, !b9
  if ((u & 4) > 0 && ((~u) & 7168) == 7168) return pa0n2 * n[2]; // face 023 : b2, !b10, !b11, !b12
  if ((u & 8) > 0 && ((~u) & 57344) == 57344) return pa1n3 * n[3]; // face 123 : b3, !b13, !b14, !b15

  // Edges
  if ((u & 589824) == 589824) return pa1 - pa1e01 * e01; //   au dessus de e01 : 012 & 013 : b16, b19
  if ((u & 4325376) == 4325376) return pa2 - pa2e02 * e02; //   au dessus de e02 : 021 & 023 : b17, b22
  if ((u & 9437184) == 9437184) return pa3 - pa3e03 * e03; //   au dessus de e03 : 031 & 032 : b20, b23
  if ((u & 33816576) == 33816576) return pa2 - pa2e12 * e12; //   au dessus de e12 : 120 & 123 : b18, b25
  if ((u & 69206016) == 69206016) return pa3 - pa3e13 * e13; //   au dessus de e13 : 130 & 132 : b21, b26
  if ((u & 150994944) == 150994944) return pa3 - pa3e23 * e23; //   au dessus de e23 : 230 & 231 : b24, b27

  // Vertex
  Vector petit;
  double npetit = pa0 * pa0;
  petit = pa0;
  if (pa1 * pa1 < npetit) { petit = pa1; npetit = pa1 * pa1; }
  if (pa2 * pa2 < npetit) { petit = pa2; npetit = pa2 * pa2; }
  if (pa2 * pa3 < npetit) { petit = pa3; }
  return petit;
}

/*!
\brief Computes the vector distance between the tetrahedron and a point.
\param x Point.
*/
Vector Tetrahedra::Normal2(const Vector& x) const
{
  Vector pa0 = x - p[0];
  Vector pa1 = x - p[1];
  Vector pa2 = x - p[2];
  Vector pa3 = x - p[3];

  // Edges: could be pre-processed in the constructor
  Vector e01 = Normalized(p[1] - p[0]);
  Vector e02 = Normalized(p[2] - p[0]);
  Vector e03 = Normalized(p[3] - p[0]);
  Vector e12 = Normalized(p[2] - p[1]);
  Vector e13 = Normalized(p[3] - p[1]);
  Vector e23 = Normalized(p[3] - p[2]);

  Vector ne012 = e01 / n[0];
  Vector ne120 = e12 / n[0];
  Vector ne201 = -(e02 / n[0]);

  // Face 031
  Vector ne013 = (n[1] / e01);
  Vector ne130 = (n[1] / e13);
  Vector ne301 = -(n[1] / e03);

  // Face 023
  Vector ne023 = -(n[2] / e02);
  Vector ne230 = -(n[2] / e23);
  Vector ne302 = (n[2] / e03);

  // Face 132
  Vector ne123 = (n[3] / e12);
  Vector ne231 = (n[3] / e23);
  Vector ne312 = -(n[3] / e13);

  // Face 012 ---------------------------------------------------
  double t = pa0 * n[0];
  if (t > 0.0)
  {
    // Offside Edge [01] Plane 012
    if (pa0 * ne012 >= 0.0)
    {
      double s = pa0 * e01;
      // Vertex [0]
      if (s <= 0.0)
      {
        if (pa0 * e02 <= 0.0)
        {
          if (pa0 * e03 <= 0.0)
          {
            return pa0;
          }
        }
      }
      // Edge [01]
      else if (pa1 * e01 < 0.0)
      {
        // Offside Edge [01] Plane 012
        if (pa0 * ne013 > 0.0)
        {
          return pa0 - s * e01;
        }
      }
      // Vertex [1]
      else
      {
        if (pa1 * e12 <= 0.0)
        {
          if (pa1 * e13 <= 0.0)
          {
            return pa1;
          }
        }
      }
    }
    // Offside Edge [12] Plane 012
    else if (pa1 * ne120 >= 0.0)
    {
      double s = pa1 * e12;
      // Vertex [1]
      if (s <= 0.0)
      {
        if (pa1 * e01 >= 0.0)
        {
          if (pa1 * e13 <= 0.0)
          {
            return pa1;
          }
        }
      }
      // Edge [12]
      else if (pa2 * e12 < 0.0)
      {
        if (pa1 * ne123 > 0.0)
        {
          return pa1 - s * e12;
        }
      }
      // Vertex [2]
      else
      {
        if (pa2 * e02 >= 0.0)
        {
          if (pa2 * e23 <= 0.0)
          {
            return pa2;
          }
        }
      }
    }
    // Offside Edge [20] Plane 012
    else if (pa2 * ne201 >= 0.0)
    {
      double s = pa2 * e02;
      // Vertex [2]
      if (s >= 0.0)
      {
        if (pa2 * e12 >= 0.0)
        {
          if (pa2 * e23 <= 0.0)
          {
            return pa2;
          }
        }
      }
      // Edge [20]
      else if (pa0 * e02 > 0.0)
      {
        if (pa2 * ne023 > 0.0)
        {
          return pa2 - s * e02;
        }
      }
      // Vertex [0]
      else
      {
        if (pa0 * e01 <= 0.0)
        {
          if (pa0 * e03 <= 0.0)
          {
            return pa0;
          }
        }
      }
    }
    else
    {
      // Face
      return t * n[0];
    }
  }

  // Face 031 ---------------------------------------------------
  t = pa0 * n[1];
  if (t > 0.0)
  {
    // Offside Edge [30] Plane 013
    if (pa0 * ne301 >= 0.0)
    {
      // Vertex [0]
      double s = pa0 * e03;
      if (s <= 0.0)
      {
        if (pa0 * e01 <= 0.0)
        {
          if (pa0 * e02 <= 0.0)
          {
            return pa0;
          }
        }
      }
      // Edge [03]
      else if (pa3 * e03 < 0.0)
      {
        if (pa0 * ne302 > 0.0)
        {
          return pa0 - s * e03;
        }
      }
      // Vertex [3]
      else
      {
        if (pa3 * e13 >= 0.0)
        {
          if (pa3 * e23 >= 0.0)
          {
            return pa3;
          }
        }
      }
    }
    else
    {
      // Offside Edge [13] Plane 013
      if (pa3 * ne130 >= 0.0)
      {
        double s = pa3 * e13;
        // Vertex [3]
        if (s >= 0.0)
        {
          if (pa3 * e03 >= 0.0)
          {
            if (pa3 * e23 >= 0.0)
            {
              return pa3;
            }
          }
        }
        // Edge [13]
        else if (pa1 * e13 > 0.0)
        {
          if (pa3 * ne312 > 0.0)
          {
            return pa3 - s * e13;
          }
        }
        // Vertex [1]
        else
        {
          if (pa1 * e01 >= 0.0)
          {
            if (pa1 * e12 <= 0.0)
            {
              return pa1;
            }
          }
        }
      }
      else
      {
        // Offside Edge [01] Plane 013
        if (pa1 * ne013 >= 0.0)
        {
          // Vertex [1]
          double s = pa1 * e01;
          if (s >= 0.0)
          {
            if (pa1 * e13 <= 0.0)
            {
              if (pa1 * e12 <= 0.0)
              {
                //std::cout<<"Vertex 1 of face 031"<<endl;
                return pa1;
              }
            }
          }
          // Edge [10]
          else if (pa0 * e01 > 0.0)
          {
            if (pa1 * ne012 > 0.0)
            {
              //std::cout<<"Edge 10"<<endl;
              return pa1 - s * e01;
            }
          }
          // Vertex [0]
          else
          {
            if (pa0 * e02 <= 0.0)
            {
              if (pa0 * e03 <= 0.0)
              {
                //std::cout<<"Vertex 0 of face 013"<<endl;
                return pa0;
              }
            }
          }
        }
        else
        {

          // Face
          return t * n[1];
        }
      }
    }
  }

  // Face 032 ---------------------------------------------------
  t = pa0 * n[2];
  if (t > 0.0)
  {
    //std::cout<<"Face 023"<<endl;
    // Offside Edge [02] Plane 023
    if (pa0 * ne023 >= 0.0)
    {
      double s = pa0 * e02;
      // Vertex [0]
      if (s <= 0.0)
      {
        if (pa0 * e03 <= 0.0)
        {
          if (pa0 * e01 <= 0.0)
          {
            //std::cout<<"Vertex 0 of face 023"<<endl;
            return pa0;
          }
        }
      }
      // Edge [02]
      else if (pa2 * e02 < 0.0)
      {
        if (pa0 * ne201 > 0.0)
        {
          //std::cout<<"Edge 02"<<endl;
          return pa0 - s * e02;
        }
      }
      // Vertex [2]
      else
      {
        if (pa2 * e12 >= 0.0)
        {
          if (pa2 * e23 <= 0.0)
          {
            //std::cout<<"Vertex 2 of face 023"<<endl;
            return pa0;
          }
        }
      }
    }
    else
    {
      // Offside Edge [23] Plane 023
      if (pa2 * ne230 >= 0.0)
      {
        double s = pa2 * e23;
        // Vertex [2]
        if (s <= 0.0)
        {
          if (pa2 * e02 >= 0.0)
          {
            if (pa2 * e12 >= 0.0)
            {
              //std::cout<<"Vertex 2 of face 023"<<endl;
              return pa2;
            }
          }
        }
        // Edge [23]
        else if (pa3 * e23 < 0.0)
        {
          if (pa2 * ne231 > 0.0)
          {
            //std::cout<<"Edge 23"<<endl;
            return pa2 - s * e23;
          }
        }
        // Vertex [3]
        else
        {
          if (pa3 * e03 >= 0.0)
          {
            if (pa3 * e13 >= 0.0)
            {
              //std::cout<<"Vertex 3 of face 023"<<endl;
              return pa0;
            }
          }
        }
      }
      else
      {
        // Offside Edge [30] Plane 023
        if (pa3 * ne302 >= 0.0)
        {
          double s = pa3 * e03;
          // Vertex [3]
          if (s >= 0.0)
          {
            if (pa3 * e23 >= 0.0)
            {
              if (pa3 * e13 >= 0.0)
              {
                //std::cout<<"Vertex 3 of face 023"<<endl;
                return pa3;
              }
            }
          }
          // Edge [30]
          else if (pa0 * e03 > 0.0)
          {
            if (pa3 * ne301 > 0.0)
            {
              //std::cout<<"Edge 30"<<endl;
              return pa3 - s * e03;
            }
          }
          // Vertex [0]
          else
          {
            if (pa0 * e01 <= 0.0)
            {
              if (pa0 * e02 <= 0.0)
              {
                //std::cout<<"Vertex 0 of face 023"<<endl;
                return pa0;
              }
            }
          }
        }
        else
        {
          // Face
          return t * n[2];
        }
      }
    }
  }

  // Face 132 ---------------------------------------------------
  t = pa1 * n[3];
  if (t > 0.0)
  {
    //std::cout<<"Face 132"<<endl;
    // Offside Edge [31] Plane 132
    if (pa1 * ne312 >= 0.0)
    {
      double s = pa1 * e13;
      // Vertex [1]
      if (s <= 0.0)
      {
        if (pa1 * e12 <= 0.0)
        {
          if (pa1 * e01 >= 0.0)
          {
            //std::cout<<"Vertex 1 of face 132"<<endl;
            return pa1;
          }
        }
      }
      // Edge [13]
      else if (pa3 * e13 < 0.0)
      {
        if (pa1 * ne130 > 0.0)
        {
          //std::cout<<"Edge 13"<<endl;
          return pa1 - s * e13;
        }
      }
      // Vertex [3]
      else
      {
        if (pa3 * e03 >= 0.0)
        {
          if (pa3 * e23 >= 0.0)
          {
            //std::cout<<"Vertex 3 of face 132"<<endl;
            return pa3;
          }
        }
      }
    }
    else
    {
      // Offside Edge [23] Plane 132
      if (pa3 * ne231 >= 0.0)
      {
        double s = pa3 * e23;
        // Vertex [3]
        if (s >= 0.0)
        {
          if (pa3 * e13 >= 0.0)
          {
            if (pa3 * e03 >= 0.0)
            {
              //std::cout<<"Vertex 3 of face 132"<<endl;
              return pa3;
            }
          }
        }
        // Edge [23]
        else if (pa2 * e23 > 0.0)
        {
          if (pa3 * ne230 > 0.0)
          {
            //std::cout<<"Edge 23"<<endl;
            return pa3 - s * e23;
          }
        }
        // Vertex [2]
        else
        {
          if (pa2 * e02 >= 0.0)
          {
            if (pa2 * e12 >= 0.0)
            {
              //std::cout<<"Vertex 2 of face 032"<<endl;
              return pa2;
            }
          }
        }
      }
      else
      {
        // Offside Edge [12] Plane 132
        if (pa2 * ne123 >= 0.0)
        {
          // Scalar product to tell side
          double s = pa2 * e12;
          // Vertex [2]
          if (s >= 0.0)
          {
            if (pa2 * e23 <= 0.0)
            {
              if (pa2 * e02 >= 0.0)
              {
                //std::cout<<"Vertex 2 of face 132"<<endl;
                return pa2;
              }
            }
          }
          // Edge [12]
          else if (pa1 * e12 > 0.0)
          {
            if (pa2 * ne120 > 0.0)
            {
              //std::cout<<"Edge 12"<<endl;
              return pa2 - s * e12;
            }
          }
          // Vertex [1]
          else
          {
            if (pa1 * e01 >= 0.0)
            {
              if (pa1 * e13 <= 0.0)
              {
                //std::cout<<"Vertex 1 of face 032"<<endl;
                return pa1;
              }
            }
          }
        }
        else
        {
          // Face
          return t * n[3];
        }
      }
    }
  }

  // Inside tetrahedron
  return Vector::Null;
}

/*!
\brief Computes the squared distance between the tetrahedron and a point.

\param x Point.
*/
double Tetrahedra::R(const Vector& x) const
{
  return SquaredNorm(Normal(x));
}

/*!
\brief Computes the signed distance between the tetrahedron and a point.

This function does not generate an Euclidean distance.

\param x Point.
*/
double Tetrahedra::Signed(const Vector& x) const
{
  double d = Math::Max(
    Math::Max(GetPlane(0).Signed(x), GetPlane(1).Signed(x)),
    Math::Max(GetPlane(2).Signed(x), GetPlane(3).Signed(x)));
  return d;
}

/*!
\brief Scale the tetrahedron.

\param s Scaling vector.
*/
Tetrahedra Tetrahedra::Scaled(const Vector& s) const
{
  return Tetrahedra(p[0].Scaled(s), p[1].Scaled(s), p[2].Scaled(s), p[3].Scaled(s));
}

/*!
\brief Uniform scaling.

\param s Scaling factor.
*/
Tetrahedra Tetrahedra::Scaled(const double& s) const
{
  return Tetrahedra(p[0] * s, p[1] * s, p[2] * s, p[3] * s);
}

/*!
\brief Translate the tetrahedron.

\param t Translation vector.
*/
Tetrahedra Tetrahedra::Translated(const Vector& t) const
{
  return Tetrahedra(p[0] + t, p[1] + t, p[2] + t, p[3] + t);
}

/*!
\brief Rotate the tetrahedron.

\param r Rotation matrix.
*/
Tetrahedra Tetrahedra::Rotated(const Matrix& r) const
{
  return Tetrahedra(r * p[0], r * p[1], r * p[2], r * p[3]);
}

/*!
\brief Compute a Poisson sphere distribution inside a tetrahedron.
\param r Radius of the sphere.
\param n Number of candidate points.
\param random %Random number generator.
*/
QVector<Vector> Tetrahedra::Poisson(const double& r, int n, Random& random) const
{
  QVector<Vector> p;

  // Collision radius
  double c = 4.0 * r * r;

  // Create instances
  for (int i = 0; i < n; i++)
  {
    Vector t = RandomInside(random);
    bool hit = false;
    for (int j = 0; j < p.size(); j++)
    {
      if (SquaredNorm(t - p.at(j)) < c)
      {
        hit = true;
        break;
      }
    }
    if (hit == false)
      p.append(t);
  }
  return p;
}