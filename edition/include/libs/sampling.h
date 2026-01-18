// Sampling  

#pragma once

#include <QtCore/QVector>

#include "libs/box.h"
#include "libs/circle.h"
#include "libs/random.h"
#include "libs/sphere.h"

class ScalarField2;

class DiscTile
{
protected:
  double s = 0.0; //!< Size.
  double r = 0.0; //!< Radius.
  double e = 0.0; //!< Squared radius with small epsilon tolerance.
  QVector<Vector2> p; //!< Set of samples.
public:
  //! Empty
  explicit DiscTile();
  explicit DiscTile(const double&, const double&, const QVector<Vector2>&);
  explicit DiscTile(const double&, const double&, int, Random & = Random::R239);
  explicit DiscTile(const double&, const double&, int, const QVector<Vector2>&, Random & = Random::R239);

  //! Empty
  ~DiscTile() {}
  void Relaxation(const double& = 1.0 / 4.0);
  void Relaxation(int, const double& = 1.0 / 4.0);

  void Scramble();

  void Draw(QGraphicsScene&) const;

  Box2 GetBox() const;

  Vector2 Vertex(int) const;
  Vector2 Vertex(int, int, int) const;
  Vector2 Vertex(const QPoint&, int) const;

  const QVector<Vector2> GetSet() const;

  int Size() const;

  QVector<Vector2> Sample(const Box2&, double(*)(const Vector2&)) const;
  QVector<Vector2> Sample(const ScalarField2&) const;
  QVector<Vector2> Sample(const Box2&) const;
  QVector<Vector2> Sample(const Circle2&) const;

  static bool Check(const Vector2&, const QVector<Vector2>&, double);
  double Packing() const;

  friend std::ostream& operator<<(std::ostream&, const DiscTile&);

protected:
  void Generate(int, Random & = Random::R239);
  bool Intersect(const Vector2&, bool = true) const;
private:
  bool Check(const Vector2&) const;
public:
  // Provide range-based for loops
  auto begin() { return p.begin(); }
  auto end() { return p.end(); }
  auto cbegin() const { return p.begin(); }
  auto cend() const { return p.end(); }
  auto begin() const { return p.begin(); }
  auto end() const { return p.end(); }
};

/*!
\brief Return the box.
*/
inline Box2 DiscTile::GetBox() const
{
  return Box2(Vector2::Null, Vector2(s));
}

/*!
\brief Getter on the i-th sample.
\param i Sample index.
*/
inline Vector2 DiscTile::Vertex(int i) const
{
  return p.at(i);
}

/*!
\brief Compute the i-th sample for a displaced tile.
\param x,y Integer coordinates of the tile.
\param i Sample index.
*/
inline Vector2 DiscTile::Vertex(int x, int y, int i) const
{
  return GetBox().Tile(x, y)[0] + p.at(i);
}

/*!
\brief Compute the i-th sample for a displaced tile.
\param t Integer coordinates of the tile.
\param i Sample index.
*/
inline Vector2 DiscTile::Vertex(const QPoint& t, int i) const
{
  return GetBox().Tile(t.x(), t.y())[0] + p.at(i);
}

/*!
\brief Return the number of samples.
*/
inline int DiscTile::Size() const
{
  return p.size();
}

/*!
\brief Return the vector containing the samples.
*/
inline const QVector<Vector2> DiscTile::GetSet() const
{
  return p;
}

class SphereTile
{
protected:
  double a = 0.0; //!< Size.
  double r = 0.0; //!< Radius.
  double e = 0.0; //!< Squared radius with small epsilon tolerance.
  QVector<Vector> p; //!< Array of samples.
  Random random; //!< %Random number generator.
public:
  //! Empty
  explicit SphereTile(const Random & = Random(239));
  explicit SphereTile(const double&, const double&, int, const Random & = Random(239));
  //! Empty
  ~SphereTile() {}
  void Relaxation(const double& = 1.0 / 4.0);

  Vector Vertex(int) const;
  Sphere At(int) const;

  const QVector<Vector> GetSet() const;

  int Size() const;
  Box GetBox() const;

  static bool Check(const Vector&, const QVector<Vector>&, double);

protected:
  void Generate(int);
};

/*!
\brief Getter on the i-th sample.
\param i Sample index.
*/
inline Vector SphereTile::Vertex(int i) const
{
  return p.at(i);
}

/*!
\brief Getter on the i-th sample as a sphere.
\param i Sample index.
*/
inline Sphere SphereTile::At(int i) const
{
  return Sphere(p.at(i), r);
}

/*!
\brief Return the number of samples.
*/
inline int SphereTile::Size() const
{
  return p.size();
}

/*!
\brief Return the vector containing the samples.
*/
inline const QVector<Vector> SphereTile::GetSet() const
{
  return p;
}

/*!
\brief Return the box.
*/
inline Box SphereTile::GetBox() const
{
  return Box(Vector::Null, Vector(a));
}


class DiscTiles {
protected:
  double s = 0.0; //!< Size.
  double r = 0.0; //!< Radius.
  double e = 0.0; //!< Squared radius with small epsilon tolerance.
  DiscTile tiles[2]; //!< Tiles.
public:
  explicit DiscTiles() {}
  explicit DiscTiles(const double&, const double&, int, Random & = Random::R239);
  //! Empty.
  ~DiscTiles() {}

  Box2 GetBox() const;

  Vector2 Vertex(int, int) const;
  Vector2 Vertex(int, int, int) const;

  int Size(int) const;

  QVector<Vector2> Sample(const ScalarField2&) const;
  QVector<Vector2> Sample(const Box2&) const;

  void DrawDebug() const;

  double Packing() const;
private:
  int Select(int, int) const;
};

class DiscTileCorner {
protected:
  double s = 0.0; //!< Size.
  double r = 0.0; //!< Radius.
  double e = 0.0; //!< Squared radius with small epsilon tolerance.
  DiscTile tiles[16]; //!< Tiles.
public:
  explicit DiscTileCorner() {}
  explicit DiscTileCorner(const double&, const double&, int, Random & = Random::R239);
  //! Empty.
  ~DiscTileCorner() {}

  Box2 GetBox() const;

  Vector2 Vertex(int, int) const;
  Vector2 Vertex(int, int, int) const;

  int Size(int) const;
  double SideLength() const;
  double Radius() const;

  QVector<Vector2> Sample(const ScalarField2&) const;
  QVector<Vector2> Sample(const Box2&) const;

  DiscTile Get(int) const;
  double Packing() const;
public:
  int Select(int, int) const;
};

/*!
\brief Return the i-th disc tile.
*/
inline DiscTile DiscTileCorner::Get(int i) const
{
  return tiles[i];
}

/*!
\brief Return the side length of the tile.
*/
inline double DiscTileCorner::SideLength() const
{
  return s;
}

/*!
\brief Return the radius of particles.
*/
inline double DiscTileCorner::Radius() const
{
  return r;
}
