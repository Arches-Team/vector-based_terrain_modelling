#pragma once

#include "libs/color.h"
#include "libs/mesh.h"

class MeshColor : public Mesh
{
protected:
  QVector<Color> colors; //!< Array of colors.
  QVector<int> carray;  //!< Indexes.
public:
  explicit MeshColor();
  explicit MeshColor(const Mesh&);
  explicit MeshColor(const Mesh&, const QVector<Color>&, const QVector<int>&);
  explicit MeshColor(const Mesh&, const QVector<Color>&);
  ~MeshColor();

  Color GetColor(int) const;
  void SetColor(int, const Color&);
  void Merge(const MeshColor&);

  int ColorIndex(int, int) const;
  int ColorIndex(int) const;

  QVector<Color> GetColors() const;
  QVector<int> ColorIndexes() const;
};

/*!
\brief Get a color.
\param i The index of the desired color.
\return The color.
*/
inline Color MeshColor::GetColor(int i) const
{
  return colors[i];
}

/*!
\brief Set the color of a vertex.
\param i The index of the desired color.
\param c The color.
*/
inline void MeshColor::SetColor(int i, const Color& c)
{
  colors[i] = c;
}

/*!
\brief Get the array of colors.
*/
inline QVector<Color> MeshColor::GetColors() const
{
  return colors;
}

/*!
\brief Return the set of color indices.
*/
inline QVector<int> MeshColor::ColorIndexes() const
{
  return carray;
}

/*!
\brief Get the color index of a given triangle.
\param t Triangle index.
\param i Normal index.
*/
inline int MeshColor::ColorIndex(int t, int i) const
{
  return carray.at(t * 3 + i);
}

/*!
\brief Get the color index.
\param i Index.
*/
inline int MeshColor::ColorIndex(int i) const
{
  return carray.at(i);
}
