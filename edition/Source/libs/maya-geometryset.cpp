// Maya

#include "libs/mayageometry.h"
#include "libs/voxel.h"

/*!
\class MayaGeometrySet maya.h
\brief The %MayaGeometrySet class provides a geometric instance linked to several frames.

\ingroup MayaCore
*/

/*!
\brief Create a set of geometric instances.

Note that no reference to the geometric object is created.

\param mg Geometric object.
*/
MayaGeometrySet::MayaGeometrySet(const MayaGeometry& mg) :MayaGeometry(mg)
{
}

/*!
\brief Create a set of geometric instances.

One reference to the geometric instance is created using the argument frame.

\param mg Geometric object.
\param mf Frame.
*/
MayaGeometrySet::MayaGeometrySet(const MayaGeometry& mg, const FrameScaled& mf) :MayaGeometry(mg)
{
  MayaGeometrySet::frames.clear();
  MayaGeometrySet::frames.append(mf);
}

/*!
\brief Create a set of geometric instances.

Several references to the geometric instance are created using the argument set of frames.

\param mg Geometric object.
\param mfs Set of translation vectors.
*/
MayaGeometrySet::MayaGeometrySet(const MayaGeometry& mg, const QVector<Vector>& mfs) :MayaGeometry(mg)
{
  MayaGeometrySet::frames.clear();
  for (int i = 0; i < mfs.size(); i++)
  {
    frames.append(FrameScaled::Translation(mfs.at(i)));
  }
}

/*!
\brief Create a set of geometric instances.

Several references to the geometric instance are created using the argument set of frames.

\param mg Geometric object.
\param mfs Set of frames.
*/
MayaGeometrySet::MayaGeometrySet(const MayaGeometry& mg, const QVector<FrameScaled>& mfs) :MayaGeometry(mg)
{
  MayaGeometrySet::frames.clear();
  MayaGeometrySet::frames = mfs;
}

/*!
\brief Rotates a set of geometrical objects.

This function simply updates the frame of the object.
\param r Rotation vector in Euler coordinates.
*/
void MayaGeometrySet::Rotate(const Vector& r)
{
  for (int i = 0; i < frames.size(); i++)
  {
    frames[i].Rotate(r);
  }
}

/*!
\brief Translates a set of geometrical objects.

This function simply updates the frame of the object.
\param t Translation vector.
*/
void MayaGeometrySet::Translate(const Vector& t)
{
  for (int i = 0; i < frames.size(); i++)
  {
    frames[i].Translate(t);
  }
}

/*!
\brief Scales a set of geometrical objects.

This function simply updates the frame of the object.
\param s Scaling vector.
*/
void MayaGeometrySet::Scale(const Vector& s)
{
  for (int i = 0; i < frames.size(); i++)
  {
    frames[i].Scale(s);
  }
}

/*!
\brief Adds a frame to the list of instances.
\param frame The frame.
*/
void MayaGeometrySet::Append(const FrameScaled& frame)
{
  frames.append(frame);
}

/*!
\brief Add a frame to the list of instances.
\param set The set of frames.
*/
void MayaGeometrySet::ApplyFrames(QVector<FrameScaled> set)
{
  QVector<FrameScaled> mfs_tmp;
  for (int i = 0; i < frames.size(); i++)
  {
    for (int j = 0; j < set.size(); j++)
    {
      FrameScaled mf = frames.at(i).Composed(set.at(j));
      mfs_tmp.append(mf);
    }
  }
  frames = mfs_tmp;
}

/*!
\brief Merge all the instances into a signel Maya geometry structure.
*/
MayaGeometry MayaGeometrySet::Collapse() const
{
  // Complete set
  MayaGeometry set;
  set.SetName(name + "_c");
  set.SetMaterial(mat);

  // Escape if there are no instances
  if (frames.size() == 0) return set;

  // Collapse all instances into the same geometry
  for (int i = 0; i < frames.size(); i++)
  {
    MayaGeometry truc(name, vertex, normal, UVmap, indexes);
    set.Merge(truc.Transform(frames.at(i)));
  }

  return set;
}

/*!
\brief Clear the instance set.
*/
void MayaGeometrySet::Clear()
{
  frames.clear();
}

/*!
\brief Compute the bounding box.
*/
Box MayaGeometrySet::GetBox() const
{
  Box b = MayaGeometry::GetBox();

  if (frames.size() == 0)
  {
    return b;
  }

  Box r(b * Math::Max(frames[0].S()[0], frames[0].S()[1], frames[0].S()[2]), frames[0]);

  for (int i = 1; i < frames.size(); i++)
  {
    r = Box(r, Box(b * Math::Max(frames[i].S()[0], frames[i].S()[1], frames[i].S()[2]), frames[i]));
  }

  return r;
}

/*!
\brief Compute the statistics of the geometry set.
\return The statistics.
*/
MayaStatistics MayaGeometrySet::GetStatistics() const
{
  return MayaStatistics(1, vertex.size(), indexes.size() / 3, frames.size(), frames.size() * vertex.size(), frames.size() * indexes.size() / 3);
}

/*!
\brief Get the text information.
\param spaces Spacing for indenting the text.
\param html Flag to specify syntax highlighting.
*/
QString MayaGeometrySet::GetText(int spaces, bool html) const
{
  return (html == true ? QString("&nbsp;") : QString(spaces, ' ')) + (html == true ? QString("<span style=\"color:#248\">") : QString("")) + QString("MayaGeometrySet") + (html == true ? QString("</span>") : QString("")) + QString(" (") + GetName() + QString(" , ") + QString("%1 )").arg(frames.size());
}

/*!
\brief Create a set of cubes representing a voxel.
\param voxel The voxel.
\param surface Only produce cubes that participate to the visible surface.
*/
MayaGeometrySet MayaGeometrySet::CreateVoxel(const Voxel& voxel, bool surface) 
{
  Box u = voxel.UnitCell();
  MayaGeometry cube = MayaGeometry("Cube", Mesh(Box(u)));

  MayaGeometrySet cubes(cube);
  cubes.SetName("VoxelCube");

  if (surface == false)
  {
    // Parse all cells
    for (int i = 0; i < voxel.CellSize(); i++)
    {
      if (voxel.At(i) != 0)
      {
        cubes.Append(FrameScaled::Translation(voxel.Cell(i).Center()));
      }
    }
  }
  // Discard cubes that do not have a frontier 
  else
  {
    // Parse all cells
    for (int i = 0; i < voxel.CellSizeX(); i++)
    {
      for (int j = 0; j < voxel.CellSizeY(); j++)
      {
        for (int k = 0; k < voxel.CellSizeZ(); k++)
        {
          // Compute configuration
          if (voxel.At(i, j, k) != 0)
          {
            bool s = false;
            if (i - 1 < 0)
            {
              s = true;
            }
            else if (voxel.At(i - 1, j, k) == 0)
            {
              s = true;
            }

            if (i + 1 > voxel.CellSizeX() - 1)
            {
              s = true;
            }
            else if (voxel.At(i + 1, j, k) == 0)
            {
              s = true;
            }
            if (j - 1 < 0)
            {
              s = true;
            }
            else if (voxel.At(i, j - 1, k) == 0)
            {
              s = true;
            }
            if (j + 1 > voxel.CellSizeY() - 1)
            {
              s = true;
            }
            else if (voxel.At(i, j + 1, k) == 0)
            {
              s = true;
            }
            if (k - 1 < 0)
            {
              s = true;
            }
            else if (voxel.At(i, j, k - 1) == 0)
            {
              s = true;
            }

            if (k + 1 > voxel.CellSizeZ() - 1)
            {
              s = true;
            }
            else if (voxel.At(i, j, k + 1) == 0)
            {
              s = true;
            }

            if (s == true)
            {
              cubes.Append(FrameScaled::Translation(voxel.Cell(i, j, k).Center()));
            }
          }
        }
      }
    }
  }
  return cubes;
}
