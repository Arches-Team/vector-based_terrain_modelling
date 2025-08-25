#include "libs/maya.h"

/*!
\class MayaCameraSet maya.h
\brief A set of cameras.

\ingroup MayaCore
*/

/*!
\brief Create set of camera.

Initialize the set as empty and camera as default camera.
*/
MayaCameraSet::MayaCameraSet()
{
}

/*!
\brief Set active camera as the next one in the set.
*/
void MayaCameraSet::Next(Camera& camera)
{
  if (cameras.size() != 0)
  {
    n++;
    n %= cameras.size();
    camera = cameras.at(n);
  }
}

/*!
\brief Add a camera to the set.
\param c %Camera.
*/
void MayaCameraSet::Push(const Camera& c)
{
  cameras.append(c);
}

/*!
\brief Remove the current camera from the set.
*/
void MayaCameraSet::Pop()
{
  if (cameras.size() != 0)
  {
    cameras.remove(n);
    if (cameras.size() == 0)
    {
      n = 0;
    }
    else
    {
      n %= cameras.size();
    }
  }
}

/*!
\brief Get current camera.
*/
Camera MayaCameraSet::Current() const
{
  return cameras.at(n);
}

/*!
\brief Get the set of cameras.
*/
QVector<Camera> MayaCameraSet::All() const
{
  return cameras;
}

