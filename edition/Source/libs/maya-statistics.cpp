// MayaMaterial
#include "libs/mayageometry.h"

/*!
\class MayaStatistics maya.h
\brief Statistics about objects and instances.

\ingroup MayaCore
*/

/*!
\brief Create statistics.
\param objects Number of unique objects.
\param tobjects, vobjects Number of triangles and vertices of unique objects.
\param instances Number of instances.
\param tinstances, vinstances Number of triangles and vertices for instances.
*/
MayaStatistics::MayaStatistics(int objects, int tobjects, int vobjects, int instances, int tinstances, int vinstances)
{
  s[0] = objects;
  s[1] = tobjects;
  s[2] = vobjects;
  s[3] = instances;
  s[4] = tinstances;
  s[5] = vinstances;
}

/*!
\brief Return a string showing the statistics.
*/
QString MayaStatistics::GetText() const
{
  return QString("Statistics: %1 %2 %3 %4 %5 %6\n").arg(s[0]).arg(s[1]).arg(s[2]).arg(s[3]).arg(s[4]).arg(s[5]);
}

/*!
\brief Overloaded operator for adding and gathering statistics.
*/
MayaStatistics operator+ (const MayaStatistics& a, const MayaStatistics& b)
{
  return MayaStatistics(a.s[0] + b.s[0], a.s[1] + b.s[1], a.s[2] + b.s[2], a.s[3] + b.s[3], a.s[4] + b.s[4], a.s[5] + b.s[5]);
}

/*
\brief Clear the statistics of all the instances.
*/
void MayaStatistics::ClearInstances()
{
  s[3] = s[4] = s[5] = 0;
}

/*
\brief Add one or more instances to the statistics.
\param n Number of instances.
*/
void MayaStatistics::AddInstances(int n)
{
  s[3] += n;
  s[4] += n*s[1];
  s[5] += n*s[2];
}