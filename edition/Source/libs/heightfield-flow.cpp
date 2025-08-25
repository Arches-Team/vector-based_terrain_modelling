// Heightfield
#include "libs/heightfield.h"

/*!
\brief Add a new direction to the flowing structure.

This function does not compute the normalized slopes.

\param point Point corresponding to the neighboring cell.
\param index Index of the neighboring cell.
\param height Relative elevation.
\param slope Slope.
*/
void FlowStruct::Add(const QPoint& point, int index, const double& height, const double& slope)
{
  q[num] = point;
  h[num] = height;
  s[num] = slope;
  //slopesum += s[n];
  mask |= 1 << index;
  i[num] = index;
  // Steepest slope
  if (num == 0)
  {
    steepest = 0;
  }
  else if (s[num] > s[steepest])
  {
    steepest = num;
  }
  num++;
}

/*!
\brief Compute the flow directions at a given point.

Compute both the flow directions using an L<SUP>p</SUP> metric, and the steepest slope.

\param p Point.
\param flow Flow information.
*/
int HeightField::CheckFlowSlope(const QPoint& p, FlowStruct& flow) const
{
  int n = 0;

  double zp = at(p);

  double slopesum = 0.0;
  flow.mask = 0;
  // By default, steepest does not exist
  flow.steepest = -1;

  for (int i = 0; i < 8; i++)
  {
    QPoint b = p + Array2::next[i];
    // Skip if point is not inside the domain
    if (!InsideVertexIndex(b))
    {
      continue;
    }

    double step = at(b) - zp;
    if (step < -HeightField::flat) // Should be 0.0, but very small negative values might crash
    {
      flow.q[n] = b;
      flow.h[n] = -step;
      flow.s[n] = -step * inverselength[i];
      slopesum += flow.s[n];
      flow.mask |= 1 << i;
      flow.i[n] = i;

      // Steepest slope
      if (n == 0)
      {
        flow.steepest = 0;
      }
      else if (flow.s[n] > flow.s[flow.steepest])
      {
        flow.steepest = n;
      }
      n++;
    }
  }

  // Relative slope
  for (int k = 0; k < n; k++)
  {
    flow.sn[k] = flow.s[k] / slopesum;
  }

  return n;
}

/*!
\brief Compute the flow directions at a given point, using an L<SUP>2</SUP> metric.

\param p Point.
\param flow Flow information.
\param power Power, set to 1.0 for diffusive flow routing, converges to steepest slope as power increases to infinity.
*/
int HeightField::CheckFlowSlopeWeighted(const QPoint& p, FlowStruct& flow, const double& power) const
{
  int n = 0;

  double zp = at(p);

  double slopesum = 0.0;
  flow.mask = 0;

  for (int i = 0; i < 8; i++)
  {
    QPoint b = p + next[i];
    // Skip if point is not inside the domain
    if (!InsideVertexIndex(b))
    {
      continue;
    }

    double step = at(b) - zp;
    if (step < -HeightField::flat) // Should be 0.0, but very small negative values might crash
    {
      flow.mask |= 1 << i;
      flow.i[n] = i;
      flow.q[n] = b;
      flow.h[n] = -step;
      flow.s[n] = -step * inverselength[i];

      // Sum of squared slopes, instead of simple weigthed average
      slopesum += Math::Pow(flow.s[n], power);
      n++;
    }
  }

  // Relative weighted squared slopes
  for (int k = 0; k < n; k++)
  {
    flow.sn[k] = Math::Pow(flow.s[k], power) / slopesum;
  }

  return n;
}

/*!
\brief Compute the flow directions at a given point.

\param p Point.
\param s Slope.
\param flow Flow information.
*/
int HeightField::CheckFlowDirectionsAngle(const QPoint& p, const double& s, FlowStruct& flow) const
{
  int n = 0;

  double zp = at(p);

  double slopesum = 0.0;
  flow.mask = 0;

  for (int i = 0; i < 8; i++)
  {
    QPoint b = p + next[i];
    // Skip if point is not inside the domain
    if (!InsideVertexIndex(b))
    {
      continue;
    }

    // Vertical step which involves slope times side length
    double step = zp - at(b) - s * celldiagonal[0] * length[i];
    if (step > 0.0)
    {
      flow.mask |= 1 << i;
      flow.i[n] = i;
      flow.q[n] = b;
      flow.h[n] = step;
      flow.s[n] = step * inverselength[i];

      slopesum += flow.s[n];
      n++;
    }
  }

  // Relative slope
  for (int k = 0; k < n; k++)
  {
    flow.sn[k] = flow.s[k] / slopesum;
  }

  return n;
}