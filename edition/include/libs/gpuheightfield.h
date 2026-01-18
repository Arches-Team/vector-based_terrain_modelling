#pragma once

#include "libs/gpu-wrapper.h"
#include "libs/array.h"
#include "libs/heightfield.h"

class GPUHeightField : public Array2
{
protected:
  GLBuffer inBuffer = 0;
  GLBuffer outBuffer = 0;
  bool owner = true;

public:
  GPUHeightField() { }
  explicit GPUHeightField(const HeightField&);
  GPUHeightField(GLuint buffer, const Box2&, int, int);
  ~GPUHeightField() { }

  void Swap();
  void Destroy();
  void GetData(HeightField&) const;
  GLBuffer GetInBuffer() const;
  GLBuffer GetOutBuffer() const;
};

/*!
\brief Get the input GPU elevation buffer.
*/
inline GLBuffer GPUHeightField::GetInBuffer() const
{
  return inBuffer;
}

/*!
\brief Get the output GPU elevation buffer.
*/
inline GLBuffer GPUHeightField::GetOutBuffer() const
{
  return outBuffer;
}

/*!
\brief
*/
inline void GPUHeightField::Swap()
{
  std::swap(inBuffer, outBuffer);
}
