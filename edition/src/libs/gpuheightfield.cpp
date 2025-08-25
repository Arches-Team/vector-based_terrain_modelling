#include "libs/gpuheightfield.h"

/*!
\class GPUHeightField gpuheightfield.h
\brief A heightfield stored on the GPU.

This class is a work in progress. There are several open questions:
  => Use double, or float ?
     Most of the time on the GPU we only use float, but sometimes we need precision, or we work on large terrains which require lots of precision.
     Ideally, this class should allow to switch between single and double precision, but then this means that the simulation shaders should be adapted accordingly.
     And what happens when we work on double precision for the rendering ? The buffers cannot be shared between the widget and the simulation/analysis.

  => Most of the time we work with multiple layer, which means that we might need a GPULayerField class ?
     With at least bedrock, sediment, and maybe water as layers.
*/

/*!
\brief Constructor from a heightfield.
*/
GPUHeightField::GPUHeightField(const HeightField& hf) : Array2(Box2(hf.GetBox()), hf.GetSizeX(), hf.GetSizeY())
{
  owner = true;

  std::vector<float> tmpData;
  tmpData.resize(hf.VertexSize());
  for (int i = 0; i < tmpData.size(); i++)
    tmpData[i] = float(hf.at(i));

  inBuffer.Generate();
  inBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * hf.VertexSize(), &tmpData.front(), GL_STREAM_READ);

  outBuffer.Generate();
  outBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * hf.VertexSize(), &tmpData.front(), GL_STREAM_READ);
}

/*!
\brief Constructor from an existing opengl buffer.
*/
GPUHeightField::GPUHeightField(GLuint buffer, const Box2& box, int nx, int ny) : Array2(box, nx, ny)
{
  owner = false;
  inBuffer = GLBuffer(buffer);
}

/*!
\brief Release the GPU data.
*/
void GPUHeightField::Destroy()
{
  if (owner)
  {
    inBuffer.Destroy();
    outBuffer.Destroy();
  }
}

/*!
\brief Get the elevation data back on the CPU into an initialized heightfield.
\param hf the heightfield, should be initialized to correct sizes.
*/
void GPUHeightField::GetData(HeightField& hf) const
{
  std::vector<float> tmpData;
  tmpData.resize(hf.VertexSize());
  outBuffer.GetData(tmpData);
  for (int i = 0; i < hf.VertexSize(); i++)
    hf[i] = double(tmpData[i]);
}
