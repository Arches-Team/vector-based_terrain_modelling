#include "libs/gputerrainsimu.h"
#include "libs/cpu.h"

/*!
\class GPUHeightFieldAnalysis gputerrainsimu.h
\brief Parallelized heightfield analysis.
*/

/*!
\brief Constructor.

Initialize everthing to zero.
*/
GPUHeightFieldAnalysis::GPUHeightFieldAnalysis()
{
  hfTextureBuffer = 0;
  nx = ny = 0;
}

/*!
\brief Destroys OpenGl buffers and shader program.
*/
GPUHeightFieldAnalysis::~GPUHeightFieldAnalysis()
{
  shader.Destroy();
  // inHfBuffer is not destroyed because we are not owner of it
  inHfBufferDouble.Destroy();
  outBuffer.Destroy();
  outVecBuffer.Destroy();
  outIntBuffer.Destroy();
  outDoubleBuffer.Destroy();
  glDeleteTextures(1, &hfTextureBuffer);
}

/*!
\brief Init function, used in the context of an already existing GPU buffer from the QOpenGL widget.

\param nx, ny %Heightfield dimensions, should match the size of the GPU buffer.
*/
void GPUHeightFieldAnalysis::Init(int nx, int ny)
{
  this->nx = nx;
  this->ny = ny;
  box = Box::Unit;
  totalBufferSize = nx * ny;
  dispatchSize = (max(nx, ny) / 8) + 1;

  QString fullPath = System::GetResource(QString::fromStdString(std::string(SOLUTION_DIR) + "/shaders/heightfield_analysis.glsl"));
  if (fullPath.isEmpty())
  {
    std::cout << "GPUHeightFieldAnalysis::Init() : error loading shader" << std::endl;
    std::cin.get();
    exit(-1);
  }
  QByteArray ba = fullPath.toLocal8Bit();
  shader.Initialize(ba.data());

  if (hfTextureBuffer == 0)
    glGenTextures(1, &hfTextureBuffer);

  outBuffer.Generate();
  outBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, 0, GL_STREAM_READ);

  outVecBuffer.Generate();
  outVecBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * 2 * totalBufferSize, 0, GL_STREAM_READ);

  outIntBuffer.Generate();
  outIntBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(int) * 1, 0, GL_STREAM_READ);

  shader.Bind();
  shader.SetUniform("a", box[0][0], box[0][1], box[0][2]);
  shader.SetUniform("b", box[1][0], box[1][1], box[1][2]);
  shader.SetUniform("nx", nx);
  shader.SetUniform("ny", ny);
  shader.Unbind();
}

/*!
\brief Init function from a given heightfield.

Initialize a GPU texture (and not a buffer).

This function is useful when the single precision is not an issue for the calling function.

\param hf %Heightfield.
*/
void GPUHeightFieldAnalysis::InitForTexture(const HeightField& hf)
{
  nx = hf.GetSizeX();
  ny = hf.GetSizeY();
  box = hf.GetBox();
  totalBufferSize = hf.VertexSize();
  dispatchSize = (max(nx, ny) / 8) + 1;

  QString fullPath = System::GetResource(QString::fromStdString(std::string(SOLUTION_DIR) + "/shaders/heightfield_analysis.glsl"));
  if (fullPath.isEmpty())
  {
    std::cout << "GPUHeightFieldAnalysis::Init() : error loading shader" << std::endl;
    std::cin.get();
    exit(-1);
  }
  QByteArray ba = fullPath.toLocal8Bit();
  shader.Initialize(ba.data());

  if (hfTextureBuffer == 0)
    glGenTextures(1, &hfTextureBuffer);
  QImage hfTexture = hf.CreateImage(false);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, hfTextureBuffer);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexImage2D(GL_TEXTURE_2D, 0, 3, hfTexture.width(), hfTexture.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, hfTexture.bits());

  outBuffer.Generate();
  outBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, 0, GL_DYNAMIC_READ);

  outVecBuffer.Generate();
  outVecBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * 2 * totalBufferSize, 0, GL_STREAM_READ);

  outIntBuffer.Generate();
  outIntBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(int) * 1, 0, GL_STREAM_READ);

  shader.Bind();
  shader.SetUniform("heightfield", 0);
  shader.SetUniform("a", box[0][0], box[0][1], box[0][2]);
  shader.SetUniform("b", box[1][0], box[1][1], box[1][2]);
  Vector2 cellDiagonal = hf.CellDiagonal();
  shader.SetUniform("cellDiagonal", float(cellDiagonal[0]), float(cellDiagonal[1]));
  shader.SetUniform("nx", nx);
  shader.SetUniform("ny", ny);
  shader.SetUniform("lipschitz", (float)hf.K());
  shader.Unbind();
}

/*!
\brief Init function from a given heightfield.

Initialize a single precision GPU buffer (and not a texture).

This function is useful when the single precision is not an issue for the calling function. It also allows
for a simple API for functions taking both an existing buffer or a heightfield.

\param hf %Heightfield.
*/
void GPUHeightFieldAnalysis::InitForSingleBuffer(const HeightField& hf)
{
  nx = hf.GetSizeX();
  ny = hf.GetSizeY();
  box = hf.GetBox();
  totalBufferSize = hf.VertexSize();
  dispatchSize = (max(nx, ny) / 8) + 1;

  QString fullPath = System::GetResource(QString::fromStdString(std::string(SOLUTION_DIR) + "/shaders/heightfield_analysis.glsl"));
  if (fullPath.isEmpty())
  {
    std::cout << "GPUHeightFieldAnalysis::Init() : error loading shader" << std::endl;
    std::cin.get();
    exit(-1);
  }
  QByteArray ba = fullPath.toLocal8Bit();
  shader.Initialize(ba.data());

  std::vector<float> data;
  data.resize(totalBufferSize);
  for (int i = 0; i < totalBufferSize; i++)
    data[i] = float(hf.at(i));
  inHfBuffer.Generate();
  inHfBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &data.front(), GL_STREAM_READ);

  outBuffer.Generate();
  outBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, 0, GL_STREAM_READ);

  outVecBuffer.Generate();
  outVecBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * 2 * totalBufferSize, 0, GL_STREAM_READ);

  outIntBuffer.Generate();
  outIntBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(int) * 1, 0, GL_STREAM_READ);

  shader.Bind();
  shader.SetUniform("a", box[0][0], box[0][1], box[0][2]);
  shader.SetUniform("b", box[1][0], box[1][1], box[1][2]);
  Vector2 cellDiagonal = hf.CellDiagonal();
  shader.SetUniform("cellDiagonal", float(cellDiagonal[0]), float(cellDiagonal[1]));
  shader.SetUniform("nx", nx);
  shader.SetUniform("ny", ny);
  shader.SetUniform("lipschitz", (float)hf.K());
  shader.Unbind();
}

/*!
\brief Init function from a given heightfield.

Initialize a double precision GPU buffer.
\param hf %Heightfield
*/
void GPUHeightFieldAnalysis::InitForDoubleBuffer(const HeightField& hf)
{
  nx = hf.GetSizeX();
  ny = hf.GetSizeY();
  box = hf.GetBox();
  totalBufferSize = hf.VertexSize();
  dispatchSize = (max(nx, ny) / 8) + 1;

  QString fullPath = QString::fromStdString(std::string(SOLUTION_DIR) + "/shaders/heightfield_analysis.glsl");
  QByteArray ba = fullPath.toLocal8Bit();
  shader.Initialize(ba.data());

  std::vector<double> data;
  data.resize(totalBufferSize);
  for (int i = 0; i < totalBufferSize; i++)
    data[i] = hf.at(i);
  inHfBufferDouble.Generate();
  inHfBufferDouble.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(double) * totalBufferSize, &data.front(), GL_STREAM_READ);

  outBuffer.Generate();
  outBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, 0, GL_STREAM_READ);

  outVecBuffer.Generate();
  outVecBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * 2 * totalBufferSize, 0, GL_STREAM_READ);

  outIntBuffer.Generate();
  outIntBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(int) * 1, 0, GL_STREAM_READ);

  outDoubleBuffer.Generate();
  outDoubleBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(double) * totalBufferSize, 0, GL_STREAM_READ);

  shader.Bind();
  shader.SetUniform("a", box[0][0], box[0][1], box[0][2]);
  shader.SetUniform("b", box[1][0], box[1][1], box[1][2]);
  Vector2 cellDiagonal = hf.CellDiagonal();
  shader.SetUniform("cellDiagonalDouble", cellDiagonal[0], cellDiagonal[1]);
  shader.SetUniform("nx", nx);
  shader.SetUniform("ny", ny);
  shader.SetUniform("lipschitz", (float)hf.K());
  shader.Unbind();
}

/*!
\brief Compute the fractional Laplacian.
\param hf %Heightfield
\param m mask size
\param s Fractional parameter, should be in [0,1].
\param N not used at the moment.
*/
ScalarField2 GPUHeightFieldAnalysis::FractionalLaplacian(const HeightField& hf, int m, double s, int /*N*/)
{
  InitForTexture(hf);

  shader.Bind();

  // Uniforms
  shader.SetUniform("func", 1);
  shader.SetUniform("m", m);
  shader.SetUniform("sFrac", (float)s);

  // Buffers
  outBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 0);
  outVecBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 1);
  outIntBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 3);

  glDispatchCompute(dispatchSize, dispatchSize, 1);
  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  // Temporary float array for OpenGL
  std::vector<float> data;
  data.resize(totalBufferSize);
  outBuffer.GetData(data);

  ScalarField2 ret = ScalarField2(Box2(box), nx, ny, 0.0);
  for (int i = 0; i < totalBufferSize; i++)
    ret[i] = double(data[i]);

  shader.Unbind();
  return ret;
}

/*!
\brief Computes the Laplacian of the heightfield in double precision on the GPU.
\param hf the heightfield
*/
ScalarField2 GPUHeightFieldAnalysis::Laplacian(const HeightField& hf)
{
  InitForDoubleBuffer(hf);

  shader.Bind();

  // Uniforms
  shader.SetUniform("func", 9);

  // Buffers
  outBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 0);
  outVecBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 1);
  outIntBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 3);
  inHfBufferDouble.BindAt(GL_SHADER_STORAGE_BUFFER, 4);
  outDoubleBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 5);

  glDispatchCompute(dispatchSize, dispatchSize, 1);
  glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

  // Get data
  std::vector<double> data;
  data.resize(totalBufferSize);
  outDoubleBuffer.GetData(data);
  ScalarField2 laplacian = ScalarField2(Box2(box), nx, ny, 0.0);
  for (int i = 0; i < totalBufferSize; i++)
    laplacian[i] = data[i];

  shader.Unbind();
  return laplacian;
}

/*!
\brief Compute fractional gradient and fractional gradient norm fields.
\param hf heightfield
\param s Fractional parameter, should be in [0,1].
\param m mask size
\param grad returned fractional gradient field
\param gradNorm returned fractional gradient norm
*/
void GPUHeightFieldAnalysis::FractionalGradient(const HeightField& hf, double s, int m, VectorField2& grad, ScalarField2& gradNorm)
{
  InitForTexture(hf);

  shader.Bind();

  // Uniforms
  shader.SetUniform("func", 2);
  shader.SetUniform("m", m);
  shader.SetUniform("sFrac", (float)s);

  // Buffers
  outBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 0);
  outVecBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 1);
  outIntBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 3);

  glDispatchCompute(dispatchSize, dispatchSize, 1);
  glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

  // Get data scalar
  std::vector<float> data;
  data.resize(totalBufferSize);
  outBuffer.GetData(data);
  gradNorm = ScalarField2(Box2(box), nx, ny, 0.0);
  for (int i = 0; i < totalBufferSize; i++)
    gradNorm[i] = double(data[i]);

  // Get data vector 2D
  data.resize(totalBufferSize * 2);
  outVecBuffer.GetData(data);
  grad = VectorField2(Box2(box), nx, ny, Vector2::Null);
  for (int i = 0; i < totalBufferSize; i++)
    grad[i] = Vector2(double(data[i * 2 + 0]), double(data[i * 2 + 1]));

  shader.Unbind();
}

/*!
\brief Computes the gradient vector field and gradient norm scalar field of a given heightfield.
\param hf the heightfield
\param grad returned gradient vector field
\param gradNorm returned gardient norm scalar field
*/
void GPUHeightFieldAnalysis::Gradient(const HeightField& hf, VectorField2& grad, ScalarField2& gradNorm)
{
  InitForDoubleBuffer(hf);

  shader.Bind();

  // Uniforms
  shader.SetUniform("func", 8);

  // Buffers
  outBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 0);
  outVecBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 1);
  outIntBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 3);
  inHfBufferDouble.BindAt(GL_SHADER_STORAGE_BUFFER, 4);
  outDoubleBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 5);

  glDispatchCompute(dispatchSize, dispatchSize, 1);
  glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

  // Get data scalar
  std::vector<double> data;
  data.resize(totalBufferSize);
  outDoubleBuffer.GetData(data);
  gradNorm = ScalarField2(Box2(box), nx, ny, 0.0);
  for (int i = 0; i < totalBufferSize; i++)
    gradNorm[i] = data[i];

  // Get data vector 2D
  std::vector<float> dataFloat;
  dataFloat.resize(totalBufferSize * 2);
  outVecBuffer.GetData(dataFloat);
  grad = VectorField2(Box2(box), nx, ny, Vector2::Null);
  for (int i = 0; i < totalBufferSize; i++)
    grad[i] = Vector2(double(dataFloat[i * 2 + 0]), double(dataFloat[i * 2 + 1]));

  shader.Unbind();
}

/*!
\brief Clear sky.
\param hf %Heightfield.
\param r maximum distance along a ray.
\param n Number of samples per cell.
*/
ScalarField2 GPUHeightFieldAnalysis::ClearSky(const HeightField& hf, double r, int n)
{
  InitForTexture(hf);

  shader.Bind();

  // Uniforms
  shader.SetUniform("func", 3);
  shader.SetUniform("r", float(r));
  shader.SetUniform("N", n);

  // Buffers
  outBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 0);
  outVecBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 1);
  outIntBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 3);

  glDispatchCompute(dispatchSize, dispatchSize, 1);
  glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

  // Temporary float array for OpenGL
  std::vector<float> data;
  data.resize(totalBufferSize);
  outBuffer.GetData(data);

  ScalarField2 ret = ScalarField2(Box2(box), nx, ny, 0.0);
  for (int i = 0; i < totalBufferSize; i++)
    ret[i] = double(data[i]);

  shader.Unbind();
  return ret;
}

/*!
\brief Compute accessibility.

\sa HeightField::Accessbility
*/
ScalarField2 GPUHeightFieldAnalysis::Accessibility(const HeightField& hf, const double& r, int n)
{
  InitForTexture(hf);

  shader.Bind();

  // Uniforms
  shader.SetUniform("func", 0);
  shader.SetUniform("r", float(r));
  shader.SetUniform("N", n);

  // Buffers
  outBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 0);
  outVecBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 1);
  outIntBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 3);

  glDispatchCompute(dispatchSize, dispatchSize, 1);
  glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

  // Temporary float array for OpenGL
  std::vector<float> data;
  data.resize(totalBufferSize);
  outBuffer.GetData(data);

  ScalarField2 ret = ScalarField2(Box2(box), nx, ny, 0.0);
  for (int i = 0; i < totalBufferSize; i++)
    ret[i] = double(data[i]);

  shader.Unbind();
  return ret;
}

/*!
\brief Compute soft shadows.
\param hf %Heightfield.
\param l Light direction, should be normalized
\param s Shadow strength
*/
ScalarField2 GPUHeightFieldAnalysis::SoftShadows(const HeightField& hf, const Vector& l, double s)
{
  InitForTexture(hf);
  shader.Bind();

  // Uniforms
  shader.SetUniform("func", 5);
  shader.SetUniform("aSoft", float(s));
  shader.SetUniform("light", l[0], l[1], l[2]);

  // Buffers
  outBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 0);
  outVecBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 1);
  outIntBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 3);

  glDispatchCompute(dispatchSize, dispatchSize, 1);
  glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

  // Temporary float array for OpenGL
  std::vector<float> data;
  data.resize(totalBufferSize);
  outBuffer.GetData(data);
  ScalarField2 ret = ScalarField2(Box2(box), nx, ny, 0.0);
  for (int i = 0; i < totalBufferSize; i++)
    ret[i] = double(data[i]);

  shader.Unbind();
  return ret;
}

/*!
\brief Compute shadows.

\param hf %Heightfield.
\param l Light direction, should be normalized
*/
ScalarField2 GPUHeightFieldAnalysis::Shadow(const HeightField& hf, const Vector& l)
{
  InitForTexture(hf);
  shader.Bind();

  // Uniforms
  shader.SetUniform("func", 13);
  shader.SetUniform("light", l[0], l[1], l[2]);

  // Buffers
  outBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 0);
  outVecBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 1);
  outIntBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 3);

  glDispatchCompute(dispatchSize, dispatchSize, 1);
  glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

  ScalarField2 ret = GetDataToScalarField();

  shader.Unbind();
  return ret;
}

/*!
\brief Compute diffuse shading.

\param hf %Heightfield.
\param l Light direction, should be normalized
\param s Shadow strength
*/
ScalarField2 GPUHeightFieldAnalysis::Diffuse(const HeightField& hf, const Vector& l, double s)
{
  InitForTexture(hf);
  shader.Bind();

  // Uniforms
  shader.SetUniform("func", 12);
  shader.SetUniform("diffusepower", float(s));
  shader.SetUniform("light", l[0], l[1], l[2]);

  // Buffers
  outBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 0);
  outVecBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 1);
  outIntBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 3);

  glDispatchCompute(dispatchSize, dispatchSize, 1);
  glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

  ScalarField2 ret = GetDataToScalarField();

  shader.Unbind();
  return ret;
}

/*!
\brief Get float buffer data from the GPU and convert it to the ScalarField2 double format.
*/
ScalarField2 GPUHeightFieldAnalysis::GetDataToScalarField() const
{
  // Temporary float array for OpenGL
  std::vector<float> data;
  data.resize(totalBufferSize);
  outBuffer.GetData(data);

  ScalarField2 ret = ScalarField2(Box2(box), nx, ny, 0.0);
  for (int i = 0; i < totalBufferSize; i++)
  {
    ret[i] = double(data[i]);
  }
  return ret;
}

/*!
\brief Compute the amount of pits in a given GPU heightfield buffer with single precision.
Results may differ from a CPU computation to the floating-point precision.

\param hfBuffer %Heightfield single precision gpu buffer, should already be initialized.
\param nx, ny %Heightfield resolution.
*/
int GPUHeightFieldAnalysis::CountPits(GLuint hfBuffer, int nx, int ny)
{
  Init(nx, ny);

  shader.Bind();

  // Uniforms
  shader.SetUniform("func", 4);
  shader.SetUniform("useDouble", 0);

  // Buffers
  outBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 0);
  outVecBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 1);
  glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, hfBuffer);
  outIntBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 3);

  glDispatchCompute(dispatchSize, dispatchSize, 1);
  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  int pitCount;
  glGetNamedBufferSubData(outIntBuffer.GetBuffer(), 0, sizeof(int) * 1, &pitCount);

  shader.Unbind();
  return pitCount;
}

/*!
\brief Compute the amount of pits in a given GPU heightfield.

This function uses a double precision buffer, so the results should be equivalent to a CPU computation.

\param hf %Heightfield.
*/
int GPUHeightFieldAnalysis::CountPits(const HeightField& hf)
{
  InitForDoubleBuffer(hf);

  shader.Bind();

  // Uniforms
  shader.SetUniform("func", 4);
  shader.SetUniform("useDouble", 1);

  // Buffers
  outBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 0);
  outVecBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 1);
  outIntBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 3);
  inHfBufferDouble.BindAt(GL_SHADER_STORAGE_BUFFER, 4);

  glDispatchCompute(dispatchSize, dispatchSize, 1);
  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  int pitCount;
  glGetNamedBufferSubData(outIntBuffer.GetBuffer(), 0, sizeof(int) * 1, &pitCount);

  shader.Unbind();
  inHfBufferDouble.Destroy();
  return pitCount;
}

/*!
\brief Computes the range of a heightfield on the GPU.
\param hf the %HeightField
*/
Vector2 GPUHeightFieldAnalysis::GetRange(const HeightField& hf)
{
  InitForSingleBuffer(hf);

  shader.Bind();

  // Uniforms
  shader.SetUniform("func", 6);
  shader.SetUniform("groupSize", 256);

  // Buffers
  outBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 0);
  outVecBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 1);
  inHfBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 2);
  outIntBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 3);

  // Part 1: map
  glDispatchCompute(64, 1, 1);
  glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

  // Part 2: reduce
  shader.SetUniform("func", 7);
  glDispatchCompute(1, 1, 1);

  // Get range as a 2D vector
  std::vector<float> data;
  data.resize(2);
  outVecBuffer.GetData(data);

  shader.Unbind();
  inHfBuffer.Destroy();
  return Vector2(data[0], data[1]);
}

/*!
\brief Computes the range of a heightfield given an already initialized GPU heightfield buffer.

The architecture is not the most efficient, but it avoids the need to get the data back on the CPU
for computing the elevation range, which is slow for large heightfields.

\param hfBuffer %Heightfield gpu buffer, should be initialized.
\param nx, ny %Heightfield dimensions.
*/
Vector2 GPUHeightFieldAnalysis::GetRange(GLuint hfBuffer, int nx, int ny)
{
  Init(nx, ny);

  shader.Bind();

  // Uniforms
  shader.SetUniform("func", 6);
  shader.SetUniform("groupSize", 256);

  // Buffers
  outBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 0);
  outVecBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 1);
  glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, hfBuffer);
  outIntBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 3);

  // Part 1: map
  glDispatchCompute(64, 1, 1);
  glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

  // Part 2: reduce
  shader.SetUniform("func", 7);
  glDispatchCompute(1, 1, 1);

  // Get range as a 2D vector
  std::vector<float> data;
  data.resize(2);
  outVecBuffer.GetData(data);

  shader.Unbind();
  return Vector2(data[0], data[1]);
}

/*!
\brief Not used in the class, but here as a reminder for the fractional gradient computation.
*/
double GPUHeightFieldAnalysis::C(int n, const double& s) const
{
  return (s * pow(4.0, s) / pow(Math::Pi, n / 2.0)) * (tgamma(s + n / 2.0) / tgamma(1.0 - s));
}


ScalarField2 GPUHeightFieldAnalysis::Viewshed(const HeightField& hf, const QPoint& p) {
  InitForTexture(hf);

  shader.Bind();

  // Uniforms
  shader.SetUniform("func", 10);
  shader.SetUniform("view_point", p.x(), p.y());

  // Buffers
  outBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 0);
  outVecBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 1);
  outIntBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 3);

  glDispatchCompute(dispatchSize, dispatchSize, 1);
  glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

  // Temporary float array for OpenGL
  std::vector<float> data;
  data.resize(totalBufferSize);
  outBuffer.GetData(data);

  ScalarField2 ret = ScalarField2(Box2(box), nx, ny, 0.0);
  for (int i = 0; i < totalBufferSize; i++)
    ret[i] = double(data[i]);

  shader.Unbind();
  return ret;
}

ScalarField2 GPUHeightFieldAnalysis::VisibilityIndex(const HeightField& hf)
{
  InitForTexture(hf);

  shader.Bind();


  // Uniforms
  shader.SetUniform("func", 11);

  for (int i = 0; i < hf.GetSizeX(); i++)
  {
    std::cout << "viz gpu " << i << std::endl;
    for (int j = 0; j < hf.GetSizeY(); j++)
    {
      shader.SetUniform("view_point", i, j);
      outBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 0);

      // Buffers

      glDispatchCompute(dispatchSize, dispatchSize, 1);
      glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT); //
    }
    glFinish();
  }

  std::cout << "fin for" << std::endl;
  std::vector<float> data(totalBufferSize, 0.);
  outBuffer.GetData(data);
  std::cout << "fin get data" << std::endl;


  // Temporary float array for OpenGL

  ScalarField2 ret = ScalarField2(Box2(box), nx, ny, 0.0);
  for (int i = 0; i < totalBufferSize; i++)
    ret[i] = double(data[i]);

  shader.Unbind();
  return ret;
}