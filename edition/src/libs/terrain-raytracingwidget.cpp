#include "libs/realtime.h"
#include "libs/scalarfield.h"
#include "libs/cubic.h"
#include "libs/plane.h"
#include "libs/cpu.h"

#include <QtGui/QKeyEvent>
#include <QtGui/QMouseEvent>

/*!
\class TerrainRaytracingWidget realtime.h
\brief %Heightfield rendering widget based on sphere-tracing.
*/

/*!
\brief Default constructor.
*/
TerrainRaytracingWidget::TerrainRaytracingWidget(QWidget* parent) : MeshWidget(parent)
{
}

/*!
\brief Destructor. Release shader program, vao and textures.
*/
TerrainRaytracingWidget::~TerrainRaytracingWidget()
{
  release_program(shaderProgram);
  glDeleteVertexArrays(1, &raytraceVAO);
  glDeleteTextures(1, &albedoTextureBuffer);

  hfBuffer.Destroy();
  shadingBuffer.Destroy();
}

/*!
\brief Initialize OpenGL, shaders and a camera centered at origin.
*/
void TerrainRaytracingWidget::initializeGL()
{
  MeshWidget::initializeGL();

  QString pPath = QString::fromStdString(std::string(SOLUTION_DIR));
  if (pPath.isEmpty())
  {
    std::cout << "MeshWidget::reloadShaders() : SOLUTION_DIR undefined" << std::endl;
    std::cin.get();
    exit(-1);
  }
  QString fullPath = pPath + "/shaders/libs/heightfield_raytrace.glsl";
  QByteArray ba = fullPath.toLocal8Bit();
  shaderProgram = read_program(ba.data());

  camera = Camera::View(Box(10.0));
  cameraAngleOfViewV = camera.GetAngleOfViewV(width(), height());

  glGenVertexArrays(1, &raytraceVAO);
  hfBuffer.Generate();
  shadingBuffer.Generate();
  glGenTextures(1, &albedoTextureBuffer);

  glBindVertexArray(0);
  glUseProgram(0);
}

/*!
\brief Renders the scene.
*/
void TerrainRaytracingWidget::paintGL()
{
  // Custom update from user
  emit _signalUpdate();

  profiler.drawCallPerFrame = 0;

  // Clear
  glClearColor(float(backgroundColor[0]), float(backgroundColor[1]), float(backgroundColor[2]), float(backgroundColor[3]));
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Move camera
  if (MoveAt)
  {
    if (stepAt < 15)
    {
      stepAt++;
      double alpha = double(stepAt) / 15.;
      camera.SetAt(currentAt * (1. - alpha) + toAt * (alpha));
    }
    else
      MoveAt = false;
  }

  // Draw
  glUseProgram(shaderProgram);

  // Uniforms - Camera (in GL coordinate system, xzy)
  glUniform3f(0, camera.Eye()[0], camera.Eye()[1], camera.Eye()[2]);
  glUniform3f(1, camera.At()[0], camera.At()[1], camera.At()[2]);
  glUniform3f(2, camera.Up()[0], camera.Up()[1], camera.Up()[2]);
  glUniform1f(3, camera.GetAngleOfViewV(width(), height()));
  glUniform2f(4, width(), height());

  // Albedo
  if (useAlbedo)
  {
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, albedoTextureBuffer);
  }

  if (hfBuffer.GetBuffer() != 0)
  {
    // Bind buffers
    hfBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 0);
    shadingBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 1);

    // Uniforms - Heightfield
    glUniform2f(5, bbox[0][0], bbox[0][1]);
    glUniform2f(6, bbox[1][0], bbox[1][1]);
    glUniform2f(7, zMin, zMax);
    glUniform1f(8, K);
    glUniform2i(9, nx, ny);
  }
  glUniform1i(10, int(useAlbedo));
  glUniform1i(11, int(useWireframe));
  glUniform1i(12, int(useCost));
  glUniform1i(13, int(useElevationShading));
  glUniform1i(14, int(useGreenBrownYellow));
  glUniform1i(16, int(useShadingBuffer));
  glUniform2f(30, float(nearAndFarPlane[0]), float(nearAndFarPlane[1]));

  // Draw
  profiler.BeginGPU();

  glBindVertexArray(raytraceVAO);
  glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
  glBindVertexArray(0);
  glUseProgram(0);
  profiler.drawCallPerFrame++;

  // Draw meshes
  RenderMeshes();
  profiler.EndGPU();

  // Draw sky
  glEnable(GL_DEPTH_TEST);
  if (useSkyShader)
  {
    RenderSky();
    profiler.drawCallPerFrame++;
    glBindVertexArray(0);
    glUseProgram(0);
  }

  // CPU Profiling
  if (profiler.enabled)
    profiler.Update();
  RenderUiPanels();

  update();
}

/*!
\brief Internal update function of the widget.

Updates the internal height buffer from the heightfield pointer, recomputes the Lipschitz constant and min/max elevations.
*/
void TerrainRaytracingWidget::SetHeightField(ScalarField2* hfPtr)
{
  hf = hfPtr;
  UpdateInternal();
}

/*!
\brief Update the GPU data from the internal CPU heightfield pointer.

This function is computationally intensive as it recomputes the Lipschitz constant of the terrain and send data to the GPU.
*/
void TerrainRaytracingWidget::UpdateInternal()
{
  makeCurrent(); // do not remove

  tmpData.resize(hf->VertexSize());

  // Min/Max elevation
  double zMinDouble, zMaxDouble;
  hf->GetRange(zMinDouble, zMaxDouble);
  zMin = float(zMinDouble);
  zMax = float(zMaxDouble);

  nx = hf->GetSizeX();
  ny = hf->GetSizeY();
  bbox = hf->GetBox();

  // This avoids problem with rendering a flat heightfield.
  if (zMin == zMax)
    zMax += 10.0f;

  // Global Lipschitz constant
  K = hf->K();

  // Heightfield data texture
  for (int i = 0; i < hf->VertexSize(); i++)
    tmpData[i] = hf->at(i);
  hfBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * hf->VertexSize(), &tmpData.front(), GL_STREAM_READ);

  for (int i = 0; i < hf->VertexSize(); i++)
    tmpData[i] = 1.0f;
  shadingBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * hf->VertexSize(), &tmpData.front(), GL_STREAM_READ);
}

/*!
\brief Reload the shaders of the widget. Useful for realtime editing and fine tuning of the rendering.
*/
void TerrainRaytracingWidget::reloadShaders()
{
  MeshWidget::reloadShaders();

  QString pPath = QString::fromStdString(std::string(SOLUTION_DIR));
  if (pPath.isEmpty())
  {
    std::cout << "MeshWidget::TerrainRaytracingWidget() : SOLUTION_DIR undefined" << std::endl;
    std::cin.get();
    exit(-1);
  }
  QString fullPath = pPath + "/shaders/libs/heightfield_raytrace.glsl";
  QByteArray ba = fullPath.toLocal8Bit();
  shaderProgram = read_program(ba.data());
}


/*!
\brief Function override from AbstractTerrainWidget.

Used to store the camera vertical angle of view, which avoids calling atan() intensively.

\param cam new camera
*/
void TerrainRaytracingWidget::SetCamera(const Camera& cam)
{
  MeshWidget::SetCamera(cam);
  cameraAngleOfViewV = camera.GetAngleOfViewV(width(), height());
}

/*!
\brief Update the albedo texture. UseAlbedo(true) should have been
called before, or no update will be performed.
\param albedo texture
*/
void TerrainRaytracingWidget::SetAlbedo(const QImage& albedo)
{
  if (hf)
  {
    makeCurrent();

    QImage albedoGL = albedo.convertedTo(QImage::Format::Format_RGBA8888);

    glUseProgram(shaderProgram);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, albedoTextureBuffer);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glTexImage2D(GL_TEXTURE_2D, 0, 3, albedoGL.width(), albedoGL.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, albedoGL.bits());
    glProgramUniform1i(shaderProgram, glGetUniformLocation(shaderProgram, "albedo"), 1);

    glUseProgram(0);
  }
}

/*!
\brief Update the internal albedo texture flag.
\param u new value for albedo flag.
*/
void TerrainRaytracingWidget::UseAlbedo(const bool& u)
{
  useAlbedo = u;
}

/*!
\brief Set the wireframe flag.
\param wire flag
*/
void TerrainRaytracingWidget::UseWireframe(const bool& wire)
{
  useWireframe = wire;
}

/*!
\brief Set the cost shading flag.
\param cost flag
*/
void TerrainRaytracingWidget::UseCost(const bool& cost)
{
  useCost = cost;
}

/*!
\brief Set the elevation shading flag.
\param u flag
*/
void TerrainRaytracingWidget::UseElevationShading(const bool& u)
{
  useElevationShading = u;
  if (useElevationShading)
    useGreenBrownYellow = false;
}

/*!
\brief Set the shading from color map flag.
\param u flag
*/
void TerrainRaytracingWidget::UseGreenBrownYellowShading(const bool& u)
{
  useGreenBrownYellow = u;
  if (useGreenBrownYellow)
    useElevationShading = false;
}

/*!
\brief Set the elevation range for the terrain.
\param min, max elevation range
*/
void TerrainRaytracingWidget::SetElevationRange(double min, double max)
{
  zMin = float(min);
  zMax = float(max);
  if (zMin == zMax)
    zMax += 10.0f;
}

/*!
\brief Set the maximum number of raymarching step.
\param s steps
*/
void TerrainRaytracingWidget::SetSteps(int s)
{
  glUseProgram(shaderProgram);
  glUniform1i(28, s);
  glUseProgram(0);
}

/*!
\brief Set the smooth shadow flag.
\param use flag
*/
void TerrainRaytracingWidget::UseSmoothShadow(const bool& use)
{
  glUseProgram(shaderProgram);
  glUniform1i(17, int(use));
  glUseProgram(0);
}

/*!
\brief Set the number of smooth shadow steps.
\param steps number of steps
*/
void TerrainRaytracingWidget::SetSmoothShadowSteps(int steps)
{
  glUseProgram(shaderProgram);
  glUniform1i(18, steps);
  glUseProgram(0);
}

/*!
\brief Set the maximum number of smooth shadow raymarching steps.
\param steps max number of steps
*/
void TerrainRaytracingWidget::SetSmoothShadowMarchingSteps(int steps)
{
  glUseProgram(shaderProgram);
  glUniform1i(19, steps);
  glUseProgram(0);
}

/*!
\brief Set the minimum marching epsilon for smooth shadow rays.
\param eps marching epsilon
*/
void TerrainRaytracingWidget::SetSmoothShadowMarchingEpsilon(float eps)
{
  glUseProgram(shaderProgram);
  glUniform1f(20, eps);
  glUseProgram(0);
}

/*!
\brief Set the smooth shadow strength.
\param strength
*/
void TerrainRaytracingWidget::SetSmoothShadowStrength(float strength)
{
  glUseProgram(shaderProgram);
  glUniform1f(21, strength);
  glUseProgram(0);
}

/*!
\brief Set the smooth shadow flag.
\param use flag
*/
void TerrainRaytracingWidget::UseSelfShadow(const bool& use)
{
  glUseProgram(shaderProgram);
  glUniform1i(22, int(use));
  glUseProgram(0);
}

/*!
\brief Set the maximum number of self shadow raymarching steps.
\param steps max number of steps
*/
void TerrainRaytracingWidget::SetSelfShadowMarchingSteps(int steps)
{
  glUseProgram(shaderProgram);
  glUniform1i(23, steps);
  glUseProgram(0);
}

/*!
\brief Set the minimum marching epsilon for self shadow rays.
\param eps marching epsilon
*/
void TerrainRaytracingWidget::SetSelfShadowMarchingEpsilon(float eps)
{
  glUseProgram(shaderProgram);
  glUniform1f(24, eps);
  glUseProgram(0);
}

/*!
\brief Set the self shadow strength.
\param strength
*/
void TerrainRaytracingWidget::SetSelfShadowStrength(float strength)
{
  glUseProgram(shaderProgram);
  glUniform1f(25, strength);
  glUseProgram(0);
}

/*!
\brief Set the minimum stepping distance for raymarching.
\param eps minimum stepping distance
*/
void TerrainRaytracingWidget::SetEpsilon(float eps)
{
  glUseProgram(shaderProgram);
  glUniform1f(29, eps);
  glUseProgram(0);
}

/*!
\brief Set the directional light.
\param light light direction
*/
void TerrainRaytracingWidget::SetLight(const Vector& light)
{
  glUseProgram(shaderProgram);
  VectorFloat l = VectorFloat(light);
  glUniform3f(26, l.x, l.y, l.z);
  glUseProgram(0);
}

/*!
\brief Set the amount of antialiasing.
\param AA antialiasing, should be between [1, 128].
*/
void TerrainRaytracingWidget::SetAntialiasing(int AA)
{
  glUseProgram(shaderProgram);
  glUniform1i(27, AA);
  glUseProgram(0);
}

/*!
\brief Manually change the underlying heightfield buffer.

This function is useful when performing erosion on the gpu - no need to retrieve the data on the CPU
and manually recreate all the buffers.

Note that the buffer must contain data whose size is equal to the base heightfield.

\param buffer initialized buffer with elevation data.
*/
void TerrainRaytracingWidget::UpdateBuffer(GLuint buffer)
{
  hfBuffer = GLBuffer(buffer);
}

/*!
\brief Set the additional shading buffer that can be used by the glsl.
\param shading 2D shading map, must be the same size as the underlying heightfield.
*/
void TerrainRaytracingWidget::SetShading(const ScalarField2& shading)
{
  for (int i = 0; i < hf->VertexSize(); i++)
    tmpData[i] = shading.at(i);
  shadingBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * hf->VertexSize(), &tmpData.front(), GL_STREAM_READ);
}

/*!
\brief Set the additional shading flag.
\param use flag
*/
void TerrainRaytracingWidget::UseShading(const bool& use)
{
  useShadingBuffer = use;
}

/*!
\brief Hot reload of the glsl program. Properly rebinds the albedo texture.
*/
void TerrainRaytracingWidget::Reload()
{
  release_program(shaderProgram);

  QString pPath = QString::fromStdString(std::string(SOLUTION_DIR));
  if (pPath.isEmpty())
  {
    std::cout << "TerrainRaytracingWidget::initializeGL() : SOLUTION_DIR undefined" << std::endl;
    std::cin.get();
    exit(-1);
  }

  QString fullPath = pPath + "/shaders/libs/heightfield_raytrace.glsl";
  QByteArray ba = fullPath.toLocal8Bit();
  shaderProgram = read_program(ba.data());

  glUseProgram(shaderProgram);
  glProgramUniform1i(shaderProgram, glGetUniformLocation(shaderProgram, "albedo"), 1);
  glUseProgram(0);
}
