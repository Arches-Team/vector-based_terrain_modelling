#include "libs/realtime.h"
#include "libs/cpu.h"
#include <QtWidgets/QApplication>
#include <fstream>

/*!
\brief Default constructor.
\author Hubert-Brierre Pierre
*/
SphereTracingWidget::SDFGL::SDFGL()
{
  enabled = true;
  SetFrame(FrameScaled::Id);
  shader = GLShader();
  shading = 0;
  sampleFieldNormal = Vector(0, 0, 1);
  sampleFieldPoint = Vector(0, 0, 0.5);
  sampleFieldFreq = 8;
  maxCost = 256;
  r = Random();
  nb_sample = 0;
  old_transformation = Matrix4Float();

  resetColorTexture(512, 512);
}


/*!
\brief constructor from a sdf glsl code
\author Hubert-Brierre Pierre
*/
SphereTracingWidget::SDFGL::SDFGL(const QString& sdf_code_glsl, const FrameScaled& fr) : SDFGL()
{
  SetFrame(fr);

  //shader = GLShader();

  std::ifstream in("../../libs/LibMaya/Shaders/sphereTracing.glsl");
  QString pPath = QString("../../libs");
  if (in.good() == false)
  {
    std::cout << "Fail to find local shader sphere tracing shader. Searching for a global shader." << std::endl;
    pPath = System::GetResource("ARCHESLIBDIR");
    if (pPath.isEmpty())
    {
      std::cout << "MeshWidget::ReloadShaders() : variable d'environnement ARCHESLIBDIR non défini" << std::endl;
      std::cin.get();
      exit(-1);
    }
  }

  QString fullPath = pPath + QString("/LibMaya/Shaders/sphereTracing.glsl");
  QByteArray ba = fullPath.toLocal8Bit();

  //shader.Initialize(ba.constData());
  shader.Initialize(ba.constData(), sdf_code_glsl.toLocal8Bit().constData());

  //std::string source = prepare_source(common_source, std::string(definitions).append("#define ").append(shader_keys[i]).append("\n"));
}


/*!
\brief Reset the color texture with the given sizes
\author Hubert-Brierre Pierre
\param with : with of the color texture
\param height : height of the color texture
*/
void SphereTracingWidget::SDFGL::resetColorTexture(int width, int height) {
    nb_sample = 0; //reset the colors;
    textureWidth = width;
    textureHeight = height;

    glDeleteTextures(1, &sexyColorTexture);
    glDeleteFramebuffers(1, &sexyFramebuffer);

    //crate texture
    glActiveTexture(GL_TEXTURE0);
    glGenTextures(1, &sexyColorTexture);
    glBindTexture(GL_TEXTURE_2D, sexyColorTexture);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0); // no mipmap

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, width, height, 0, GL_RGBA, GL_FLOAT, nullptr);

    //create FBO to attach the texture and render to it
    glGenFramebuffers(1, &sexyFramebuffer);
    glBindFramebuffer(GL_READ_FRAMEBUFFER, sexyFramebuffer);
    glFramebufferTexture(GL_READ_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, sexyColorTexture, 0);
    glReadBuffer(GL_COLOR_ATTACHMENT0);
    glActiveTexture(GL_TEXTURE0);


    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glBindTexture(GL_TEXTURE_2D, 0);
}

/*!
\brief Set the frame for the sdf.
\author Hubert-Brierre Pierre
\param fr new frame
*/
void SphereTracingWidget::SDFGL::SetFrame(const FrameScaled& fr)
{
  fr.GetMatrix4().Float(TRSMatrix);
}

/*!
\brief set the shading
\param 0 : normal shading
\param 1 : nb steps shading
\author Hubert-Brierre Pierre
*/
void SphereTracingWidget::SDFGL::setShading(int newShading)
{
  shading = newShading;
}

/*!
\brief Delete the sdf shader
*/
void SphereTracingWidget::SDFGL::Delete()
{
  shader.Destroy();
}




/*!
\brief Default constructor.
*/
SphereTracingWidget::PointGL::PointGL()
{
  vao = 0;
  fullBuffer = GLBuffer();
  indexBuffer = GLBuffer();
  pointCount = 0;
  SetFrame(FrameScaled::Id);
}

/*!
\brief Constructor from a List of points and a frame scaled.
\author Hubert-Brierre Pierre
*/
SphereTracingWidget::PointGL::PointGL(const std::vector<Vector> points, const std::vector<float> nbRecCalls, const FrameScaled& fr) : PointGL()
{
  SetFrame(fr);

  int nbVertex = points.size();
  int singleBufferSize = nbVertex * 3;
  float* vertices = new float[singleBufferSize];
  float* colors = new float[nbVertex];
  for (int i = 0; i < nbVertex; i++)
  {
    //int indexVertex = vertexIndexes[i];

    Vector vertex = points[i];
    vertices[i * 3 + 0] = float(vertex[0]);
    vertices[i * 3 + 1] = float(vertex[1]);
    vertices[i * 3 + 2] = float(vertex[2]);

    colors[i] = nbRecCalls[i];

  }
  pointCount = nbVertex;

  // Generate vao & buffers
  if (vao == 0)
    glGenVertexArrays(1, &vao);
  fullBuffer.Generate();

  glBindVertexArray(vao);
  size_t fullSize = sizeof(float) * singleBufferSize + sizeof(float) * nbVertex;
  fullBuffer.SetData(GL_ARRAY_BUFFER, fullSize, nullptr);

  // Vertices(0)
  size_t size = 0;
  size_t offset = 0;
  size = sizeof(float) * singleBufferSize;
  fullBuffer.SetSubData(GL_ARRAY_BUFFER, offset, size, vertices);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (const void*)offset);
  glEnableVertexAttribArray(0);



  // Color(1)
  offset = offset + size;
  size = sizeof(float) * nbVertex;
  fullBuffer.SetSubData(GL_ARRAY_BUFFER, offset, size, colors);
  glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 0, (const void*)offset);
  glEnableVertexAttribArray(1);

  // Free data
  delete[] vertices;
  delete[] colors;
}




/*!
\brief Set the frame for the mesh.
\author Hubert-Brierre Pierre
\param fr New frame.
*/
void SphereTracingWidget::PointGL::SetFrame(const FrameScaled& fr)
{
  fr.GetMatrix4().Float(TRSMatrix);
}

/*!
\brief Delete all opengl buffers and vao.
\author Hubert-Brierre Pierre
*/
void SphereTracingWidget::PointGL::Delete()
{
  glDeleteVertexArrays(1, &vao);
  fullBuffer.Destroy();
  indexBuffer.Destroy();
}





/*!
\class SphereTracingWidget realtime.h
\brief SDF rendering widget based on sphere-tracing.
\author Hubert-Brierre Pierre
*/

/*!
\brief Default constructor.
*/
SphereTracingWidget::SphereTracingWidget(QWidget* parent) : MeshWidget(parent)
{
}

/*!
\brief Destructor. Release shader program.
*/
SphereTracingWidget::~SphereTracingWidget()
{
  ClearAll();
  pointShader.Destroy();
}

/*!
\brief Add a new sdf in the scene.
\param sdf_code_glsl new sdf glsl code
\param frame sdf frame, identity by default.
*/
void SphereTracingWidget::AddSDF(const QString& name, const QString& sdf_code_glsl)
{
  makeCurrent();
  sdfShaders.insert(name, new SDFGL(sdf_code_glsl));
}



/*!
\brief Check if a sdf exists in the widget.
\param name sdf name
*/
bool SphereTracingWidget::HasSDF(const QString& name) const
{
  return sdfShaders.contains(name);
}

/*!
\brief Delete a sdf in the scene from its name.
\param sdf mesh sdf
*/
void SphereTracingWidget::DeleteSDF(const QString& name)
{
  makeCurrent();
  if (sdfShaders.contains(name))
  {
    sdfShaders[name]->Delete();
    sdfShaders.remove(name);
  }
}

/*!
\brief Updates the transform of a sdf given its name.
\param name sdf name
\param frame new frame
*/
void SphereTracingWidget::UpdateSDF(const QString& name, const FrameScaled& frame)
{
  makeCurrent();
  if (sdfShaders.contains(name))
  {
    sdfShaders[name]->SetFrame(frame);
  }
}


/*!
\brief Enable a sdf given its name.
\param name sdf name
*/
void SphereTracingWidget::EnableSDF(const QString& name)
{
  if (sdfShaders.contains(name))
    sdfShaders[name]->enabled = true;
}

/*!
\brief Disable a sdf given its name. The sdf will not be rendered until it's enabled again.
\param name sdf name
*/
void SphereTracingWidget::DisableSDF(const QString& name)
{
  if (sdfShaders.contains(name))
    sdfShaders[name]->enabled = false;
}

/*!
\brief Destroys all sdf objects in the scene.
*/
void SphereTracingWidget::ClearAll()
{
  makeCurrent();
  for (SDFGL* it : sdfShaders)
  {
    it->Delete();
    delete it;
  }
  sdfShaders.clear();
}


/*!
\brief Renders the scene.
*/
void SphereTracingWidget::paintGL()
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

  // Draw SDF
  RenderSDF();

  // Draw meshes
  RenderMeshes();

  // Draw points
  RenderPoints();

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
\brief Reload the shaders of the widget. Useful for realtime editing and fine tuning of the rendering.
*/
void SphereTracingWidget::ReloadShaders()
{
  //search for a local shader first, then in the ARCHE LIBS.
  std::ifstream in("../../libs/LibMaya/Shaders/point.glsl");
  QString pPath = QString("../../libs");
  if (in.good() == false)
  {
    std::cout << "Fail to find local shader point rendering. Searching for a global shader." << std::endl;
    pPath = System::GetResource("ARCHESLIBDIR");
    if (pPath.isEmpty())
    {
      std::cout << "MeshWidget::ReloadShaders() : variable d'environnement ARCHESLIBDIR non défini" << std::endl;
      std::cin.get();
      exit(-1);
    }
  }
  // Point
  QString fullPath = pPath + QString("/LibMaya/Shaders/point.glsl");
  QByteArray ba = fullPath.toLocal8Bit();
  pointShader.Initialize(ba.constData());

  // Shader/Camera/Profiler
  pPath = System::GetResource("ARCHESLIBDIR");
  if (pPath.isEmpty())
  {
    std::cout << "MeshWidget::ReloadShaders() : variable d'environnement ARCHESLIBDIR non défini" << std::endl;
    std::cin.get();
    exit(-1);
  }

  // Sky
  fullPath = pPath + QString("/LibMaya/Shaders/skybox.glsl");
  ba = fullPath.toLocal8Bit();
  skyShader.Initialize(ba.constData());

  // Mesh
  fullPath = pPath + QString("/LibMaya/Shaders/mesh.glsl");
  ba = fullPath.toLocal8Bit();
  meshShader.Initialize(ba.constData());

  // Box
  fullPath = pPath + QString("/LibMaya/Shaders/boundingBox.glsl");
  ba = fullPath.toLocal8Bit();
  boxShader.Initialize(ba.constData());
}

/*!
\brief Render all the evaluation points.
\author Hubert-Brierre Pierre
*/
void SphereTracingWidget::RenderPoints()
{
  // Camera matrices
  Matrix4Float modelView = Matrix4Float::LookAt(camera.Eye(), camera.At(), camera.Up());
  Matrix4Float projection;
  if (perspectiveProjection)
    projection = Matrix4Float::Perspective(Math::RadianToDegree(camera.GetAngleOfViewV(width(), height())), nearAndFarPlane[0], nearAndFarPlane[1], width(), height());
  else
    projection = Matrix4Float::Orthographic(-cameraOrthoSize, cameraOrthoSize, -cameraOrthoSize, cameraOrthoSize, nearAndFarPlane[0], nearAndFarPlane[1]);


  // Shared uniforms
  pointShader.Bind();
  glEnable(GL_PROGRAM_POINT_SIZE);
  pointShader.SetUniform("ModelViewMatrix", &modelView(0, 0));
  pointShader.SetUniform("ProjectionMatrix", &projection(0, 0));
  pointShader.SetUniform("TRSMatrix", &evaluationsPoints.TRSMatrix[0]);

  pointShader.SetUniform("Width", width());
  pointShader.SetUniform("Height", height());


  // Draw
  glBindVertexArray(evaluationsPoints.vao);
  glDrawArrays(GL_POINTS, 0, (GLsizei)evaluationsPoints.pointCount);

  profiler.drawCallPerFrame++;

  profiler.EndGPU();

}

/*!
\brief Render all the sdf.
\author Hubert-Brierre Pierre
*/
void SphereTracingWidget::RenderSDF()
{
  // Camera matrices
  Matrix4Float modelView = Matrix4Float::LookAt(camera.Eye(), camera.At(), camera.Up());
  Matrix4Float projection;
  if (perspectiveProjection)
    projection = Matrix4Float::Perspective(Math::RadianToDegree(camera.GetAngleOfViewV(width(), height())), nearAndFarPlane[0], nearAndFarPlane[1], width(), height());
  else
    projection = Matrix4Float::Orthographic(-cameraOrthoSize, cameraOrthoSize, -cameraOrthoSize, cameraOrthoSize, nearAndFarPlane[0], nearAndFarPlane[1]);


  // Iterate over all sdf and draw
  for (SDFGL* it : sdfShaders)
  {
    if (!it->enabled)
      continue;

    Matrix4Float TSR = Matrix4Float(it->TRSMatrix);
    Matrix4Float matrix = projection * modelView * TSR;
    Matrix4Float invMatrix = Matrix4Float::inverse(matrix);

    if (matrix != it->old_transformation) {
        it->nb_sample = 0;
        std::cout << "mooved" << std::endl;
    }
    it->old_transformation = matrix;
    if (width() != it->textureWidth || height() != it->textureHeight)
        it->resetColorTexture(width(), height());

    // Shared uniforms
    it->shader.Bind();


    //loading uniforms:
    it->shader.SetUniform("invMatrix", &invMatrix(0, 0));
    it->shader.SetUniform("iResolution", width(), height());
    it->shader.SetUniform("zfarnear", (float)nearAndFarPlane[1], (float)nearAndFarPlane[0]);
    it->shader.SetUniform("cameraPos", camera.Eye());
    it->shader.SetUniform("maxSteps", nbSteps); //Max step to find an intersection. Can be link in a QT structure to be modified in live.
    it->shader.SetUniform("shading", it->shading); //used shading
    it->shader.SetUniform("sampleFieldNormal", it->sampleFieldNormal);
    it->shader.SetUniform("sampleFieldPoint", it->sampleFieldPoint);
    it->shader.SetUniform("sampleFieldFreq", it->sampleFieldFreq);
    it->shader.SetUniform("maxCost", it->maxCost);
    it->shader.SetUniform("kg", it->kg);
    it->shader.SetUniform("kn", it->kn);
    it->shader.SetUniform("aa", aa);
    it->shader.SetUniform("nb_bounce", 4);

    if (it->shading == it->SEXY_SHADING) {
        glBindFramebuffer(GL_FRAMEBUFFER, it->sexyFramebuffer);
        it->shader.SetUniform("light_dir", Normalized(Vector(1., 1., 1.)));
        /*it->shader.SetUniform("light_color", Vector(0.65, 0.6, 0.55) * 1.5*5);
        it->shader.SetUniform("ambiant_color", Vector(0.2, 0.21, 0.22));*/
        it->shader.SetUniform("light_color", Vector(0.65, 0.6, 0.55) * 1.2*2);
        it->shader.SetUniform("ambiant_color", Vector(0.65, 0.6, 0.55) *1.2);
        it->shader.SetUniform("nb_sample", it->nb_sample);
        it->shader.SetUniform("seed", float(it->r.Uniform()));
        
        it->shader.SetUniform("readColor", 0);
        it->shader.SetUniform("texture_size", width(), height());
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, it->sexyColorTexture);

        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT | GL_SHADER_IMAGE_ACCESS_BARRIER_BIT | GL_FRAMEBUFFER_BARRIER_BIT);

        glBindFramebuffer(GL_FRAMEBUFFER, defaultFramebufferObject());
    }

    it->shader.SetUniform("nb_sample", -1);
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);



    it->nb_sample++;
    if (it->shading != it->SEXY_SHADING)
        it->nb_sample = 0;

    profiler.drawCallPerFrame++;
  }
  profiler.EndGPU();

}

/*!
\brief Set all the evaluations points.
\author Hubert-Brierre Pierre
*/
void SphereTracingWidget::setEvaluationPoints(std::vector<Vector> points, std::vector<float> nbRecCalls)
{
  evaluationsPoints = PointGL(points, nbRecCalls);
}


/*!
\brief Set the number of step of the sphere tracing.
\author Hubert-Brierre Pierre
*/
void SphereTracingWidget::setNbSteps(int newNbSteps) {
  nbSteps = newNbSteps;
}

/*!
\brief Set the shading to use.
\author Hubert-Brierre Pierre
*/
void SphereTracingWidget::setShading(int newShading) {
  for (SDFGL* it : sdfShaders)
  {
    it->setShading(newShading);
  }
}

/*!
\brief Set the normal of the plane showing the sdf field.
\author Hubert-Brierre Pierre
*/
void SphereTracingWidget::setSampleFieldNormal(const Vector& sampleFieldNormal) {
  for (SDFGL* it : sdfShaders)
  {
    it->sampleFieldNormal = sampleFieldNormal;
  }
}

/*!
\brief Set the offset of the plane showing the sdf field.
\author Hubert-Brierre Pierre
*/
void SphereTracingWidget::setSampleFieldPoint(const Vector& sampleFieldPoint) {
  for (SDFGL* it : sdfShaders)
  {
    it->sampleFieldPoint = sampleFieldPoint;
  }
}

/*!
\brief Set the frequence of the sdf field.
\author Hubert-Brierre Pierre
*/
void SphereTracingWidget::setSampleFieldFreq(float sampleFieldFreq) {
  for (SDFGL* it : sdfShaders)
  {
    it->sampleFieldFreq = sampleFieldFreq;
  }
}

/*!
\brief Set the maximum cost to show on the cost shader.
\author Hubert-Brierre Pierre
*/
void SphereTracingWidget::setMaxCost(int maxCost)
{
  for (SDFGL* it : sdfShaders)
  {
    it->maxCost = maxCost;
  }
}

/*!
\brief Set the geometric interpolation coeficient.
\author Hubert-Brierre Pierre
*/
void SphereTracingWidget::setKG(float kg)
{
  for (SDFGL* it : sdfShaders)
  {
    it->kg = kg;
  }
}

/*!
\brief Set the normal interpolation coeficient.
\author Hubert-Brierre Pierre
*/
void SphereTracingWidget::setKN(float kn)
{
  for (SDFGL* it : sdfShaders)
  {
    it->kn = kn;
  }
}

/*!
\brief Set the anti-aliasing level.
\param a Level.

\author Hubert-Brierre Pierre
*/
void SphereTracingWidget::setAA(int a)
{
  aa = a;
}
