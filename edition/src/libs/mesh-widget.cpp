#include "libs/realtime.h"
#include "libs/mesh.h"
#include "libs/meshcolor.h"
#include "libs/cpu.h"
#include "libs/plane.h"

#include <QtWidgets/QApplication>
#include <QtGui/QKeyEvent>
#include <QtGui/QMouseEvent>
#include <QtGui/QPainter>

/*!
\brief Default constructor.
*/
MeshWidget::MeshGL::MeshGL()
{
  enabled = true;
  useWireframe = false;
  shading = MeshShading::Triangles;
  material = MeshMaterial::Normal;

  vao = 0;
  fullBuffer = GLBuffer();
  indexBuffer = GLBuffer();
  triangleCount = 0;
  SetFrame(FrameScaled::Id);
}

/*!
\brief Constructor from a Mesh and a frame scaled.
*/
MeshWidget::MeshGL::MeshGL(const Mesh& mesh, const FrameScaled& fr) : MeshGL()
{
  SetFrame(fr);

  // Compute plain arrays of sorted vertices & normals
  QVector<int> vertexIndexes = mesh.VertexIndexes();
  QVector<int> normalIndexes = mesh.NormalIndexes();
  assert(vertexIndexes.size() == normalIndexes.size());

  int nbVertex = vertexIndexes.size();
  int singleBufferSize = nbVertex * 3;
  float* vertices = new float[singleBufferSize];
  float* normals = new float[singleBufferSize];
  for (int i = 0; i < nbVertex; i++)
  {
    int indexVertex = vertexIndexes[i];
    int indexNormal = normalIndexes[i];

    Vector vertex = mesh.Vertex(indexVertex);
    vertices[i * 3 + 0] = float(vertex[0]);
    vertices[i * 3 + 1] = float(vertex[1]);
    vertices[i * 3 + 2] = float(vertex[2]);

    Vector normal = mesh.Normal(indexNormal);
    normals[i * 3 + 0] = float(normal[0]);
    normals[i * 3 + 1] = float(normal[1]);
    normals[i * 3 + 2] = float(normal[2]);
  }
  // Indices are now sorted
  int* indices = new int[nbVertex];
  for (int i = 0; i < nbVertex; i++)
    indices[i] = i;
  triangleCount = nbVertex;

  // Generate vao & buffers
  if (vao == 0)
    glGenVertexArrays(1, &vao);
  fullBuffer.Generate();
  indexBuffer.Generate();

  glBindVertexArray(vao);
  size_t fullSize = sizeof(float) * singleBufferSize + sizeof(float) * singleBufferSize;
  fullBuffer.SetData(GL_ARRAY_BUFFER, fullSize, nullptr);

  // Vertices(0)
  size_t size = 0;
  size_t offset = 0;
  size = sizeof(float) * singleBufferSize;
  fullBuffer.SetSubData(GL_ARRAY_BUFFER, offset, size, vertices);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (const void*)offset);
  glEnableVertexAttribArray(0);

  // Normals(1)
  offset = offset + size;
  size = sizeof(float) * singleBufferSize;
  fullBuffer.SetSubData(GL_ARRAY_BUFFER, offset, size, normals);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (const void*)offset);
  glEnableVertexAttribArray(1);

  // Triangles
  indexBuffer.SetData(GL_ELEMENT_ARRAY_BUFFER, sizeof(int) * nbVertex, indices);

  // Free data
  delete[] vertices;
  delete[] normals;
  delete[] indices;
}

/*!
\brief Constructor from a MeshColor and a frame scaled.
*/
MeshWidget::MeshGL::MeshGL(const MeshColor& mesh, const FrameScaled& fr) : MeshGL()
{
  SetFrame(fr);

  // Compute plain arrays of sorted vertices & normals

  /*
  QVector<int> vertexIndexes = mesh.VertexIndexes();
   QVector<int> normalIndexes = mesh.NormalIndexes();
   QVector<int> colorIndexes = mesh.ColorIndexes();
    assert(vertexIndexes.size() == normalIndexes.size());
    */
  assert(mesh.Vertexes() == mesh.Normals());

  int nbVertex = mesh.Vertexes();
  //	int nbVertex = vertexIndexes.size();
  int singleBufferSize = nbVertex * 3;
  int colorBufferSize = nbVertex * 4;

  float* vertices = new float[singleBufferSize];
  float* normals = new float[singleBufferSize];
  float* colors = new float[colorBufferSize];

  for (int i = 0; i < nbVertex; i++)
  {
    /*
    int indexVertex = vertexIndexes[i];
    int indexNormal = normalIndexes[i];
    int indexColor = colorIndexes[i];
    */

    int indexVertex = mesh.VertexIndex(i);
    int indexNormal = mesh.NormalIndex(i);
    int indexColor = mesh.ColorIndex(i);

    Vector vertex = mesh.Vertex(indexVertex);
    vertices[i * 3 + 0] = float(vertex[0]);
    vertices[i * 3 + 1] = float(vertex[1]);
    vertices[i * 3 + 2] = float(vertex[2]);

    Vector normal = mesh.Normal(indexNormal);
    normals[i * 3 + 0] = float(normal[0]);
    normals[i * 3 + 1] = float(normal[1]);
    normals[i * 3 + 2] = float(normal[2]);

    Color color = mesh.GetColor(indexColor);
    colors[i * 4 + 0] = float(color[0]);
    colors[i * 4 + 1] = float(color[1]);
    colors[i * 4 + 2] = float(color[2]);
    colors[i * 4 + 3] = float(color[3]);
  }
  // Indices are now sorted
  int* indices = new int[nbVertex];
  for (int i = 0; i < nbVertex; i++)
    indices[i] = i;
  triangleCount = nbVertex;

  // Generate vao & buffers
  if (vao == 0)
    glGenVertexArrays(1, &vao);
  fullBuffer.Generate();
  indexBuffer.Generate();

  glBindVertexArray(vao);
  size_t fullSize =
    sizeof(float) * singleBufferSize	// Vertices
    + sizeof(float) * singleBufferSize	// Normals
    + sizeof(float) * colorBufferSize;	// Colors
  fullBuffer.SetData(GL_ARRAY_BUFFER, fullSize, nullptr);

  // Vertices(0)
  size_t size = 0;
  size_t offset = 0;
  size = sizeof(float) * singleBufferSize;
  fullBuffer.SetSubData(GL_ARRAY_BUFFER, offset, size, vertices);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (const void*)offset);
  glEnableVertexAttribArray(0);

  // Normals(1)
  offset = offset + size;
  size = sizeof(float) * singleBufferSize;
  fullBuffer.SetSubData(GL_ARRAY_BUFFER, offset, size, normals);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (const void*)offset);
  glEnableVertexAttribArray(1);

  // Colors(2)
  offset = offset + size;
  size = sizeof(float) * colorBufferSize;
  fullBuffer.SetSubData(GL_ARRAY_BUFFER, offset, size, colors);
  glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, 0, (const void*)offset);
  glEnableVertexAttribArray(2);

  // Triangles
  indexBuffer.SetData(GL_ELEMENT_ARRAY_BUFFER, sizeof(int) * nbVertex, indices);

  // Free data
  delete[] vertices;
  delete[] normals;
  delete[] colors;
  delete[] indices;
}

/*!
\brief Constructor from a bounding box. Creates an indexed mesh suited for wireframe rendering only.
\param box bounding box
*/
MeshWidget::MeshGL::MeshGL(const Box& box, const FrameScaled& fr) : MeshGL()
{
  // Generate vao & buffers
  if (vao == 0)
    glGenVertexArrays(1, &vao);
  fullBuffer.Generate();
  indexBuffer.Generate();
  glBindVertexArray(vao);

  // Vertices
  std::vector<float> vertices;
  vertices.resize(3 * 8);
  for (int i = 0; i < 8; i++)
  {
    Vector v = box.Vertex(i);
    vertices[(i * 3) + 0] = float(v[0]);
    vertices[(i * 3) + 1] = float(v[1]);
    vertices[(i * 3) + 2] = float(v[2]);
  }
  // Vertices (0)
  fullBuffer.SetData(GL_ARRAY_BUFFER, sizeof(float) * vertices.size(), vertices.data());
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
  glEnableVertexAttribArray(0);

  // Triangles (indexed just as LibCore box vertex indices)
  int elements[] = {
    0, 1, 3, 2,
    4, 5, 7, 6,
    0, 4, 5, 1,
    2, 6, 7, 3
  };
  indexBuffer.SetData(GL_ELEMENT_ARRAY_BUFFER, sizeof(elements), elements);

  // Apply frame to cube
  SetFrame(fr);

  // Bounding box mesh is disabled by default
  enabled = false;
}

/*!
\brief Delete all opengl buffers and vao.
*/
void MeshWidget::MeshGL::Delete()
{
  glDeleteVertexArrays(1, &vao);
  fullBuffer.Destroy();
  indexBuffer.Destroy();
}

/*!
\brief Set the frame for the mesh.
\param fr new frame
*/
void MeshWidget::MeshGL::SetFrame(const FrameScaled& fr)
{
  fr.GetMatrix4().Float(TRSMatrix);
}


/*!
\brief Default constructor.
*/
MeshWidget::MeshWidget(QWidget* parent) : QOpenGLWidget(parent)
{
  perspectiveProjection = true;
  cameraOrthoSize = 1000.0f;
  backgroundColor = Color::White;

  setMouseTracking(true);
  setFocusPolicy(Qt::StrongFocus);
}

/*!
\brief Destructor.
*/
MeshWidget::~MeshWidget()
{
  // Destroy all meshes
  ClearAll();

  // Release shader
  meshShader.Destroy();
  boxShader.Destroy();
  skyShader.Destroy();
}

/*!
\brief Initialize OpenGL, shaders and a camera centered at origin.
*/
void MeshWidget::initializeGL()
{
  GLenum err = glewInit();
  if (err != GLEW_OK)
  {
    std::cout << "GLEW Error: " << glewGetErrorString(err) << std::endl;
    std::cin.get();
    exit(-1);
  }

  glEnable(GL_LINE_SMOOTH);
  glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);

  glEnable(GL_DEPTH_TEST);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  // Disable vsync for proper performance measurements
  // QT6
  //QGLFormat fmt;
  //fmt.setSwapInterval(0);
  //QGLFormat::setDefaultFormat(fmt);

  // Shader/Camera/Profiler
  ReloadShaders();

  glGenVertexArrays(1, &skyboxVAO);

  camera = Camera::View(Box(10.0));
  SetNearAndFarPlane(nearAndFarPlane);

  profiler.Init();
}

/*!
\brief Renders the scene.
*/
void MeshWidget::paintGL()
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

  // Sky
  glEnable(GL_DEPTH_TEST);
  if (useSkyShader)
  {
    RenderSky();
    profiler.drawCallPerFrame++;
    glBindVertexArray(0);
    glUseProgram(0);
  }

  // Draw meshes
  profiler.BeginGPU();
  RenderMeshes();

  // CPU Profiling
  if (profiler.enabled)
    profiler.Update();

  RenderUiPanels();

  // Schedule next draw
  update();
}

/*!
\brief Resize window.
\param w, h Width and height.
*/
void MeshWidget::resizeGL(int w, int h)
{
  glViewport(0, 0, (GLint)w, (GLint)h);
}

/*!
\brief Render the sky.
*/
void MeshWidget::RenderSky()
{
  glDepthFunc(GL_LEQUAL);
  skyShader.Bind();
  glBindVertexArray(skyboxVAO);
  skyShader.SetUniform(0, camera.Eye());
  skyShader.SetUniform(1, camera.At());
  skyShader.SetUniform(2, camera.Up());
  skyShader.SetUniform(3, float(width()), float(height()));
  glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
}

/*!
\brief Render all the meshes.
*/
void MeshWidget::RenderMeshes()
{
  // Camera matrices
  Matrix4Float modelView = Matrix4Float::LookAt(camera.Eye(), camera.At(), camera.Up());
  Matrix4Float projection;
  if (perspectiveProjection)
    projection = Matrix4Float::Perspective(Math::RadianToDegree(camera.GetAngleOfViewV(width(), height())), nearAndFarPlane[0], nearAndFarPlane[1], width(), height());
  else
    projection = Matrix4Float::Orthographic(-cameraOrthoSize, cameraOrthoSize, -cameraOrthoSize, cameraOrthoSize, nearAndFarPlane[0], nearAndFarPlane[1]);

  // Shared uniforms
  meshShader.Bind();
  meshShader.SetUniform("ModelViewMatrix", &modelView(0, 0));
  meshShader.SetUniform("ProjectionMatrix", &projection(0, 0));
  meshShader.SetUniform("WIN_SCALE", width() / 2.0f, height() / 2.0f);
  meshShader.SetUniform("viewDir", Normalized(camera.At() - camera.Eye()));

  // Iterate over all meshes and draw
  for (MeshGL* it : objects)
  {
    if (!it->enabled)
      continue;
    meshShader.SetUniform("TRSMatrix", &it->TRSMatrix[0]);
    meshShader.SetUniform("useWireframe", it->useWireframe ? 1 : 0);
    meshShader.SetUniform("material", (int)it->material);
    meshShader.SetUniform("shading", (int)it->shading);

    // Draw
    glBindVertexArray(it->vao);
    glDrawElements(GL_TRIANGLES, (GLsizei)it->triangleCount, GL_UNSIGNED_INT, 0);

    profiler.drawCallPerFrame++;
  }
  profiler.EndGPU();

  // Draw bounding box
  boxShader.Bind();
  boxShader.SetUniform("ModelViewMatrix", &modelView(0, 0));
  boxShader.SetUniform("ProjectionMatrix", &projection(0, 0));
  glLineWidth(2.0f);
  for (MeshGL* it : boxObjects)
  {
    if (!it->enabled)
      continue;
    boxShader.SetUniform("TRSMatrix", &it->TRSMatrix[0]);

    glBindVertexArray(it->vao);
    glDrawElements(GL_LINE_LOOP, 4, GL_UNSIGNED_INT, 0);
    glDrawElements(GL_LINE_LOOP, 4, GL_UNSIGNED_INT, (GLvoid*)(4 * sizeof(int)));
    glDrawElements(GL_LINES, 8, GL_UNSIGNED_INT, (GLvoid*)(8 * sizeof(int)));

    profiler.drawCallPerFrame++;
  }
  glLineWidth(1.0f);
}

/*!
\brief Render the UI panels of the widget.
*/
void MeshWidget::RenderUiPanels()
{
  // We need to unbind VAO and program for now because it causes problem with below command.
  glBindVertexArray(0);
  glUseProgram(0);
  glBindVertexArray(0);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  QPainter painter(this);
  painter.setRenderHint(QPainter::Antialiasing);
  QPen penLineGrey(QColor(50, 50, 50));
  QPen penLineWhite(QColor(250, 250, 250));

  // Stats panel
  int offsetY = 0;
  if (profiler.enabled)
  {
    const int bX = 10;
    const int bY = 10;
    const int sizeX = 200;
    const int sizeY = 80;

    offsetY += sizeY + 10;

    // Background
    painter.setPen(penLineGrey);
    painter.fillRect(QRect(bX, bY, sizeX, sizeY), QColor(0, 0, 255, 25));
    painter.drawRect(bX, bY, sizeX, sizeY);

    // Text
    QFont f1, f2;
    f1.setBold(true);
    painter.setFont(f1);
    painter.setPen(penLineWhite);
    painter.drawText(10 + 5, bY + 10 + 5, "Statistics");
    painter.setFont(f2);
    painter.drawText(10 + 5, bY + 10 + 20, "CPU FPS:\t" + QString::number(profiler.framePerSecond));
    painter.drawText(10 + 5, bY + 10 + 35, "CPU Frame:\t" + QString::number(profiler.msPerFrame) + "ms");
    painter.drawText(10 + 5, bY + 10 + 50, "GPU:\t" + QString::number(profiler.elapsedTimeGPU / 1000000.0) + "ms");
    painter.drawText(10 + 5, bY + 10 + 65, "Draw calls:\t" + QString::number(profiler.drawCallPerFrame));
  }

  // Camera panel
  if (renderCameraPanel)
  {
    const int bX = 10;
    const int bY = offsetY + 10;
    const int sizeX = 400;
    const int sizeY = 30;

    painter.setPen(penLineGrey);
    painter.fillRect(QRect(bX, (bY), sizeX, sizeY), QColor(0, 0, 255, 25));
    painter.drawRect(bX, (bY), sizeX, sizeY);

    painter.setPen(penLineWhite);
    painter.drawText(bX + 10, bY + 12 + 1 * 15, camera.ToString());

    QFont f1; f1.setBold(true);
    painter.setFont(f1);
    painter.drawText(bX + bX, bY + 12, QString("Camera Position"));
  }

  painter.end();

  // Reset GL depth test
  glEnable(GL_DEPTH_TEST);
}

/*!
\brief Reload the shaders of the widget. Useful for realtime editing and fine tuning of the rendering.
*/
void MeshWidget::ReloadShaders()
{
  // Shader/Camera/Profiler
  QString pPath = System::GetResource("ARCHESLIBDIR");
  if (pPath.isEmpty())
  {
    std::cout << "MeshWidget::ReloadShaders() : variable d'environnement ARCHESLIBDIR non défini" << std::endl;
    std::cin.get();
    exit(-1);
  }

  // Sky
  QString fullPath = pPath + QString("/LibMaya/Shaders/skybox.glsl");
  QByteArray ba = fullPath.toLocal8Bit();
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
\brief Add a new mesh in the scene.
\param mesh new mesh
\param frame mesh frame, identity by default.
*/
void MeshWidget::AddMesh(const QString& name, const Mesh& mesh, const FrameScaled& frame)
{
  makeCurrent();
  objects.insert(name, new MeshGL(mesh, frame));
  boxObjects.insert(name, new MeshGL(mesh.GetBox(), frame));
}

/*!
\brief Add a new colored mesh in the scene.
\param mesh new colored mesh
\param frame mesh frame, identity by default.
*/
void MeshWidget::AddMesh(const QString& name, const MeshColor& mesh, const FrameScaled& frame)
{
  makeCurrent();
  objects.insert(name, new MeshGL(mesh, frame));
  boxObjects.insert(name, new MeshGL(mesh.GetBox(), frame));
}

/*!
\brief Check if a mesh exists in the widget.
\param name mesh name
*/
bool MeshWidget::HasMesh(const QString& name) const
{
  return objects.contains(name);
}

/*!
\brief Delete a mesh in the scene from its name.
\param name mesh name
*/
void MeshWidget::DeleteMesh(const QString& name)
{
  makeCurrent();
  if (objects.contains(name))
  {
    objects[name]->Delete();
    boxObjects[name]->Delete();
    objects.remove(name);
    boxObjects.remove(name);
  }
}

/*!
\brief Updates the transform of a mesh given its name.
\param name mesh name
\param frame new frame
*/
void MeshWidget::UpdateMesh(const QString& name, const FrameScaled& frame)
{
  makeCurrent();
  if (objects.contains(name))
  {
    objects[name]->SetFrame(frame);
    boxObjects[name]->SetFrame(frame);
  }
}

/*!
\brief Set the sky rendering mode.
\param useShader flag for using the sky shader or not.
\param clearColor background color used in case the shader is not enabled.
*/
void MeshWidget::SetSky(bool useShader, const Color& clearColor)
{
  useSkyShader = useShader;
  backgroundColor = clearColor;
}

/*!
\brief Enable a mesh given its name.
\param name mesh name
*/
void MeshWidget::EnableMesh(const QString& name)
{
  if (objects.contains(name))
    objects[name]->enabled = true;
}

/*!
\brief Disable a mesh given its name. The mesh will not be rendered until it's enabled again.
\param name mesh name
*/
void MeshWidget::DisableMesh(const QString& name)
{
  if (objects.contains(name))
    objects[name]->enabled = false;
}

/*!
\brief Set the bounding box rendering flag for a given mesh.
\param name mesh name
\param use render flag
*/
void MeshWidget::UseBoundingBox(const QString& name, bool use)
{
  if (boxObjects.contains(name))
    boxObjects[name]->enabled = use;
}

/*!
\brief Destroys all mesh objects in the scene.
*/
void MeshWidget::ClearAll()
{
  makeCurrent();
  for (MeshGL* it : objects)
  {
    it->Delete();
    delete it;
  }
  objects.clear();

  for (MeshGL* it : boxObjects)
  {
    it->Delete();
    delete it;
  }
  boxObjects.clear();
}

/*!
\brief Compute a ray starting from the camera position at a given pixel.
\param pix the pixel
*/
Ray MeshWidget::ConvertPixelToRay(const QPoint& pix) const
{
  return camera.PixelToRay(pix.x(), pix.y(), width(), height());
}

/*!
\brief Capture the rendering viewport and save it to disk (as png and jpg).
\param w, h Width and height of the image.
*/
void MeshWidget::SaveScreen(int w, int h)
{
  QSize qs = size();
  resize(w, h);
  QImage image(size(), QImage::Format_ARGB32);

  // Set name according to date and time
  QString name = "../screen-" + System::DateTime();

  // Get image
  image = grabFramebuffer();

  // Save
  image.save(name + ".png", "PNG");
  image.save(name + ".jpg", "JPG", 96);

  // Restore buffer size
  resize(qs);
}

/*!
\brief Capture the rendering viewport and save it to disk.
\param w, h Width and height of the image.
\param url Name of the exported file
*/
void MeshWidget::SaveScreen(int w, int h, const QString& url)
{
  QSize qs = size(); // Save buffersize

  resize(w, h);
  QImage image(size(), QImage::Format_ARGB32);

  // Get image
  image = grabFramebuffer();

  // Save
  image.save(url, "PNG");

  // Restore buffer size
  resize(qs);
}

/*!
\brief Set the current camera.
\param cam new camera
*/
void MeshWidget::SetCamera(const Camera& cam)
{
  camera = cam;
  camera.SetPlanes(nearAndFarPlane[0], nearAndFarPlane[1]);
}

/*!
\brief Set the near and far plane of the camera.
\param planes near and far planes
*/
void MeshWidget::SetNearAndFarPlane(const Vector2& planes)
{
  nearAndFarPlane = planes;
  camera.SetPlanes(nearAndFarPlane[0], nearAndFarPlane[1]);
}

/*!
\brief Returns the current camera.
*/
Camera MeshWidget::GetCamera() const
{
  return camera;
}

/*!
\brief Capture and return the rendering viewport.
\param w, h dimensions of the capture
*/
QImage MeshWidget::GrabScreen(int w, int h)
{
  QSize qs = size();
  resize(w, h);
  QImage image(size(), QImage::Format_ARGB32);

  // Set name according to date and time
  QString name = "../screen-" + System::DateTime();

  // Get image
  image = grabFramebuffer();

  resize(qs);

  return image;
}


/*!
\brief Set the camera mode (perspective or orthographic).
\param perspective
*/
void MeshWidget::SetCameraMode(bool perspective)
{
  perspectiveProjection = perspective;
}

/*!
\brief Changes the material for a mesh given its name.
\param name mesh name
\param mat the new material
*/
void MeshWidget::SetMaterial(const QString& name, MeshMaterial mat)
{
  if (objects.contains(name))
    objects[name]->material = mat;
}

/*!
\brief Changes the material for all meshes.
\param mat the new material
*/
void MeshWidget::SetMaterialGlobal(MeshMaterial mat)
{
  for (MeshIterator i = objects.begin(); i != objects.end(); i++)
    i.value()->material = mat;
}

/*!
\brief Changes the wireframe render flag for a mesh given its name.
\param name mesh name
\param wireframe new wireframe flag value
*/
void MeshWidget::UseWireframe(const QString& name, bool wireframe)
{
  if (objects.contains(name))
    objects[name]->useWireframe = wireframe;
}

/*!
\brief Changes the wireframe render flag for all meshes.
\param u new value for wireframe flag.
*/
void MeshWidget::UseWireframeGlobal(bool wireframe)
{
  for (MeshIterator i = objects.begin(); i != objects.end(); i++)
    i.value()->useWireframe = wireframe;
}

/*!
\brief Changes the shading mode for a given mesh.
\param shading new shading mode
*/
void MeshWidget::SetShading(const QString& name, MeshShading shading)
{
  if (objects.contains(name))
    objects[name]->shading = shading;
}

/*!
\brief Changes the shading mode for all meshes.
\param shading new shading mode
*/
void MeshWidget::SetShadingGlobal(MeshShading shading)
{
  for (MeshIterator i = objects.begin(); i != objects.end(); i++)
    i.value()->shading = shading;
}

/*!
\brief Returns the current mouse position.
*/
QPoint MeshWidget::GetMousePosition() const
{
  return QPoint(x0, y0);
}

/*!
\brief Process mouse click events.
\param e Events.
*/
void MeshWidget::mousePressEvent(QMouseEvent* e)
{
  x0 = e->globalPosition().x();
  y0 = e->globalPosition().y();
  if (e->modifiers() & Qt::ControlModifier)
  {
    if (e->buttons() == Qt::LeftButton)
      emit _signalEditSceneLeft(ConvertPixelToRay(e->pos()));
    else if (e->buttons() == Qt::RightButton)
      emit _signalEditSceneRight(ConvertPixelToRay(e->pos()));
    emit _signalMouseMove();
  }
}

/*!
\brief Process the mouse release events.
\param e Events.
*/
void MeshWidget::mouseReleaseEvent(QMouseEvent*)
{
  QApplication::setOverrideCursor(QCursor(Qt::ArrowCursor));
  emit _signalMouseRelease();
  update();
}

/*!
\brief Process the mouse double click events.
\param e Events.
*/
void MeshWidget::mouseDoubleClickEvent(QMouseEvent* e)
{
  double t;
  QPoint pixel = e->pos();
  Ray ray = camera.PixelToRay(pixel.x(), pixel.y() - 1, width(), height());
  MoveAt = false;

  // Intersection with plane
  Plane plane(Vector(0., 0., 1.0), 1.);
  if (plane.Intersect(ray, t))
  {
    MoveAt = true;
    currentAt = camera.At();
    toAt = ray(t);
    stepAt = 0;
  }
  emit _signalMouseMove();
  update();
}

/*!
\brief Mouse move event with camera movements.
*/
void MeshWidget::mouseMoveEvent(QMouseEvent* e)
{
  int x = e->globalPosition().x();
  int y = e->globalPosition().y();
  if ((e->modifiers() & Qt::AltModifier))
  {
    // Displacement mode 
    double MoveScale = Norm(camera.View()) * 0.015 * 0.05;
    if (e->buttons() & Qt::LeftButton)
    {
      // Alt + Left Mouse Move    : Rotation 
      camera.LeftRightRound((x0 - x) * 0.01);
      camera.UpDownRound((y0 - y) * 0.005);
    }
    else if (e->buttons() & Qt::RightButton)
    {
      // Alt + Right Mouse Move   : Forward and Backward
      camera.BackForth((y - y0) * MoveScale);
      QApplication::setOverrideCursor(QCursor(Qt::SplitVCursor));
    }
    else if (e->buttons() & Qt::MiddleButton)
    {
      // Alt + Left Mouse Move    : Plan displacement
      camera.LeftRightHorizontal((x - x0) * MoveScale);
      camera.UpDownVertical((y - y0) * MoveScale);
      QApplication::setOverrideCursor(QCursor(Qt::SizeAllCursor));
    }
    x0 = e->globalPosition().x();
    y0 = e->globalPosition().y();

    emit _signalMouseMove();
  }
  if (e->modifiers() & Qt::ControlModifier)
  {
    emit _signalMouseMoveEdit(ConvertPixelToRay(e->pos()));
  }
}

/*!
\brief Mouse wheel event with zoom.
*/
void MeshWidget::wheelEvent(QWheelEvent* e)
{
  if (!(e->modifiers() & Qt::ControlModifier) && !(e->modifiers() & Qt::ShiftModifier))
  {
    double MoveScale = Norm(camera.View()) * 0.025;
    if (perspectiveProjection)
    {
      if (e->angleDelta().y() > 0)
        camera.BackForth(MoveScale);
      else
        camera.BackForth(-MoveScale);
    }
    else
    {
      if (e->angleDelta().y() > 0)
        cameraOrthoSize -= MoveScale * 0.25;
      else
        cameraOrthoSize += MoveScale * 0.25;
    }
    emit _signalMouseMove();
  }
  update();
}

/*!
\brief Key press event.
*/
void MeshWidget::keyPressEvent(QKeyEvent* e)
{
  switch (e->key())
  {
    // F1: 720p screenshot in the app folder
  case Qt::Key_F1:
    SaveScreen();
    break;
    // F2: 4k screen in the app folder
  case Qt::Key_F2:
    SaveScreen(3840, 2160);
    break;
    // F3: hot reload shader
  case Qt::Key_F3:
    ReloadShaders();
    break;
  case Qt::Key_S:
    // Alt + S: Statistics
    if (e->modifiers() & Qt::AltModifier)
      profiler.enabled = !profiler.enabled;
    break;
  case Qt::Key_C:
    // Alt + C: Camera
    if (e->modifiers() & Qt::AltModifier)
      renderCameraPanel = !renderCameraPanel;
    break;
  case Qt::Key_B:
    // Alt + B: Bounding box
    if (e->modifiers() & Qt::AltModifier)
    {
      static bool useBboxGlobal = false;
      useBboxGlobal = !useBboxGlobal;
      for (MeshGL* it : boxObjects)
        it->enabled = useBboxGlobal;
    }
    break;
  default:
    QOpenGLWidget::keyPressEvent(e);
  }
}

/*!
\brief Key release event.
*/
void MeshWidget::keyReleaseEvent(QKeyEvent* e)
{
  QOpenGLWidget::keyReleaseEvent(e);
  update();
}
