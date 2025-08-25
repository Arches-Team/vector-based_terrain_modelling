// Maya
#include "libs/maya.h"
#include "libs/cpu.h"

/*!
\class MayaWidget maya.h
\brief A core widget for displaying objects.

This class implements a widget for displaying geometric objects.

\ingroup MayaCore
*/

/*!
\brief Create a widget for real time rendering.
\param parent Parent widget.
*/
MayaWidget::MayaWidget(QWidget* parent) :QOpenGLWidget(parent)
{
  // Interface
  setMouseTracking(true);
  setFocusPolicy(Qt::StrongFocus);

  // Init Camera Position
  camera = Camera(Vector(-10.0, -10.0, 10.0), Vector::Null);

  // Light
  light = Vector(10.0, 10.0, 10.0);

  usePlane = true;
  useWorld = true;
  useBackG = true;
  useAncre = false;
  useLight = false;
  useBBox = false;
  useStats = false;
  useVAxis = false;
  useTime = false;
  useCamera = false;
  useCamRatio = 0;
  gpuall = MayaGpuAll();

  gpubackground = nullptr;

  QSurfaceFormat format = QSurfaceFormat::defaultFormat();
  format.setDepthBufferSize(32);
  setFormat(format);
}

/*!
\brief Release allocated resources.
*/
MayaWidget::~MayaWidget()
{
  if (gpubackground != nullptr)
  {
    glDeleteVertexArrays(1, &backgroundVAO);
    delete gpubackground;
    gpubackground = nullptr;
  }
}

/*!
\brief Set up the OpenGL view port and matrix mode according to the size of the window.
\param w, h Width and height of the window.
*/
void MayaWidget::resizeGL(int w, int h)
{
  glViewport(0, 0, (GLint)w, (GLint)h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  gluPerspective(Math::RadianToDegree(camera.GetAngleOfViewV(w, h)), (GLdouble)w / (GLdouble)h, camera.GetNear(), camera.GetFar());
}

/*!
\brief Renders the scene.

The
s for drawing the scene are performed here.
*/
void MayaWidget::paintGL()
{
  QElapsedTimer t;
  t.start();

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  gluLookAt(camera.Eye()[0], camera.Eye()[1], camera.Eye()[2], camera.At()[0], camera.At()[1], camera.At()[2], camera.Up()[0], camera.Up()[1], camera.Up()[2]);

  updateGPU();
  if (useBackG) RenderBackGround();
  if (useLight) RenderLight();
  if (useAncre) RenderAnchor();
  if (useWorld) gpuall.Render(gpuparam);
  if (useBBox) gpuall.RenderBBox();
  if (usePlane) plane.Render();
  RenderForeGround(t);
}

/*!
\brief Change the parameter of the plane.

It updates the MayaWidget view.

\param box Rectangle.
\param h Height of the plane.
\param l Space between lines.
*/
void MayaWidget::SetPlane(const Box2& box, const double& h, const double& l)
{
  plane = MayaPlane(box, h, l);
  update();
}

/*!
\brief Set the near and far planes.

This function modifies the current camera. It updates the MayaWidget view.

\param n, f Near and far planes distance.
*/
void MayaWidget::SetPlanes(const double& n, const double& f)
{
  camera.SetPlanes(n, f);
  update();
}

/*!
\brief Change the light position.
\param p Position.
*/
void MayaWidget::SetLight(const Vector& p)
{
  light = p;
}

/*!
\brief Set the geometry instances to be rendered.
\param mga Set of geometry instances.
*/
void MayaWidget::SetWorld(const MayaGeometryAll& mga)
{
  gpuall.ClearFrames();
  gpuall.Append(mga);
  useWorld = true;
  update();
}

/*!
\brief Append objects or instances to be rendered.
\param mga Set of objects.
*/
void MayaWidget::AppendWorld(const MayaGeometryAll& mga)
{
  gpuall.Append(mga);
  update();
}

/*!
\brief Replace instances to be rendered.
\param mga Set of geometry instances.
*/
void MayaWidget::ReplaceInWorld(const MayaGeometryAll& mga)
{
  gpuall.Replace(mga);
  update();
}

/*!
\brief Replace material of all instances.
\param mga Set of geometry instances.
*/
void MayaWidget::ReplaceMaterialInWorld(MayaGeometryAll& mga)
{
  gpuall.UpdateMaterial(mga);
  update();
}

/*!
\brief Clear all objects in the scene.
*/
void MayaWidget::ClearWorld()
{
  gpuall.Clear();
  useWorld = true;
  update();
}

/*!
\brief Render the focus point.
*/
void MayaWidget::RenderAnchor()
{
  anchor.ClearFrames();
  anchor.Append(FrameScaled::Translation(ancre));
  anchor.RefreshFrames();
  anchor.Render(gpuparam);
}

/*!
\brief Render the Light position.
*/
void MayaWidget::RenderLight()
{
  glPushMatrix();
  glTranslated(light[0], light[1], light[2]);
  lightshape.Render(gpuparam);
  glPopMatrix();
}

/*!
\brief Render the Background.
*/
void MayaWidget::RenderBackGround()
{
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);

  //Bind program and VAO
  glUseProgram(gpubackground->GetProgram());
  glBindVertexArray(backgroundVAO);

  //Update camera data
  glUniform3f(0, camera.Eye()[0], camera.Eye()[1], camera.Eye()[2]);
  glUniform3f(1, camera.At()[0], camera.At()[1], camera.At()[2]);
  glUniform3f(2, camera.Up()[0], camera.Up()[1], camera.Up()[2]);
  glUniform2f(3, width(), height());

  glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
}

/*!
\brief Render the foreground.
\param t Time.
*/
void MayaWidget::RenderForeGround(QElapsedTimer t)
{
  if (useStats || useTime || useCamera || useCamRatio != 0)
  {
    //We need to unbind VAO and program for now because it causes problem with below command.
    glBindVertexArray(0);
    glUseProgram(0);
    glBindVertexArray(0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    GLint viewport[4];

    // Init Matrix
    glGetIntegerv(GL_VIEWPORT, viewport);

    QPainter painter;

    painter.begin(this);
    painter.setRenderHint(QPainter::Antialiasing);

    QPen penLineGrey(QColor(50, 50, 50));
    QPen penLineWhite(QColor(250, 250, 250));

    if (useStats)
    {
      int bX = 10;
      int bY = 10;
      int sizeX = 200;
      int sizeY = 85;
      MayaStatistics ms = gpuall.GetStatistics();

      painter.setPen(penLineGrey);
      painter.fillRect(QRect(bX, bY, sizeX, sizeY), QColor(0, 0, 255, 25));
      painter.drawRect(bX, bY, sizeX, sizeY);

      painter.setPen(penLineWhite);

      painter.drawText(10 + 1 * 70, bY + 12 + 2 * 15, System::LongInteger(ms.GetObjects()));
      painter.drawText(10 + 2 * 70, bY + 12 + 2 * 15, System::LongInteger(ms.GetInstances()));

      painter.drawText(10 + 1 * 70, bY + 12 + 3 * 15, System::LongInteger(ms.GetObjectsVertices()));
      painter.drawText(10 + 2 * 70, bY + 12 + 3 * 15, System::LongInteger(ms.GetInstancesVertices()));

      painter.drawText(10 + 1 * 70, bY + 12 + 4 * 15, System::LongInteger(ms.GetObjectsTriangles()));
      painter.drawText(10 + 2 * 70, bY + 12 + 4 * 15, System::LongInteger(ms.GetInstancesTriangles()));

      QFont f1; f1.setBold(true);
      painter.setFont(f1);
      painter.drawText(bX + 5, bY + 12, QString("Statistics"));
      painter.drawText(10 + 1 * 70, bY + 12 + 1 * 15, QString("Objects"));
      painter.drawText(10 + 2 * 70, bY + 12 + 1 * 15, QString("Instances"));
      painter.drawText(bX + 5, bY + 12 + 2 * 15, QString("- Number "));
      painter.drawText(bX + 5, bY + 12 + 3 * 15, QString("- Vertices "));
      painter.drawText(bX + 5, bY + 12 + 4 * 15, QString("- Triangles "));
    }

    if (useTime)
    {
      int bX = 10;
      int bY = 10;
      int sizeX = 110;
      int sizeY = 55;

      int dX = viewport[2] - sizeX - bX;

      painter.setPen(penLineGrey);
      painter.fillRect(QRect(dX, (bY), sizeX, sizeY), QColor(0, 0, 255, 25));
      painter.drawRect(dX, (bY), sizeX, sizeY);

      painter.setPen(penLineWhite);
      painter.drawText(dX + 10, bY + 12 + 1 * 15, QString().setNum(t.elapsed()).append(" ms/frames"));
      painter.drawText(dX + 10, bY + 12 + 2 * 15, QString().setNum(1000.0 / t.elapsed()).append(" frames/s"));

      QFont f1; f1.setBold(true);
      painter.setFont(f1);
      painter.drawText(dX + bX, bY + 12, QString("Times"));
    }

    if (useCamRatio != 0)
    {
      QPen penDashLine(QColor(50, 50, 50));
      penDashLine.setStyle(Qt::DashLine);

      double s_ratio = 9. / 16.; // Format 16/9
      if (useCamRatio == 2) s_ratio = 3. / 4.;  // Format 4/3

      QPoint p1(0 + 5, (viewport[3] - viewport[2] * s_ratio) / 2.);
      QPoint p2(0 + 5, (viewport[3] + viewport[2] * s_ratio) / 2.);
      QPoint p3(viewport[2] - 5, (viewport[3] - viewport[2] * s_ratio) / 2.);
      QPoint p4(viewport[2] - 5, (viewport[3] + viewport[2] * s_ratio) / 2.);

      painter.setPen(penDashLine);
      painter.drawLine(p1, p2);
      painter.drawLine(p1, p3);
      painter.drawLine(p2, p4);
      painter.drawLine(p3, p4);
    }

    if (useCamera)
    {
      int bX = 10;
      int bY = 10;
      int sizeX = 400;
      int sizeY = 30;

      int dX = viewport[2] - sizeX - bX;

      painter.setPen(penLineGrey);
      painter.fillRect(QRect(dX, (bY), sizeX, sizeY), QColor(0, 0, 255, 25));
      painter.drawRect(dX, (bY), sizeX, sizeY);

      painter.setPen(penLineWhite);
      painter.drawText(dX + 10, bY + 12 + 1 * 15, camera.ToString());

      QFont f1; f1.setBold(true);
      painter.setFont(f1);
      painter.drawText(dX + bX, bY + 12, QString("Camera Position"));
    }
    painter.end();
  }
}

/*!
\brief Key events.

Left and right keys: move sideway.

Up and down keys: move back and forth.

L : show light
A : show anchor
P : show plane

\param e Key events.
*/
void MayaWidget::keyPressEvent(QKeyEvent* e)
{
  switch (e->key())
  {
  case Qt::Key_Left:
    camera.SideWay(1.0); update();
    break;
  case Qt::Key_Right:
    camera.SideWay(-1.0); update();
    break;
  case Qt::Key_Up:
    camera.BackForth(1.0); update();
    break;
  case Qt::Key_Down:
    camera.BackForth(-1.0); update();
    break;

  case Qt::Key_A:
    if (e->modifiers() & Qt::AltModifier) showAncre(!useAncre);
    break;
  case Qt::Key_L:
    if (e->modifiers() & Qt::AltModifier) showLight(!useLight);
    break;
  case Qt::Key_P:
    if (e->modifiers() & Qt::AltModifier) showPlane(!usePlane);
    break;
  case Qt::Key_W:
    if (e->modifiers() & Qt::AltModifier) showWorld(!useWorld);
    break;
  case Qt::Key_B:
    if (e->modifiers() & Qt::AltModifier) showBBox(!useBBox);
    break;
  case Qt::Key_G:
    if (e->modifiers() & Qt::AltModifier) showBackG(!useBackG);
    break;

  case Qt::Key_R:
    if (e->modifiers() & Qt::AltModifier) showCameraRatio((++useCamRatio) % 3);// Alt + R : Camera ratio
    break;

  case Qt::Key_S:
    if (e->modifiers() & Qt::AltModifier) showStats(!useStats);// Alt + S : Statistics
    break;
  case Qt::Key_V:
    if (e->modifiers() & Qt::AltModifier) showVAxis(!useVAxis);
    break;

  case Qt::Key_T:
    if (e->modifiers() & Qt::AltModifier) showTime(!useTime);// Alt + T : Time
    break;
  case Qt::Key_C:
    // Shift + C : add camera
    if (e->modifiers() & Qt::ShiftModifier)
    {
      AddCamera();
    }
    // Alt + C : show camera
    else if (e->modifiers() & Qt::AltModifier)
    {
      showCamera(!useCamera);

    }
    // Delete camera
    else
    {
      DeleteCamera();
    }
    break;
  case Qt::Key_D:
    ChangeCamera();
    break;
  case Qt::Key_F1:
    SaveScreen();
    break;
  default:
    QOpenGLWidget::keyPressEvent(e);
  }
}

/*!
\brief Capture the rendering viewport and save it to disk.
*/
void MayaWidget::SaveScreen(int w, int h)
{
  QSize qs = size();
  resize(w, h);
  QImage image(size(), QImage::Format_ARGB32);

  // Set name according to date and time
  QString name = "../screen-" + System::DateTime();

  // Get image
  image = grabFramebuffer();

  // Save
  std::cout << int(image.save(name + ".png", "PNG")) << std::endl;
  std::cout << int(image.save(name + "-96.jpg", "JPG", 96)) << std::endl;
  std::cout << int(image.save(name + "-92.jpg", "JPG", 92)) << std::endl;

  resize(qs);
}

/*!
\brief Capture the rendering viewport and save it to disk.
\param w width of the image
\param h height of the image
\param url url and name of the exported file
*/
void MayaWidget::SaveScreen(int w, int h, const QString& url)
{
  QSize qs = size(); // Save buffersize

  resize(w, h);
  QImage image(size(), QImage::Format_ARGB32);

  // Get image
  image = grabFramebuffer();

  QUrl urlQT(url);
  QDir().mkpath(urlQT.adjusted(QUrl::RemoveFilename).path());

  // Save
  image.save(url, "PNG");

  resize(qs); // Restore buffersize
}


/*!
\brief This function has been overloaded mainly to enable and disable the display of the tool on the screen.
\param e Key events.
*/
void MayaWidget::keyReleaseEvent(QKeyEvent* e)
{
  // Update view
  QOpenGLWidget::keyReleaseEvent(e);
  update();
}

/*!
\brief Process mouse click events.
\param e Events.
*/
void MayaWidget::mousePressEvent(QMouseEvent* e)
{
  x0 = e->globalPosition().x();
  y0 = e->globalPosition().y();

  if (e->modifiers() & Qt::ControlModifier)
  {
    if (e->buttons() & Qt::LeftButton)
    {
      // Ctrl + Left Mouse Click 
      emit _signalEditSceneLeft(ancre);
    }
    else if (e->buttons() & Qt::RightButton)
    {
      // Ctrl + Right Mouse Click 
      emit _signalEditSceneRight(ancre);
    }
  }
  update();
}

/*!
\brief Overloaded function to capture the mouse events.

This function implements controls.
\param e Mouse move events.
*/
void MayaWidget::mouseMoveEvent(QMouseEvent* e)
{
  int x = e->globalPosition().x();
  int y = e->globalPosition().y();

  if (e->modifiers() & Qt::AltModifier)
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
    }
    else if (e->buttons() & Qt::MiddleButton)
    {
      // Alt + Left Mouse Move    : Plan displacement
      camera.LeftRightHorizontal((x - x0) * MoveScale);
      camera.UpDownVertical((y - y0) * MoveScale);
    }

    x0 = e->globalPosition().x();
    y0 = e->globalPosition().y();
    emit _signalMouseMove();
  }

  UpdateAnchor(e);
  update();
}

/*!
\brief Process the mouse double click events.
\param e Events.
*/
#pragma warning(push)
#pragma warning(disable: 4100)  
void MayaWidget::mouseDoubleClickEvent(QMouseEvent* e)
{
  camera.SetAt(ancre);
  update();
}
#pragma warning(pop)

/*!
\brief Overloaded function to capture the wheel events.
\param e Wheel events.
*/
void MayaWidget::wheelEvent(QWheelEvent* e)
{
  int numDegrees = e->angleDelta().y() / 8;
  int numSteps = numDegrees / 15;

  if (e->modifiers() & Qt::ControlModifier)
  {
    emit _signalEditHeight(ancre, numSteps);
  }
  else
  {
    // Zoom
    double MoveScale = Norm(camera.View()) * 0.025;
    if (e->angleDelta().y() > 0)  camera.BackForth(MoveScale);
    else                      camera.BackForth(-MoveScale);
    emit _signalMouseMove();
  }

  update();
}

/*!
\brief Process the mouse release events.
\param e Events.
*/
#pragma warning(push)
#pragma warning(disable: 4100)  
void MayaWidget::mouseReleaseEvent(QMouseEvent* e)
{
  update();
}
#pragma warning(pop)

/*!
\brief Update the position of the Anchor on the intersection between the ray (mouse direction and the plane).
\param e Mouse event.
*/
void MayaWidget::UpdateAnchor(QMouseEvent* e)
{
  if (!hasFocus()) setFocus(Qt::MouseFocusReason);

  double t;
  Ray ray = ConvertPixelToRay(e->pos());

  // Intersection with plane
  if (plane.Intersect(ray, t))
  {
    SetAnchor(ray(t));
  }
}

/*!
\brief Modify the anchor position.

Tha anchor is not rendered if is out of the plane.
\param v Position of the anchor
*/
void MayaWidget::SetAnchor(const Vector& v)
{
  ancre = v;
}

/*!
\brief Set the viewing camera.
\param c %Camera.
*/
void MayaWidget::SetCamera(const Camera& c)
{
  makeCurrent();
  camera = c;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(Math::RadianToDegree(camera.GetAngleOfViewV(width(), height())), (GLdouble)width() / (GLdouble)height(), camera.GetNear(), camera.GetFar());
  update();
}

/*!
\brief Get the camera.
\return Camera.
*/
Camera MayaWidget::GetCamera() const
{
  return camera;
}

/*!
\brief Return the set of cameras.
*/
QVector<Camera> MayaWidget::GetCameras() const
{
  return cameras.All();
}

/*!
\brief Add the current camera to the set of reference cameras.
*/
void MayaWidget::AddCamera()
{
  cameras.Push(camera);
}

/*!
\brief Remove the current camera from the set.
*/
void MayaWidget::DeleteCamera()
{
  cameras.Pop();
}

/*!
\brief Charge une nouvelle position de la camera
*/
void MayaWidget::ChangeCamera()
{
  cameras.Next(camera);
}

/*!
\brief Update the GPU parameters information (light, camera).
*/
void MayaWidget::updateGPU()
{
  glGetFloatv(GL_MODELVIEW_MATRIX, gpuparam.transform.ModelViewMatrix);
  glGetFloatv(GL_PROJECTION_MATRIX, gpuparam.transform.ProjectionMatrix);

  glBindBuffer(GL_UNIFORM_BUFFER, gpuparam.buffer_transform);
  glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(gpuparam.transform), &gpuparam.transform);
  glBindBufferBase(GL_UNIFORM_BUFFER, 1, gpuparam.buffer_transform);

  gpuparam.light = light;
  gpuparam.camera = camera.Eye();

  // Activation du shader, draw
  glUseProgram(gpuparam.shader_program->GetProgram());

  // Set shader parameters
  glUniform3f(glGetUniformLocation(gpuparam.shader_program->GetProgram(), "lightPos_vert"), gpuparam.light[0], gpuparam.light[1], gpuparam.light[2]);
  glUniform3f(glGetUniformLocation(gpuparam.shader_program->GetProgram(), "light_dir_axel"), gpuparam.light[0], gpuparam.light[1], gpuparam.light[2]);
  glUniform3f(glGetUniformLocation(gpuparam.shader_program->GetProgram(), "cameraPos_vert"), gpuparam.camera[0], gpuparam.camera[1], gpuparam.camera[2]);
  // Bind GeometryShader
  glUniform2f(glGetUniformLocation(gpuparam.shader_program->GetProgram(), "WIN_SCALE"), this->size().width() / 2., this->size().height() / 2.);

  glUseProgram(0);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

/*!
\brief Convert a pixel on the screen of the viewport into a ray.

This function relies on the camera model to compute the ray.

\param pixel The pixel coordinates.
*/
Ray MayaWidget::ConvertPixelToRay(const QPoint& pixel) const
{
  return camera.PixelToRay(pixel.x(), pixel.y() - 1, width(), height());
}

/*!
\brief Convert a pixel on the screen of the viewport into a ray, with anti-aliasing.

\param pixel The pixel coordinates
\param a Anti-aliasing.
\param x,y Integer coordinates of the sub-pixel.
*/
Ray MayaWidget::ConvertPixelToRay(const QPoint& pixel, int a, int x, int y) const
{
  return camera.PixelToRay(pixel.x(), pixel.y() - 1, width(), height(), a, x, y);
}
