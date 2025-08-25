// Maya  

#pragma once

#ifdef _MSC_VER
#include "libs/GL/glew.h"
#else
#include "GL/glew.h"
#endif

#include <QtCore/QVector> 
#include <QtGui/QKeyEvent>
#include <QtGui/QMouseEvent>
#include <QtWidgets/QApplication>
#include <QtOpenGLWidgets/QOpenGLWidget>
#include <QtGui/QPainter>

#include "libs/mayashader.h"
#include "libs/mayarender.h"
#include "libs/mayageometry.h"
#include "libs/camera.h"

struct Transform {
  GLfloat ModelViewMatrix[16];			//!< modelview matrix of the transformation
  GLfloat ProjectionMatrix[16];			//!< projection matrix of the transformation
};

struct GpuUniform
{
  //Uniforms
  GLint loc_wireframe;
  GLint loc_texture;
  GLint loc_applyNormalMap;
  GLint loc_normalMapStrength;
  GLint loc_textureScaling;

  //Atributes
  GLint loc_position;
  GLint loc_normal;
  GLint loc_color;
  GLint loc_UV;
  GLint loc_inst_pos;
  GLint loc_inst_scl;
  GLint loc_inst_rot;

  //Subroutine indexes
  GLuint phongIndex;
  GLuint phongVCIndex;
  GLuint normalIndex;
  GLuint classiqueIndex;
  GLuint gridIndex;
  GLuint ignIndex;
  GLuint goochIndex;
  GLuint terrain3DIndex;
  GLuint UVTextureModel;
  GLuint bumpNormalIndex;
  GLuint ThinSnowIndex;
  GLuint GreyCoolWarmIndex;
  GLuint AxelTerrainIndex;
  GLuint WireframeIndex;
  GLuint WireframeAspectIndex;
  GLuint TriPlanarIndex;
};

struct GpuParameter
{
  MayaShader* shader_program;
  GpuUniform uniforms;
  GLuint buffer_textureT[8];	//!< Buffer of texture : Terrain(4) + AxelTerrain(4)

  // ENV
  GLuint buffer_textureE[6];	//!< Buffer of texture : env

  Vector light; //!< Light
  Vector camera; //!< %Camera

  GLuint buffer_transform;	//!< Uniform buffer transform
  Transform transform;		//!< Projection and modelview matrix 
};

struct MayaGpuMaterial
{
  GLfloat ambient[4];
  GLfloat diffuse[4];
  GLfloat specular[4];
  GLfloat shininess;
};

class MayaTimeFrame {
public:
  FrameScaled mf;//!< Frame.
  int time;//!< Time of the frame
  bool acc;//!< Acceleration
};

// Frame class for instances
class MayaAnimatedFrame : public FrameScaled
{
protected:
  QVector<MayaTimeFrame> mafs; //!< List of key-frames
  int currentTime;
public:
  MayaAnimatedFrame();
  //! Empty
  ~MayaAnimatedFrame() {}

  void Append(const FrameScaled&, int, bool = false);
  void Delete(int);

  int start();
  int end();

  int size() { return mafs.size(); }

  bool setTime(int);
  bool setTimeSmooth(int);


  MayaTimeFrame operator[](int) const;
};

inline MayaTimeFrame MayaAnimatedFrame::operator[](int num) const
{
  return mafs[num];
}


class MayaGpu;
class MayaGpuSet;
class MayaGpuAll;

class MayaResources;


// Maya instance
class MayaInstance {
protected:
  QString name;      //!< Signature of the object, implemented as the name of the corresponding mesh file
  QString address;   //!< Location of the object
public:
  MayaInstance();
  explicit MayaInstance(const QString&);
  explicit MayaInstance(const QString&, const QString&);
  //! Empty 
  ~MayaInstance() {}

  // Setters
  void SetName(const QString&);
  void SetAddress(const QString&);

  // Getters
  QString GetName() const;
  QString GetAddress() const;
};

/*!
\brief Get name.
*/
inline QString MayaInstance::GetName() const
{
  return name;
}

/*!
\brief Get the address.
*/
inline QString MayaInstance::GetAddress() const
{
  return address;
}

/*!
\brief Set the address.
\param n The name.
*/
inline void MayaInstance::SetName(const QString& n)
{
  name = n;
}

/*!
\brief Set the address.
\param a The address.
*/
inline void MayaInstance::SetAddress(const QString& a)
{
  address = a;
}

// Set of Maya instances
class MayaInstanceSet : public MayaInstance
{
protected:
  QVector<FrameScaled> frames; //!< Array of frames
public:
  //! Empty.
  MayaInstanceSet() :MayaInstance() {}
  MayaInstanceSet(const QString&);
  MayaInstanceSet(const MayaInstance&);
  MayaInstanceSet(const QString&, const QString&);
  MayaInstanceSet(const QString&, const FrameScaled&);
  MayaInstanceSet(const QString&, const QString&, const FrameScaled&);

  //! Empty
  ~MayaInstanceSet() {}

  // Setters
  void Append(const FrameScaled&);

  void Rotate(const Vector&);
  void Translate(const Vector&);
  void Scale(const Vector&);

  void Clear();
  void Remove(int);

  MayaInstanceSet Crop(const Box2& box) const;

  FrameScaled GetFrameScaled(int) const;
  QVector<FrameScaled> GetFrameScaleds() const;
  int count() const;

  QString GetText(int = 0, bool = false) const;

  QString ExportVueScript(int, int);

  void exportSceneJson(QTextStream& out) const;
};

/*!
\brief Create an instance set.
\param name The name of the set.
*/
inline MayaInstanceSet::MayaInstanceSet(const QString& name) :MayaInstance(name)
{
}

/*!
\brief Create an instance set.
\param i Input instance.
*/
inline MayaInstanceSet::MayaInstanceSet(const MayaInstance& i) :MayaInstance(i.GetName(), i.GetAddress())
{
}

/*!
\brief Create an instance set.
\param name The name of the set.
\param address Directory.
*/
inline MayaInstanceSet::MayaInstanceSet(const QString& name, const QString& address) :MayaInstance(name, address)
{
}

/*!
\brief Create an instance set.
\param name The name of the set.
\param frame The frame.
*/
inline MayaInstanceSet::MayaInstanceSet(const QString& name, const FrameScaled& frame) :MayaInstance(name)
{
  frames.append(frame);
}

/*!
\brief Create an instance set.
\param name The name of the set.
\param address Directory.
\param frame The frame.
*/
inline MayaInstanceSet::MayaInstanceSet(const QString& name, const QString& address, const FrameScaled& frame) :MayaInstance(name, address)
{
  frames.append(frame);
}

/*!
\brief Get the i-th frame.
\param i Index.
*/
inline FrameScaled MayaInstanceSet::GetFrameScaled(int i) const
{
  return frames[i];
}

/*!
\brief Get the set of frames.
*/
inline QVector<FrameScaled> MayaInstanceSet::GetFrameScaleds() const
{
  return frames;
}

/*!
\brief Get the number of frames.
*/
inline int MayaInstanceSet::count() const
{
  return frames.count();
}

/*!
\brief Remove the i-th element.
\param i Element.
*/
inline void MayaInstanceSet::Remove(int i)
{
  frames.remove(i);
}

// Set of instances along with their frames
class MayaInstanceAll {
protected:
  QMap<QString, MayaInstanceSet> instances; //!< Map between names and set of instances.
public:
  //! Empty 
  MayaInstanceAll() {}
  //! Empty 
  ~MayaInstanceAll() {}

  // Setters
  void Append(const MayaInstance&);
  void Append(const MayaInstanceSet&);
  void Append(const MayaInstanceAll&);
  void Append(const QString&, const FrameScaled&);

  void AppendNew(const MayaInstanceSet&);

  // Getters
  QList<QString> GetNames() const;
  MayaInstanceSet operator[](const QString&) const;

  void Rotate(const Vector&);
  void Translate(const Vector&);
  void Scale(const Vector&);

  void Remove(const QString&);

  void Clear();
  void ClearFrames();

  MayaInstanceAll Crop(const Box2&) const;

  // Getters
  QVector<FrameScaled> GetFrameScaleds(const QString&) const;

  // Load and save functions
  bool LoadInstances(const QString&);
  bool Save(const QString&, const char* = nullptr) const;
  QString ExportVueScript();

  void exportSceneJson(QTextStream& out) const;

  QString GetText(int = 0, bool = false) const;

  // Friends
  friend class MayaGeometryAll;
};

/*!
\brief Get the list of names from the set of instances.
*/
inline QList<QString> MayaInstanceAll::GetNames() const
{
  return instances.keys();
}

/*!
\brief Get one set from the list of all instances.
\param name The name.
*/
inline MayaInstanceSet MayaInstanceAll::operator[](const QString& name) const
{
  return instances[name];
}


/*!
\brief Get the set of frames from the list of all instances.
\param n The name.
*/
inline QVector<FrameScaled> MayaInstanceAll::GetFrameScaleds(const QString& n) const
{
  return instances[n].GetFrameScaleds();
}

// Maya Gpu instance
class MayaGpu {
protected:
  QString name;              //!< Signature of the object, implemented as the name of the corresponding mesh file
  MayaStatistics statistics; //!< Statistics of the object, which are derived from the arguments of the constructor
  MayaMaterial mat;          //!< Material structure of the Object.

  GLuint buffer_vertex; //!< Buffer of vertex
  GLuint buffer_normal; //!< Buffer of normal
  GLuint buffer_color;  //!< Buffer of color
  GLuint buffer_UVmap;  //!< Buffer of UVMaps

  GLuint buffer_index;  //!< Buffer of indexes

  MayaGpuMaterial material;	//!< Material data
  GLuint buffer_material;	//!< Uniform material buffer 
  GLuint buffer_texture;    //!< Buffer of texture

  int nv;               //!< Number of vertices
  int nt;               //!< Number of triangles
  Box bbox;             //!< Bounding box

  MayaSimpleRenderer* boxRenderer; //!< Renderer used to draw the bounding box

  void InitRenderer();
  void CreateGpu(const MayaGeometry&);
  void DeleteBuffers();

public:
  //! Empty
  MayaGpu();
  MayaGpu(const MayaGeometry&);
  ~MayaGpu();

  // Renderer
  void RenderBBox();

  // Setters
  void SetName(const QString&);
  void SetMaterial(const MayaMaterial&);

  void InitMaterial();
  void InitTexture();

  // Getters
  QString GetName() const;
  Box getBox();

  MayaStatistics GetStatistics() const;
};

/*!
\brief Set the name of the object
\param n Name.
*/
inline void MayaGpu::SetName(const QString& n)
{
  name = n;
}

/*!
\brief Get the name of the object.
*/
inline QString MayaGpu::GetName() const
{
  return name;
}

/*!
\brief Get the bounding box of the object.
*/
inline Box MayaGpu::getBox()
{
  return bbox;
}

// Set of Maya instances
class MayaGpuSet :public MayaGpu
{
protected:
  QVector<FrameScaled> frames;        //!< Array of frames
  GLuint buffer_instances_position = 0; //!< Buffer of array of positions
  GLuint buffer_instances_rotation = 0; //!< Buffer of array of rotations
  GLuint buffer_instances_scale = 0;    //!< Buffer of array of scales

  GLuint buffer_Vao = 0;
public:
  //! Empty
  MayaGpuSet() { }
  MayaGpuSet(const MayaGeometrySet&);
  ~MayaGpuSet() {}

  void InitInstances();
  void InitVAO(GpuParameter&);

  void RefreshFrames();

  bool isWireframe() const;

  // Renderer
  void Render(GpuParameter&);
  void RenderBBox();
  void SetSubroutine(GpuParameter&) const;

  // Setters
  void Append(const FrameScaled&);

  void ClearFrames();

  // Getters
  int count() const;
  FrameScaled GetFrameScaled(int) const;
  MayaStatistics GetStatistics() const;
public:
  void DeleteBuffers();
};

/*!
\brief Add an instance to the set.
\param frame The frame of the instance.
*/
inline void MayaGpuSet::Append(const FrameScaled& frame)
{
  // Statistics
  statistics.AddInstances();

  frames.append(frame);
}

/*!
\brief Get the number of frames in the set.
*/
inline int MayaGpuSet::count() const
{
  return frames.size();
}

/*!
\brief Get the number i-th frame in the set.
\param i Index.
*/
inline FrameScaled MayaGpuSet::GetFrameScaled(int i) const
{
  return frames[i];
}

// Set of instances along with their frames
class MayaGpuAll {
protected:
  QMap<QString, MayaGpuSet> instances; //!< Map of instances with their set of frames
public:
  //! Empty
  MayaGpuAll() {}
  MayaGpuAll(const MayaGeometryAll&);
  //! Empty
  ~MayaGpuAll() {}

  // Renderer
  void Render(GpuParameter&);
  void RenderBBox();

  // Setters
  void Append(const MayaGeometrySet&);
  void Append(const MayaGeometryAll&);
  void Replace(const MayaGeometryAll&);

  void UpdateMaterial(MayaGeometryAll&);
  void Remove(const QString&);
  void Clear();
  void ClearFrames();

  // Getters
  Box getBox();

  // Statistics
  MayaStatistics GetStatistics() const;
};


class MayaGeometryStack :public MayaGeometry
{
protected:
  FrameScaled theframe; //!< Current transformation.
  QVector<FrameScaled> stackframe; //!< Stack of frames.
public:
  MayaGeometryStack();
  MayaGeometryStack(const MayaGeometry&);
  //! Empty
  ~MayaGeometryStack() {}

  // Stack management
  void Push();
  void Pop();
  void PopAll();

  void Transform(const FrameScaled&);
  void Translate(const Vector&);
  void Rotate(const Vector&);
  void Scale(const Vector&);

  // Add simple primitives
  void AddTriangle(const Vector&, const Vector&, const Vector&);

  // Add a complete geometry
  void AddGeometry(const MayaGeometryStack&);

  FrameScaled GetFrame() const;
protected:
  void Reserve(int, int, int);
  void AddTriangleOffset(int, int, int, int);
};

class MayaSceneStack :public MayaGeometryStack
{
protected:
  MayaResources atlas; //!< Set of resources.
  MayaInstanceAll scene; //!< The scene.
public:
  MayaSceneStack();
  MayaSceneStack(const MayaGeometry&);
  //! Empty
  ~MayaSceneStack() {}

  // Scene creation members
  void AddGeometryToAtlas(const QString&);
  void AddGeometryToAtlas(const QString&, const MayaGeometry&);
  void AddInstanceToAtlas(const QString&, const QString&);

  void AddInstanceToScene(const QString&);
  void AddInstanceToScene(const QString&, const FrameScaled&);
  void AddGeometryToScene(const QString&);

  // Create the scene
  MayaGeometryAll GetScene() const;

  QString GetText(bool = false) const;
};

class MayaPlane {
protected:
  Box2 area = Box2::Unit; //!< Area of the grid.
  double height = 0.0; //!< Elevation.
  double line = 1.0; //!< Line interval.
  MayaSimpleRenderer* planeRendererWhiteLine = nullptr;//!< Renderer used to draw the white line of the plane
  MayaSimpleRenderer* planeRendererQuad = nullptr;//!< Renderer used to draw the transparent quad of the plane
public:
  explicit MayaPlane(const Box2 & = Box2(10.0), const double& = 0.0, const double& = 1.0);
  ~MayaPlane();
  double Z() const;
  Box2 GetArea() const;
  void SetArea(const Box2&);
  void Translate(const double&);
  bool Intersect(const Ray&, double&) const;
  void Render();
private:
  static const double epsilon; //!< Epsilon value for offseting lines when rendering the plane.
protected:
  void InitRenderer();
};

class MayaCameraSet
{
protected:
  QVector<Camera> cameras; //!< Set of cameras.
  int n; // Current camera index.
public:
  MayaCameraSet();
  void Push(const Camera&);
  void Pop();
  Camera Current() const;
  void Next(Camera&);
  QVector<Camera> All() const;
};

class MayaWidget :public QOpenGLWidget
{
  Q_OBJECT
protected:
  MayaGpuAll gpuall;  //!< Instances for storing the objects on the GPU.
  MayaGpuSet anchor; //!< Internal instance for drawing the anchor point.
  MayaGpuSet lightshape; //!< Light shape.

  MayaGpuSet sphereshape; //!< Transparent sphere.

  Camera camera;//!< Current camera.

  int x0, y0; //!< Reference mouse coordinates.
  bool useAncre, usePlane, useWorld, useLight, useBBox, useBackG, useStats, useVAxis, useTime, useCamera; //!< Boolean flags.
  short int useCamRatio;
  MayaPlane plane; //!< Horizontal plane.
  Vector ancre = Vector::Null; //!< Focus point.
  Vector light; //!< Light point.
  GpuParameter gpuparam; //!< Set of parameters defining the way objects are rendered.
  MayaShader* gpubackground = nullptr; //!< Program used to draw background;
  GLuint backgroundVAO; //!< VAO used when drawing the background.

  MayaCameraSet cameras; //!< Set of cameras.
public:
  MayaWidget(QWidget* = nullptr);
  ~MayaWidget();

  // Camera
  void SetCamera(const Camera&);
  QVector<Camera> GetCameras() const;
  Camera GetCamera() const;
  void AddCamera();
  void DeleteCamera();
  void ChangeCamera();

  void SetPlanes(const double&, const double&);

  // Plane
  void SetPlane(const Box2&, const double&, const double&);

  // Lighting
  void SetLight(const Vector&);

  // Anchor
  void SetAnchor(const Vector&);
  Ray ConvertPixelToRay(const QPoint&) const;
  Ray ConvertPixelToRay(const QPoint&, int, int, int) const;

  // World Operators
  void SetWorld(const MayaGeometryAll&);
  void AppendWorld(const MayaGeometryAll&);
  void ReplaceInWorld(const MayaGeometryAll&);
  void ReplaceMaterialInWorld(MayaGeometryAll&);
  void ClearWorld();

  void showAncre(bool b) { useAncre = b; }
  void showLight(bool b) { useLight = b; }
  void showPlane(bool b) { usePlane = b; }
  void showWorld(bool b) { useWorld = b; }
  void showBBox(bool b) { useBBox = b; }
  void showBackG(bool b) { useBackG = b; }
  void showStats(bool b) { useStats = b; }
  void showVAxis(bool b) { useVAxis = b; }
  void showTime(bool b) { useTime = b; }
  void showCamera(bool b) { useCamera = b; }
  void showCameraRatio(int b) { useCamRatio = b; }

  void SaveScreen(int = 1280, int = 720);
  void SaveScreen(int, int, const QString&);
protected:
  void initializeGL();
  void resizeGL(int, int);

  virtual void UpdateAnchor(QMouseEvent*);

  virtual void RenderAnchor();
  virtual void RenderLight();
  void RenderBackGround();
  void RenderForeGround(QElapsedTimer);

  virtual void InitGPU();
  void InitProgram();
  void InitCamera();

  virtual void InitAnchorShape();
  virtual void InitLightShape();
  virtual void InitSphereShape();
  void InitTerrainTextures();
  void InitEnvTextures();

  void InitLight();
  void updateGPU();

public:
  static bool Verbose;			//!< Static flag for verbose output in the console (useful for OpenGL debugging)
  static GLenum GLVerboseLevel;	//!< OpenGL verbosity level. Everything below is ignored.
  static void GLAPIENTRY OpenGLDebugMessageCallback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar* message, const void* = nullptr);

public slots:
  void paintGL();

  virtual void mousePressEvent(QMouseEvent*);
  virtual void mouseReleaseEvent(QMouseEvent*);
  virtual void mouseDoubleClickEvent(QMouseEvent*);
  virtual void mouseMoveEvent(QMouseEvent*);
  virtual void wheelEvent(QWheelEvent*);
  virtual void keyPressEvent(QKeyEvent*);
  virtual void keyReleaseEvent(QKeyEvent*);

signals:
  void _signalEditHeight(const Vector&, int);
  void _signalEditSceneLeft(const Vector&);
  void _signalEditSceneRight(const Vector&);
  void _signalMouseMove();
};

