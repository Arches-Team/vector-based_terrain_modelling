#pragma once

#include "libs/gpu-shader.h"
#include "libs/gpu-wrapper.h"
#include "libs/camera.h"
#include "libs/evectorfloat.h"
#include "libs/color.h"
#include <QtCore/QMap> 

class ScalarField2;
class MeshColor;
class Mesh;

// Utility class for profiling CPU & GPU
typedef std::chrono::time_point<std::chrono::high_resolution_clock> MyChrono;

class RenderingProfiler
{
public:
  bool enabled = false;			//!< Flag linked to UI.

  GLuint query;					//!< GL Query for stats
  GLuint64 elapsedTimeGPU;		//!< GPU rendering time for a frame.

  int nbframes = 0;				//!< CPU Frame counter.
  MyChrono start;					//!< CPU profiler.
  double msPerFrame = 0;			//!< Recorded info.
  double framePerSecond = 0;		//!< Recorded info. 
  int drawCallPerFrame = 0;		//!< Recorded info.

  /*!
  \brief Init the profiler. Only has to be done once in the program.
  */
  inline void Init()
  {
    glGenQueries(1, &query);
    start = std::chrono::high_resolution_clock::now();
  }

  /*!
  \brief Starts profiling the GPU if enabled.
  */
  inline void BeginGPU()
  {
    if (enabled)
      glBeginQuery(GL_TIME_ELAPSED, query);
  }

  /*!
  \brief Ends the GPU profiling if enabled.
  */
  inline void EndGPU()
  {
    if (enabled)
    {
      glEndQuery(GL_TIME_ELAPSED);
      int done = 0;
      while (!done)
        glGetQueryObjectiv(query, GL_QUERY_RESULT_AVAILABLE, &done);
      glGetQueryObjectui64v(query, GL_QUERY_RESULT, &elapsedTimeGPU);
    }
  }

  /*!
  \brief Update the CPU profiling.
  */
  inline void Update()
  {
    nbframes++;
    auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count();
    double seconds = double(microseconds) / 1000000.0;
    if (seconds >= 1.0)
    {
      msPerFrame = seconds * 1000.0 / nbframes;
      framePerSecond = nbframes / seconds;
      nbframes = 0;
      start = std::chrono::high_resolution_clock::now();
    }
  }
};

enum class MeshMaterial
{
  Normal = 0,
  Color = 1,
  Aspect = 2,
};

enum class MeshShading
{
  Triangles = 0,
  Lines = 1,
};

class MeshWidget : public QOpenGLWidget
{
  // Must include this if you use Qt signals/slots
  Q_OBJECT

protected:
  // Internal definition of renderable mesh
  class MeshGL
  {
  public:
    bool enabled;				//!< Render flag. Mesh is not rendered if enabled equals false.
    GLuint vao;					//!< Mesh VAO.
    GLBuffer fullBuffer;		//!< Mesh buffer. Contains 3D vertices, 3D normals, and possibly colors.
    GLBuffer indexBuffer;		//!< Mesh index buffer.
    int triangleCount;			//!< Number of triangles.
    float TRSMatrix[16];		//!< Translation-Rotation-Scale matrix (computed from a FrameScaled).

    MeshShading shading;		//!< Render flag.
    MeshMaterial material;		//!< Render flag.
    bool useWireframe;			//!< Render flag.

  public:
    MeshGL();
    MeshGL(const Mesh&, const FrameScaled&);
    MeshGL(const MeshColor&, const FrameScaled&);
    MeshGL(const Box&, const FrameScaled & = FrameScaled::Id);

    void Delete();
    void SetFrame(const FrameScaled& fr);
  };
  typedef QMap<QString, MeshGL*>::iterator MeshIterator;

protected:
  // Scene
  int x0 = 0, y0 = 0;  //!< Reference mouse coordinates.
  Camera camera; //!< %Camera.
  bool perspectiveProjection;
  float cameraOrthoSize;
  bool MoveAt = false;
  Vector currentAt = Vector::Null;
  Vector toAt = Vector::Null;
  int stepAt = 0;
  Vector2 nearAndFarPlane = Vector2(1.0, 50000.0); //!< Near and far distance.

  // Skybox
  bool useSkyShader = true;
  Color backgroundColor;
  GLShader skyShader;
  GLuint skyboxVAO = 0;

  // Profiler & panel flags
  RenderingProfiler profiler;
  bool renderCameraPanel = false;

  // Meshes
  GLShader meshShader;
  QMap<QString, MeshGL*> objects;

  // Boxes
  GLShader boxShader;
  QMap<QString, MeshGL*> boxObjects;

public:
  MeshWidget(QWidget* parent = nullptr);
  ~MeshWidget();

  // Meshes
  void AddMesh(const QString& name, const Mesh& mesh, const FrameScaled & = FrameScaled::Id);
  void AddMesh(const QString& name, const MeshColor& mesh, const FrameScaled & = FrameScaled::Id);
  void DeleteMesh(const QString& name);
  bool HasMesh(const QString& name) const;
  void ClearAll();
  void UpdateMesh(const QString& name, const FrameScaled& frame);
  void EnableMesh(const QString& name);
  void DisableMesh(const QString& name);

  // Material, shading type, bounding box
  void SetMaterial(const QString& name, MeshMaterial mat);
  void SetMaterialGlobal(MeshMaterial mat);
  void UseWireframe(const QString& name, bool wireframe);
  void UseWireframeGlobal(bool wireframe);
  void SetShading(const QString& name, MeshShading shading);
  void SetShadingGlobal(MeshShading shading);
  void UseBoundingBox(const QString& name, bool use);

  // Camera, sky, screen
  QPoint GetMousePosition() const;
  void SetCamera(const Camera& cam);
  void SetCameraMode(bool);
  void SetNearAndFarPlane(const Vector2&);
  Camera GetCamera() const;
  Ray ConvertPixelToRay(const QPoint&) const;
  void SetSky(bool useShader, const Color& clearColor = Color::White);
  QImage GrabScreen(int = 1920, int = 1080);
  void SaveScreen(int = 1920, int = 1080);
  void SaveScreen(int, int, const QString&);

protected:
  virtual void initializeGL();
  virtual void resizeGL(int, int);
  virtual void paintGL();
  virtual void RenderSky();
  virtual void RenderUiPanels();
  virtual void RenderMeshes();
  virtual void reloadShaders();

private:
signals:
  void _signalUpdate();
  void _signalMouseMove();
  void _signalMouseRelease();
  void _signalMouseMoveEdit(const Ray&);
  void _signalEditSceneLeft(const Ray&);
  void _signalEditSceneRight(const Ray&);

public slots:
  virtual void mousePressEvent(QMouseEvent*);
  virtual void mouseReleaseEvent(QMouseEvent*);
  virtual void mouseDoubleClickEvent(QMouseEvent*);
  virtual void mouseMoveEvent(QMouseEvent*);
  virtual void wheelEvent(QWheelEvent*);
  virtual void keyPressEvent(QKeyEvent*);
  virtual void keyReleaseEvent(QKeyEvent*);
};

class TerrainRaytracingWidget : public MeshWidget
{
  // Must include this if you use Qt signals/slots
  Q_OBJECT

protected:
  // CPU Data
  ScalarField2* hf = nullptr;				//!< Pointer to heightfield.
  Box2 bbox;						//!< Bounding box of the heightfield.
  float zMin = 0.0, zMax = 0.0;				//!< Min/max elevation of the heighfield.
  float K = 1.0;						//!< Global Lipschitz constant of the heightfield.
  int nx = 0, ny = 0;						//!< Grid sizes.

  // GPU Data
  GLuint shaderProgram = 0;			//!< GL program for shader.
  GLuint raytraceVAO = 0;				//!< Raytracer VAO.
  GLBuffer hfBuffer;				//!< Heightfield elevation buffer.
  GLBuffer shadingBuffer;			//!< Shading buffer.
  std::vector<float> tmpData;		//!< Temporary float vector.

  GLuint albedoTextureBuffer = 0;			//!< Albedo texture buffer.
  bool useAlbedo = false;				//!< Albedo texture flag.
  bool useWireframe = false;			//!< Wireframe flag.
  bool useCost = false;				//!< Cost shading flag.
  bool useElevationShading = false;	//!< Elevation shading flag.
  bool useGreenBrownYellow = false;	//!< Shading from color map flag.
  bool useShadingBuffer = false;		//!< Additional shading buffer use flag.
  float cameraAngleOfViewV = 0.0;

public:
  TerrainRaytracingWidget(QWidget* parent = nullptr);
  ~TerrainRaytracingWidget();

  virtual void SetHeightField(ScalarField2* hfPtr);
  virtual void UpdateInternal();

  void UpdateBuffer(GLuint);
  void SetShading(const ScalarField2&);
  void UseShading(const bool& = true);
  void SetAlbedo(const QImage&);
  void UseAlbedo(const bool& = true);
  void UseWireframe(const bool& = true);
  void UseCost(const bool& = true);
  void UseElevationShading(const bool&);
  void UseGreenBrownYellowShading(const bool&);
  void SetCamera(const Camera&);
  void Reload();

  void SetElevationRange(double, double);
  void SetSteps(int);
  void SetEpsilon(float);
  void SetLight(const Vector&);
  void SetAntialiasing(int);

  void UseSmoothShadow(const bool& = true);
  void SetSmoothShadowSteps(int);
  void SetSmoothShadowMarchingSteps(int);
  void SetSmoothShadowMarchingEpsilon(float);
  void SetSmoothShadowStrength(float);

  void UseSelfShadow(const bool& = true);
  void SetSelfShadowMarchingSteps(int);
  void SetSelfShadowMarchingEpsilon(float);
  void SetSelfShadowStrength(float);

protected:
  void paintGL() override;
  void initializeGL() override;
  void reloadShaders() override;
};



class SphereTracingWidget : public MeshWidget
{
  // Must include this if you use Qt signals/slots
  Q_OBJECT

protected:

  class PointGL
  {
  public:
    GLuint vao;					//!< Mesh VAO.
    GLBuffer fullBuffer;		//!< Mesh buffer. Contains 3D vertices, 3D normals, and possibly colors.
    GLBuffer indexBuffer;		//!< Mesh index buffer.
    int pointCount;			//!< Number of triangles.
    float TRSMatrix[16];		//!< Translation-Rotation-Scale matrix (computed from a FrameScaled).

  public:
    PointGL();
    PointGL(const std::vector<Vector>, const std::vector<float>, const FrameScaled& fr = FrameScaled::Id);

    void Delete();
    void SetFrame(const FrameScaled& fr);
  };

  class SDFGL
  {
  public:
    bool enabled;				//!< Render flag. sdf is not rendered if enabled equals false.
    float TRSMatrix[16];		//!< Translation-Rotation-Scale matrix (computed from a FrameScaled).
    GLShader shader;			//!< compiled shader to render the sdf using sphere tracing
    int shading;
    Vector sampleFieldNormal;				//!< normal of the plane to visualize the value of the sdf
    Vector sampleFieldPoint;				//! point on the plane to visualize the value of the sdf
    float sampleFieldFreq;		//!< frequence of the sampleField
    int maxCost;				//!< maxCost for cost shading
    float kg = 1.0;					//! Geometric coeficient for interpolation node
    float kn = 0.0;					//! Normal coeficient for interpolation node

    Random r;                       //!random generator
    int nb_sample;                  //! nb of sample currently in the color texture
    Matrix4Float old_transformation;
    int textureWidth;               //! width of the texture
    int textureHeight;              //! height of the texture
    GLuint sexyFramebuffer;            //! FBO to render the accumulated sample of the sexy shading
    GLuint sexyColorTexture;       //! buffer to store all samples of the sexy shading


    static const int SEXY_SHADING = 5;

  public:
    SDFGL();
    SDFGL(const QString& sdf_code_glsl, const FrameScaled& fr = FrameScaled::Id);

    void Delete();
    void SetFrame(const FrameScaled& fr);
    void setShading(int newShading);
    void resetColorTexture(int width, int height);

  };

protected:
  QMap<QString, SDFGL*> sdfShaders;
  int nbSteps = 256; // Ray-tracing steps.
  int aa = 0; //!< Anti aliasing size.

  //evaluation points for BID
  GLShader pointShader;
  PointGL evaluationsPoints; //!< 3D points where we count the number of reccursive call done when evaluating the tree.


public:

  SphereTracingWidget(QWidget* parent = nullptr);
  ~SphereTracingWidget();

  virtual void paintGL();
  virtual void reloadShaders();

  void AddSDF(const QString&, const QString&);
  void DeleteSDF(const QString&);
  bool HasSDF(const QString&) const;
  void ClearAll();
  void UpdateSDF(const QString&, const FrameScaled&);
  void EnableSDF(const QString&);
  void DisableSDF(const QString&);
  void setShading(int);
  void setSampleFieldNormal(const Vector&);
  void setSampleFieldPoint(const Vector&);
  void setSampleFieldFreq(float);
  void setMaxCost(int);
  void setKG(float);
  void setKN(float);
  void setAA(int);

  void setEvaluationPoints(std::vector<Vector>, std::vector<float>);

  void setNbSteps(int);

protected:
  void RenderSDF();
  void RenderPoints();
};
