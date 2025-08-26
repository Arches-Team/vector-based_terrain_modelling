// Maya

#include <QtCore/QString>

#include "libs/cpu.h"
#include "libs/maya.h"

//!< Static flag for verbose output in the console (useful for OpenGL debugging)
bool MayaWidget::Verbose = false;

//!< OpenGL verbosity level. Everything below is ignored.
GLenum MayaWidget::GLVerboseLevel = GL_DEBUG_SEVERITY_MEDIUM;

/*!
\brief Custom callback for OpenGL message.
\param source
\param type
\param id
\param severity
\param length
\param message
\param userParam
*/
void GLAPIENTRY MayaWidget::OpenGLDebugMessageCallback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar* message, const void* userParam)
{
  // Filter below a given level
  if (severity < GLVerboseLevel)
    return;

  // Type of message
  std::cout << "GL CALLBACK: ";
  switch (type) {
  case GL_DEBUG_TYPE_ERROR:
    std::cout << "ERROR, ";
    break;
  default:
    std::cout << "INFO, ";
    break;
  }

  // Severity
  switch (severity) {
  case GL_DEBUG_SEVERITY_HIGH:
    std::cout << "severity = HIGH, ";
    break;
  case GL_DEBUG_SEVERITY_MEDIUM:
    std::cout << "severity = MEDIUM, ";
    break;
  case GL_DEBUG_SEVERITY_LOW:
    std::cout << "severity = LOW, ";
    break;
  case GL_DEBUG_SEVERITY_NOTIFICATION:
    std::cout << "severity = NOTIFICATION, ";
    break;
  }

  // Message
  std::cout << message << std::endl;
}

/*!
\brief Set up the rendering state.
*/
void MayaWidget::initializeGL()
{
  GLenum err = glewInit();
  if (err != GLEW_OK)
  {
    std::cout << "GLEW Error: " << glewGetErrorString(err) << std::endl;
    return;
  }

  if (Verbose)
  {
    std::cout << "OpenGL device info - device: " << (char*)glGetString(GL_VENDOR) << std::endl;
    std::cout << "OpenGL device info - version: " << (char*)glGetString(GL_VERSION) << std::endl;
    std::cout << "OpenGL device info - renderer: " << (char*)glGetString(GL_RENDERER) << std::endl;
    std::cout << "OpenGL device info - GLSL: " << (char*)glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;

    glEnable(GL_DEBUG_OUTPUT);
    glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
    glDebugMessageCallback(&OpenGLDebugMessageCallback, 0);
  }

  glEnable(GL_LINE_SMOOTH);
  glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_BLEND);

  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  InitLight();

  InitGPU();
}

/*!
\brief Set up the OpenGL Light rendering state.
*/
void MayaWidget::InitLight()
{
  //
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  //
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1);

  // Lumière globale
  GLfloat global_ambient[] = { 0.05f, 0.05f, 0.05f, 1.0 };
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);

  // Lumière 0
  GLfloat light_ambient0[] = { 0.0f, 0.0f, 0.0f, 1.0f };
  GLfloat light_diffuse0[] = { 0.6f, 0.6f, 0.6f, 1.0f };
  GLfloat light_specular0[] = { 0.2f, 0.2f, 0.2f, 1.0f };
  GLfloat light_position0[] = { 20.0f, 10.0f, 20.0f, 1.0f };

  glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient0);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
  glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular0);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position0);

  //
  GLfloat no_mat[] = { 0.0f, 0.0f, 0.0f, 1.0f };
  GLfloat mat_specular[] = { 0.2f, 0.2f, 0.2f, 1.0f };
  GLfloat mat_shininess[] = { 80.0 };
  GLfloat mat_ambient[] = { 0.5f, 0.5f, 0.5f, 1.0f };
  GLfloat mat_back[] = { 0.0f, 0.0f, 0.0f, 1.0f };

  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
  glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
  glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
  glMaterialfv(GL_BACK, GL_AMBIENT_AND_DIFFUSE, mat_back);
}

/*!
\brief Initialize GPU information (camera, textures, shaders).
*/
void MayaWidget::InitGPU()
{
  // Camera
  InitCamera();
  InitTerrainTextures();
  InitEnvTextures();
  InitProgram();
  // Anchor shapes
  InitAnchorShape();
  // Light shapes
  InitLightShape();
  // Transparent spherical shape
  InitSphereShape();
}

/*!
\brief Initialize the program shader.
*/
void MayaWidget::InitProgram()
{
  QString pPath = QString::fromStdString(std::string(SOLUTION_DIR));
  if (pPath.isEmpty()) std::cout << "MayaWidget::initProgram() : SOLUTION_DIR not defined" << std::endl;

  gpuparam.shader_program = new MayaShader(pPath + QString("/shaders/libs/default"));

  // Bind uniform texture
  GLuint loc_tex1 = glGetUniformLocation(gpuparam.shader_program->GetProgram(), "tex1");
  GLuint loc_tex2 = glGetUniformLocation(gpuparam.shader_program->GetProgram(), "tex2");
  GLuint loc_tex3 = glGetUniformLocation(gpuparam.shader_program->GetProgram(), "tex3");
  GLuint loc_tex4 = glGetUniformLocation(gpuparam.shader_program->GetProgram(), "tex4");
  GLuint loc_tex5 = glGetUniformLocation(gpuparam.shader_program->GetProgram(), "tex5");
  GLuint loc_tex6 = glGetUniformLocation(gpuparam.shader_program->GetProgram(), "tex6");
  GLuint loc_tex7 = glGetUniformLocation(gpuparam.shader_program->GetProgram(), "tex7");
  GLuint loc_tex8 = glGetUniformLocation(gpuparam.shader_program->GetProgram(), "tex8");

  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_2D, gpuparam.buffer_textureT[0]);
  glProgramUniform1i(gpuparam.shader_program->GetProgram(), loc_tex1, 1);

  glActiveTexture(GL_TEXTURE2);
  glBindTexture(GL_TEXTURE_2D, gpuparam.buffer_textureT[1]);
  glProgramUniform1i(gpuparam.shader_program->GetProgram(), loc_tex2, 2);

  glActiveTexture(GL_TEXTURE3);
  glBindTexture(GL_TEXTURE_2D, gpuparam.buffer_textureT[2]);
  glProgramUniform1i(gpuparam.shader_program->GetProgram(), loc_tex3, 3);

  glActiveTexture(GL_TEXTURE4);
  glBindTexture(GL_TEXTURE_2D, gpuparam.buffer_textureT[3]);
  glProgramUniform1i(gpuparam.shader_program->GetProgram(), loc_tex4, 4);

  glActiveTexture(GL_TEXTURE5);
  glBindTexture(GL_TEXTURE_2D, gpuparam.buffer_textureT[4]);
  glProgramUniform1i(gpuparam.shader_program->GetProgram(), loc_tex5, 5);

  glActiveTexture(GL_TEXTURE6);
  glBindTexture(GL_TEXTURE_2D, gpuparam.buffer_textureT[5]);
  glProgramUniform1i(gpuparam.shader_program->GetProgram(), loc_tex6, 6);

  glActiveTexture(GL_TEXTURE7);
  glBindTexture(GL_TEXTURE_2D, gpuparam.buffer_textureT[6]);
  glProgramUniform1i(gpuparam.shader_program->GetProgram(), loc_tex7, 7);

  glActiveTexture(GL_TEXTURE8);
  glBindTexture(GL_TEXTURE_2D, gpuparam.buffer_textureT[7]);
  glProgramUniform1i(gpuparam.shader_program->GetProgram(), loc_tex8, 8);


  // Bind uniform texture
  // ENV
  GLuint loc_down = glGetUniformLocation(gpuparam.shader_program->GetProgram(), "environment_down");
  GLuint loc_up = glGetUniformLocation(gpuparam.shader_program->GetProgram(), "environment_up");
  GLuint loc_left = glGetUniformLocation(gpuparam.shader_program->GetProgram(), "environment_left");
  GLuint loc_right = glGetUniformLocation(gpuparam.shader_program->GetProgram(), "environment_right");
  GLuint loc_front = glGetUniformLocation(gpuparam.shader_program->GetProgram(), "environment_front");
  GLuint loc_back = glGetUniformLocation(gpuparam.shader_program->GetProgram(), "environment_back");

  glActiveTexture(GL_TEXTURE9);
  glBindTexture(GL_TEXTURE_2D, gpuparam.buffer_textureE[0]);
  glProgramUniform1i(gpuparam.shader_program->GetProgram(), loc_down, 9);

  glActiveTexture(GL_TEXTURE10);
  glBindTexture(GL_TEXTURE_2D, gpuparam.buffer_textureE[1]);
  glProgramUniform1i(gpuparam.shader_program->GetProgram(), loc_up, 10);

  glActiveTexture(GL_TEXTURE11);
  glBindTexture(GL_TEXTURE_2D, gpuparam.buffer_textureE[2]);
  glProgramUniform1i(gpuparam.shader_program->GetProgram(), loc_left, 11);

  glActiveTexture(GL_TEXTURE12);
  glBindTexture(GL_TEXTURE_2D, gpuparam.buffer_textureE[3]);
  glProgramUniform1i(gpuparam.shader_program->GetProgram(), loc_right, 12);

  glActiveTexture(GL_TEXTURE13);
  glBindTexture(GL_TEXTURE_2D, gpuparam.buffer_textureE[4]);
  glProgramUniform1i(gpuparam.shader_program->GetProgram(), loc_front, 13);

  glActiveTexture(GL_TEXTURE14);
  glBindTexture(GL_TEXTURE_2D, gpuparam.buffer_textureE[5]);
  glProgramUniform1i(gpuparam.shader_program->GetProgram(), loc_back, 14);

  // Set shader parameters
  glUniformBlockBinding(gpuparam.shader_program->GetProgram(), glGetUniformBlockIndex(gpuparam.shader_program->GetProgram(), "transform"), 1);
  glUniformBlockBinding(gpuparam.shader_program->GetProgram(), glGetUniformBlockIndex(gpuparam.shader_program->GetProgram(), "material"), 0);

  //Get Uniform location
  gpuparam.uniforms.loc_wireframe = glGetUniformLocation(gpuparam.shader_program->GetProgram(), "wireframe");
  gpuparam.uniforms.loc_applyNormalMap = glGetUniformLocation(gpuparam.shader_program->GetProgram(), "applyNormalMapping");
  gpuparam.uniforms.loc_normalMapStrength = glGetUniformLocation(gpuparam.shader_program->GetProgram(), "normalMappingStrength");
  gpuparam.uniforms.loc_textureScaling = glGetUniformLocation(gpuparam.shader_program->GetProgram(), "textureScaling");

  // On recupere la position des attributs generiques de sommets 
  gpuparam.uniforms.loc_position = glGetAttribLocation(gpuparam.shader_program->GetProgram(), "position_vert");
  gpuparam.uniforms.loc_normal = glGetAttribLocation(gpuparam.shader_program->GetProgram(), "normal_vert");
  gpuparam.uniforms.loc_color = glGetAttribLocation(gpuparam.shader_program->GetProgram(), "color_vert");
  gpuparam.uniforms.loc_UV = glGetAttribLocation(gpuparam.shader_program->GetProgram(), "UV_vert");
  gpuparam.uniforms.loc_inst_pos = glGetAttribLocation(gpuparam.shader_program->GetProgram(), "inst_pos");
  gpuparam.uniforms.loc_inst_scl = glGetAttribLocation(gpuparam.shader_program->GetProgram(), "inst_scl");
  gpuparam.uniforms.loc_inst_rot = glGetAttribLocation(gpuparam.shader_program->GetProgram(), "inst_rot");

  gpuparam.uniforms.loc_texture = glGetUniformLocation(gpuparam.shader_program->GetProgram(), "textmat");

  // Get Subroutine Index
  gpuparam.uniforms.phongIndex = glGetSubroutineIndex(gpuparam.shader_program->GetProgram(), GL_FRAGMENT_SHADER, "phongModel");
  gpuparam.uniforms.phongVCIndex = glGetSubroutineIndex(gpuparam.shader_program->GetProgram(), GL_FRAGMENT_SHADER, "phongVertexColorModel");
  gpuparam.uniforms.normalIndex = glGetSubroutineIndex(gpuparam.shader_program->GetProgram(), GL_FRAGMENT_SHADER, "normalColorModel");
  gpuparam.uniforms.gridIndex = glGetSubroutineIndex(gpuparam.shader_program->GetProgram(), GL_FRAGMENT_SHADER, "gridModel");
  gpuparam.uniforms.ignIndex = glGetSubroutineIndex(gpuparam.shader_program->GetProgram(), GL_FRAGMENT_SHADER, "ignModel");
  gpuparam.uniforms.goochIndex = glGetSubroutineIndex(gpuparam.shader_program->GetProgram(), GL_FRAGMENT_SHADER, "goochModel");
  gpuparam.uniforms.terrain3DIndex = glGetSubroutineIndex(gpuparam.shader_program->GetProgram(), GL_FRAGMENT_SHADER, "terrain3DModel");
  gpuparam.uniforms.UVTextureModel = glGetSubroutineIndex(gpuparam.shader_program->GetProgram(), GL_FRAGMENT_SHADER, "UVTextureModel");
  gpuparam.uniforms.bumpNormalIndex = glGetSubroutineIndex(gpuparam.shader_program->GetProgram(), GL_FRAGMENT_SHADER, "bumpNormal");
  gpuparam.uniforms.ThinSnowIndex = glGetSubroutineIndex(gpuparam.shader_program->GetProgram(), GL_FRAGMENT_SHADER, "ThinSnow");
  gpuparam.uniforms.GreyCoolWarmIndex = glGetSubroutineIndex(gpuparam.shader_program->GetProgram(), GL_FRAGMENT_SHADER, "GreyCoolWarm");
  gpuparam.uniforms.AxelTerrainIndex = glGetSubroutineIndex(gpuparam.shader_program->GetProgram(), GL_FRAGMENT_SHADER, "AxelTerrain");
  gpuparam.uniforms.WireframeIndex = glGetSubroutineIndex(gpuparam.shader_program->GetProgram(), GL_FRAGMENT_SHADER, "Wireframe");
  gpuparam.uniforms.WireframeAspectIndex = glGetSubroutineIndex(gpuparam.shader_program->GetProgram(), GL_FRAGMENT_SHADER, "WireframeAspect");
  gpuparam.uniforms.TriPlanarIndex = glGetSubroutineIndex(gpuparam.shader_program->GetProgram(), GL_FRAGMENT_SHADER, "EricTriPlanar");

  // Create the program to draw the background
  gpubackground = new MayaShader(pPath + QString("/shaders/libs/background"));

  // Create the VAO used when drawing the background
  glGenVertexArrays(1, &backgroundVAO);
}

/*!
\brief Initialize camera transforms.
*/
void MayaWidget::InitCamera()
{
  // Create uniform buffer and store camera data
  glGenBuffers(1, &gpuparam.buffer_transform);
  glBindBuffer(GL_UNIFORM_BUFFER, gpuparam.buffer_transform);
  glBufferData(GL_UNIFORM_BUFFER, sizeof(gpuparam.transform), &gpuparam.transform, GL_DYNAMIC_DRAW);
  glBindBufferBase(GL_UNIFORM_BUFFER, 0, gpuparam.buffer_transform);
}

/*!
\brief Initialize the terrain texture.

Basically load four texture images corresponding to bedrock, sand, snow, and grass, plus some additional textures with normal maps.
*/
void MayaWidget::InitTerrainTextures()
{
  QString pPath = QString::fromStdString(std::string(SOLUTION_DIR));

  if (pPath.isEmpty()) std::cout << "MayaWidget::InitTerrainTextures() : SOLUTION_DIR undefined" << std::endl;

  // Load textures
  QImage rock(pPath + QString("/shaders/rock.png"));
  QImage stone(pPath + QString("/shaders/grass.png"));
  QImage tex3(pPath + QString("/shaders/sand.png"));
  QImage tex4(pPath + QString("/shaders/snow.png"));
  if (rock.isNull() || stone.isNull() || tex3.isNull() || tex4.isNull())
    std::cout << "Errors: loading images ! #1" << std::endl;

  // Axel textures (with normal maps)
  QImage niceRock(pPath + QString("/shaders/limestone-rock-albedo.png"));
  QImage niceRockNormal(pPath + QString("/shaders/limestone-rock-normal.png"));
  QImage niceRock2(pPath + QString("/shaders/StonesBeach-albedo.png"));
  QImage niceRock2Normal(pPath + QString("/shaders/StonesBeach-normal.png"));
  if (niceRock.isNull() || niceRockNormal.isNull() || niceRock2.isNull() || niceRock2Normal.isNull())
    std::cout << "Errors: loading images ! #2" << std::endl;

  // Convert to GPU format
  rock.convertTo(QImage::Format::Format_RGBA8888);
  stone.convertTo(QImage::Format::Format_RGBA8888);
  tex3.convertTo(QImage::Format::Format_RGBA8888);
  tex4.convertTo(QImage::Format::Format_RGBA8888);

  niceRock.convertTo(QImage::Format::Format_RGBA8888);
  niceRockNormal.convertTo(QImage::Format::Format_RGBA8888);
  niceRock2.convertTo(QImage::Format::Format_RGBA8888);
  niceRock2Normal.convertTo(QImage::Format::Format_RGBA8888);

  // Create the Textures
  glGenTextures(8, gpuparam.buffer_textureT);

  // Typical Texture Generation Using Data From The Bitmap
  glBindTexture(GL_TEXTURE_2D, gpuparam.buffer_textureT[0]);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexImage2D(GL_TEXTURE_2D, 0, 3, rock.width(), rock.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, rock.bits());

  glBindTexture(GL_TEXTURE_2D, gpuparam.buffer_textureT[1]);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexImage2D(GL_TEXTURE_2D, 0, 3, stone.width(), stone.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, stone.bits());

  glBindTexture(GL_TEXTURE_2D, gpuparam.buffer_textureT[2]);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexImage2D(GL_TEXTURE_2D, 0, 3, tex3.width(), tex3.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, tex3.bits());

  glBindTexture(GL_TEXTURE_2D, gpuparam.buffer_textureT[3]);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexImage2D(GL_TEXTURE_2D, 0, 3, tex4.width(), tex4.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, tex4.bits());

  // Textures Axel
  // Mipmap generation is required in order to avoid pixel sparkling in the distance.
  // OpenGL automatically choose which mipmap to use for each fragment based on the fragment ratio on the screen.
  glBindTexture(GL_TEXTURE_2D, gpuparam.buffer_textureT[4]);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexImage2D(GL_TEXTURE_2D, 0, 3, niceRock.width(), niceRock.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, niceRock.bits());
  glGenerateMipmap(GL_TEXTURE_2D);

  glBindTexture(GL_TEXTURE_2D, gpuparam.buffer_textureT[5]);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexImage2D(GL_TEXTURE_2D, 0, 3, niceRockNormal.width(), niceRockNormal.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, niceRockNormal.bits());
  glGenerateMipmap(GL_TEXTURE_2D);

  glBindTexture(GL_TEXTURE_2D, gpuparam.buffer_textureT[6]);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexImage2D(GL_TEXTURE_2D, 0, 3, niceRock2.width(), niceRock2.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, niceRock2.bits());
  glGenerateMipmap(GL_TEXTURE_2D);

  glBindTexture(GL_TEXTURE_2D, gpuparam.buffer_textureT[7]);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexImage2D(GL_TEXTURE_2D, 0, 3, niceRock2Normal.width(), niceRock2Normal.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, niceRock2Normal.bits());
  glGenerateMipmap(GL_TEXTURE_2D);

  glBindTexture(GL_TEXTURE_2D, 0);
}

/*!
\brief Initialize skybox (envmap) textures.
*/
void MayaWidget::InitEnvTextures()
{
  QString pPath = System::GetResource("ARCHESLIBDIR");
  if (pPath.isEmpty())
    std::cout << "MayaWidget::InitTerrainTextures() : environment variable ARCHESLIBDIR undefined" << std::endl;

  // Load textures
  QImage left(pPath + QString("/LibMaya/Shaders/storm-left.png"));
  QImage right(pPath + QString("/LibMaya/Shaders/storm-right.png"));
  QImage up(pPath + QString("/LibMaya/Shaders/storm-up.png"));
  QImage down(pPath + QString("/LibMaya/Shaders/storm-down.png"));
  QImage front(pPath + QString("/LibMaya/Shaders/storm-front.png"));
  QImage back(pPath + QString("/LibMaya/Shaders/storm-back.png"));

  if (left.isNull() || up.isNull() || right.isNull() || down.isNull() || front.isNull() || back.isNull())
    std::cout << "Errors : loading env images !" << std::endl;

  // Convert to GPU format
  left.convertTo(QImage::Format::Format_RGBA8888);
  right.convertTo(QImage::Format::Format_RGBA8888);
  up.convertTo(QImage::Format::Format_RGBA8888);
  down.convertTo(QImage::Format::Format_RGBA8888);
  front.convertTo(QImage::Format::Format_RGBA8888);
  back.convertTo(QImage::Format::Format_RGBA8888);

  // Create the Textures
  glGenTextures(6, gpuparam.buffer_textureE);
  // ENV
  // Typical Texture Generation Using Data From The Bitmap
  glBindTexture(GL_TEXTURE_2D, gpuparam.buffer_textureE[0]);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexImage2D(GL_TEXTURE_2D, 0, 3, left.width(), left.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, left.bits());

  glBindTexture(GL_TEXTURE_2D, gpuparam.buffer_textureE[1]);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexImage2D(GL_TEXTURE_2D, 0, 3, right.width(), right.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, right.bits());

  glBindTexture(GL_TEXTURE_2D, gpuparam.buffer_textureE[2]);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexImage2D(GL_TEXTURE_2D, 0, 3, up.width(), up.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, up.bits());

  glBindTexture(GL_TEXTURE_2D, gpuparam.buffer_textureE[3]);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexImage2D(GL_TEXTURE_2D, 0, 3, down.width(), down.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, down.bits());

  glBindTexture(GL_TEXTURE_2D, gpuparam.buffer_textureE[4]);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexImage2D(GL_TEXTURE_2D, 0, 3, front.width(), front.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, front.bits());

  glBindTexture(GL_TEXTURE_2D, gpuparam.buffer_textureE[5]);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexImage2D(GL_TEXTURE_2D, 0, 3, back.width(), back.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, back.bits());

  glBindTexture(GL_TEXTURE_2D, 0);
}

/*!
\brief Initialize anchor shape.
*/
void MayaWidget::InitAnchorShape()
{
  MayaGeometry mg = MayaGeometry("Anchor", Mesh(Sphere(2.0), 14), MayaMaterial::None);
  MayaGeometrySet mgs(mg, FrameScaled::Id);
  anchor = MayaGpuSet(mgs);
}

/*!
\brief Initialize light shape.
*/
void MayaWidget::InitLightShape()
{
  MayaGeometry mg("Light", Mesh(Box(1.0)), MayaMaterial(ShaderPhong, Color(0.8, 0.7, 0.2), Color(0.1, 0.1, 0.1), Color(0.1, 0.1, 0.1)));

  MayaGeometrySet mgs(mg, FrameScaled::Id);
  lightshape = MayaGpuSet(mgs);
}

/*!
\brief Initialize transparent sphere.
*/
void MayaWidget::InitSphereShape()
{
  MayaGeometry mg = MayaGeometry("Transparent Sphere", Mesh(Sphere(1.0), 14), MayaMaterial(ShaderPhong, Color(0.8, 0.7, 0.2, 0.5), Color(0.1, 0.1, 0.1), Color(0.1, 0.1, 0.1)));

  MayaGeometrySet mgs(mg, FrameScaled::Id);
  sphereshape = MayaGpuSet(mgs);
}
