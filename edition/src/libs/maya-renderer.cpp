#include "libs/GL/glew.h"
#include "libs/mayarender.h"

#include "libs/cpu.h"

/*!
\brief Build a simple renderer.
*/
MayaSimpleRenderer::MayaSimpleRenderer() :lineType(lineType), color(color), program(QString::fromStdString(std::string(SOLUTION_DIR) + "/shaders/libs/simple"))// Load the program for the simple renderer
{
  QVector<Vector> lines;
  lines.append(Vector(0., 0., 0.));
  lines.append(Vector::Z);
  Color color(1., 0., 0.);
  GLenum lineType = GL_LINES;

  initSizeLines = lines.size();
  currentSizeLines = lines.size();

  //Build the VBO and map it
  glGenBuffers(1, &VBO);
  glBindBuffer(GL_ARRAY_BUFFER, VBO);
  glBufferStorage(GL_ARRAY_BUFFER, lines.size() * sizeof(VectorFloat), nullptr, GL_MAP_WRITE_BIT | GL_MAP_READ_BIT | GL_MAP_PERSISTENT_BIT);
  map = (VectorFloat*)glMapBufferRange(GL_ARRAY_BUFFER, 0, lines.size() * sizeof(VectorFloat), GL_MAP_WRITE_BIT | GL_MAP_READ_BIT | GL_MAP_PERSISTENT_BIT);

  //Upload data and signal OpenGL
  for (int i = 0; i < lines.size(); i++)
  {
    map[i] = VectorFloat(lines[i][0], lines[i][1], lines[i][2]);
  }

  glMemoryBarrier(GL_CLIENT_MAPPED_BUFFER_BARRIER_BIT);

  //Create the VAO and enable the position attrib.
  glGenVertexArrays(1, &VAO);
  glBindVertexArray(VAO);

  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(VectorFloat), (void*)offsetof(VectorFloat, x));
  glEnableVertexAttribArray(0);

  glBindBuffer(GL_ARRAY_BUFFER, 0);
}

/*!
\brief Build a simple renderer.
\param lines All points on the line(s).
\param color Color of the lines.
\param lineType Type of line to draw.
*/
MayaSimpleRenderer::MayaSimpleRenderer(const QVector<Vector>& lines, const Color& color, GLenum lineType) : lineType(lineType), color(color), program(QString::fromStdString(std::string(SOLUTION_DIR) + "/shaders/libs/simple")) // Load the program for the simple renderer
{
  initSizeLines = lines.size();
  currentSizeLines = lines.size();

  //Build the VBO and map it
  glGenBuffers(1, &VBO);
  glBindBuffer(GL_ARRAY_BUFFER, VBO);
  glBufferStorage(GL_ARRAY_BUFFER, lines.size() * sizeof(VectorFloat), nullptr, GL_MAP_WRITE_BIT | GL_MAP_READ_BIT | GL_MAP_PERSISTENT_BIT);
  map = (VectorFloat*)glMapBufferRange(GL_ARRAY_BUFFER, 0, lines.size() * sizeof(VectorFloat), GL_MAP_WRITE_BIT | GL_MAP_READ_BIT | GL_MAP_PERSISTENT_BIT);

  //Upload data and signal OpenGL
  for (int i = 0; i < lines.size(); i++)
  {
    map[i] = VectorFloat(lines[i]);
  }

  glMemoryBarrier(GL_CLIENT_MAPPED_BUFFER_BARRIER_BIT);

  //Create the VAO and enable the position attrib.
  glGenVertexArrays(1, &VAO);
  glBindVertexArray(VAO);

  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(VectorFloat), (void*)offsetof(VectorFloat, x));
  glEnableVertexAttribArray(0);

  glBindBuffer(GL_ARRAY_BUFFER, 0);
}

/*!
\brief Build a simple renderer.
\param lines All points on the line(s).
\param n Number of points.
\param color Color of the lines.
\param lineType Type of line to draw.
*/
MayaSimpleRenderer::MayaSimpleRenderer(Vector* lines, int n, const Color& color, GLenum lineType) : lineType(lineType), color(color), program(QString::fromStdString(std::string(SOLUTION_DIR) + "/shaders/libs/simple"))// Load the program for the simple renderer
{
  initSizeLines = n;
  currentSizeLines = n;

  //Build the VBO and map it
  glGenBuffers(1, &VBO);
  glBindBuffer(GL_ARRAY_BUFFER, VBO);
  glBufferStorage(GL_ARRAY_BUFFER, n * sizeof(VectorFloat), nullptr, GL_MAP_WRITE_BIT | GL_MAP_READ_BIT | GL_MAP_PERSISTENT_BIT);
  map = (VectorFloat*)glMapBufferRange(GL_ARRAY_BUFFER, 0, n * sizeof(VectorFloat), GL_MAP_WRITE_BIT | GL_MAP_READ_BIT | GL_MAP_PERSISTENT_BIT);

  //Upload data and signal OpenGL
  for (int i = 0; i < n; i++)
  {
    map[i] = VectorFloat(lines[i]);
  }

  glMemoryBarrier(GL_CLIENT_MAPPED_BUFFER_BARRIER_BIT);

  //Create the VAO and enable the position attrib.
  glGenVertexArrays(1, &VAO);
  glBindVertexArray(VAO);

  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(VectorFloat), (void*)offsetof(VectorFloat, x));
  glEnableVertexAttribArray(0);

  glBindBuffer(GL_ARRAY_BUFFER, 0);
}
/*!
\brief Build a simple renderer.
\param lines All points on the line(s).
\param n Number of points.
\param color Color of the lines.
\param lineType Type of line to draw.
*/
MayaSimpleRenderer::MayaSimpleRenderer(Segment* lines, int n, const Color& color, GLenum lineType) : lineType(lineType), color(color), program(QString::fromStdString(std::string(SOLUTION_DIR) + "/shaders/libs/simple"))// Load the program for the simple renderer
{
  initSizeLines = 2 * n;
  currentSizeLines = 2 * n;

  //Build the VBO and map it
  glGenBuffers(1, &VBO);
  glBindBuffer(GL_ARRAY_BUFFER, VBO);
  glBufferStorage(GL_ARRAY_BUFFER, 2 * n * sizeof(VectorFloat), nullptr, GL_MAP_WRITE_BIT | GL_MAP_READ_BIT | GL_MAP_PERSISTENT_BIT);
  map = (VectorFloat*)glMapBufferRange(GL_ARRAY_BUFFER, 0, 2 * n * sizeof(VectorFloat), GL_MAP_WRITE_BIT | GL_MAP_READ_BIT | GL_MAP_PERSISTENT_BIT);

  //Upload data and signal OpenGL
  for (int i = 0; i < n; i++)
  {
    map[2 * i] = VectorFloat(lines[i].Vertex(0));
    map[2 * i + 1] = VectorFloat(lines[i].Vertex(1));
  }

  glMemoryBarrier(GL_CLIENT_MAPPED_BUFFER_BARRIER_BIT);

  //Create the VAO and enable the position attrib.
  glGenVertexArrays(1, &VAO);
  glBindVertexArray(VAO);

  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(VectorFloat), (void*)offsetof(VectorFloat, x));
  glEnableVertexAttribArray(0);

  glBindBuffer(GL_ARRAY_BUFFER, 0);
}
/*!
\brief Update the system with a set of lines.
\param lines Set of lines.
*/
void MayaSimpleRenderer::Update(QVector<Vector>& lines)
{
  //We check that lines.size if not > to our buffer (let's avoid segfault)
  currentSizeLines = min((unsigned int)lines.size(), initSizeLines);

  //Upload data and signal OpenGL
  for (int i = 0; i < lines.size(); i++)
  {
    map[i] = VectorFloat(lines[i]);
  }

  glMemoryBarrier(GL_CLIENT_MAPPED_BUFFER_BARRIER_BIT);
}

/*!
\brief Draw the renderer.
*/
void MayaSimpleRenderer::Draw()
{
  if (depthTest)
  {
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
  }
  else
    glDisable(GL_DEPTH_TEST);

  glLineWidth(lineWidth);

  //Bind program and VAO
  glUseProgram(program.GetProgram());
  glBindVertexArray(VAO);
  glEnableVertexAttribArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, VBO);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(VectorFloat), (void*)offsetof(VectorFloat, x));

  // Should update when camera can create it's own matrix

  //Get matrices and send it to the gpu
  float viewMatTemp[16];

  glGetFloatv(GL_MODELVIEW_MATRIX, viewMatTemp);
  glUniformMatrix4fv(0, 1, GL_FALSE, viewMatTemp);

  float projMatTemp[16];

  glGetFloatv(GL_PROJECTION_MATRIX, projMatTemp);
  glUniformMatrix4fv(1, 1, GL_FALSE, projMatTemp);

  //Update color and draw
  glUniform4f(2, color[0], color[1], color[2], color[3]);

  glDrawArrays(lineType, 0, currentSizeLines);
  glBindVertexArray(0);
  glUseProgram(0);
  glDisable(GL_DEPTH_TEST);
}

/*!
\brief Delete all buffers
*/
void MayaSimpleRenderer::DeleteBuffers()
{
  if (VBO != 0) {
    glDeleteBuffers(1, &VBO); VBO = 0;
  }
  if (VAO != 0) {
    glDeleteVertexArrays(1, &VAO); VAO = 0;
  }
  // Should delete program here (possible bug)
}
