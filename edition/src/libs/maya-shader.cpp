#include <QtCore/QDebug>
#include <ostream>

#include "libs/mayashader.h"

/*!
\class MayaShader mayashader.h
\brief A helper class for loading and creating program shaders.
*/

/*!
\brief Create a MayaShader.

Shaders must be name : [baseFilename]_[x]s.glsl where x = v || f || g

\param baseFilename Base filename of the shaders.
*/
MayaShader::MayaShader(const QString& baseFilename)
{
  //Load the shaders
  GLuint VS = LoadShaderFromFile(baseFilename + "_vs.glsl", GL_VERTEX_SHADER);
  GLuint FS = LoadShaderFromFile(baseFilename + "_fs.glsl", GL_FRAGMENT_SHADER);

  GLuint GS = 0;
  if (QFile::exists(baseFilename + "_gs.glsl"))
    GS = LoadShaderFromFile(baseFilename + "_gs.glsl", GL_GEOMETRY_SHADER);

  name = glCreateProgram();

  //Attach and link
  glAttachShader(name, VS);
  glAttachShader(name, FS);

  if (GS != 0)
    glAttachShader(name, GS);

  glLinkProgram(name);

  GLint status;

  //Check that the program was linked successfully
  glGetProgramiv(name, GL_LINK_STATUS, &status);

  if (status != GL_TRUE)
  {
    GLint logSize;

    //Get the log size then the log
    glGetProgramiv(name, GL_INFO_LOG_LENGTH, &logSize);

    GLchar* log = new GLchar[logSize];
    glGetProgramInfoLog(name, logSize, nullptr, log);
    qDebug() << "Program : " << baseFilename << " has link error.";
    qDebug() << qPrintable(QString(log));

    delete[] log;
  }

  // Delete shaders because there are not useful anymore
  glDeleteShader(VS);
  glDeleteShader(FS);

  if (GS != 0)
    glDeleteShader(GS);
}

/*!
\brief Destructor.

Removes the program on the gpu.
*/
MayaShader::~MayaShader()
{
  glDeleteProgram(name);
}

/*!
\brief Load a shader from a file and optionnaly append some code to it
\param fileName Path to the shader file on disk
\param shaderType Type of the shader to load
\return Name of the newly created shader
*/
GLuint MayaShader::LoadShaderFromFile(const QString& fileName, GLenum shaderType)
{
  //Open the shader file
  QFile shaderFile(fileName);
  if (!shaderFile.open(QFile::ReadOnly))
  {
    qDebug() << "Can't open shader " << qPrintable(fileName);
    return 0;
  }

  QByteArray shaderByteArray = shaderFile.readAll();

  const char* shaderChar = shaderByteArray.data();

  //Create the shader here
  GLuint shader = glCreateShader(shaderType);

  glShaderSource(shader, 1, &shaderChar, nullptr);

  glCompileShader(shader);

  //CHeck the compilation status of the shader
  GLint compileStatus = 0;
  glGetShaderiv(shader, GL_COMPILE_STATUS, &compileStatus);

  if (compileStatus != GL_TRUE)
  {
    GLint infoLogSize = 0;
    //display the error log
    glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infoLogSize);
    GLchar* log = new GLchar[infoLogSize];
    glGetShaderInfoLog(shader, infoLogSize, nullptr, log);
    qDebug() << "Shader : " << fileName << "\n";

    qDebug() << "Error : " << "\n";
    qDebug() << log << "\n";

    delete[] log;
  }

  return shader;
}