#pragma once

#ifdef _MSC_VER
#include "libs/GL/glew.h"
#else
#include "GL/glew.h"
#endif

#include <QtCore/QtCore>

class MayaShader
{
protected:
  GLuint name = 0; //!< Program name on the gpu
public:
  MayaShader(const QString&);
  ~MayaShader();

  //! Get the name of the program on the GPU
  inline GLuint GetProgram() { return name; }

  static GLuint LoadShaderFromFile(const QString& fileName, GLenum shaderType);

};

