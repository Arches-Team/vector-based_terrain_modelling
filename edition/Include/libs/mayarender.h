#pragma once

#include <QtCore/QtCore>
#include "libs/evector.h"
#include "libs/segment.h"
#include "libs/mayashader.h"
#include "libs/color.h"

#include "libs/evectorfloat.h"

class MayaSimpleRenderer
{
protected:
  MayaShader program; //!< Program used to draw line.
  GLuint VBO; //!< Where we store our vertex.
  GLuint VAO; //!< VAO.
  Color color = Color::White; //!< Color of our lines.
  VectorFloat* map = nullptr; //!< CPU map of our VBO.

  GLuint initSizeLines; //!< Initial size of our lines vector.
  GLuint currentSizeLines; //!< How many vertex there is in our VBO.

  GLenum lineType; //!< type of line to draw: lines or line loops.

  bool depthTest = true;
  float lineWidth = 1.0; //!< Width.
public:
  MayaSimpleRenderer(); 
  MayaSimpleRenderer(const QVector<Vector>&, const Color&, GLenum);
  MayaSimpleRenderer(Segment*, int, const Color&, GLenum);
  MayaSimpleRenderer(Vector*, int, const Color&, GLenum);

  void DeleteBuffers();
  void Update(QVector<Vector>&);
  void Draw();

  void setDepthTest(bool d) { depthTest = d; }
  void setLineWidth(float l) { lineWidth = l; }
};

