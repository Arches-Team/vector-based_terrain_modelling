#pragma once

#include <QtCore/QtCore>

#include "libs/mayashader.h"
#include "libs/color.h"
#include "libs/evectorfloat.h"

class ColorFloat
{
protected:
    float c[4] = { 0.0,0.0,0.0,1.0 };
public:
    ColorFloat() { ColorFloat(0.f, 0.f, 0.f, 0.f); }
    ColorFloat(float r, float g, float b, float a = 1.f)
    {
        c[0] = r;
        c[1] = g;
        c[2] = b;
        c[3] = a;
    }

    explicit ColorFloat(const Color& color)
    {
        for (int i = 0; i < 4; ++i)
            c[i] = static_cast<float>(color[i]);
    }

    float& operator[](int i) { return c[i]; }
};

class MayaSimpleRendererColors
{
protected:
    MayaShader program; //!< Program used to draw line.
    GLuint VBO; //!< Where we store our vertex.
    GLuint VAO; //!< VAO.
    GLuint CBO; //!< Where we store our colors.
    VectorFloat* map; //!< CPU map of our VBO.
    ColorFloat* colorMap; //!< CPU map of our CBO.

    GLuint initSizeLines; //!< Initial size of our lines vector.
    GLuint currentSizeLines; //!< how many vertex there is in our VBO

    GLenum lineType; //!< type of line to draw : GL_LINES, GL_LINE_LOOP

    const int sizeOfColorFloat = 4 * sizeof(float);

    bool depthTest = true;
    float lineWidth = 1.;
public:
    MayaSimpleRendererColors();
    MayaSimpleRendererColors(int, int, GLenum);
    MayaSimpleRendererColors(QVector<VectorFloat>, QVector<ColorFloat>, GLenum);
    ~MayaSimpleRendererColors()
    {
        DeleteBuffers();
    }

    void DeleteBuffers();
    void Update(QVector<VectorFloat> lines, QVector<ColorFloat> colors);
    void Draw();

    void setDepthTest(bool d) { depthTest = d; }
    void setLineWidth(float l) { lineWidth = l; }
};
