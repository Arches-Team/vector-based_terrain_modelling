#include "MayaSimpleRendererColors.h"

#include "cpu.h"

/*!
\brief Build a simple renderer.
*/
MayaSimpleRendererColors::MayaSimpleRendererColors() : program(QString::fromStdString(std::string(SOLUTION_DIR) + "/Shaders/simple_colors"))// Load the program for the simple renderer
{
    QVector<Vector> lines;
    lines.append(Vector(0., 0., 0.));
    lines.append(Vector::Z);
    QVector<Color> colors;
    colors.append(Color(1., 0., 0., 1.));
    colors.append(Color(0., 1., 0., 1.));
    GLenum lineType = GL_LINES;

    initSizeLines = lines.size();
    currentSizeLines = lines.size();
    this->lineType = lineType;

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

    // Create the CBO and map it
    glGenBuffers(1, &CBO);
    glBindBuffer(GL_ARRAY_BUFFER, CBO);
    glBufferStorage(GL_ARRAY_BUFFER, colors.size() * sizeOfColorFloat, nullptr, GL_MAP_WRITE_BIT | GL_MAP_READ_BIT | GL_MAP_PERSISTENT_BIT);
    colorMap = (ColorFloat*)glMapBufferRange(GL_ARRAY_BUFFER, 0, colors.size() * sizeOfColorFloat, GL_MAP_WRITE_BIT | GL_MAP_READ_BIT | GL_MAP_PERSISTENT_BIT);

    // Upload color data and signal OpenGL
    for (int i = 0; i < colors.size(); i++)
    {
        colorMap[i] = ColorFloat(colors[i]);
    }

    glMemoryBarrier(GL_CLIENT_MAPPED_BUFFER_BARRIER_BIT);

    //Create the VAO and enable the position and color attribs.
    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(VectorFloat), (void*)offsetof(VectorFloat, x));
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, CBO);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeOfColorFloat, (void*)0);
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

MayaSimpleRendererColors::MayaSimpleRendererColors(int linesSize, int colorsSize, GLenum lineType) : program(QString::fromStdString(std::string(SOLUTION_DIR) + "/Shaders/simple_colors"))// Load the program for the simple renderer
{
    this->lineType = lineType;

    //Build the VBO and map it
    glGenBuffers(1, &VBO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferStorage(GL_ARRAY_BUFFER, linesSize * sizeof(VectorFloat), nullptr, GL_MAP_WRITE_BIT | GL_MAP_READ_BIT | GL_MAP_PERSISTENT_BIT);
    map = (VectorFloat*)glMapBufferRange(GL_ARRAY_BUFFER, 0, linesSize * sizeof(VectorFloat), GL_MAP_WRITE_BIT | GL_MAP_READ_BIT | GL_MAP_PERSISTENT_BIT);

    glMemoryBarrier(GL_CLIENT_MAPPED_BUFFER_BARRIER_BIT);

    // Create the CBO and map it
    glGenBuffers(1, &CBO);
    glBindBuffer(GL_ARRAY_BUFFER, CBO);
    glBufferStorage(GL_ARRAY_BUFFER, colorsSize * sizeof(ColorFloat), nullptr, GL_MAP_WRITE_BIT | GL_MAP_READ_BIT | GL_MAP_PERSISTENT_BIT);
    colorMap = (ColorFloat*)glMapBufferRange(GL_ARRAY_BUFFER, 0, colorsSize * sizeof(ColorFloat), GL_MAP_WRITE_BIT | GL_MAP_READ_BIT | GL_MAP_PERSISTENT_BIT);

    glMemoryBarrier(GL_CLIENT_MAPPED_BUFFER_BARRIER_BIT);

    //Create the VAO and enable the position and color attribs.
    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(VectorFloat), (void*)offsetof(VectorFloat, x));
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, CBO);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeOfColorFloat, (void*)0);
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

/*!
\brief Build a simple renderer.
\param lines All points on the line(s).
\param colors Colors of the lines.
\param lineType Type of line to draw.
*/
MayaSimpleRendererColors::MayaSimpleRendererColors(QVector<VectorFloat> lines, QVector<ColorFloat> colors, GLenum lineType) : program(QString::fromStdString(std::string(SOLUTION_DIR) + "/Shaders/simple_colors"))// Load the program for the simple renderer
{
    initSizeLines = lines.size();
    currentSizeLines = lines.size();
    this->lineType = lineType;

    //Build the VBO and map it
    glGenBuffers(1, &VBO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferStorage(GL_ARRAY_BUFFER, lines.size() * sizeof(VectorFloat), nullptr, GL_MAP_WRITE_BIT | GL_MAP_READ_BIT | GL_MAP_PERSISTENT_BIT);
    map = (VectorFloat*)glMapBufferRange(GL_ARRAY_BUFFER, 0, lines.size() * sizeof(VectorFloat), GL_MAP_WRITE_BIT | GL_MAP_READ_BIT | GL_MAP_PERSISTENT_BIT);

    //Upload data and signal OpenGL
    for (int i = 0; i < lines.size(); i++)
    {
        map[i] = lines[i];
    }

    glMemoryBarrier(GL_CLIENT_MAPPED_BUFFER_BARRIER_BIT);

    // Create the CBO and map it
    glGenBuffers(1, &CBO);
    glBindBuffer(GL_ARRAY_BUFFER, CBO);
    glBufferStorage(GL_ARRAY_BUFFER, colors.size() * sizeof(ColorFloat), nullptr, GL_MAP_WRITE_BIT | GL_MAP_READ_BIT | GL_MAP_PERSISTENT_BIT);
    colorMap = (ColorFloat*)glMapBufferRange(GL_ARRAY_BUFFER, 0, colors.size() * sizeof(ColorFloat), GL_MAP_WRITE_BIT | GL_MAP_READ_BIT | GL_MAP_PERSISTENT_BIT);

    // Upload color data and signal OpenGL
    for (int i = 0; i < colors.size(); i++)
    {
        colorMap[i] = colors[i];
    }

    glMemoryBarrier(GL_CLIENT_MAPPED_BUFFER_BARRIER_BIT);

    //Create the VAO and enable the position and color attribs.
    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(VectorFloat), (void*)offsetof(VectorFloat, x));
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, CBO);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeOfColorFloat, (void*)0);
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

/*!
\brief Update the system with a set of lines and colors.
\param lines Set of lines.
\param colors Set of colors.
*/
void MayaSimpleRendererColors::Update(QVector<VectorFloat> lines, QVector<ColorFloat> colors)
{
    /*glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferSubData(GL_ARRAY_BUFFER, 0, lines.size() * sizeof(VectorFloat), lines.constData());

    glBindBuffer(GL_ARRAY_BUFFER, CBO);
    glBufferSubData(GL_ARRAY_BUFFER, 0, colors.size() * sizeof(ColorFloat), colors.constData());

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    GLenum err;
    while ((err = glGetError()) != GL_NO_ERROR) {
        std::cerr << err << std::endl;
    }*/

    /*map = lines.data();
    colorMap = colors.data();*/

    //Upload data and signal OpenGL
    //for (int i = 0; i < lines.size(); i++)
    //{
    //    map[i] = lines[i];
    //}

    //// Upload color data and signal OpenGL
    //for (int i = 0; i < colors.size(); i++)
    //{
    //    colorMap[i] = colors[i];
    //}

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    map = (VectorFloat*)glMapBufferRange(GL_ARRAY_BUFFER, 0, lines.size() * sizeof(VectorFloat), GL_MAP_WRITE_BIT | GL_MAP_READ_BIT | GL_MAP_PERSISTENT_BIT);

    //Upload data and signal OpenGL
    for (int i = 0; i < lines.size(); i++)
    {
        map[i] = lines[i];
    }

    glBindBuffer(GL_ARRAY_BUFFER, CBO);
    colorMap = (ColorFloat*)glMapBufferRange(GL_ARRAY_BUFFER, 0, colors.size() * sizeof(ColorFloat), GL_MAP_WRITE_BIT | GL_MAP_READ_BIT | GL_MAP_PERSISTENT_BIT);

    // Upload color data and signal OpenGL
    for (int i = 0; i < colors.size(); i++)
    {
        colorMap[i] = colors[i];
    }

    glMemoryBarrier(GL_CLIENT_MAPPED_BUFFER_BARRIER_BIT);
}

/*!
\brief Draw the renderer.
*/
void MayaSimpleRendererColors::Draw()
{
    if (depthTest)
    {
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LEQUAL);
    }
    else
        glDisable(GL_DEPTH_TEST);

    glEnable(GL_BLEND);

    glLineWidth(lineWidth);

    //Bind program and VAO
    glUseProgram(program.GetProgram());
    glBindVertexArray(VAO);
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(VectorFloat), (void*)offsetof(VectorFloat, x));

    glBindBuffer(GL_ARRAY_BUFFER, CBO);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(ColorFloat), (void*)0);

    // Should update when camera can create it's own matrix

    //Get matrices and send it to the gpu
    float viewMatTemp[16];

    glGetFloatv(GL_MODELVIEW_MATRIX, viewMatTemp);
    glUniformMatrix4fv(2, 1, GL_FALSE, viewMatTemp);

    float projMatTemp[16];

    glGetFloatv(GL_PROJECTION_MATRIX, projMatTemp);
    glUniformMatrix4fv(3, 1, GL_FALSE, projMatTemp);

    // Draw
    glDrawArrays(lineType, 0, currentSizeLines);
    glBindVertexArray(0);
    glUseProgram(0);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);
}

/*!
\brief Delete all buffers
*/
void MayaSimpleRendererColors::DeleteBuffers()
{
    if (VBO != 0) {
        glDeleteBuffers(1, &VBO); VBO = 0;
    }
    if (VAO != 0) {
        glDeleteVertexArrays(1, &VAO); VAO = 0;
    }
    if (CBO != 0) {
        glDeleteBuffers(1, &CBO); CBO = 0;
    }
    // Should delete program here (possible bug)
}