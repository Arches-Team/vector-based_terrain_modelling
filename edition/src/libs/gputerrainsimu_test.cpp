#include "libs/gputerrainsimu.h"


#include <QtCore/QFileInfo>

/*!
\class GPUHydraulicErosion gputerrainsimu.h
\brief Hydraulic erosion simulation on the GPU.

The algorithm uses a grid-based approach with no race condition as compared to the other one.
The implementation is based on https://www.shadertoy.com/view/tt2Szh

Results are less good than the other approach, see GPUHydraulicErosionParticle.
*/

/*!
\brief Default constructor. No call to OpenGL because context is not initialized.
*/
GPUHydraulicErosionGrid_Test::GPUHydraulicErosionGrid_Test()
{
    bedrockBuffer = sedDepoBuffer = waterBuffer = sedimentBuffer = 0;
    inDeltaH = outDeltaH = 0;
    tempBedrockBuffer = tempSedDepoBuffer = tempWaterBuffer = tempSedimentBuffer = 0;
    outWaterSpeedBuffer = outFullHeightBuffer = 0;
    outHardnessBuffer = 0;
    simulationShader = 0;
    hardnessFunctionIndex = 0;
    totalBufferSize = dispatchSize = 0;
    outIntDebug = outFloatDebug = 0;
}

/*!
\brief Destructor. Destroys the OpenGl buffers and shader program.
*/
GPUHydraulicErosionGrid_Test::~GPUHydraulicErosionGrid_Test()
{
    glDeleteBuffers(1, &waterBuffer);
    glDeleteBuffers(1, &tempWaterBuffer);
    glDeleteBuffers(1, &bedrockBuffer);
    glDeleteBuffers(1, &inDeltaH);
    glDeleteBuffers(1, &outDeltaH);
    glDeleteBuffers(1, &tempBedrockBuffer);
    glDeleteBuffers(1, &sedimentBuffer);
    glDeleteBuffers(1, &tempSedimentBuffer);
    glDeleteBuffers(1, &outWaterSpeedBuffer);
    glDeleteBuffers(1, &outHardnessBuffer);
    glDeleteBuffers(1, &outFullHeightBuffer);
    glDeleteBuffers(1, &outIntDebug);
    glDeleteBuffers(1, &outFloatDebug);
    glDeleteBuffers(1, &sedDepoBuffer);
    glDeleteBuffers(1, &tempSedDepoBuffer);

    release_program(simulationShader);
}

/*!
\brief Init function, must be called if the underlying heightfield changed.
Initialize a default water and sediment layer to zero and send it to the GPU;

OpenGL buffers are also initialized (just once) here.
*/
void GPUHydraulicErosionGrid_Test::Init(const HeightField& hf)
{
    // Prepare data for first step
    int nx = hf.GetSizeX();
    int ny = hf.GetSizeY();
    totalBufferSize = hf.VertexSize();
    dispatchSize = (max(nx, ny) / 8) + 1;

    std::vector<float> bedrock, water, sediment;
    bedrock.resize(totalBufferSize);
    water.resize(totalBufferSize, 0.0f);
    sediment.resize(totalBufferSize, 0.0f);
    for (int i = 0; i < totalBufferSize; i++)
        bedrock[i] = hf.at(i);

    // Prepare shader & Init buffer - Just done once
    if (simulationShader == 0)
    {
        QString fullPath = System::GetResource("ARCHESLIBDIR", "/LibHeightField/LibHeightField/Shaders/heightfield_hydraulic_test.glsl");

        if (fullPath.isEmpty())
        {
            std::cout << "GPUHydraulicErosion::Init() : variable d'environnement ARCHESLIBDIR non définie" << std::endl;
            std::cin.get();
            exit(-1);
        }
        QByteArray ba = fullPath.toLocal8Bit();
        simulationShader = read_program(ba.data());
    }
    if (bedrockBuffer == 0)
        glGenBuffers(1, &bedrockBuffer);
    if (sedDepoBuffer == 0)
        glGenBuffers(1, &sedDepoBuffer);
    if (waterBuffer == 0)
        glGenBuffers(1, &waterBuffer);
    if (sedimentBuffer == 0)
        glGenBuffers(1, &sedimentBuffer);
    if (tempBedrockBuffer == 0)
        glGenBuffers(1, &tempBedrockBuffer);
    if (tempSedDepoBuffer == 0)
        glGenBuffers(1, &tempSedDepoBuffer);
    if (tempWaterBuffer == 0)
        glGenBuffers(1, &tempWaterBuffer);
    if (tempSedimentBuffer == 0)
        glGenBuffers(1, &tempSedimentBuffer);
    if (inDeltaH == 0)
        glGenBuffers(1, &inDeltaH);
    if (outDeltaH == 0)
        glGenBuffers(1, &outDeltaH);
    if (outWaterSpeedBuffer == 0)
        glGenBuffers(1, &outWaterSpeedBuffer);
    if (hardnessBuffer.buffer == 0) {
        glGenBuffers(1, &hardnessBuffer.buffer);
        //std::cout << "test\n";
    }
    if (outHardnessBuffer == 0)
        glGenBuffers(1, &outHardnessBuffer);
    if (outFullHeightBuffer == 0)
        glGenBuffers(1, &outFullHeightBuffer);
    if (outIntDebug == 0)
        glGenBuffers(1, &outIntDebug);
    if (outFloatDebug == 0)
        glGenBuffers(1, &outFloatDebug);

    // Storage buffer
    glUseProgram(simulationShader);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, bedrockBuffer);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &bedrock.front(), GL_STREAM_READ);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, sedDepoBuffer);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &sediment.front(), GL_STREAM_READ);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, waterBuffer);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &water.front(), GL_STREAM_READ);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, sedimentBuffer);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &sediment.front(), GL_STREAM_READ);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, tempBedrockBuffer);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &bedrock.front(), GL_STREAM_READ);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, tempSedDepoBuffer);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &sediment.front(), GL_STREAM_READ);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, tempWaterBuffer);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &water.front(), GL_STREAM_READ);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, tempSedimentBuffer);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &sediment.front(), GL_STREAM_READ);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, inDeltaH);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &sediment.front(), GL_STREAM_READ);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, outDeltaH);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &sediment.front(), GL_STREAM_READ);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, outWaterSpeedBuffer);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &sediment.front(), GL_STREAM_READ);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, hardnessBuffer.buffer);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &sediment.front(), GL_STREAM_READ);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, outHardnessBuffer);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &sediment.front(), GL_STREAM_READ);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, outFullHeightBuffer);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &sediment.front(), GL_STREAM_READ);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, outIntDebug);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &sediment.front(), GL_STREAM_READ);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, outFloatDebug);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &sediment.front(), GL_STREAM_READ);

    // Uniforms - just once
    Box box = hf.GetBox();
    glUniform1i(glGetUniformLocation(simulationShader, "nx"), nx);
    glUniform1i(glGetUniformLocation(simulationShader, "ny"), ny);
    glUniform3f(glGetUniformLocation(simulationShader, "a"), box[0][0], box[0][1], box[0][2]);
    glUniform3f(glGetUniformLocation(simulationShader, "b"), box[1][0], box[1][1], box[1][2]);

    glUseProgram(0);

    Vector2 cellDiag = hf.CellDiagonal();
    processingShader.Init(nx, ny, cellDiag[0]);
}

/*!
\brief Swap input/output buffer of the lastest simulation step.
*/
void GPUHydraulicErosionGrid_Test::SwapGPUBuffers() const
{
    // Swap in/out buffers
    /*glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, pair ? outBedrockBuffer : bedrockBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, pair ? bedrockBuffer : outBedrockBuffer);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, pair ? outWaterBuffer : waterBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, pair ? waterBuffer : outWaterBuffer);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, pair ? outSedimentBuffer : sedimentBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, pair ? sedimentBuffer : outSedimentBuffer);

    pair = !pair;*/
}

/*!
\brief Simulate n step of hydraulic erosion.

This function doesn't get any data back on the CPU.
\param n number of step
*/
void GPUHydraulicErosionGrid_Test::Step(int n)
{
    
    for (int i = 0; i < n; i++) {

        processingShader.Step(smooth_steps, inDeltaH, outDeltaH);

        glUseProgram(simulationShader);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, tempBedrockBuffer);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, tempWaterBuffer);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, tempSedimentBuffer);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, bedrockBuffer);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, waterBuffer);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, sedimentBuffer);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 6, outWaterSpeedBuffer);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 7, hardnessBuffer.buffer);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 8, outHardnessBuffer);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 9, outFullHeightBuffer);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 10, outIntDebug);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 11, outFloatDebug);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 12, inDeltaH);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 13, outDeltaH);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 14, sedDepoBuffer);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 15, tempSedDepoBuffer);

        glUniform1i(glGetUniformLocation(simulationShader, "useHardness"), hardnessBuffer.used);
        glUniform1i(glGetUniformLocation(simulationShader, "hardnessFunctionIndex"), hardnessFunctionIndex);

        glDispatchCompute(dispatchSize, dispatchSize, 1);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

        // dual buffering
        std::swap(bedrockBuffer,  tempBedrockBuffer);
        std::swap(sedDepoBuffer,  tempSedDepoBuffer);
        std::swap(waterBuffer,    tempWaterBuffer);
        std::swap(sedimentBuffer, tempSedimentBuffer);
    }
    glUseProgram(0);
}

/*!
\brief Set the uniforms for the shader.
*/
void GPUHydraulicErosionGrid_Test::SetUniforms(float erosionSpeed, float depositionSpeed, float rain, float evaporation)
{
    glUseProgram(simulationShader);
    glUniform1f(glGetUniformLocation(simulationShader, "erosionSpeed"), erosionSpeed);
    glUniform1f(glGetUniformLocation(simulationShader, "depositionSpeed"), depositionSpeed);
    glUniform1f(glGetUniformLocation(simulationShader, "evaporation"), evaporation);
    glUniform1f(glGetUniformLocation(simulationShader, "rain"), rain);
    glUseProgram(0);
}

/*!
\brief Set the global erosion speed uniform which affects both erosion and deposition.
\param globalSpeed speed, should be between [0, 1]
*/
void GPUHydraulicErosionGrid_Test::SetGlobalSpeed(float globalSpeed)
{
    glUseProgram(simulationShader);
    glUniform1f(glGetUniformLocation(simulationShader, "globalSpeed"), globalSpeed);
    glUseProgram(0);
}

/*!
\brief Set maxWater, water will be clamped between [0, maxWater].
\param maxWater, should be between [0, +inf]
*/
void GPUHydraulicErosionGrid_Test::SetMaxWater(float maxWater)
{
    glUseProgram(simulationShader);
    //glUniform1f(glGetUniformLocation(simulationShader, "maxWater"), maxWater);
    glUseProgram(0);
}

/*!
\brief Set flowRate, ie fraction of water & sediment that are transported at each iteration.
\param flowRate, should be between [0, 1]
*/
void GPUHydraulicErosionGrid_Test::SetFlowRate(float flowRate)
{
    glUseProgram(simulationShader);
    glUniform1f(glGetUniformLocation(simulationShader, "flowRate"), flowRate);
    glUseProgram(0);
}


/*!
\brief Set flow_p, ie exponant in flow computation, same for water and sediments.
\param flow_p, should be between [0, +inf]
*/
void GPUHydraulicErosionGrid_Test::SetFlowP(float flow_p)
{
    glUseProgram(simulationShader);
    glUniform1f(glGetUniformLocation(simulationShader, "flow_p"), flow_p);
    glUseProgram(0);
}


/*!
\brief Get the latest full elevation data buffer with no transfer back on the CPU.
*/
GLuint GPUHydraulicErosionGrid_Test::GetData() const
{
    return outFullHeightBuffer;
}

/*!
\brief Reads back newest bedrock and sediment elevation from the GPU, into initialized scalarfields.
*/
void GPUHydraulicErosionGrid_Test::GetData(ScalarField2& hf) const
{
    // Temporary float array for OpenGL.
    std::vector<float> data;
    data.resize(hf.VertexSize());

    // Initialize to bedrock elevation
    glGetNamedBufferSubData(bedrockBuffer, 0, sizeof(float) * totalBufferSize, data.data());
    for (int i = 0; i < hf.VertexSize(); i++)
        hf[i] = double(data[i]);

    // Add sediments
    glGetNamedBufferSubData(sedDepoBuffer, 0, sizeof(float) * totalBufferSize, data.data());
    for (int i = 0; i < hf.VertexSize(); i++)
        hf[i] += double(data[i]);
}

/*!
\brief Reads back newest bedrock and sediment elevation from the GPU, into initialized scalarfields.
*/
void GPUHydraulicErosionGrid_Test::GetSediments(ScalarField2& sed) const
{
    // Temporary float array for OpenGL.
    std::vector<float> data;
    data.resize(sed.VertexSize());

    // Get sediments
    glGetNamedBufferSubData(sedDepoBuffer, 0, sizeof(float) * totalBufferSize, data.data());
    for (int i = 0; i < sed.VertexSize(); i++)
        sed[i] = double(data[i]);
}

/*!
\brief Reads back newest bedrock and sediment elevation from the GPU, into initialized scalarfields.
*/
void GPUHydraulicErosionGrid_Test::GetData(ScalarField2& bedrock, ScalarField2& sediments) const
{
    // Temporary float array for OpenGL
    const int size = bedrock.VertexSize();
    std::vector<float> data;
    data.resize(size);

    // Get bedrock
    glGetNamedBufferSubData(bedrockBuffer, 0, sizeof(float) * totalBufferSize, data.data());
    for (int i = 0; i < size; i++)
        bedrock[i] = double(data[i]);

    // Get sediment
    glGetNamedBufferSubData(sedDepoBuffer, 0, sizeof(float) * totalBufferSize, data.data());
    for (int i = 0; i < size; i++)
        sediments[i] = double(data[i]);
}

/*!
\brief Reads back newest bedrock, water and sediment elevation from the GPU, into initialized scalarfields.
*/
void GPUHydraulicErosionGrid_Test::GetData(ScalarField2& bedrock, ScalarField2& sediments, ScalarField2& water) const
{
    // Temporary float array for OpenGL.
    const int size = bedrock.VertexSize();
    std::vector<float> data;
    data.resize(size);

    std::cout << bedrock.VertexSize() << " " << sediments.VertexSize() << " " << water.VertexSize() << std::endl;

    // Get bedrock
    glGetNamedBufferSubData(bedrockBuffer, 0, sizeof(float) * totalBufferSize, data.data());
    for (int i = 0; i < size; i++)
        bedrock[i] = double(data[i]);

    // Get Sediments
    glGetNamedBufferSubData(sedDepoBuffer, 0, sizeof(float) * totalBufferSize, data.data());
    for (int i = 0; i < size; i++)
        sediments[i] = double(data[i]);    

    // Get water
    glGetNamedBufferSubData(waterBuffer, 0, sizeof(float) * totalBufferSize, data.data());
    for (int i = 0; i < size; i++)
        water[i] = double(data[i]);
}

void GPUHydraulicErosionGrid_Test::GetDataDebug(ScalarField2& intMap, ScalarField2& floatMap) const {
    // Temporary arrays for OpenGL.
    const int size = totalBufferSize;
    std::vector<int> int_data;
    int_data.resize(size);
    std::vector<float> float_data;
    float_data.resize(size);

    // Get intDebug
    glGetNamedBufferSubData(outIntDebug, 0, sizeof(float) * totalBufferSize, int_data.data());
    for (int i = 0; i < size; i++)
        intMap[i] = int_data[i];

    // Get floatDebug
    glGetNamedBufferSubData(outFloatDebug, 0, sizeof(float) * totalBufferSize, float_data.data());
    for (int i = 0; i < size; i++)
        floatMap[i] = float_data[i];
}

/*!
\brief Flag for using the hardness alpha map or not.
\param use
*/
void GPUHydraulicErosionGrid_Test::UseHardness(bool use)
{
    hardnessBuffer.used = use;
}

/*!
\brief Set the hardness function to use in the shader.

See function GetHardness in the shader for more details.
\param index
0 = Simple strata
1 = Warped strata
2 = 3D noise
3 = Faulting
The hardness function can be combined with an alpha map
*/
void GPUHydraulicErosionGrid_Test::SetHardnessFunction(int index)
{
    hardnessFunctionIndex = index;
}

/*!
\brief Set the hardness alpha map in a buffer.
\param hardness Hardness alpha.
*/
void GPUHydraulicErosionGrid_Test::SetHardnessAlpha(const ScalarField2& hardness)
{
    glUseProgram(simulationShader);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, hardnessBuffer.buffer);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &hardness.GetAsFloats()[0], GL_STREAM_READ);
    glUseProgram(0);
}

/*!
\brief Get hardness data back from the GPU and store it in a scalar field.
\param hardness the returned hardness data
*/
void GPUHydraulicErosionGrid_Test::GetDataHardness(ScalarField2& hardness) const
{
    //if (hardnessBuffer.used == false) {
    //    //std::cout << "false\n";
    //    return;
    //}
    std::vector<float> data;
    data.resize(totalBufferSize);
    glGetNamedBufferSubData(outHardnessBuffer, 0, sizeof(float) * totalBufferSize, data.data());
    for (int i = 0; i < totalBufferSize; i++) {
        hardness[i] = data[i];
        //if (i % 1000 == 0) std::cout << data[i] << std::endl;
    }
}


// GPUMapProcessing ----------------------------------------------------------
GPUMapProcessing::GPUMapProcessing(const QString& _shader_program_name) {
    shader_program_name = _shader_program_name;
}

GPUMapProcessing::~GPUMapProcessing() {

    /*glDeleteBuffers(1, &outBuffer);
    glDeleteBuffers(1, &tempBuffer);*/

    release_program(shader_program_id);
}

void GPUMapProcessing::Init(int _buffer_size_x, int _buffer_size_y, float _cellSize) {
    buffer_size_x = _buffer_size_x;
    buffer_size_y = _buffer_size_y;
    cellSize = _cellSize;
    total_buffer_size = buffer_size_x * buffer_size_y;
    dispatch_size = (max(buffer_size_x, buffer_size_x) / 8) + 1;

    //std::vector<float> tmpZeroes(total_buffer_size, 0);

    // Prepare shader & Init buffer - Just done once
    if (shader_program_id == 0) {
        QString fullPath = shader_program_name;
        QByteArray ba = fullPath.toLocal8Bit();
        shader_program_id = read_program(ba.data());
    }

    /*glGenBuffers(1, &outBuffer);
    glGenBuffers(1, &tempBuffer);*/

    // Storage buffer
    glUseProgram(shader_program_id);

    /*glBindBuffer(GL_SHADER_STORAGE_BUFFER, outBuffer);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * total_buffer_size, &tmpZeroes.front(), GL_STREAM_READ);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, tempBuffer);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * total_buffer_size, &tmpZeroes.front(), GL_STREAM_READ);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);*/

    // Uniforms - just once
    glUniform1i(glGetUniformLocation(shader_program_id, "buffer_size_x"), buffer_size_x);
    glUniform1i(glGetUniformLocation(shader_program_id, "buffer_size_y"), buffer_size_y);

    glUseProgram(0);
}

void GPUMapProcessing::Step(int smooth_steps, GLuint& in, GLuint& out) {
    glUseProgram(shader_program_id);

    for (int i = 0; i < smooth_steps; i++) {

        //glUniform1f(glGetUniformLocation(shader_program_id, "threshold"), threshold);

        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, out);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, in);

        glDispatchCompute(dispatch_size, dispatch_size, 1);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

        std::swap(in, out);
    }

    std::swap(in, out);

    glUseProgram(0);
}