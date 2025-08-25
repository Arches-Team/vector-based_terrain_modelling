#include "libs/gputerrainsimu.h"
#include "libs/cpu.h"

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
GPUHydraulicErosionGrid::GPUHydraulicErosionGrid()
{
	hardnessFunctionIndex = 0;
	totalBufferSize = dispatchSize = 0;

}

/*!
\brief Destructor. Destroys the OpenGl buffers and shader program.
*/
GPUHydraulicErosionGrid::~GPUHydraulicErosionGrid()
{
	bedrockBuffer.Destroy();
	outBedrockBuffer.Destroy();
	waterBuffer.Destroy();
	outWaterBuffer.Destroy();
	sedimentBuffer.Destroy();
	outSedimentBuffer.Destroy();
	outWaterSpeedBuffer.Destroy();
	outHardnessBuffer.Destroy();
	outFullHeightBuffer.Destroy();
	hardnessBuffer.Destroy();

	simulationShader.Destroy();
}

/*!
\brief Init function, must be called if the underlying heightfield changed.
Initialize a default water and sediment layer to zero and send it to the GPU;

OpenGL buffers are also initialized (just once) here.
*/
void GPUHydraulicErosionGrid::Init(const HeightField& hf)
{
	// Prepare data for first step
	const int nx = hf.GetSizeX();
	const int ny = hf.GetSizeY();
	totalBufferSize = hf.VertexSize();
	dispatchSize = (max(nx, ny) / 8) + 1;

	std::vector<float> bedrock;
	std::vector<float> zeros(totalBufferSize, 0.);
	double h_min, h_max;
	hf.GetRange(h_min, h_max);
	bedrock.resize(totalBufferSize);
	for (int i = 0; i < totalBufferSize; i++)
		bedrock[i] = (hf.at(i) - h_min) / (h_max - h_min) * scale ;

	// Prepare shader
	QString fullPath = System::GetResource("ARCHESLIBDIR", "/LibHeightField/LibHeightField/Shaders/heightfield_hydraulic.glsl");
	if (fullPath.isEmpty())
	{
		std::cout << "GPUHydraulicErosion::Init(): variable d'environnement ARCHESLIBDIR non définie" << std::endl;
		std::cin.get();
		exit(-1);
	}
	QByteArray ba = fullPath.toLocal8Bit();
	simulationShader.Initialize(ba.data());

	// Storage buffers
	bedrockBuffer.Generate();
	bedrockBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &bedrock.front(), GL_STREAM_READ);

	outBedrockBuffer.Generate();
	outBedrockBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &bedrock.front(), GL_STREAM_READ);

	waterBuffer.Generate();
	waterBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &zeros.front(), GL_STREAM_READ);

	outWaterBuffer.Generate();
	outWaterBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &zeros.front(), GL_STREAM_READ);

	sedimentBuffer.Generate();
	sedimentBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &zeros.front(), GL_STREAM_READ);

	outSedimentBuffer.Generate();
	outSedimentBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &zeros.front(), GL_STREAM_READ);

	outWaterSpeedBuffer.Generate();
	outWaterSpeedBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &zeros.front(), GL_STREAM_READ);

	outFullHeightBuffer.Generate();
	outFullHeightBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &zeros.front(), GL_STREAM_READ);

	outHardnessBuffer.Generate();
	outHardnessBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &zeros.front(), GL_STREAM_READ);

	hardnessBuffer.Generate();
	hardnessBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &zeros.front(), GL_STREAM_READ);

	// Uniforms - just once
	simulationShader.Bind();
	Box2 box = hf.Array2::GetBox();
	Vector2 cellDiag = hf.CellDiagonal();
	simulationShader.SetUniform("nx", nx);
	simulationShader.SetUniform("ny", ny);
	simulationShader.SetUniform("a", float(box[0][0]), float(box[0][1]));
	simulationShader.SetUniform("b", float(box[1][0]), float(box[1][1]));
	simulationShader.SetUniform("cellDiag", float(cellDiag[0]), float(cellDiag[1]));
	simulationShader.Unbind();
}

/*!
\brief Simulate n step of hydraulic erosion.

This function doesn't get any data back on the CPU.
\param n number of step
*/
void GPUHydraulicErosionGrid::Step(int n)
{
	simulationShader.Bind();

	outWaterSpeedBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 6);
	hardnessBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 7);
	outHardnessBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 8);
	outFullHeightBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 9);
	for (int i = 0; i < n; i++)
	{
		bedrockBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 0);
		waterBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 1);
		sedimentBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 2);
		outBedrockBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 3);
		outWaterBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 4);
		outSedimentBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 5);

		simulationShader.SetUniform("local", int(false));
		simulationShader.SetUniform("useHardness", int(hardnessBuffer.used));
		simulationShader.SetUniform("hardnessFunctionIndex", hardnessFunctionIndex);

		glDispatchCompute(dispatchSize, dispatchSize, 1);
		glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

		std::swap(outBedrockBuffer, bedrockBuffer);
		std::swap(outWaterBuffer, waterBuffer);
		std::swap(outSedimentBuffer, sedimentBuffer);
	}
	simulationShader.Unbind();
}

/*!
\brief Simulate n step of hydraulic erosion in a circular domain.

This function doesn't get any data back on the CPU.
\param n number of step
\param c center of the domain
\param r radius of the domain
*/
void GPUHydraulicErosionGrid::Step(int n, const Vector2& c, double r)
{
	simulationShader.Bind();

	outWaterSpeedBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 6);
	hardnessBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 7);
	outHardnessBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 8);
	outFullHeightBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 9);
	for (int i = 0; i < n; i++)
	{
		bedrockBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 0);
		waterBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 1);
		sedimentBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 2);
		outBedrockBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 3);
		outWaterBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 4);
		outSedimentBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 5);

		simulationShader.SetUniform("local", int(true));
		simulationShader.SetUniform("center", float(c[0]), float(c[1]));
		simulationShader.SetUniform("radius", float(r));

		simulationShader.SetUniform("useHardness", int(hardnessBuffer.used));
		simulationShader.SetUniform("hardnessFunctionIndex", hardnessFunctionIndex);

		glDispatchCompute(dispatchSize, dispatchSize, 1);
		glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

		std::swap(outBedrockBuffer, bedrockBuffer);
		std::swap(outWaterBuffer, waterBuffer);
		std::swap(outSedimentBuffer, sedimentBuffer);
	}
	simulationShader.Unbind();
}

/*!
\brief Set the uniforms for the shader.
*/
void GPUHydraulicErosionGrid::SetUniforms(float erosionSpeed, float depositionSpeed, float rain, float evaporation)
{
	simulationShader.Bind();
	simulationShader.SetUniform("erosionSpeed", erosionSpeed);
	simulationShader.SetUniform("depositionSpeed", depositionSpeed);
	simulationShader.SetUniform("rain", rain);
	simulationShader.SetUniform("evaporation", evaporation);
	simulationShader.Unbind();
}

/*!
\brief Set the global erosion speed uniform which affects both erosion and deposition.
\param globalSpeed speed, should be between [0, 1]
*/
void GPUHydraulicErosionGrid::SetGlobalSpeed(float globalSpeed)
{
	simulationShader.Bind();
	simulationShader.SetUniform("globalSpeed", globalSpeed);
	simulationShader.Unbind();
}

/*!
\brief Get the latest full elevation data buffer with no transfer back on the CPU.
*/
GLuint GPUHydraulicErosionGrid::GetData() const
{
	return outFullHeightBuffer.GetBuffer();
}

/*!
\brief Get the latest water speed data buffer with no transfer back on the CPU.
*/
void GPUHydraulicErosionGrid::GetSpeeds(ScalarField2& speeds) const {
	std::vector<float> data;
	data.resize(speeds.VertexSize());

	outWaterSpeedBuffer.GetData(data);
	for (int i = 0; i < speeds.VertexSize(); i++)
		speeds[i] = double(data[i]);
}


/*!
\brief Reads back newest bedrock and sediment elevation from the GPU, into initialized scalarfields.
*/
void GPUHydraulicErosionGrid::GetData(ScalarField2& hf) const
{
	// Temporary float array for OpenGL.
	std::vector<float> data;
	data.resize(hf.VertexSize());

	outBedrockBuffer.GetData(data);
	for (int i = 0; i < hf.VertexSize(); i++)
		hf[i] = 3000. / scale * double(data[i]);
	/*outSedimentBuffer.GetData(data);
	for (int i = 0; i < hf.VertexSize(); i++)
		hf[i] += double(data[i]);*/
}

/*!
\brief Reads back newest bedrock and sediment elevation from the GPU, into initialized scalarfields.
*/
void GPUHydraulicErosionGrid::GetData(ScalarField2& bedrock, ScalarField2& sediments) const
{
	// Temporary float array for OpenGL
	const int size = bedrock.VertexSize();
	std::vector<float> data;
	data.resize(size);

	outBedrockBuffer.GetData(data);
	for (int i = 0; i < size; i++)
		bedrock[i] = double(data[i]);
	outSedimentBuffer.GetData(data);
	for (int i = 0; i < size; i++)
		sediments[i] = double(data[i]);
}

/*!
\brief Reads back newest bedrock, water and sediment elevation from the GPU, into initialized scalarfields.
*/
void GPUHydraulicErosionGrid::GetData(ScalarField2& bedrock, ScalarField2& sediments, ScalarField2& water) const
{
	// Temporary float array for OpenGL.
	const int size = bedrock.VertexSize();
	std::vector<float> data;
	data.resize(size);

	outBedrockBuffer.GetData(data);
	for (int i = 0; i < size; i++)
		bedrock[i] = double(data[i]);
	outSedimentBuffer.GetData(data);
	for (int i = 0; i < size; i++)
		sediments[i] = double(data[i]);
	outWaterBuffer.GetData(data);
	for (int i = 0; i < size; i++)
		water[i] = double(data[i]);
}

/*!
\brief Flag for using the hardness alpha map or not.
\param use
*/
void GPUHydraulicErosionGrid::UseHardness(bool use)
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
void GPUHydraulicErosionGrid::SetHardnessFunction(int index)
{
	hardnessFunctionIndex = index;
}

/*!
\brief Set the hardness alpha map in a buffer.
\param hardness Hardness alpha.
*/
void GPUHydraulicErosionGrid::SetHardnessAlpha(const ScalarField2& hardness)
{
	hardnessBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &hardness.GetAsFloats()[0], GL_STREAM_READ);
}

/*!
\brief Get hardness data back from the GPU and store it in a scalar field.
\param hardness the returned hardness data
*/
void GPUHydraulicErosionGrid::GetDataHardness(ScalarField2& hardness) const
{
	if (hardnessBuffer.used == false)
		return;
	std::vector<float> data;
	data.resize(size_t(hardness.VertexSize()));
	outHardnessBuffer.GetData(data);
	for (int i = 0; i < hardness.VertexSize(); i++)
		hardness[i] = double(data[i]);
}
