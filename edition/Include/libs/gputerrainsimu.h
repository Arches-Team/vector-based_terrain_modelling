#pragma once

#include "libs/gpu-wrapper.h"
#include "libs/gpuheightfield.h"
#include "libs/heightfield.h"
#include "libs/vectorfield.h"
#include "libs/cpu.h"

class GLOptionalStorageBuffer : public GLBuffer
{
public:
  bool used = false;
public:
  inline GLOptionalStorageBuffer() { }
  inline ~GLOptionalStorageBuffer() { }
};

class GLOptionalBuffer
{
public:
  GLuint buffer;
  bool used = false;
};

class GPUHydraulicErosionGrid
{
protected:
  GLShader simulationShader;	//!< Compute shader.

  GLBuffer bedrockBuffer;		//!< Bedrock elevation buffer
  GLBuffer outBedrockBuffer;	//!< Out bedrock elevation buffer

  GLBuffer waterBuffer;			//!< Water elevation buffer
  GLBuffer outWaterBuffer;		//!< Out water elevation buffer

  GLBuffer sedimentBuffer;		//!< Sediment elevation buffer
  GLBuffer outSedimentBuffer;	//!< Out Sediment elevation buffer

  GLOptionalStorageBuffer hardnessBuffer;	//!< Optional Hardness buffer provided by user
  int hardnessFunctionIndex;
  GLBuffer outHardnessBuffer;		//!< Output hardness buffer (computed on GPU)

  GLBuffer outWaterSpeedBuffer;		//!< Output water speed buffer
  GLBuffer outFullHeightBuffer;		//!< Output latest full height (bedrock + sediment) buffer.

  int totalBufferSize;				//!< Total buffer size defined as nx * ny
  int dispatchSize;					//!< Single dispatch size

  double scale = 0.2;

public:
  GPUHydraulicErosionGrid();
  ~GPUHydraulicErosionGrid();

  void Init(const HeightField&);
  void Step(int n);
  void Step(int n, const Vector2& p, double r);
  void SetUniforms(float erosionSpeed, float depositionSpeed, float rain, float evaporation);
  void SetGlobalSpeed(float globalSpeed);
  void UseHardness(bool use = true);
  void SetHardnessFunction(int index);
  void SetHardnessAlpha(const ScalarField2& hardness);

  GLuint GetData() const;
  void GetSpeeds(ScalarField2& speeds) const;
  void GetData(ScalarField2& hf) const;
  void GetData(ScalarField2& bedrock, ScalarField2& sediments) const;
  void GetData(ScalarField2& bedrock, ScalarField2& sediments, ScalarField2& water) const;
  void GetDataHardness(ScalarField2& hardness) const;
};

class GPUMapProcessing {
private:

  QString shader_program_name;
  GLuint shader_program_id = 0;

  int buffer_size_x = 0;
  int buffer_size_y = 0;
  float cellSize;
  float kw;
  int total_buffer_size = 0;
  int dispatch_size = 0;

public:
  GPUMapProcessing(const QString& _shader_program_name);
  ~GPUMapProcessing();

  void Init(int _buffer_size_x, int _buffer_size_y, float _cellSize);
  void Step(int smooth_steps, GLuint& in, GLuint& out);

  //GLuint GetData() const { return outBuffer; };

  void SetUniforms(float _kw) { kw = _kw; };
};

class GPUHydraulicErosionGrid_Test
{
protected:
  GPUMapProcessing processingShader = GPUMapProcessing(System::GetResource(QString::fromStdString(std::string(SOLUTION_DIR) + "/shaders/libs/map_processing.glsl")));

  GLuint simulationShader;			//!< Compute shader.

  GLuint bedrockBuffer;				//!< Bedrock elevation buffer
  GLuint sedDepoBuffer;             //!< Deposited sediment elevation buffer
  GLuint waterBuffer;				//!< Water elevation buffer
  GLuint sedimentBuffer;			//!< Sediment elevation buffer

  GLOptionalBuffer hardnessBuffer;	//!< Optional Hardness buffer provided by user
  int hardnessFunctionIndex;
  GLuint outHardnessBuffer;			//!< Output hardness buffer (computed on GPU or not)

  GLuint tempBedrockBuffer;			//!< Temporary bedrock elevation buffer
  GLuint tempSedDepoBuffer;         //!< Temporarty deposited sediment elevation buffer
  GLuint tempWaterBuffer;			//!< Temporary water elevation buffer
  GLuint tempSedimentBuffer;		//!< Temporary sediment elevation buffer

  GLuint outWaterSpeedBuffer;		//!< Output water speed buffer
  GLuint outFullHeightBuffer;		//!< Output latest full height (bedrock + sediment) buffer.

  GLuint inDeltaH;
  GLuint outDeltaH;
  int smooth_steps = 0;

  GLuint outIntDebug;
  GLuint outFloatDebug;

  int totalBufferSize;				//!< Total buffer size defined as nx * ny
  int dispatchSize;					//!< Single dispatch size

  void SwapGPUBuffers() const;

public:
  GPUHydraulicErosionGrid_Test();
  ~GPUHydraulicErosionGrid_Test();

  void Init(const HeightField&);
  void Step(int n);
  void SetUniforms(float erosionSpeed, float depositionSpeed, float rain, float evaporation);
  void SetGlobalSpeed(float globalSpeed);
  void SetMaxWater(float maxWater);
  void SetFlowRate(float flowRate);
  void SetFlowP(float flow_p);
  void SetSmoothSteps(int _smooth_steps) { smooth_steps = _smooth_steps; };

  void UseHardness(bool use = true);
  void SetHardnessFunction(int index);
  void SetHardnessAlpha(const ScalarField2& hardness);

  GLuint GetData() const;
  void GetData(ScalarField2& hf) const;
  void GetSediments(ScalarField2& sed) const;
  void GetData(ScalarField2& bedrock, ScalarField2& sediments) const;
  void GetData(ScalarField2& bedrock, ScalarField2& sediments, ScalarField2& water) const;
  void GetDataDebug(ScalarField2& intMap, ScalarField2& floatMap) const;
  void GetDataHardness(ScalarField2& hardness) const;
};

class GPUHydraulicErosionParticle
{
protected:
  GLShader simulationShader;		//!< Compute shader

  GLBuffer bedrockBuffer;			//!< Bedrock buffer
  GLBuffer sedimentBuffer;			//!< Sediment buffer
  GLBuffer randomIndicesBuffer;		//!< Random number buffer. Easier than generating random number on the GPU
  GLBuffer brushIndicesBuffer;		//!< Compute buffer
  GLBuffer brushWeightsBuffer;		//!< Compute buffer

  GLOptionalStorageBuffer hardnessBuffer;	//!< Optional Hardness buffer
  int hardnessFunctionIndex;
  GLBuffer outHardnessBuffer;		//!< Output hardness buffer (useful if we use the GPU function to generate hardness on the fly)
  GLOptionalStorageBuffer precipitationBuffer;	//!< Precipitation buffer

  int totalBufferSize;						//!< Total buffer size defined as nx * ny
  int dispatchSize;							//!< Single dispatch size
  int nx;									//!< Heightfield dimension
  int ny;									//!< Heightfield dimension

  double scaleZa;							//!< Heightfield min elevation
  double scaleZb;							//!< Heightfield max elevation
  std::vector<int> randomIndices;			//!< Random indices for the algorithm
  std::vector<float> tmpData;			//!< Temporary array for retreiving GPU data.

public:
  GPUHydraulicErosionParticle();
  ~GPUHydraulicErosionParticle();

  void Init(const HeightField& hf);
  void Step(int n);

  void SetUniforms(float inertia, float sedimentCapacity, float depositSpeed, float evaporateSpeed, float erodeSpeed, int maxLifetime, int erosionBrushRadius);
  void SetUniforms(float depositSpeed, float erosionSpeed, float evaporationSpeed, float randomSedimentation, float randomSedFreq, float randomSedAmp);
  void SetUniforms(float depositSpeed, float erosionSpeed, float evaporationSpeed, int erosionBrushRadius, int maxLifeTime);
  void SetGlobalSpeed(float globalSpeed);

  void SetHardnessAlpha(const ScalarField2& hardness);
  void UseHardness(bool use = true);
  void SetHardnessFunction(int index);
  void UsePrecipitation(bool use = true);
  void SetPrecipitation(const ScalarField2& precipitation);

  void GetData(ScalarField2& hf);
  void GetData(ScalarField2& bedrock, ScalarField2& sediments);
  void GetDataHardness(ScalarField2& hardness);
};

class GPUStreamPowerOptimized
{
protected:
  GLShader simulationShader;			//!< Compute shader

  GLBuffer bedrockBuffer;				//!< Bedrock elevation buffer
  GLBuffer outBedrockBuffer;			//!< Output bedrock elevation buffer

  GLBuffer streamBuffer;				//!< Water elevation buffer
  GLBuffer outStreamBuffer;				//!< Output water elevation buffer

  GLBuffer outHardnessBuffer;			//!< Output hardness buffer (to read back from the GPU).

  GLOptionalStorageBuffer upliftBuffer;		//!< Optional uplift buffer
  GLOptionalStorageBuffer boundaryBuffer;	//!< Optional boundary buffer
  GLOptionalStorageBuffer mParamBuffer;		//!< Optional exponent map buffer
  GLOptionalStorageBuffer nParamBuffer;		//!< Optional exponent map buffer
  GLOptionalStorageBuffer kParamBuffer;		//!< Optional erosion factor map buffer
  GLOptionalStorageBuffer retargetBuffer;	//!< Optional retargeting alpha buffer
  GLOptionalStorageBuffer laplacianBuffer;	//!< Optional laplacian factor buffer
  bool dynamicHardness;						//!< Flag for dynamic hardness		

  int nx, ny;
  int totalBufferSize;					//!< Total buffer size defined as nx * ny
  int dispatchSize;						//!< Single dispatch size
  std::vector<float> tmpData;			//!< Temporary array for retreiving GPU data

public:
  GPUStreamPowerOptimized();
  ~GPUStreamPowerOptimized();

  void Init(const HeightField&);
  void Step(int n);
  void Step(int, const Vector2&, double);
  void SetUniforms(int erosionMode, float uplift, float k, float p_sa, float p_sl, float k_h, float k_d);
  void SetDt(float dt) const;

  void SetKParameter(ScalarField2 kParam);
  void UseKParameter(bool use = true);
  void SetMParameter(ScalarField2 mParam);
  void UseMParameter(bool use = true);
  void SetNParameter(ScalarField2 nParam);
  void UseNParameter(bool use = true);
  void SetLaplacianParameter(ScalarField2 laplacianParam);
  void UseLaplacianParameter(bool use = true);
  void SetFlowMode(int flowMode);
  void SetFlowExponent(float flow_p);
  void SetUplift(const ScalarField2& uplift);
  void UseUplift(bool use = true);
  void SetBoundary(Array2I boundary);
  void SetRetargeting(ScalarField2 retargetAlpha);
  void UseRetargeting(bool use = true);
  void UseHardness(bool use = true);
  void SetHardnessFunc(int hardnessFunc);
  void UseDynamicHardness(bool use = true);

  GLuint GetData() const;
  void GetData(ScalarField2&);
  void GetData(ScalarField2&, ScalarField2&);
  void GetData(ScalarField2&, ScalarField2&, ScalarField2&);
  void GetHardness(ScalarField2&);
  double GetLastStepElevationDifference();
};

class GPUThermalErosion
{
protected:
  GLShader simulationShader;	//!< Compute shader.

  GLBuffer bedrockBuffer;		//!< Bedrock elevation buffer
  GLBuffer sedimentBuffer;		//!< Sediment elevation buffer
  GLBuffer outSedimentBuffer;	//!< Output sediment elevation buffer

  int nx, ny;
  int totalBufferSize;			//!< Total buffer size defined as nx * ny
  int dispatchSize;				//!< Single dispatch size
  std::vector<float> tmpData;	//!< Temporary array for retreiving GPU data

public:
  GPUThermalErosion();
  ~GPUThermalErosion();

  void Init(const HeightField&);
  void Init(const ScalarField2&, const ScalarField2&);
  void Step(int);
  void Step(int, const Vector2&, double);
  void SetUniforms(float eps, float tanThresholdAngle, bool noisifiedAngle, float noiseWavelength);
  void GetData(ScalarField2& hf);
  void GetData(ScalarField2& bedrock, ScalarField2& sediments);
  GLuint GetData() const;
  double GetLastStepElevationDifference();
};

class GPUThermalErosionSimple
{
protected:
  GLuint simulationShader;	//!< Compute shader.

  GLuint terrain_buffer;
  GLuint temp_terrain_buffer;

  GLuint in_threshold_buffer;
  GLuint out_threshold_buffer;

  int nx, ny;
  int totalBufferSize;			//!< Total buffer size defined as nx * ny
  int dispatchSize;				//!< Single dispatch size

public:
  GPUThermalErosionSimple();
  ~GPUThermalErosionSimple();

  void Init(const HeightField&);
  void Step(int);
  void SetUniforms(float eps, float tanThresholdAngle, bool noisifiedAngle, float noiseWavelength);
  void SetNoise(float noise_min, float noise_max);
  void UseThreshold(bool use);
  void SetThreshold(const ScalarField2& threshold);
  void GetHeight(ScalarField2& hf);
  void GetThreshold(ScalarField2& sf);
};

class GPUThermalErosionDouble
{
protected:
  GLuint simulationShader;	//!< Compute shader.

  GLuint terrain_buffer;
  GLuint sediment_buffer;
  GLuint temp_sediment_buffer;

  int nx, ny;
  int totalBufferSize;			//!< Total buffer size defined as nx * ny
  int dispatchSize;				//!< Single dispatch size

public:
  GPUThermalErosionDouble();
  ~GPUThermalErosionDouble();

  void Init(const HeightField& terrain, const HeightField& sediment);
  void Step(int);
  void SetUniforms(float eps, float tanThresholdAngle, bool noisifiedAngle, float noiseWavelength);
  void GetData(ScalarField2& hf) const;
  void GetSediment(ScalarField2& sed) const;
};

class GPUHeightFieldAnalysis
{
protected:
  GLShader shader;				//!< Glsl program.
  GLuint hfTextureBuffer;		//!< %Heightfield texture.
  GLBuffer inHfBuffer;			//!< %Heightfield buffer, used by CountPits(), with single precision.
  GLBuffer inHfBufferDouble;	//!< %Heightfield buffer, with double precision.
  GLBuffer outBuffer;			//!< Output scalar buffer.
  GLBuffer outVecBuffer;		//!< Output 2D vector buffer.
  GLBuffer outIntBuffer;		//!< Output integer buffer.
  GLBuffer outDoubleBuffer;		//!< Output double buffer.

  Box box;
  int nx, ny;
  int totalBufferSize;
  int dispatchSize;

public:
  GPUHeightFieldAnalysis();
  ~GPUHeightFieldAnalysis();

  ScalarField2 FractionalLaplacian(const HeightField&, int, double, int);
  ScalarField2 Laplacian(const HeightField&);

  void FractionalGradient(const HeightField&, double, int, VectorField2&, ScalarField2&);
  void Gradient(const HeightField&, VectorField2&, ScalarField2&);

  ScalarField2 Accessibility(const HeightField&, const double&, int);
  ScalarField2 ClearSky(const HeightField&, double, int);
  ScalarField2 Shadow(const HeightField&, const Vector&);
  ScalarField2 SoftShadows(const HeightField&, const Vector&, double);
  ScalarField2 Diffuse(const HeightField&, const Vector&, double);

  int CountPits(const HeightField&);
  int CountPits(GLuint hfBuffer, int nx, int ny);
  Vector2 GetRange(const HeightField&);
  Vector2 GetRange(GLuint hfBuffer, int nx, int ny);

  ScalarField2 Viewshed(const HeightField&, const QPoint&);
  ScalarField2 VisibilityIndex(const HeightField&);

protected:
  void Init(int nx, int ny);
  void InitForSingleBuffer(const HeightField&);
  void InitForDoubleBuffer(const HeightField&);
  void InitForTexture(const HeightField&);

protected:
  double C(int, const double&) const;
  ScalarField2 GetDataToScalarField() const;
};
