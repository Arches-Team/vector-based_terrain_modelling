#version 450 core

#ifdef COMPUTE_SHADER

// In
layout(binding = 0, std430)  readonly buffer InElevation       { float hf[]; };
layout(binding = 14, std430) readonly buffer InSedElevation    { float sed_hf[]; };
layout(binding = 1, std430)  readonly buffer InWater           { float water[]; };
layout(binding = 2, std430)  readonly buffer InSediment        { float sediments[]; };
layout(binding = 12, std430) readonly buffer InDeltaH          { float inDeltaH[]; };

// Out
layout(binding = 3, std430)  writeonly buffer OutElevation     { float out_hf[]; };
layout(binding = 15, std430) writeonly buffer OutSedElevation  { float out_sed_hf[]; };
layout(binding = 4, std430)  writeonly buffer OutWater         { float out_water[]; };
layout(binding = 5, std430)  writeonly buffer OutSediment      { float out_sediments[]; };
layout(binding = 6, std430)  writeonly buffer OutWaterSpeed    { float out_speeds[]; };
layout(binding = 13, std430) writeonly buffer OutDeltaH        { float outDeltaH[]; };


// Options
layout(binding = 7, std430) readonly   buffer HardnessBuffer     { float in_hardness[]; };
layout(binding = 8, std430) writeonly  buffer OutHardnessBuffer  { float out_hardness[]; };
layout(binding = 9, std430) writeonly  buffer OutFullHeight      { float out_fullHeight[]; };
layout(binding = 10, std430)           buffer IntDebug           { int out_intDebug[]; };
layout(binding = 11, std430) writeonly buffer FloatDebug         { float out_floatDebug[]; };

uniform int nx;
uniform int ny;
uniform vec3 a;
uniform vec3 b;

// Default values for uniforms
uniform float globalSpeed = 0.1f;
uniform float erosionSpeed = 1.f;
uniform float depositionSpeed = 1.f;
uniform float evaporation = .81f;
uniform float rain = 2.6f;
uniform float flowRate = .32f;
uniform float flow_p = 1.3;

uniform int useHardness = 0;
uniform int hardnessFunctionIndex = 0;

const float cellArea = (b.x - a.x) * (b.y - a.y) / float((nx - 1) * (ny - 1)) * 0.00001;
const vec2  cellSize = (b.xy - a.xy) / vec2(nx - 1, ny - 1);
const float cellDiag = length(cellSize);
const ivec2 next8[8] = ivec2[8]( ivec2(0,  1), ivec2( 1,  1), ivec2( 1, 0), ivec2( 1, -1),
								 ivec2(0, -1), ivec2(-1, -1), ivec2(-1, 0), ivec2(-1,  1) );




float GetHardness(in float z, int index) {
	//float h = 1.f;
	//return h;

	// Strata function
	/*if (hardnessFunctionIndex == 0) {
		z *= .01f;
		h = sin(z) * .5f + .5f;
	}*/

	if (useHardness != 0)
		return in_hardness[index];
	else
		return 1.f;
}

int ToIndex1D(int i, int j) { return i + nx * j; }

float Slope(ivec2 p, ivec2 q) {
    if (p.x < 0 || p.x >= nx || p.y < 0 || p.y >= ny) return 0.0f;
    if (q.x < 0 || q.x >= nx || q.y < 0 || q.y >= ny) return 0.0f;
    if (p == q) return 0.0f;
    int index_p = ToIndex1D(p.x, p.y);
    int index_q = ToIndex1D(q.x, q.y);
    float d = length((q - p) * cellSize);
    return ((hf[index_q] + sed_hf[index_q]) - (hf[index_p] + sed_hf[index_p])) / d;
}

ivec2 GetFlowSteepest(ivec2 p) {
    ivec2 d = ivec2(0, 0);
    float maxSlope = 0.0f;
    for (int i = 0; i < 8; i++) {
        float ss = Slope(p + next8[i], p);
        if (ss > maxSlope) {
            maxSlope = ss;
            d = next8[i];
        }
    }
    return d;
}

vec2 ComputeIncomingFlowSteepest(ivec2 p) {
    vec2 inFlow = vec2(0., 0.);
    for (int i = 0; i < 8; i++) {
        ivec2 q = p + next8[i];
        ivec2 fd = GetFlowSteepest(q);
        if (q + fd == p) {
			int qid = ToIndex1D(q.x, q.y);
            inFlow += vec2(water[qid], sediments[qid]);
		}
    }
    return inFlow;
}

vec2 GetFlowWeighted(ivec2 p, ivec2 q) { // return what flows from q to p
	vec2 flow = vec2(0., 0.);

	int index_p = ToIndex1D(p.x, p.y);
	int index_q = ToIndex1D(q.x, q.y);
	if (hf[index_q] + sed_hf[index_q] <= hf[index_p] + sed_hf[index_p]) return flow;
	
	float slope_sum = 0.;
	float p_pow_slope = 1.;
	for (int i = 0; i < 8; i++) {
		ivec2 nei = q + next8[i];
		float pow_slope = pow(Slope(nei, q), flow_p);

		if (pow_slope > 0.)
			slope_sum += pow_slope;
		
		if (nei == p)
			p_pow_slope = pow_slope;
	}

	slope_sum = (slope_sum == 0.0f) ? 1.0f : slope_sum;
	
	return vec2(water[index_q], sediments[index_q]) * p_pow_slope / slope_sum;
}

vec2 ComputeIncomingFlowWeighted(ivec2 p) {

	vec2 flow = vec2(0., 0.);

	for (int i = 0; i < 8; i++) {
		ivec2 q = p + next8[i];

		if (p.x < 0 || p.x >= nx || p.y < 0 || p.y >= ny) continue;
		if (q.x < 0 || q.x >= nx || q.y < 0 || q.y >= ny) continue;

		flow += GetFlowWeighted(p, q);
	}

	return flow;
}


vec4 GetAllData(int id) {
	vec4 ret;
	ret.w = hf[id];// Bedrock elevation
	ret.x = water[id];// Amount of water
	ret.y = sediments[id];// Amount of sediment
	ret.z = 0.f;// Flow speed - which is just write-only
	return ret;
}

void WriteAllData(vec4 frag, float hardness, int id) {
	out_hf[id] = frag.w;
	out_water[id] = frag.x;
	out_sediments[id] = frag.y;
	out_speeds[id] = frag.z;
	//out_hardness[id] = hardness;
	out_fullHeight[id] = frag.w + frag.y;
}


layout(local_size_x = 8, local_size_y = 8, local_size_z = 1)in;
void main() {

	int x = int(gl_GlobalInvocationID.x);
	int y = int(gl_GlobalInvocationID.y);
	if (x < 0)return;
	if (y < 0)return;
	if (x >= nx)return;
	if (y >= ny)return;

	ivec2 p = ivec2(x, y);
	int id = ToIndex1D(x, y);
	vec4 data = GetAllData(id);

	//data.x += rain * cellArea;

	float outSlope = Slope(p + GetFlowSteepest(p), p);
	//vec2 inFlow = ComputeIncomingFlowSteepest(p);
	vec2 inFlow = ComputeIncomingFlowWeighted(p);

	// Modify water & sediment values
	if (outSlope > 0. || true) {
		//data.x *= 1. - flowRate;
		data.y *= 1. - flowRate;
	}

	// Add sediment and water
	//data.x += inFlow * flowRate;
	data.x = rain * cellArea + inFlow.x;
	data.y += inFlow.y * flowRate;

	// Prevent waterfalls causing infinitely powerful erosion
	// This gives really nice control over how slopey (big values) or cliffy (small values) the terrain is.
	const float maxWaterSpeed = 1.f;

	// Compute speed // was 1.0/0.01
	const float slopeToSpeed = 1.f / 1.f;// denominator is height difference between input and output points at which we get a waterfall
	float speed = clamp(outSlope * outSlope * slopeToSpeed, 0., maxWaterSpeed);

	
	data.z = speed;

	// Deposit/pick up sediment
	float streamPower = pow(data.x, 0.8) * speed;
	//outDeltaH[id] = data.x;
	//float streamPower = pow(inDeltaH[id], 0.8) * speed;


	const float sedimentToHeight = 1.;

	float hardnessFactor = GetHardness(data.w, id);

	int erosion_criteria = 0;
	//float delta = 0.;
	//data.y = 0.;
	// Deposit
	if (data.y > streamPower) {

		float deposit = (data.y - streamPower) * globalSpeed * depositionSpeed;
		outDeltaH[id] = deposit * sedimentToHeight;
		//delta = deposit * sedimentToHeight;
		data.y -= deposit;
	}
	// Erode
	else {
		erosion_criteria = 1;
		outDeltaH[id] = -(streamPower - data.y) * erosionSpeed * globalSpeed;// *hardnessFactor;
		//delta = -(streamPower - data.y) * erosionSpeed * globalSpeed;// *hardnessFactor;
		data.y += (streamPower - data.y) * depositionSpeed * globalSpeed;//  *hardnessFactor;
		
	}

	// Evaporate
	//data.x -= evaporation * cellArea;

	// Limit amount of water
	//data.x = max(data.x, 0.f);
	

	// Apply preprocessed deltaH
	float delta = inDeltaH[id];

	// erosion
	if (delta < 0) {
		float new_sed_h = sed_hf[id] + delta;  // erode sediment layer
		out_sed_hf[id] = new_sed_h;
		if (new_sed_h < 0) {
			data.w += new_sed_h * hardnessFactor;  // continue erosion on bedrock layer
			out_sed_hf[id] = 0.;
			data.w = max(a.z, data.w);
		}
	}
	// deposit
	else {
		out_sed_hf[id] = sed_hf[id] + delta;
	}

	// Write data to output buffers
	WriteAllData(data, streamPower, id);
	out_hardness[id] = hardnessFactor;
	out_intDebug[id] += erosion_criteria;
	out_floatDebug[id] = inDeltaH[id];
}

#endif

// x : water
// y : sediment
// z : speed
// w : height
