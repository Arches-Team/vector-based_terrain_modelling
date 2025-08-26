#version 450 core

#ifdef COMPUTE_SHADER

// In
layout(binding = 0, std430) readonly buffer InElevation
{
	float hf[];
};
layout(binding = 1, std430) readonly buffer InWater
{
	float water[];
};
layout(binding = 2, std430) readonly buffer InSediment
{
	float sediments[];
};

// Out
layout(binding = 3, std430) writeonly buffer OutElevation
{
	float out_hf[];
};
layout(binding = 4, std430) writeonly buffer OutWater
{
	float out_water[];
};
layout(binding = 5, std430) writeonly buffer OutSediment
{
	float out_sediments[];
};
layout(binding = 6, std430) writeonly buffer OutWaterSpeed
{
	float out_speeds[];
};

// Options
layout(binding = 7, std430) readonly buffer HardnessBuffer
{
	float hardness[];
};
layout(binding = 8, std430) writeonly buffer OutHardnessBuffer
{
	float out_hardness[];
};
layout(binding = 9, std430) writeonly buffer OutFullHeight
{
	float out_fullHeight[];
};

// Heightfield data
uniform int nx;
uniform int ny;
uniform vec2 a;
uniform vec2 b;
uniform vec2 cellDiag;

// Hydraulic erosion parameters
uniform float globalSpeed = 1.f;
uniform float erosionSpeed = 1.f;
uniform float depositionSpeed = 1.f;
uniform float evaporation = .05f;
uniform float rain = .16f;

// Hardness
uniform int useHardness = 0;
uniform int hardnessFunctionIndex = 0;

// Local editing
uniform bool local = false;
uniform vec2 center;
uniform float radius;

// Constants
const float flowRate = .32f;
const float diagonalDist = 1.f / sqrt(2.);
float scalexy;

float random(vec2 st) { return fract(sin(dot(st.xy, vec2(12.9898, 78.233))) * 43758.5453123); }
vec3 mod289(vec3 x) { return x - floor(x * (1. / 289.)) * 289.; }
vec4 mod289(vec4 x) { return x - floor(x * (1. / 289.)) * 289.; }
vec4 permute(vec4 x) { return mod289(((x * 34.) + 1.) * x); }
vec4 taylorInvSqrt(vec4 r) { return 1.79284291400159 - .85373472095314 * r; }
float snoise(vec3 v) {
	const vec2 C = vec2(1. / 6., 1. / 3.); const vec4 D = vec4(0., .5, 1., 2.);
	vec3 i = floor(v + dot(v, C.yyy)); vec3 x0 = v - i + dot(i, C.xxx); vec3 g = step(x0.yzx, x0.xyz); vec3 l = 1. - g; vec3 i1 = min(g.xyz, l.zxy); vec3 i2 = max(g.xyz, l.zxy);
	vec3 x1 = x0 - i1 + C.xxx; vec3 x2 = x0 - i2 + C.yyy; vec3 x3 = x0 - D.yyy;
	i = mod289(i); vec4 p = permute(permute(permute(i.z + vec4(0., i1.z, i2.z, 1.)) + i.y + vec4(0., i1.y, i2.y, 1.)) + i.x + vec4(0., i1.x, i2.x, 1.));
	float n_ = .142857142857; vec3 ns = n_ * D.wyz - D.xzx; vec4 j = p - 49. * floor(p * ns.z * ns.z);//  mod(p,7*7)
	vec4 x_ = floor(j * ns.z); vec4 y_ = floor(j - 7. * x_); vec4 x = x_ * ns.x + ns.yyyy; vec4 y = y_ * ns.x + ns.yyyy; vec4 h = 1. - abs(x) - abs(y);
	vec4 b0 = vec4(x.xy, y.xy); vec4 b1 = vec4(x.zw, y.zw); vec4 s0 = floor(b0) * 2. + 1.; vec4 s1 = floor(b1) * 2. + 1.; vec4 sh = -step(h, vec4(0.));
	vec4 a0 = b0.xzyw + s0.xzyw * sh.xxyy; vec4 a1 = b1.xzyw + s1.xzyw * sh.zzww; vec3 p0 = vec3(a0.xy, h.x); vec3 p1 = vec3(a0.zw, h.y); vec3 p2 = vec3(a1.xy, h.z); vec3 p3 = vec3(a1.zw, h.w);
	vec4 norm = taylorInvSqrt(vec4(dot(p0, p0), dot(p1, p1), dot(p2, p2), dot(p3, p3))); p0 *= norm.x; p1 *= norm.y; p2 *= norm.z; p3 *= norm.w;
	vec4 m = max(.6 - vec4(dot(x0, x0), dot(x1, x1), dot(x2, x2), dot(x3, x3)), 0.); m = m * m;
	return 42. * dot(m * m, vec4(dot(p0, x0), dot(p1, x1), dot(p2, x2), dot(p3, x3)));
}

float ridge(vec3 v) {
	return 1.f - abs(snoise(v) * .5 + .5);
}

float Height(int index) {
	return hf[index];
}

int Side(vec3 p, vec3 pp, vec3 n) {
	float c = dot(pp, n);
	float r = dot(n, p) - c;
	if (r > .0000001f)
		return 1;
	else if (r < -.0000001f)
		return-1;
	else
		return 0;
}

float GetHardness(in vec3 p, int index) {
	float h = 1.f;

	// Strata function
	if (hardnessFunctionIndex == 0)
	{
		p *= .025f;
		h = sin(p.z) * .5f + .5f;
	}

	// Warped strata function
	if (hardnessFunctionIndex == 1)
	{
		float wa = 3. * snoise(p * .000004f);
		p *= .025f;
		h = sin(p.z + wa) * .5f + .5f;
	}

	// 3D Noise
	if (hardnessFunctionIndex == 2)
	{
		h = (snoise(p * .000002f) * .5f + .5f);
		h *= 2.;
	}

	// Then modulate the effect with a user map, if enabled
	if (useHardness == 1)
		return h;
	else
		return 1.f;
}

float GetHardness(float height, vec2 uv)
{
	//return 0.5;
	float offset = -.2 * sin(6.283185 * uv.x) + .2 * sin(6.283185 * uv.y);

	float a = (height - offset) * 6.283185 * 12.5;
	return mix(mix(sin(a) * .5 + .5, sin(a * 1.618 + dot(uv, vec2(4, 2)) * 6.283185) * .5 + .5, .5),
		sin(a * 6.18) * .5 + .5, .2);
}

int ToIndex1D(int i, int j) {
	return i + nx * j;
}

vec2 ArrayPoint(ivec2 p) {
	return a.xy + vec2(p) * cellDiag;
}

vec3 Point(ivec2 p)
{
	return vec3(ArrayPoint(p), hf[ToIndex1D(p.x, p.y)]);
}

bool IsInDomain(vec2 p) {
	return (length(p - center) < radius);
}

vec4 GetAllData(int id)
{
	vec4 ret;
	ret.w = hf[id];// Bedrock elevation
	ret.x = water[id];// Amount of water
	ret.y = sediments[id];// Amount of sediment
	ret.z = 0.f;// Flow speed - which is just write-only
	return ret;
}

void WriteAllData(vec4 frag, float hardness, int id)
{
	out_hf[id] = frag.w;
	out_water[id] = frag.x;
	out_sediments[id] = frag.y;
	out_speeds[id] = frag.z;
	out_hardness[id] = hardness;
	out_fullHeight[id] = frag.w + frag.y;
}

void TestOutFlow(inout float lowestSlope, inout int lowestIndex, in int indexMe, in int indexNeighbour, in float dist, in vec4 samples[25])
{
	float slope = (samples[indexNeighbour].w - samples[indexMe].w) * scalexy * dist;
	if (slope < lowestSlope)
	{
		lowestSlope = slope;
		lowestIndex = indexNeighbour;
	}
}

void ProcessNeighbour(inout float outSlope, inout float inSlope, inout vec2 inFlow, in ivec2 index, in float dist, in vec4 samples[25])
{
	int idx = index.y * 5 + index.x;
	int meIdx = 2 * 5 + 2;

	vec4 me = samples[meIdx];
	vec4 neighbour = samples[idx];

	// Flux sortant
	// On calcule juste que ca part, on s'en fout de ou, ce n'est pas le but car on n'updatera pas les voisins
	float slope = (neighbour.w - me.w) * scalexy * dist;
	outSlope = min(outSlope, slope);

	// Flux entrant : on calcule a la fois si ca rentre, et d'ou cela viendrait
	float neighbourOutSlope = 0.;
	int neighbourOutSlopeTarget = -1;

	TestOutFlow(neighbourOutSlope, neighbourOutSlopeTarget, idx, idx - 5 - 1, diagonalDist, samples);
	TestOutFlow(neighbourOutSlope, neighbourOutSlopeTarget, idx, idx - 5 + 1, diagonalDist, samples);
	TestOutFlow(neighbourOutSlope, neighbourOutSlopeTarget, idx, idx + 5 - 1, diagonalDist, samples);
	TestOutFlow(neighbourOutSlope, neighbourOutSlopeTarget, idx, idx + 5 + 1, diagonalDist, samples);
	TestOutFlow(neighbourOutSlope, neighbourOutSlopeTarget, idx, idx - 5, 1., samples);
	TestOutFlow(neighbourOutSlope, neighbourOutSlopeTarget, idx, idx + 5, 1., samples);
	TestOutFlow(neighbourOutSlope, neighbourOutSlopeTarget, idx, idx - 1, 1., samples);
	TestOutFlow(neighbourOutSlope, neighbourOutSlopeTarget, idx, idx + 1, 1., samples);

	if (neighbourOutSlopeTarget == meIdx)
	{
		inSlope += max(-slope, 0.);// do I want to count this for ones flowing away from me?
		inFlow += neighbour.xy;
	}
}

layout(local_size_x = 8, local_size_y = 8, local_size_z = 1)in;
void main()
{
	scalexy = 1. / ((b - a).x / float(nx - 1));

	int x = int(gl_GlobalInvocationID.x);
	int y = int(gl_GlobalInvocationID.y);
	if (x < 0) return;
	if (y < 0) return;
	if (x >= nx) return;
	if (y >= ny) return;

	int id = ToIndex1D(x, y);
	vec4 fragColour = GetAllData(id);

	// Sample a 5x5 grid around the pixel
	vec4 samples[25];
	for (int j = 0; j < 5; j++)
	{
		for (int i = 0; i < 5; i++)
		{
			ivec2 tapUV = (ivec2(x, y) + ivec2(i, j) - ivec2(2, 2) + ivec2(nx, ny)) % ivec2(nx, ny);
			samples[j * 5 + i] = GetAllData(ToIndex1D(tapUV.x, tapUV.y));
		}
	}

	float outSlope = 0.;		// height difference to lowest neighbour (negative)
	float inSlope = 0.;			// sum of height differences to higher neighbours (positive)
	vec2 inFlow = vec2(0.);		// sum of (water, sediment) from higher neighbours
	ProcessNeighbour(outSlope, inSlope, inFlow, ivec2(1, 1), diagonalDist, samples);
	ProcessNeighbour(outSlope, inSlope, inFlow, ivec2(1, 3), diagonalDist, samples);
	ProcessNeighbour(outSlope, inSlope, inFlow, ivec2(3, 1), diagonalDist, samples);
	ProcessNeighbour(outSlope, inSlope, inFlow, ivec2(3, 3), diagonalDist, samples);
	ProcessNeighbour(outSlope, inSlope, inFlow, ivec2(1, 2), 1., samples);
	ProcessNeighbour(outSlope, inSlope, inFlow, ivec2(2, 1), 1., samples);
	ProcessNeighbour(outSlope, inSlope, inFlow, ivec2(2, 3), 1., samples);
	ProcessNeighbour(outSlope, inSlope, inFlow, ivec2(3, 2), 1., samples);

	// Modify water & sediment values
	if (outSlope < 0.)
	{
		// reduce my water and sediment - they'll be passed to a neighbour
		fragColour.xy *= 1. - flowRate;
	}

	// Add sediment and water
	fragColour.xy += inFlow * flowRate;

	const float maxWater = 10.f;

	// Prevent waterfalls causing infinitely powerful erosion
	// This gives really nice control over how slopey (big values) or cliffy (small values) the terrain is.
	const float maxWaterSpeed = 3.f;

	// Compute speed // was 1.0/0.01
	const float slopeToSpeed = 1.f / 1.f;	// denominator is height difference between input and output points at which we get a waterfall
	float speed = clamp(-outSlope * slopeToSpeed, 0., maxWaterSpeed);	// USE ONLY OUT - so water in dips is static! oops!

	// Evaporate
	fragColour.x -= evaporation;

	// Limit amount of water
	fragColour.x = min(fragColour.x, maxWater);
	fragColour.x = max(fragColour.x, 0.f);
	fragColour.z = speed;

	// Deposit/pick up sediment
	float maxSediment = fragColour.x * speed;
	const float sedimentToHeight = 0.1;

	float hardnessFactor = GetHardness(Point(ivec2(x, y)), id);
	hardnessFactor = GetHardness(fragColour.w, vec2(x, y) / vec2(nx, ny));

	// Deposit
	if (fragColour.y > maxSediment)
	{
		float deposit = (fragColour.y - maxSediment) * globalSpeed * depositionSpeed;
		fragColour.w += deposit;
		fragColour.y -= deposit * sedimentToHeight;
	}
	// Erode
	else
	{
		const float sedimentPickUp = .1 * globalSpeed * erosionSpeed;
		float pickupRate = sedimentPickUp * hardnessFactor;

		float pickup = min(min(maxSediment - fragColour.y, fragColour.w / sedimentToHeight), pickupRate);
		fragColour.w -= pickup;
		fragColour.y += pickup * sedimentToHeight;
	}

	// Rain
	if (!local || (local && IsInDomain(ArrayPoint(ivec2(x, y)))))
	{
		fragColour.x += rain;
	}

	// Debug to show value
	//fragColour.x = speed;

	// Write data to output buffers
	WriteAllData(fragColour, hardnessFactor, id);
}

#endif

// x : water
// y : sediment
// z : speed
// w : height
