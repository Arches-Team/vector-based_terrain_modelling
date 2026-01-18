#version 450 core

#ifdef COMPUTE_SHADER

layout(binding = 0, std430) coherent buffer Out
{
	float outBuffer[];
};
layout(binding = 1, std430) coherent buffer InOutVec
{
	vec2 outBufferVec[];
};
layout(binding = 2, std430) readonly buffer InHeightBuffer
{
	float hf[];
};
layout(binding = 3, std430) coherent buffer OutIntBuffer
{
	int outInt;
};
layout(binding = 4, std430) coherent buffer InHeightBufferDouble
{
	double hfDouble[];
};
layout(binding = 5, std430) coherent buffer OutHeightBufferDouble
{
	double outBufferDouble[];
};

uniform vec3 a;
uniform vec3 b;
uniform vec2 cellDiagonal;
uniform dvec2 cellDiagonalDouble;
uniform int nx;
uniform int ny;
uniform float lipschitz; // Lipschitz constant of the terrain

uniform float r;
uniform int N;

uniform float epsilon = 0.005f; // 5 millimeters
uniform sampler2D heightfield;
uniform int useDouble;
uniform ivec2 view_point = ivec2(0, 0);

uniform int func;
uniform float sFrac;
uniform int m = 50; // Mask size for fractional computation

uniform float s;
uniform float aSoft = 0.1f;
uniform float diffusepower = 4.0f; // Gundy
uniform vec3 light;
uniform int softShadowSteps = 256;

uniform int groupSize = 64;

const float M_PI = 3.14159265359f;

// This blue noise in disk array is used for sampling the sun disk.
// these were generated using a modified mitchell's best candidate algorithm.
// 1) It was not calculated on a torus (no wrap around distance for points)
// 2) Candidates were forced to be in the unit circle (through rejection sampling)
const vec2 BlueNoiseInDisk[64] = vec2[64](
	vec2(0.478712, 0.875764),
	vec2(-0.337956, -0.793959),
	vec2(-0.955259, -0.028164),
	vec2(0.864527, 0.325689),
	vec2(0.209342, -0.395657),
	vec2(-0.106779, 0.672585),
	vec2(0.156213, 0.235113),
	vec2(-0.413644, -0.082856),
	vec2(-0.415667, 0.323909),
	vec2(0.141896, -0.939980),
	vec2(0.954932, -0.182516),
	vec2(-0.766184, 0.410799),
	vec2(-0.434912, -0.458845),
	vec2(0.415242, -0.078724),
	vec2(0.728335, -0.491777),
	vec2(-0.058086, -0.066401),
	vec2(0.202990, 0.686837),
	vec2(-0.808362, -0.556402),
	vec2(0.507386, -0.640839),
	vec2(-0.723494, -0.229240),
	vec2(0.489740, 0.317826),
	vec2(-0.622663, 0.765301),
	vec2(-0.010640, 0.929347),
	vec2(0.663146, 0.647618),
	vec2(-0.096674, -0.413835),
	vec2(0.525945, -0.321063),
	vec2(-0.122533, 0.366019),
	vec2(0.195235, -0.687983),
	vec2(-0.563203, 0.098748),
	vec2(0.418563, 0.561335),
	vec2(-0.378595, 0.800367),
	vec2(0.826922, 0.001024),
	vec2(-0.085372, -0.766651),
	vec2(-0.921920, 0.183673),
	vec2(-0.590008, -0.721799),
	vec2(0.167751, -0.164393),
	vec2(0.032961, -0.562530),
	vec2(0.632900, -0.107059),
	vec2(-0.464080, 0.569669),
	vec2(-0.173676, -0.958758),
	vec2(-0.242648, -0.234303),
	vec2(-0.275362, 0.157163),
	vec2(0.382295, -0.795131),
	vec2(0.562955, 0.115562),
	vec2(0.190586, 0.470121),
	vec2(0.770764, -0.297576),
	vec2(0.237281, 0.931050),
	vec2(-0.666642, -0.455871),
	vec2(-0.905649, -0.298379),
	vec2(0.339520, 0.157829),
	vec2(0.701438, -0.704100),
	vec2(-0.062758, 0.160346),
	vec2(-0.220674, 0.957141),
	vec2(0.642692, 0.432706),
	vec2(-0.773390, -0.015272),
	vec2(-0.671467, 0.246880),
	vec2(0.158051, 0.062859),
	vec2(0.806009, 0.527232),
	vec2(-0.057620, -0.247071),
	vec2(0.333436, -0.516710),
	vec2(-0.550658, -0.315773),
	vec2(-0.652078, 0.589846),
	vec2(0.008818, 0.530556),
	vec2(-0.210004, 0.519896)
);

float Remap(float x, float oldMin, float oldMax, float newMin, float newMax) {
	return newMin + (newMax - newMin) * ((x - oldMin) / (oldMax - oldMin));
}

int ToIndex1D(int i, int j) {
	return i + nx * j;
}

// Box signed distance field
// p : Point
// va, vb : Vertices of the box
float Box(vec3 p, vec3 va, vec3 vb) {
	vec3 c = .5 * (va + vb);
	vec3 r = .5 * (vb - va);
	vec3 q = abs(p - c) - r;
	float d = length(max(q, 0.)) + min(max(q.x, max(q.y, q.z)), 0.);
	return d;
}

// Box signed distance field
// p : Point
// va, vb : Vertices of the box
float Box(vec2 p, vec2 va, vec2 vb) {
	vec2 c = .5 * (va + vb);
	vec2 r = .5 * (vb - va);
	vec2 q = abs(p - c) - r;
	float d = length(max(q, 0.)) + min(max(q.x, q.y), 0.);
	return d;
}

float Intersection(float va, float vb) {
	return max(va, vb);
}

double at(int i, int j) {
	return hfDouble[ToIndex1D(i, j)];
}

float Height(vec2 p) {
	vec2 q = p - a.xy;
	vec2 d = b.xy - a.xy;
	vec2 uv = q / d;
	vec4 col = texture(heightfield, uv);
	float zu = col.b + col.g / 255.f; // Gundy ? Was 256.0
	return Remap(zu, 0.f, 1.f, a.z, b.z);
}

float Height(int i, int j) {
	vec2 uv = vec2(float(i) / (nx - 1), float(j) / (ny - 1));
	vec4 col = texture(heightfield, uv);
	float zu = col.b + col.g / 255.f; // Gundy ? Was 256.0
	return Remap(zu, 0.f, 1.f, a.z, b.z);
}

vec2 ArrayPoint(int i, int j) {
	return a.xy + vec2(i, j) * cellDiagonal;
}

vec3 Vertex(int i, int j) {
	return vec3(ArrayPoint(i, j), Height(i, j));
}

// Normal at point.
// i, j : Integer coordinates.
vec3 Normal(int i, int j) 
{
	vec2 p = ArrayPoint(i, j);
	const vec3 e = vec3(cellDiagonal, 0.f);
	return normalize(vec3(
			-(Height(i+1,j) - Height(i-1,j)), // Gundy ? Was below, check if ok
			-(Height(i,j+1) - Height(i,j-1)),
			length(e.xy)
		)
	);

/*
return normalize(vec3(
			-(Height(p + e.xz) - Height(p - e.xz)),
			-(Height(p + e.zy) - Height(p - e.zy)),
			length(e.xy)
		)
	);
*/
}

// Implicit elevation.
// p : Point.
float Implicit(vec3 p) 
{
	return p.z - Height(p.xy);
}

bool IntersectBox(vec3 ro, vec3 rd, out float tN, out float tF)
{
	vec3 rinvDir = 1. / rd;
	// Increase vertical size of the bounding box of the terrain
	float delta = .05 * abs(b.z - a.z);
	vec3 tbot = rinvDir * (vec3(a.x, a.y, a.z - delta) - ro);
	vec3 ttop = rinvDir * (vec3(b.x, b.y, b.z + delta) - ro);
	vec3 tmin = min(ttop, tbot);
	vec3 tmax = max(ttop, tbot);
	vec2 t = max(tmin.xx, tmin.yz);
	float t0 = max(t.x, t.y);
	t = min(tmax.xx, tmax.yz);
	float t1 = min(t.x, t.y);
	tN = t0;
	tF = t1;
	return t1 > max(t0, 0.);
}

// Sphere tracing 
// Returns if intersection occurs
// ro : Ray origin
// rd : Direction
//  r : Mximum marching distance
//  t : Intersection depth
bool SphereTrace(vec3 ro, vec3 rd, float r, out float t) 
{
	t = 0.f;

	// Compute maximum slope from Lipschitz constant
	float lambda = sqrt(lipschitz * lipschitz + 1.0);

	//float uz = abs(rd.z);
	//float kr = uz + lipschitz * sqrt(1.0f - (uz * uz));

	// Test if ray is more vertical that the maximum slope, as in Genevaux 2015
	//if (rd.z > lambda * length(rd.xy))
	//	return false;

	float ta = 0.0, tb = 0.0;
	if (!IntersectBox(ro, rd, ta, tb))
		return false;

	tb = min(tb, r);
	t = max(ta + epsilon, 0.0f);
	for (int i = 0; i < 2256; i++)
	{
		if (t > tb)
			break;
		vec3 p = ro + rd * t;
		float d = Implicit(p);
		if (d < 0.0)
			return true;
		t += max(d / lambda, epsilon);
	}
	return false;
}

// Hashing 
// Returns a random number in [-1,1]
// s : Random seed
float Hash(float s) {
	return fract(sin(s) * 43758.5453);
}

// Uniform sphere sampling
// fi : Sample identifier
vec3 SphereSampling(float fi) 
{
	float theta = 2.0f * M_PI * Hash(fi * 2.18622f);
	float phi = acos(1.0f - 2.0f * Hash(fi * .94312f));
	float x = sin(phi) * cos(theta);
	float y = sin(phi) * sin(theta);
	float z = cos(phi);
	vec3 d = vec3(x, y, z);
	return d;
}

// Uniform hemisphere sampling
// fi : Sample identifier
//  n : Normal
vec3 HemisphereSampling(float fi, vec3 n) 
{
	vec3 d = SphereSampling(fi);
	if (dot(d, n) < 0.0)
		d = -d;
	return d;
}

// Clear sky
// Hemisphere sampling with up direction
// i, j : Integer coordinates
float ClearSky(int i, int j) 
{
	// Shift point vertically
	vec3 p = Vertex(i, j) + vec3(0.0, 0.0, epsilon * 10.f);

	// Contrary to ambient occlusion or accessibility, use vertical axis
	vec3 n = vec3(0.0, 0.0, 1.0);
	float cosine = 0.0;
	float sum = 0.0;
	for (int k = 0; k < N; k++)
	{
		vec3 d = HemisphereSampling(float(k), n);

		float c = dot(n, d);
		float t;
		// Accumulate cosine if not intersection occurs
		if (!SphereTrace(p, d, r, t))
			cosine += c;
		sum += c;
	}
	return cosine / sum;
}

// Accessibility 
// Uniform hemisphere sampling oriented along normal.
// i, j: Integer coordinates
float Accessibility(int i, int j) 
{
	
	// Retrieve normal
	vec3 n = Normal(i, j);

	// Shift point in the direction of the normal
	vec3 p = Vertex(i, j) + epsilon * 10.f*n;

	int hit = 0;
	for (int k = 0; k < N; k++)
	{
		vec3 d = HemisphereSampling(float(k), n);
		float t;
		// Accumulate hit
		if (SphereTrace(p, d, r, t))
			hit++;
	}
	return 1.0f - (float(hit) / float(N));
}

float Viewshed(int i, int j) {
	vec3 p = Vertex(view_point.x, view_point.y);// +vec3(0., 0., 10000.);
	vec3 q = Vertex(i, j) + vec3(0., 0., 1.8);
	vec3 d = p - q;
	float t;
	SphereTrace(q, normalize(d), 1000000000, t);
	float error = length(p - (q + normalize(d) * t));
	if (error < 250.) return 1.;
	return 0.;
}

// Fractional Laplacian
// x, y : Integer coordinates
float LaplacianFractional(int x, int y) {
	float sum = 0.f;

	// Factorize value at center
	float hxy = 2.f * Height(x, y);

	// Area of integrand : not used as we often work on scaled Fractional Laplacian
	// double dxy=cellDiagonal.x*cellDiagonal.y;

	// Reference implementation
	/*
	for(int i=-size;i<=size;i++)
	{
		for(int j=-size;j<=size;j++)
		{
			// Skip central point for which value is undefined (distribution)
			if((i==0)&&(j==0))
			continue;
			sum+=(Height(x+i,y+j)-hxy+Height(x-i,y-j))/pow(length(vec2(float(i),float(j))),2.+2.*sLapl);
		}
	}
	// The constant C(2,s) is not computed in the shader
	return -0.5f*sum;
	*/

	// Integration on half of the domain RxR-* allows to remove the test inside the loop
	for (int i = -m; i <= m; i++)
	{
		for (int j = -m; j < 0; j++)
		{
			sum += (Height(x + i, y + j) - hxy + Height(x - i, y - j)) / pow(length(vec2(float(i), float(j))), 2. + 2. * sFrac);
		}
	}
	// Complete integration on R- x {0}
	for (int i = -m; i < 0; i++)
		sum += (Height(x + i, y + 0) - hxy + Height(x - i, y - 0)) / pow(length(vec2(float(i), float(0))), 2. + 2. * sFrac);

	// The constant C(2,s) is not computed in the shader
	return -sum;
}

// Fractional Gradient
// x, y : Integer coordinates
vec2 GradientFractional(int x, int y) 
{
	vec2 sum = vec2(0.0, 0.0);

	// Factorize value at center
	float hxy = Height(x, y);

	// Reference implementation
	for (int i = -m; i <= m; i++)
	{
		for (int j = -m; j <= m; j++)
		{
			// Skip central point for which value is undefined (distribution)
			if ((i == 0) && (j == 0))
				continue;
			// Vector
			vec2 vij = vec2(float(i), float(j));
			// Its length
			float lvij = length(vij);
			sum += ((hxy - Height(x + i, y + j)) / pow(lvij, 2. + sFrac)) * vij / lvij;
		}
	}
	// The constant C(2, s) is not computed in the shader
	return sum;
}

// Gradient
// i, j: Integer coordinates
dvec2 Gradient(int i, int j) {
	dvec2 n;

	// Gradient along x axis
	if (i == 0)
	{
		n.x = (at(i + 1, j) - at(i, j)) / cellDiagonalDouble[0];
	}
	else if (i == nx - 1)
	{
		n.x = (at(i, j) - at(i - 1, j)) / cellDiagonalDouble[0];
	}
	else
	{
		n.x = (at(i + 1, j) - at(i - 1, j)) / (2.0 * cellDiagonalDouble[0]);
	}

	// Gradient along y axis
	if (j == 0)
	{
		n.y = (at(i, j + 1) - at(i, j)) / cellDiagonalDouble[1];
	}
	else if (j == ny - 1)
	{
		n.y = (at(i, j) - at(i, j - 1)) / cellDiagonalDouble[1];
	}
	else
	{
		n.y = (at(i, j + 1) - at(i, j - 1)) / (2.0 * cellDiagonalDouble[1]);
	}

	return n;
}

// Laplacian
// i, j: Integer coordinates
double Laplacian(int i, int j) {
	double laplacian = 0.0;

	// Divergence along x axis
	if (i == 0)
	{
		laplacian += (at(i, j) - 2.0 * at(i + 1, j) + at(i + 2, j)) / (cellDiagonalDouble[0] * cellDiagonalDouble[0]);
	}
	else if (i == nx - 1)
	{
		laplacian += (at(i, j) - 2.0 * at(i - 1, j) + at(i - 2, j)) / (cellDiagonalDouble[0] * cellDiagonalDouble[0]);
	}
	else
	{
		laplacian += (at(i + 1, j) - 2.0 * at(i, j) + at(i - 1, j)) / (cellDiagonalDouble[0] * cellDiagonalDouble[0]);
	}

	// Divergence along y axis
	if (j == 0)
	{
		laplacian += (at(i, j) - 2.0 * at(i, j + 1) + at(i, j + 2)) / (cellDiagonalDouble[1] * cellDiagonalDouble[1]);
	}
	else if (j == ny - 1)
	{
		laplacian += (at(i, j) - 2.0 * at(i, j - 1) + at(i, j - 2)) / (cellDiagonalDouble[1] * cellDiagonalDouble[1]);
	}
	else
	{
		laplacian += (at(i, j + 1) - 2.0 * at(i, j) + at(i, j - 1)) / (cellDiagonalDouble[1] * cellDiagonalDouble[1]);
	}

	return laplacian;
}

// Count the amount of pits of a given cell
// x, y : Integer coordinates
int CountPitsDouble(int x, int y)  {
	const int DIRECTIONS = 8;
	const int dx[DIRECTIONS] = { -1, -1, 0, 1, 1, 1, 0, -1 };
	const int dy[DIRECTIONS] = { 0, -1, -1, -1, 0, 1, 1, 1 };
	if (x == 0 || x == (nx - 1) || y == 0 || y == (ny - 1)) 
		return 0;
	double lowest_neighbour = 1e8f;
	for (int n = 0; n < DIRECTIONS; n++) {
		int px = x + dx[n];
		int py = y + dy[n];
		lowest_neighbour = min(hfDouble[ToIndex1D(px, py)], lowest_neighbour);
	}
	if (hfDouble[ToIndex1D(x, y)] <= lowest_neighbour)
		return 1;
	return 0;
}

// Count the amount of pits of a given cell
// x, y : Integer coordinates
int CountPitsFloat(int x, int y)  {
	const int DIRECTIONS = 8;
	const int dx[DIRECTIONS] = { -1, -1, 0, 1, 1, 1, 0, -1 };
	const int dy[DIRECTIONS] = { 0, -1, -1, -1, 0, 1, 1, 1 };
	if (x == 0 || x == (nx - 1) || y == 0 || y == (ny - 1)) 
		return 0;
	float lowest_neighbour = 1e8f;
	for (int n = 0; n < DIRECTIONS; n++) {
		int px = x + dx[n];
		int py = y + dy[n];
		lowest_neighbour = min(hf[ToIndex1D(px, py)], lowest_neighbour);
	}
	if (hf[ToIndex1D(x, y)] <= lowest_neighbour)
		return 1;
	return 0;
}

// Compute a perpendicular direction
// u : Direction
vec3 perp(vec3 u)
{
	vec3 a = abs(u);
	vec3 v;
	if (a.x <= a.y && a.x <= a.z)
		v = vec3(0, -u.z, u.y);
	else if (a.y <= a.x && a.y <= a.z)
		v = vec3(-u.z, 0, u.x);
	else
		v = vec3(-u.y, u.x, 0);
	return v;
}

// Soft shadows
// x, y : Integer coordinates
float SoftShadows(int x, int y) 
{
	// Retrieve normal
	vec3 no = Normal(x, y);

	// Shift point in the direction of the normal
	vec3 ro = Vertex(x, y) + epsilon * 10.f*no;

	//vec3 ro = Vertex(x, y) + vec3(0.0, 0.0, 5.0 * epsilon);
	vec3 rd = -normalize(light);
	float t;
	vec3 p;
	vec3 u = normalize(perp(rd));
	vec3 v = normalize(cross(rd, u));
	int n = 0;
	for (int i = 0; i < softShadowSteps; i++)
	{
		vec3 dir = normalize(rd + aSoft * (BlueNoiseInDisk[i % 64].x * u + BlueNoiseInDisk[i % 64].y * v));
		bool hit = SphereTrace(ro, dir, 1000000000, t);
		if (hit)
			n++;
	}
	return float(n) / float(softShadowSteps);
}

// Shadow
// x, y : Integer coordinates
float Shadow(int x, int y) 
{
	// Retrieve normal
	vec3 n = Normal(x, y);

	// Shift point in the direction of the normal
	vec3 ro = Vertex(x, y) + epsilon * 10.f*n;
	
	//vec3 ro = Vertex(x, y) + vec3(0.0, 0.0, 10.0 * epsilon);
	vec3 rd = -normalize(light);
	float t;
	vec3 p;

		vec3 dir = normalize(rd);
		bool hit = SphereTrace(ro, dir, 1000000000,t);
		if (hit)
			return 1.0;
	
	return 0.0;
}

float Diffuse(int i,int j)
{
  // Normal
      vec3 n = Normal(i, j);
      // Cosine
	  float s = dot(n, -normalize(light));
      s = 0.5 * (1.0 + s);
      s = pow(s,diffusepower);

     return s;
 }

// Get range
// start, end: indexes
vec2 GetRange(int start, int end) 
{
	vec2 range = vec2(100000.0f, -100000.0f);
	for (int i = start; i < end; i++) {
		if (i > (nx * ny))
			continue;
		float z = hf[i];
		range.x = min(range.x, z);
		range.y = max(range.y, z);
	}
	return range;
}


layout(local_size_x = 8, local_size_y = 8, local_size_z = 1)in;
void main()
{
	int i = int(gl_GlobalInvocationID.x);
	int j = int(gl_GlobalInvocationID.y);
	if (i < 0) return;
	if (j < 0) return;
	if (i >= nx) return;
	if (j >= ny) return;

	// Ambient occlusion
	if (func == 0)
	{
		outBuffer[ToIndex1D(i, j)] = Accessibility(i, j);
	}
	// Fractionnal Laplacian
	else if (func == 1)
	{
		outBuffer[ToIndex1D(i, j)] = LaplacianFractional(i, j);
	}
	// Fractionnal gradient
	else if (func == 2)
	{
		vec2 g = GradientFractional(i, j);
		outBuffer[ToIndex1D(i, j)] = length(g);
		outBufferVec[ToIndex1D(i, j)] = g;
	}
	// Clear sky
	else if (func == 3)
	{
		outBuffer[ToIndex1D(i, j)] = ClearSky(i, j);
	}
	// Count pits
	else if (func == 4)
	{
		int pits;
		if (useDouble == 1)
			pits = CountPitsDouble(i, j);
		else
			pits = CountPitsFloat(i, j);
		atomicAdd(outInt, pits);
	}
	// Soft shadows
	else if (func == 5)
	{
		outBuffer[ToIndex1D(i, j)] = 1.0f - SoftShadows(i, j); 
	}
	// Get range (part 1)
	else if (func == 6)
	{
		if (i >= groupSize)
			return;

		// Each thread computes its part
		int part = ((nx * ny) / groupSize);
		int start = i * part;
		int end = start + part;
		outBufferVec[i] = GetRange(start, end);
	}
	// Get range (part 2)
	else if (func == 7)
	{
		if (i > 0)
			return;
		for (int k = 0; k < groupSize; k++) {
			outBufferVec[0] = vec2(
				min(outBufferVec[0].x, outBufferVec[k].x),
				max(outBufferVec[0].y, outBufferVec[k].y)
			);
		}
	}
	// Gradient & gradient norm
	else if (func == 8)
	{
		dvec2 g = Gradient(i, j);
		outBufferDouble[ToIndex1D(i, j)] = length(g);
		outBufferVec[ToIndex1D(i, j)] = vec2(float(g.x), float(g.y));
	}
	else if (func == 9)
	{
		outBufferDouble[ToIndex1D(i, j)] = Laplacian(i, j);
	}
	else if (func == 10) {
		outBuffer[ToIndex1D(i, j)] = Viewshed(i, j);
	}
	else if (func == 11) {
		outBuffer[ToIndex1D(i, j)] += Viewshed(i, j);
	}
	else if (func == 12) {
		outBuffer[ToIndex1D(i, j)] = Diffuse(i, j);
	}
	else if (func == 13) {
		outBuffer[ToIndex1D(i, j)] = Shadow(i, j);
	}
}

#endif
