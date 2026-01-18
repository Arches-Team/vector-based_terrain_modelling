#version 430 core

#define M_PI 3.1415926535897932384626433832795

#ifdef VERTEX_SHADER
void main(void)
{
	vec4 vertices[4] = vec4[4](vec4(-1.0, -1.0, 1.0, 1.0),
                               vec4( 1.0, -1.0, 1.0, 1.0),
                               vec4(-1.0,  1.0, 1.0, 1.0),
                               vec4( 1.0,  1.0, 1.0, 1.0));
    vec4 pos = vertices[gl_VertexID];
    gl_Position = pos;
} 
#endif

#ifdef FRAGMENT_SHADER
// Camera data
layout(location = 0) uniform vec3 CamPos;
layout(location = 1) uniform vec3 CamLookAt;
layout(location = 2) uniform vec3 CamUp;
layout(location = 3) uniform float CamAngleOfViewV;
layout(location = 4) uniform vec2 iResolution;
layout(location = 30) uniform vec2 cameraPlanes;
layout(location = 31) uniform vec3 skyColor;

// Heightfield data
layout(location = 5) uniform vec2 a;
layout(location = 6) uniform vec2 b;
layout(location = 7) uniform vec2 zRange;
layout(location = 8) uniform float K;
layout(location = 9) uniform ivec2 texSize;

// Shading data
layout(location = 10) uniform int useAlbedo;
layout(location = 11) uniform int useWireframe;
layout(location = 12) uniform int useCost;
layout(location = 13) uniform int useElevation;
layout(location = 14) uniform int useGreenBrownYellow;
// uniform 15 is not used atm
layout(location = 16) uniform int useShadingBuffer;

// Smooth shadows
layout(location = 17) uniform int useSmoothShadow = 0;
layout(location = 18) uniform int smoothShadowNbStep = 1;
layout(location = 19) uniform int smoothShadowMarchingSteps = 256;
layout(location = 20) uniform float smoothShadowMarchingEps = 1.0;
layout(location = 21) uniform float smoothShadowStrength = 0.25;

// Self shadows
layout(location = 22) uniform int useSelfShadow = 0;
layout(location = 23) uniform int selfShadowMarchingSteps = 256;
layout(location = 24) uniform float selfShadowMarchingEps = 1.0;
layout(location = 25) uniform float selfShadowStrength = 0.55;

// Texture
uniform sampler2D albedo;

// Buffers: elevation + shading (optional)
layout (binding = 0, std430) buffer ElevationBuffer {
	float hf[];
};
layout (binding = 1, std430) buffer ShadingBuffer {
	float shading[];
};

// Raymarching
layout(location = 26) uniform vec3 lightDir = vec3(-1.0f, -2.0f, -0.9f);
layout(location = 27) uniform int AA = 1;
layout(location = 28) uniform int STEPS = 1024;
layout(location = 29) uniform float epsilon = 0.05f;

out vec4 Fragment;

// This "blue noise in disk" array is blue noise in a circle and is used for sampling the
// sun disk for the blue noise.
// these were generated using a modified mitchell's best candidate algorithm.
// 1) It was not calculated on a torus (no wrap around distance for points)
// 2) Candidates were forced to be in the unit circle (through rejection sampling)
const vec2 BlueNoiseInDisk[64] = vec2[64](
    vec2(0.478712,0.875764),
    vec2(-0.337956,-0.793959),
    vec2(-0.955259,-0.028164),
    vec2(0.864527,0.325689),
    vec2(0.209342,-0.395657),
    vec2(-0.106779,0.672585),
    vec2(0.156213,0.235113),
    vec2(-0.413644,-0.082856),
    vec2(-0.415667,0.323909),
    vec2(0.141896,-0.939980),
    vec2(0.954932,-0.182516),
    vec2(-0.766184,0.410799),
    vec2(-0.434912,-0.458845),
    vec2(0.415242,-0.078724),
    vec2(0.728335,-0.491777),
    vec2(-0.058086,-0.066401),
    vec2(0.202990,0.686837),
    vec2(-0.808362,-0.556402),
    vec2(0.507386,-0.640839),
    vec2(-0.723494,-0.229240),
    vec2(0.489740,0.317826),
    vec2(-0.622663,0.765301),
    vec2(-0.010640,0.929347),
    vec2(0.663146,0.647618),
    vec2(-0.096674,-0.413835),
    vec2(0.525945,-0.321063),
    vec2(-0.122533,0.366019),
    vec2(0.195235,-0.687983),
    vec2(-0.563203,0.098748),
    vec2(0.418563,0.561335),
    vec2(-0.378595,0.800367),
    vec2(0.826922,0.001024),
    vec2(-0.085372,-0.766651),
    vec2(-0.921920,0.183673),
    vec2(-0.590008,-0.721799),
    vec2(0.167751,-0.164393),
    vec2(0.032961,-0.562530),
    vec2(0.632900,-0.107059),
    vec2(-0.464080,0.569669),
    vec2(-0.173676,-0.958758),
    vec2(-0.242648,-0.234303),
    vec2(-0.275362,0.157163),
    vec2(0.382295,-0.795131),
    vec2(0.562955,0.115562),
    vec2(0.190586,0.470121),
    vec2(0.770764,-0.297576),
    vec2(0.237281,0.931050),
    vec2(-0.666642,-0.455871),
    vec2(-0.905649,-0.298379),
    vec2(0.339520,0.157829),
    vec2(0.701438,-0.704100),
    vec2(-0.062758,0.160346),
    vec2(-0.220674,0.957141),
    vec2(0.642692,0.432706),
    vec2(-0.773390,-0.015272),
    vec2(-0.671467,0.246880),
    vec2(0.158051,0.062859),
    vec2(0.806009,0.527232),
    vec2(-0.057620,-0.247071),
    vec2(0.333436,-0.516710),
    vec2(-0.550658,-0.315773),
    vec2(-0.652078,0.589846),
    vec2(0.008818,0.530556),
    vec2(-0.210004,0.519896) 
);

float random(vec2 st) { return fract(sin(dot(st.xy, vec2(12.9898, 78.233))) * 43758.5453123); }
vec3 mod289(vec3 x) { return x - floor(x * (1.0 / 289.0)) * 289.0; }
vec4 mod289(vec4 x) { return x - floor(x * (1.0 / 289.0)) * 289.0; }
vec4 permute(vec4 x) { return mod289(((x * 34.0) + 1.0) * x); }
vec4 taylorInvSqrt(vec4 r) { return 1.79284291400159 - 0.85373472095314 * r; }
float snoise(vec3 v) {
    const vec2  C = vec2(1.0 / 6.0, 1.0 / 3.0);  const vec4  D = vec4(0.0, 0.5, 1.0, 2.0);
    vec3 i = floor(v + dot(v, C.yyy));  vec3 x0 = v - i + dot(i, C.xxx);  vec3 g = step(x0.yzx, x0.xyz);  vec3 l = 1.0 - g;  vec3 i1 = min(g.xyz, l.zxy);  vec3 i2 = max(g.xyz, l.zxy);
    vec3 x1 = x0 - i1 + C.xxx;  vec3 x2 = x0 - i2 + C.yyy;   vec3 x3 = x0 - D.yyy;
    i = mod289(i);  vec4 p = permute(permute(permute(i.z + vec4(0.0, i1.z, i2.z, 1.0)) + i.y + vec4(0.0, i1.y, i2.y, 1.0)) + i.x + vec4(0.0, i1.x, i2.x, 1.0));
    float n_ = 0.142857142857;   vec3  ns = n_ * D.wyz - D.xzx;  vec4 j = p - 49.0 * floor(p * ns.z * ns.z);  //  mod(p,7*7)
    vec4 x_ = floor(j * ns.z);  vec4 y_ = floor(j - 7.0 * x_);      vec4 x = x_ * ns.x + ns.yyyy;  vec4 y = y_ * ns.x + ns.yyyy;  vec4 h = 1.0 - abs(x) - abs(y);
    vec4 b0 = vec4(x.xy, y.xy);  vec4 b1 = vec4(x.zw, y.zw);  vec4 s0 = floor(b0) * 2.0 + 1.0;  vec4 s1 = floor(b1) * 2.0 + 1.0;  vec4 sh = -step(h, vec4(0.0));
    vec4 a0 = b0.xzyw + s0.xzyw * sh.xxyy;  vec4 a1 = b1.xzyw + s1.xzyw * sh.zzww;  vec3 p0 = vec3(a0.xy, h.x);  vec3 p1 = vec3(a0.zw, h.y);  vec3 p2 = vec3(a1.xy, h.z);  vec3 p3 = vec3(a1.zw, h.w);
    vec4 norm = taylorInvSqrt(vec4(dot(p0, p0), dot(p1, p1), dot(p2, p2), dot(p3, p3)));   p0 *= norm.x;  p1 *= norm.y;  p2 *= norm.z;  p3 *= norm.w;
    vec4 m = max(0.6 - vec4(dot(x0, x0), dot(x1, x1), dot(x2, x2), dot(x3, x3)), 0.0);  m = m * m;
    return 42.0 * dot(m * m, vec4(dot(p0, x0), dot(p1, x1), dot(p2, x2), dot(p3, x3)));
}

float Remap(float x, float oldMin, float oldMax, float newMin, float newMax) {
	return newMin + (newMax - newMin) * ((x - oldMin) / (oldMax - oldMin));
}

float Bilinear(float a00, float a10, float a11, float a01, float u, float v) {
	return (1 - u) * (1 - v) * a00 + (1 - u) * (v)*a01 + (u) * (1 - v) * a10 + (u) * (v)*a11;
}

// Heighfield elevation
// i, j: Integer coordinates
float at(int i, int j) {
	return hf[j * texSize.x + i];
}

// Heighfield color
// i, j: Integer coordinates
float atShading(int i, int j) {
	return shading[j * texSize.x + i];
}

void GetUV(vec2 p, out vec2 uv, out int i, out int j) {
	p = clamp(p, a, b);
	
	vec2 q = p - a;
	vec2 d = b - a;

	uv = q / d;
	uv = uv * vec2(texSize.x - 1, texSize.y - 1);
		
	i = int(uv.x);
	j = int(uv.y);

	uv = vec2(uv.x - i, uv.y - j);
}

mat2 matmul(mat2 A, mat2 B) {
	return mat2(A[0][0]*B[0][0]+A[0][1]*B[1][0], A[0][0]*B[0][1]+A[0][1]*B[1][1],
				A[1][0]*B[0][0]+A[1][1]*B[1][0], A[1][0]*B[0][1]+A[1][1]*B[1][1]);
}

// Heighfield elevation
// Bilinear interpolation of elevations
// p: Point
float Height(vec2 p) {
	vec2 uv;
	int i, j;
	GetUV(p, uv, i, j);

	return Bilinear(at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1), uv.x, uv.y);
}

// Heighfield color
// Bilinear interpolation of colors
// p: Point
float Shading(vec2 p) {
	vec2 uv;
	int i, j;
	GetUV(p, uv, i, j);

	return Bilinear(atShading(i, j), atShading(i + 1, j), atShading(i + 1, j + 1), atShading(i, j + 1), uv.x, uv.y);
}

// Normal using central difference
// p: point
// eps: epsilon for the computation
vec3 Normal(vec3 p, vec2 eps) {
    const vec3 e = vec3(eps.x, eps.y, 0.0);
    return normalize(vec3(
        -((Height(p.xy + e.xz) - Height(p.xy - e.xz)) / (2.0f * eps.x)),
        -((Height(p.xy + e.zy) - Height(p.xy - e.zy)) / (2.0f * eps.y)),
		1.0f
    ));
}

// Build the ray direction from the fragment position
vec3 BuildRd(int xa, int ya) {	
	vec3 view = normalize(CamLookAt - CamPos);
	vec3 horizontal = normalize(cross(view, CamUp));
	vec3 vertical = normalize(cross(horizontal, view));
	
	float length = 1.0f;
	float rad = CamAngleOfViewV;
	
	float vLength = tan(rad / 2.0f) * length;
	float hLength = vLength * iResolution.x / iResolution.y;
	
	vertical *= vLength;
	horizontal *= hLength;

	vec2 p = gl_FragCoord.xy;

	p.x = p.x - iResolution.x / 2.0f;
	p.y = p.y - iResolution.y / 2.0f;

	p.x += float(xa) / float(AA) - 0.5f;
	p.y += float(ya) / float(AA) - 0.5f;

	p.x /= iResolution.x / 2.0f;
	p.y /= iResolution.y / 2.0f;
	
	return normalize(view * length + horizontal * p.x + vertical * p.y);
}

// Intersection
// va, vb : Signed distance field values
float Intersection(float va, float vb) {
	return max(va, vb);
}

// Box signed distance field 
// p : Point
// va, vb : Vertices of the box
float Box(vec3 p, vec3 va, vec3 vb) {
	vec3 c = 0.5 * (va + vb);
	vec3 r = 0.5 * (vb - va);
	vec3 q = abs(p - c) - r;
	float d = length(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)),0.0);
	return d;
}

// Box signed distance field 
// p : Point
// va, vb : Vertexes of the box
float Box(vec2 p, vec2 va, vec2 vb) {
	vec2 c = 0.5 * (va + vb);
	vec2 r = 0.5 * (vb - va);
	vec2 q = abs(p - c) - r;
	float d = length(max(q, 0.0)) + min(max(q.x, q.y), 0.0);
	return d;
}

// Signed distance field object
// p : Point
float Map(vec3 p) {
	float t = p.z - Height(p.xy);
	float delta = 0.1f * (zRange.y - zRange.x);
	return Intersection(Box(p, vec3(a.x, a.y, zRange.x - delta), vec3(b.x, b.y, zRange.y + delta)), t);
}

// Ray-box intersection
// ro, rd : Ray
bool IntersectBox(vec3 ro, vec3 rd, out float tN, out float tF) {
	vec3 rinvDir = 1.0 / rd;
	float delta = 0.1 * (zRange.y - zRange.x);
	vec3 tbot = rinvDir * (vec3(a.x, a.y, zRange.x - delta) - ro);
	vec3 ttop = rinvDir * (vec3(b.x, b.y, zRange.y + delta) - ro);
	vec3 tmin = min(ttop, tbot);
	vec3 tmax = max(ttop, tbot);
	vec2 t = max(tmin.xx, tmin.yz);
	float t0 = max(t.x, t.y);
	t = min(tmax.xx, tmax.yz);
	float t1 = min(t.x, t.y);
	tN = t0;
	tF = t1;
	return t1 > max(t0, 0.0);
}

// Linear step function
// x value
// a, b interval
float LinearStep(float x, float a, float b) {
	if (x < a)
		return 0.0f;
	else if (x > b)
		return 1.0f;
	else
		return (x - a) / (b - a);
}

// Apply a green-brown-yellow color palette
// t interpolation parameter in [0, 1]
vec3 GetColorElevation(float t) {
	//vec3 c[4] = { vec3(0.627, 0.863, 0.411), vec3(1.00, 0.90, 0.45),vec3(0.659, 0.607, 0.541), vec3(0.95, 0.95, 0.95) };
	vec3 c[4] = { vec3(0.627, 0.763, 0.411), vec3(1.00, 0.90, 0.45),vec3(0.659, 0.607, 0.541), vec3(0.95, 0.95, 0.95) };
	//float a[4] = { 0.0, 150.0, 250.0, 450.0 };
	float a[4] = { 0.0, 100.0, 150.0, 200.0 };
	if (t < a[0])
		return c[0];
	if (t > a[3])
		return c[3];
	for (int i = 0; i < 3; i++)
	{
		if (t < a[i + 1])
		{
			float s = LinearStep(t, a[i], a[i + 1]);
			return mix(c[i], c[i + 1], s);
		}
	}
	return vec3(0);
}

// Apply a white to brown color palette, computed from elevation, slope and fractional laplacian.
// p point
// n normal at p
vec3 ShadeElevation(vec3 p, vec3 n) {
	// base color: elevation
	// Note(Axel): played with alpha depending on the figure
	float alpha = 2.0;
	float t = alpha * LinearStep(p.z, zRange.x, zRange.y);

	// Note(Axel): played with aa/bb depending on the figure
	float aa = 0.0;
	float bb = 1.0;
	t = clamp(t, aa, bb);
	
	// White to brown
	vec3 col = mix(vec3(1.0), vec3(125.0/255.0, 70.0/255.0, 45.0/255.0), t); 

	// Modulate with slope
	t = 1.0f - n.z;
	//t = pow(t, 2.0); // Note(Axel): you can play with this
	col = mix(col, vec3(0.0), t);

	// Modulate with fractional laplacian, stored in a buffer
	t = Shading(p.xy);
	col = mix(col, vec3(1.0), t);

	// Uncomment for the figure
	// For endoreic lake/sea shading
	if (p.z <= 25.0f) // Note(Axel): you can play with this
	{
		// modify color
		col = vec3(107.0 / 255.0, 163.0 / 255.0, 219 / 255.0);

		// mix with shores
		col = mix(col, vec3(1.0), smoothstep(-0.1f, 22.0f, p.z));

		// Modify normal with noise to 'simulate' a water surface
		float no = 2.5f * snoise(p * 0.00004f) * 0.5f + 0.5f;
		no = clamp(no, 0.0, no);
		n.x += no; // no normalization, we don't care
	}

	// Simple lighting
	float s = dot(n, -normalize(lightDir));
	s = 0.5 * (1.0 + s);
	s *= s;
	col = col * s;

	return col;
}

float cosine(float x)
{
	return pow(0.5+0.5*x,4.0);
}
vec3 ShadeGreenBrownYellow(vec3 p, vec3 n) {
	// Color from colormap
	
	vec3 col = GetColorElevation(p.z);
	
	// Simple lighting
	float s = dot(n, -normalize(lightDir));
	s = 0.5 * (1.0 + s);
	s *= s;
	col = col * s;

	return col;
	
	float x=cosine(dot(n.xy,vec2(1.0,0.0)));
	float y=cosine(dot(n.xy,vec2(-1.0,0.0)));
	float z=cosine(dot(n.xyz,vec3(0.0,0.0,1.0)));
	vec3 c = 0.6*vec3(1.0) + 0.3*x*vec3(1.0,0.8,0.1)+0.3*y*vec3(0.5,0.3,0.9)+0.1*z*vec3(1.0,1.0,1.0);

    float d1 = pow(0.5 * (1.0 + dot(n.xy, normalize(lightDir.xy) )),2.0);
    //float d2 = pow(0.5 * (1.0 + dot(n.xy, normalize(lightDir.xy) )),2.0);
	return c*(0.8+0.2*d1);

}

    // Not even Phong shading, use weighted cosine instead for smooth transitions
vec3 Shade(vec2 p, vec3 n) {
    float d1 = pow(0.5 * (1.0 + dot(n.xy, normalize(lightDir.xy) )),2.0);
	vec2 q = p - a;
	vec2 d = b - a;
	vec2 uv = q / d;
	vec3 col = texture(albedo, uv).rgb;

	// Sides
	if (abs(Box(p,a,b))<epsilon) { col=vec3(0.3,0.29,0.31); }
    return 0.925*col+0.075*d1;   // Gundy: 12/11/2023 : need a parameter
    //return 0.79*col + 0.21 *d1; // Original
}

vec4 ShadeWire(vec3 p, vec4 c) {
	vec2 q = p.xy - a;
	vec2 d = b - a;
	vec2 uv = q / d;
	uv *= texSize;
	uv = uv - vec2(int(uv.x), int(uv.y));
	if (uv.x < 0.1f || uv.x > 0.9f || uv.y < 0.1f || uv.y > 0.9f)
		c = vec4(0, 0, 0, 1);
	return c;
}

vec4 ShadeCost(int s) {
	float t = float(s) / (float(STEPS - 1)); 
	return 0.5 + mix(vec4(0.05, 0.05, 0.5, 1.0), vec4(0.65, 0.39, 0.65, 1.0), t);
}

bool Raymarch(vec3 ro, vec3 rd, out vec3 p, out float t, out int s, int maxStep, float eps) {
    t = 0.0f;
	s = 0;

	// Check if ray intersects the bounding box of the heightfield
	float ta, tb;
	if (!IntersectBox(ro, rd, ta, tb))
	{
		t = cameraPlanes.y;
		return false;
	}
	
	// Raymarch
	t = max(ta, 0.0f);
	float d = 0.0f;
	float uz = abs(rd.z);
	float kr = uz + K * sqrt(1.0f - (uz * uz));
    for(s = 0; s < maxStep; s++) 
	{
		if (t > tb)
		{
			t = cameraPlanes.y;
		    return false;
		}
        p = ro + rd * t;
        d = Map(p);
        if(d <= 0.0)
			return true;
        t += max(d / kr, eps);
    }
	t = cameraPlanes.y;
    return false;
}

// Compute a perpendicular direction to a vector
// u : direction
vec3 perp(vec3 u) {
     vec3 a = abs(u);
     vec3 v;
     if (a.x <= a.y && a.x <= a.z)
         v = vec3(0, -u.z, u.y);
     else if (a.y <= a.x && a.y <= a.z)
         v = vec3(-u.z, 0.0, u.x);
     else
        v = vec3(-u.y, u.x, 0.0);
    return v;
}

// Compute Soft shadows
// o : Ray origin
// a : ray direction offset strength
float SmoothShadow(vec3 o, float a) {
	vec3 rd = -normalize(lightDir);
	vec3 p;
	float t;
	int s;
	vec3 u = normalize(perp(rd));
	vec3 v = normalize(cross(rd,u));
	int n = 0;
	for (int i = 0; i < smoothShadowNbStep; i++)
	{
		vec3 dir = normalize(rd + a * (BlueNoiseInDisk[i % 64].x * u + BlueNoiseInDisk[i % 64].y * v));
		bool hit = Raymarch(o, dir, p, t, s, smoothShadowMarchingSteps, smoothShadowMarchingEps);
		if (hit)
			n++;
	}
	return float(n) / float(smoothShadowNbStep);
}

// Compute sky color 
// d : Ray direction
vec3 SkyShadeBlue(in vec3 d) {
  	// light direction
	vec3 lig = normalize(vec3(0.35,0.55, 0.65));
	float sun = 0.5*(dot(lig,d)+1.0);
	vec3 color = vec3(0.65,0.85,1.0) * (0.75-0.25*d.z);
	color += vec3(0.95,0.90,0.95)*pow( sun, 18.0 );
	return color;
}

vec4 Pixel(int xa, int ya) {
	// Compute ray
	vec3 ro = CamPos;
	vec3 rd = BuildRd(xa, ya);

	// Compute Intersection with terrain
	vec4 c;
	vec3 p;
	float t;
	int s;
	bool hit = Raymarch(ro, rd, p, t, s, STEPS, epsilon);

	// Write to z-buffer (which is not linear, see https://learnopengl.com/Advanced-OpenGL/Depth-testing)
	// This allows us to render classical meshes with the raytracing widget
	if (hit)
		gl_FragDepth = ((1.0f / t) - (1.0f / cameraPlanes.x)) / ((1.0f / cameraPlanes.y) - (1.0f / cameraPlanes.x));
	else
		gl_FragDepth = 1.0f;

	if (hit)
	{
		vec3 n = Normal(p, (b - a) / vec2(texSize - ivec2(1, 1)));

		// Albedo + simple diffuse lighting
		if (useAlbedo == 1)
			c = vec4(Shade(p.xy, n), 1.0);
		// Elevation shading + diffuse lighting
		else if (useElevation == 1)
			c = vec4(ShadeElevation(p, n), 1.0);
		else if (useGreenBrownYellow == 1)
			c = vec4(ShadeGreenBrownYellow(p, n), 1.0);
		// Normal
		else
			c = vec4(0.2 * (vec3(3.0) + 2.0 * n.xyz), 1.0);	

		// Add wireframe
		if (useWireframe == 1 && abs(Box(p.xy, a, b)) > epsilon)
			c = ShadeWire(p, c);

		// Heightfield sides
		if (abs(Box(p.xy, a, b)) < epsilon || abs(p.z - zRange.x + 0.1f * (zRange.y - zRange.x)) < epsilon)
			c = vec4(0.3, 0.29, 0.31, 1.0);

		// Self shadows
		vec3 oq = p + vec3(0.0, 0.0, 3.0 * epsilon);
		vec3 rq = normalize(-lightDir);
		if (useSelfShadow == 1)
		{
			vec3 q;
			float qt;
			int sq;
			bool hit2 = Raymarch(oq, rq, q, qt, sq, selfShadowMarchingSteps, selfShadowMarchingEps);
			if (hit2)
				c.xyz *= selfShadowStrength; // Darken
		}

		// Smooth shadows
		if (useSmoothShadow == 1)
		{
			float ss = SmoothShadow(oq, 0.1);
			ss = 1.0 - smoothShadowStrength * ss;
			c.xyz *= ss;
		}
	}
	else
	{
		c = vec4(SkyShadeBlue(rd), 1.0);
		c = mix(c, vec4(0.85, 0.95, 1.0, 1.0), gl_FragCoord.x * 0.0001 + 0.95 * gl_FragCoord.y / iResolution.y);
	}

	// Shade raymarching steps
	if (useCost == 1)
		c = ShadeCost(s);
	return c;
}

void main() {
	if (AA == 1)
		Fragment = Pixel(0, 0);
	else
	{
		vec4 tot = vec4(0.0f);
		for (int m = 0; m < AA; m++)
		{
			for (int n = 0; n < AA; n++)
				tot += Pixel(m, n);
		}
		tot = tot / float(AA * AA);
		tot.a = 1.0;
		Fragment = tot;
	}
}

#endif
