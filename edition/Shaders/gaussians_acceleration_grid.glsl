#version 430 core

#ifdef COMPUTE_SHADER

layout(binding = 0, std430) coherent buffer GaussiansBuffer
{
	float gaussians[];
};

layout(std430, binding = 1) coherent buffer GridCellBuffer {
    int gridCellCounts[];
};

layout(std430, binding = 2) coherent buffer GridMappingBuffer {
    int gridCellMappings[];
};

uniform ivec2 gridResolution;
uniform int maxPerCell;
uniform int gaussianOffset;

struct bound2
{
    vec2 mMin;
    vec2 mMax;
};

// bounding box for a ellipse (https://iquilezles.org/articles/ellipses)
bound2 ellipseAABB( in vec2 c, in vec2 u, in vec2 v )  // disk: center, 1st axis, 2nd axis
{
    vec2 e = sqrt( u*u + v*v );
    return bound2( c-e, c+e );
}

vec2 ellipsePoint(float t, float sigma_x, float sigma_y, float theta)
{
    return vec2(
        sigma_x*cos(t)*cos(theta) - sigma_y*sin(t)*sin(theta),
        sigma_x*cos(t)*sin(theta) + sigma_y*sin(t)*cos(theta)
    );
}

layout(local_size_x = 64) in;
void main() {
    uint i = gl_GlobalInvocationID.x;

    int id = int(gaussians[i*gaussianOffset+0]);

    float theta = gaussians[i*gaussianOffset+3];
    
    float sigma_x = gaussians[i * gaussianOffset + 1];
    float sigma_y = gaussians[i * gaussianOffset + 2];

    vec2 c = vec2((gaussians[i*gaussianOffset+6]+1.)/2., (gaussians[i*gaussianOffset+5]+1.)/2.);
    // c = -c;

    // float p = smoothstep(0.2, 0.7, (int(i)/17077.));
    // float a = 1.;
    // if((c.y < p &&(c.x + c.y) < a) || (c.x < p && (c.x + c.y) < a) || (c.x + c.y) < 0.3 )
    //     return;

    mat2 rot = mat2(cos(theta),sin(theta),-sin(theta),cos(theta));

    vec2 urot;
    vec2 vrot;

    float sigma_influence = 3.;
    sigma_x *= sigma_influence;
    sigma_y *= sigma_influence;
    urot = vec2(0., sigma_y)*rot;
    vrot = vec2(sigma_x, 0.)*rot;

    bound2 bb = ellipseAABB(c, urot/2., vrot/2.);

    // float t0 = atan(-(sigma_y*tan(theta))/sigma_x);
    // float t1 = atan((sigma_y * cos(theta))/(sigma_x*sin(theta)));
    // vec2 p0 = ellipsePoint(t0, sigma_x, sigma_y, theta);
    // vec2 p1 = ellipsePoint(t1, sigma_x, sigma_y, theta);
    // bound2 bb;
    // bb.mMin.x = c.x - p0.x;
    // bb.mMax.x = c.x + p0.x;
    // bb.mMin.y = c.y - p1.y;
    // bb.mMax.y = c.y + p1.y;
    
    // bb.mMin = (bb.mMin+1.)/2.;
    // bb.mMax = (bb.mMax+1.)/2.;

    ivec2 minCell = ivec2(floor(bb.mMin * gridResolution));
    ivec2 maxCell = ivec2(ceil(bb.mMax * gridResolution));

    minCell = clamp(minCell, ivec2(0), gridResolution);
    maxCell = clamp(maxCell, ivec2(0), gridResolution);

    for (int y = minCell.y; y < maxCell.y; ++y) {
        for (int x = minCell.x; x < maxCell.x; ++x) {
            int cellIndex = y * gridResolution.x + x;
			int count = atomicAdd(gridCellCounts[cellIndex], 1);

            count = min(count, maxPerCell);
            gridCellMappings[cellIndex * maxPerCell + count] = int(i);
        }
    }
}

#endif
