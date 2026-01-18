#version 430 core

#ifdef COMPUTE_SHADER

layout(binding = 0, std430) coherent buffer PrimitivesBuffer
{
	float primitives[];
};

layout(std430, binding = 1) coherent buffer GridCellBuffer {
    int gridCellCounts[];
};

layout(std430, binding = 2) coherent buffer GridMappingBuffer {
    int gridCellMappings[];
};

uniform ivec2 gridResolution;
uniform int maxPerCell;
uniform int primitivesOffset;

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

    int id = int(primitives[i*primitivesOffset+0]);

    float theta = primitives[i*primitivesOffset+3];
    
    float sigma_x = primitives[i * primitivesOffset + 1];
    float sigma_y = primitives[i * primitivesOffset + 2];

    vec2 c = vec2((primitives[i*primitivesOffset+6]+1.)/2., (primitives[i*primitivesOffset+5]+1.)/2.);

    mat2 rot = mat2(cos(theta),sin(theta),-sin(theta),cos(theta));

    vec2 urot;
    vec2 vrot;

    float sigma_influence = 3.;
    sigma_x *= sigma_influence;
    sigma_y *= sigma_influence;
    urot = vec2(0., sigma_y)*rot;
    vrot = vec2(sigma_x, 0.)*rot;

    bound2 bb = ellipseAABB(c, urot/2., vrot/2.);

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
