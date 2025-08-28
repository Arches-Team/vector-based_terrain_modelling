#version 430 core

#define M_PI 3.1415926535897932384626433832795

#ifdef COMPUTE_SHADER

layout (binding = 0, std430) coherent buffer HeightfieldBuffer
{
    float hf[];
};

layout (binding = 1, std430) coherent buffer PrimitivesBuffer
{
    float primitives[];
};

layout (std430, binding = 2) buffer GridCellBuffer {
    int gridCellCounts[];
};

layout (std430, binding = 3) buffer GridMappingBuffer {
    int gridCellMappings[];
};

layout (std430, binding = 4) buffer DetailsBuffer {
    float details[];
};


uniform ivec2 gridResolution;
uniform int maxPerCell;
uniform ivec2 nxy;
uniform vec2 zRange;
uniform float noiseLevel;
uniform int primitiveOffset;
uniform int showNbPrimitives;
uniform int gaussianID;
uniform int detailsID;
uniform int detailsSize;
int nx = nxy.x;
int ny = nxy.y;


float atDetails(int i, int j) {
	return details[i * detailsSize + j];
}

float Bilinear(float a00, float a10, float a11, float a01, float u, float v) {
	return (1 - u) * (1 - v) * a00 + (1 - u) * (v)*a01 + (u) * (1 - v) * a10 + (u) * (v)*a11;
}

int ToIndex1D(int i, int j) {
    return i + nx * j;
}

mat2 matmul(mat2 A, mat2 B) {
    return mat2(A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1],
                A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]);
}

// p must be in [-1, 1]
float Height(vec2 p)
{
    float ret = 0;

    p.xy = p.yx;

    ivec2 cell = ivec2(floor(((p + 1.) / 2.) * gridResolution));
    int cellIndex = cell.y * gridResolution.x + cell.x;
    int count = gridCellCounts[cellIndex];

    for (int j = 0; j < min(count, maxPerCell); ++j)
    {
        int i = gridCellMappings[cellIndex * maxPerCell + j];
        if (i >= showNbPrimitives) continue;

        int id = int(primitives[i * primitiveOffset + 0]);
        float sigma_x = primitives[i * primitiveOffset + 1];
        float sigma_y = primitives[i * primitiveOffset + 2];
        float theta = primitives[i * primitiveOffset + 3];
        float amplitude = primitives[i * primitiveOffset + 4];
        vec2 mu = vec2(primitives[i * primitiveOffset + 5], primitives[i * primitiveOffset + 6]);

        float x = p.x;
        float y = p.y;

        x = mu.y - x;
        y = mu.x - y;
        
        float dist;

        mat2 R = mat2(cos(theta), -sin(theta),
                    sin(theta), cos(theta));

        // If the pritimives does not influence the current point, do not compute
        {
            vec2 tmpP = R * vec2(x, y);

            float sigma_influence = 6.;
            float axisX = sigma_x * sigma_influence;
            float axisY = sigma_y * sigma_influence;
            tmpP.x *= axisY / axisX;
            dist = length(tmpP) - axisY;
            if (dist > 0) continue;
        }

        mat2 S = mat2(sigma_x, 0.,
                    0., sigma_y);

        mat2 M = matmul(S, R);
        mat2 cov = matmul(transpose(M), M);
        vec3 cov_compact = vec3(cov[0][0], cov[0][1], cov[1][1]);

        float det = cov_compact[0] * cov_compact[2] - cov_compact[1] * cov_compact[1];
        float inv_det = 1. / det;

        vec3 inv_cov = vec3(cov_compact[2] * inv_det, cov_compact[1] * inv_det, cov_compact[0] * inv_det);
        float z = 0.5 * (inv_cov[0] * x * x + inv_cov[2] * y * y) + inv_cov[1] * x * y;

        if(id == gaussianID)
        {
            float beta = primitives[i * primitiveOffset + 7];
            float val = exp(-pow(z, beta)) * amplitude;
            ret += val;
        }
        else if(id == detailsID)
        {
            float val = exp(-z) * amplitude;

            // [-1, 1]
            vec2 pNorm = vec2(-x, -y);
            //theta = abs(theta);
            theta = (M_PI/2.) - theta;
            mat2 minusR = mat2(cos(-theta), -sin(-theta),
                               sin(-theta), cos(-theta));
            
            float origSigmaX = primitives[i * primitiveOffset + 9];
            float origSigmaY = primitives[i * primitiveOffset + 10];
            float origTheta = primitives[i * primitiveOffset + 11];
            float thetaNoDeform = primitives[i * primitiveOffset + 12];

            mat2 noDeformR = mat2(cos(thetaNoDeform), -sin(thetaNoDeform),
                               sin(thetaNoDeform), cos(thetaNoDeform));


            mat2 origR = mat2(cos(origTheta), -sin(origTheta),
                              sin(origTheta), cos(origTheta));
            mat2 origS = mat2(origSigmaY, 0.,
                              0., origSigmaX);
            float sx, sy;
            sx = min(sigma_y, sigma_x);
            sy = max(sigma_y, sigma_x);
            if(origSigmaY < sigma_y || origSigmaX < sigma_x)
            {
                sx = max(sigma_y, sigma_x);
                sy = min(sigma_y, sigma_x);
                
            }
                
            mat2 inverseS = mat2(1./sx, 0.,
                                 0., 1./sy);

            vec2 t = vec2(primitives[i * primitiveOffset + 7], primitives[i * primitiveOffset + 8]);
            //[ -1, 1]
            t = (t * 2.) - 1.;

            pNorm = noDeformR * pNorm;
            pNorm = inverseS * pNorm;
            pNorm = origS * pNorm;
            pNorm = origR * pNorm;

            t *= -1;
            
            t = pNorm - t;
            t = clamp(t, -1., 1.);
            t = ((t+1.)/2.)*ivec2(detailsSize, detailsSize);
            ivec2 tInt = ivec2(t);

            int m = tInt.x;
            int n = tInt.y;

            ret += Bilinear(atDetails(m, n), atDetails(m + 1, n), atDetails(m + 1, n + 1), atDetails(m, n + 1), t.x - m, t.y - n) * val * noiseLevel;
        }
    }

    ret *= zRange.y;
    ret = max(0, ret);

    return ret;
}

layout (local_size_x = 8, local_size_y = 8, local_size_z = 1) in;
void main() {
    int i = int(gl_GlobalInvocationID.x);
    int j = int(gl_GlobalInvocationID.y);
    if (i < 0) return;
    if (j < 0) return;
    if (i >= nx) return;
    if (j >= ny) return;

    hf[ToIndex1D(i, j)] = Height(vec2((float(i) / nx) * 2. - 1., (float(j) / ny) * 2. - 1.));
}

#endif
