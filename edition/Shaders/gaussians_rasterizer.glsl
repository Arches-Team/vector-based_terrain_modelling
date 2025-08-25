#version 430 core

#define M_PI 3.1415926535897932384626433832795

#ifdef COMPUTE_SHADER

layout (binding = 0, std430) coherent buffer HeightfieldBuffer
{
    float hf[];
};

layout (binding = 1, std430) coherent buffer GaussiansBuffer
{
    float gaussians[];
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
uniform int gaussianOffset;
uniform int showNbGaussians;
uniform int gaussianID;
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

    // if(gaussianID == 1)
    //     return 10.;
    p.xy = p.yx;

    ivec2 cell = ivec2(floor(((p + 1.) / 2.) * gridResolution));
    int cellIndex = cell.y * gridResolution.x + cell.x;
    int count = gridCellCounts[cellIndex];
    // return count;
    // return cellIndex;
    // return gridCellMappings[cellIndex * maxPerCell + count/2];

    // [-5 ; 5]
    // p = p * 5.;

    // Debug to check range
    // return (length(p)/sqrt(2))*zRange.y/5.;

    // return showNbGaussians;

    //bool has_beta = gaussianOffset >= 8;//mod(gaussiansSize, 8) == 0;

    //int cpt = 0;
    for (int j = 0; j < min(count, maxPerCell); ++j)
    {
        int i = gridCellMappings[cellIndex * maxPerCell + j];
        if (i >= showNbGaussians) continue;

        int id = int(gaussians[i * gaussianOffset + 0]);
        float sigma_x = gaussians[i * gaussianOffset + 1];
        float sigma_y = gaussians[i * gaussianOffset + 2];
        float theta = gaussians[i * gaussianOffset + 3];
        float amplitude = gaussians[i * gaussianOffset + 4];
        vec2 mu = vec2(gaussians[i * gaussianOffset + 5], gaussians[i * gaussianOffset + 6]);

        float x = p.x;
        float y = p.y;

        // x -= mu.y * 5.;
        // y -= mu.x * 5.;
        x = mu.y - x;
        y = mu.x - y;
        
        float dist;

        // theta = 3*M_PI/2 - theta;
        // sigma_x *= 2;
        // sigma_y *= 2;

        mat2 R = mat2(cos(theta), -sin(theta),
                    sin(theta), cos(theta));

        // If the gaussian does not influence the current point, do not compute
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
        //float z = inv_cov[0] * x * x + 2 * (inv_cov[1] * x * y) + inv_cov[2] * y * y;

        // TODO: gaussianId is not working if set as uniform. I suspect an overflow of uniforms even if I should be able to have thousands of them...
        if(id == gaussianID)
        {
            float beta = gaussians[i * gaussianOffset + 7];
            float val = exp(-pow(z, beta)) * amplitude;
            ret += val;
        }
        else if(id == 1)
        {
            float val = exp(-z) * amplitude;

            // [-1, 1]
            vec2 pNorm = vec2(-x, -y);
            //theta = abs(theta);
            theta = (M_PI/2.) - theta;
            mat2 minusR = mat2(cos(-theta), -sin(-theta),
                               sin(-theta), cos(-theta));
            
            float origSigmaX = gaussians[i * gaussianOffset + 9];
            float origSigmaY = gaussians[i * gaussianOffset + 10];
            float origTheta = gaussians[i * gaussianOffset + 11];
            float thetaNoDeform = gaussians[i * gaussianOffset + 12];

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

            vec2 t = vec2(gaussians[i * gaussianOffset + 7], gaussians[i * gaussianOffset + 8]);
            //[ -1, 1]
            t = (t * 2.) - 1.;

            // pNorm.x = mu.y - pNorm.x;
            // pNorm.y = mu.x - pNorm.y;
            // pNorm = (mu - t) + pNorm;
            pNorm = noDeformR * pNorm;
            // pNorm.x *= sigma_y/sigma_x;
            // pNorm.x /= sigma_x;
            // pNorm.x *= origSigmaX;
            // pNorm.y *= origSigmaY;
            // pNorm.x = (pNorm.x/sigma_x)*origSigmaX;
            // pNorm.y = (pNorm.y/sigma_y)*origSigmaY;
            pNorm = inverseS * pNorm;
            pNorm = origS * pNorm;
            pNorm = origR * pNorm;

            t *= -1;

            //t = abs(t);
            
            t = pNorm - t;
            t = clamp(t, -1., 1.);
            t = ((t+1.)/2.)*ivec2(detailsSize, detailsSize);
            ivec2 tInt = ivec2(t);

            int m = tInt.x;
            int n = tInt.y;

            ret += Bilinear(atDetails(m, n), atDetails(m + 1, n), atDetails(m + 1, n + 1), atDetails(m, n + 1), t.x - m, t.y - n) * val * noiseLevel;
            // ret += details[tInt.x * detailsSize.x + tInt.y] * val * noiseLevel;
        }

        //cpt++;
        
    }
    // return cpt*1.;
    // ret = ((clamp(ret, -1., 10.)+1.)/2.)*100.;
    // ret += perlinNoise(((p / 5.) + 1.) / 2., 100, 6, 0.6, 2.0, 0x578437adU) * noiseLevel / 10.;
    ret *= zRange.y;
    // ret += 50.;
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
