#version 430 core

#ifdef COMPUTE_SHADER

layout (binding = 0, std430) coherent buffer InitialHeightfieldBuffer {
    float initialHf[];
};

layout (binding = 1, std430) coherent buffer UpdatedHeightfieldBuffer {
    float updatedHf[];
};

layout (binding = 2, std430) coherent buffer DifferenceBuffer {
    float difference[];
};

uniform int nx;
uniform int ny;

int ToIndex1D(int i, int j) {
    return i + nx * j;
}

layout (local_size_x = 8, local_size_y = 8, local_size_z = 1) in;
void main() {
    int i = int(gl_GlobalInvocationID.x);
    int j = int(gl_GlobalInvocationID.y);
    if (i < 0 || j < 0 || i >= nx || j >= ny) return;

    int index = ToIndex1D(i, j);
    float initialVal = initialHf[index];
    float updatedVal = updatedHf[index];
    float diff = updatedVal - initialVal;
    float eps = 1e-6;
    float threshold = 10.;

    difference[index] = 0.;
    if (abs(diff) > eps && abs(updatedVal) < threshold)
    {
        difference[index] = 1.;
    }
}

#endif
