#version 450 core

#ifdef COMPUTE_SHADER

// In
layout(binding = 0, std430) readonly buffer InMap { float in_map[]; };

// Out
layout(binding = 1, std430) writeonly buffer OutMap { float out_map[]; };

uniform int buffer_size_x;
uniform int buffer_size_y;
uniform float threshold;

const ivec2 next8[8] = ivec2[8](ivec2(0, 1),  ivec2(1, 1),   ivec2(1, 0),  ivec2(1, -1),
                                ivec2(0, -1), ivec2(-1, -1), ivec2(-1, 0), ivec2(-1, 1));

int ToIndex1D(ivec2 p) { return p.x + buffer_size_x * p.y; }


layout(local_size_x = 8, local_size_y = 8, local_size_z = 1) in;
void main() {
    int x = int(gl_GlobalInvocationID.x);
    int y = int(gl_GlobalInvocationID.y);
    if (x < 0)   return;
    if (y < 0)   return;
    if (x >= buffer_size_x) return;
    if (y >= buffer_size_y) return;

    ivec2 pos = ivec2(x, y);
    int id = ToIndex1D(pos);
   

    // mean filter
    float mean = in_map[id];
    for (int i = 0; i < 8; i++) {
        ivec2 nei_pos = pos + next8[i];
        if (nei_pos.x < 0 || nei_pos.x >= buffer_size_x || nei_pos.y < 0 || nei_pos.y >= buffer_size_y) continue;
        mean += in_map[ToIndex1D(nei_pos)];
    }

    // apply max only if superior to threshold
    out_map[id] = mean / 9.;
}

#endif
