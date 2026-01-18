#version 430

layout(location = 0) in vec3 Position;
layout(location = 1) in vec4 VertexColor;

layout(location = 2) uniform mat4 ViewMatrix;
layout(location = 3) uniform mat4 ProjMatrix;

out vec4 FragColor;

void main()
{
    gl_Position = ProjMatrix * ViewMatrix * vec4(Position, 1.0);
    FragColor = VertexColor;
}