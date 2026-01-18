#version 430

in vec4 Position;

layout(location = 0) uniform mat4 ViewMatrix;
layout(location = 1) uniform mat4 ProjMatrix;

void main()
{
    gl_Position = ProjMatrix * ViewMatrix * Position;
}