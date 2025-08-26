#version 430

out vec4 color;

layout (location = 2) uniform vec4 FragColor;

void main()
{
	color = FragColor;
}