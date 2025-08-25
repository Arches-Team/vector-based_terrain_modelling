#version 430 core

#ifdef VERTEX_SHADER
layout(location = 0) in vec3 vertex;

uniform mat4 ModelViewMatrix;
uniform mat4 ProjectionMatrix;
uniform mat4 TRSMatrix;

void main(void)
{
	mat4 MVP      = ProjectionMatrix * ModelViewMatrix;
	gl_Position   = MVP * TRSMatrix * (vec4(vertex, 1.0)); 
} 
#endif

#ifdef FRAGMENT_SHADER
out vec4 fragment;

void main()
{
	fragment = vec4(0.0, 0.0, 0.0, 1.0);
}

#endif
