#version 430 core

#ifdef VERTEX_SHADER
layout(location = 0) in vec3 vertex;
layout(location = 1) in float color;


uniform mat4 ModelViewMatrix;
uniform mat4 ProjectionMatrix;
uniform mat4 TRSMatrix;

uniform int Width;
uniform int Height;


out float outColor;
out vec2 coord;
out float radius;



void main(void)
{
	mat4 MVP      = ProjectionMatrix * ModelViewMatrix;
	gl_Position   = MVP * TRSMatrix * (vec4(vertex, 1.0)); 
	vec4 tmp = ModelViewMatrix * TRSMatrix * (vec4(vertex, 1.0));
	gl_PointSize = -100/tmp.z; 
	outColor = color;
	coord = ((gl_Position.xy/gl_Position.w)+1)*vec2(Width, Height)/2;
	radius = gl_PointSize/2;
} 
#endif

#ifdef FRAGMENT_SHADER


in float outColor;
in vec2 coord;
in float radius;

out vec4 fragment;

vec4 colorRamp(float color)
{
	float white = (0.5-abs(color-0.5))*2;
	if(color < 0.5)
		return vec4(white, white, 1, 0.5); 
	return vec4(1, white, white, 0.5); 
}


void main()
{

	vec2 toCenter = gl_FragCoord.xy - coord;



	if(length(toCenter) > radius)
		discard;

	fragment = colorRamp(outColor);
}

#endif
