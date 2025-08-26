#version 430 core

uniform vec2 WIN_SCALE;

layout(triangles) in;
layout(triangle_strip, max_vertices=3) out;

in vec3 position_geom[]; 
in vec3 color_geom[];     
in vec3 normal_geom[];   
in vec3 eyeVec_geom[];    
in vec3 lightDir_geom[];  
in vec2 UV_geom[];  
  
out vec3 position_frag ;  // Point
out vec3 color_frag;      // Color
out vec3 normal_frag;     // Normal
out vec3 eyeVec_frag;     // Eve position
out vec3 lightDir_frag;   // Light direction
out vec2 UV_frag;  

out vec3 dist;
out float ratio;

void main()
{	
  vec2 p0 = WIN_SCALE * gl_in[0].gl_Position.xy / gl_in[0].gl_Position.w;
  vec2 p1 = WIN_SCALE * gl_in[1].gl_Position.xy / gl_in[1].gl_Position.w;
  vec2 p2 = WIN_SCALE * gl_in[2].gl_Position.xy / gl_in[2].gl_Position.w;

  vec2 v0 = p2-p1;
  vec2 v1 = p2-p0;
  vec2 v2 = p1-p0;

  float area = abs(v1.x*v2.y - v1.y * v2.x);

  // Triangle aspect ratio
  float ab=length(position_geom[0]-position_geom[1]);
  float bc=length(position_geom[1]-position_geom[2]);
  float ca=length(position_geom[2]-position_geom[0]);

  float s = 0.5*(ab+bc+ca);
  float u = (s - ab)*(s - bc)*(s - ca);

  ratio=8.0*u / (ab*bc*ca);

  dist = vec3(area/length(v0),0,0);
  gl_Position = gl_in[0].gl_Position;
  position_frag = position_geom[0]; color_frag =  color_geom[0];  normal_frag = normal_geom[0];   eyeVec_frag = eyeVec_geom[0];  lightDir_frag = lightDir_geom[0];  UV_frag = UV_geom[0]; 
  EmitVertex();

  dist = vec3(0,area/length(v1),0);
  gl_Position = gl_in[1].gl_Position;
  position_frag = position_geom[1]; color_frag =  color_geom[1];  normal_frag = normal_geom[1];   eyeVec_frag = eyeVec_geom[1];  lightDir_frag = lightDir_geom[1];  UV_frag = UV_geom[1]; 
  EmitVertex();

  dist = vec3(0,0,area/length(v2));
  gl_Position = gl_in[2].gl_Position;
  position_frag = position_geom[2]; color_frag =  color_geom[2];  normal_frag = normal_geom[2];   eyeVec_frag = eyeVec_geom[2];  lightDir_frag = lightDir_geom[2];  UV_frag = UV_geom[2]; 
  EmitVertex();

  EndPrimitive();
  
}  
