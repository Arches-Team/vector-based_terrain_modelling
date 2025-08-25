#version 430 core

layout(std140) uniform transform {		
	uniform	mat4 ModelViewMatrix;		
	uniform	mat4 ProjectionMatrix;		
} Transform;
uniform vec3 lightPos_vert;	
uniform vec3 cameraPos_vert;				

in  vec3 position_vert ; // Positions
in  vec3 color_vert;     // Colors : use per instance color or per vertex
in  vec3 normal_vert ;   // Normals
in  vec2 UV_vert ;  	 // UV
in  vec3 inst_pos;       // Attrib pour les instances - Translate
in  vec3 inst_scl;       // Attrib pour les instances - Scale
in  mat3 inst_rot;       // Attrib pour les instances - Rotate

out vec3 position_geom ; 
out vec3 color_geom;     
out vec3 normal_geom;   

out vec3 eyeVec_geom;    
out vec3 lightDir_geom;  
 
out vec2 UV_geom;  
  

void main()
{
	// Matrice de transformation
	mat4 MVP      = Transform.ProjectionMatrix * Transform.ModelViewMatrix ;
	// Position des vertex
	position_geom = position_vert*inst_scl*inst_rot+inst_pos;
	gl_Position   = MVP*(vec4(position_geom,1.0));                // La position n'est pas la meme en fonction de l'instance
	vec4 vVertex  = Transform.ModelViewMatrix * vec4(position_geom,1.0);
	
	color_geom    = color_vert;
	normal_geom   = normalize(normal_vert*inst_rot);
	
	eyeVec_geom   = normalize(vVertex.xyz-cameraPos_vert);
	lightDir_geom = normalize(lightPos_vert - position_geom);
	
	UV_geom		  = UV_vert;
}