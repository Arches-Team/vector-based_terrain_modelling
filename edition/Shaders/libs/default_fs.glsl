#version 430 core

layout(std140) uniform material {
	uniform	vec4 ambient;				
	uniform	vec4 diffuse;				
	uniform	vec4 specular;				
	float shininess;					
} Material;

uniform int wireframe;

// Terrain textures
uniform sampler2D tex1;
uniform sampler2D tex2;
uniform sampler2D tex3;
uniform sampler2D tex4;

// Axel Terrain Textures
uniform sampler2D tex5;
uniform sampler2D tex6;
uniform sampler2D tex7;
uniform sampler2D tex8;

uniform vec3 light_dir_axel;
uniform int applyNormalMapping;
uniform float normalMappingStrength;
uniform float textureScaling;

// Texture
uniform sampler2D textmat;

// Environment
uniform sampler2D environment_front;
uniform sampler2D environment_back;
uniform sampler2D environment_left;
uniform sampler2D environment_right;
uniform sampler2D environment_down;
uniform sampler2D environment_up;

in vec3 position_frag;
in vec3 color_frag; 
in vec3 normal_frag;
in vec2 UV_frag;
in vec3 lightDir_frag;
in vec3 eyeVec_frag;

out vec4 frag_color;

in vec3 dist;
in float ratio;

const vec4 WIRE_COL = vec4(0.1,0.1,0.1,1);
const vec4 FILL_COL = vec4(0.9,0.9,0.9,0.0);


subroutine vec4 shadeModelType( vec3 color, vec3 normal);
subroutine uniform shadeModelType shadeModel;

// Compute smooth diffuse color
// normal : Normal vector
// lighting : Lighting vector
float Diffuse(in vec3 normal,in vec3 lighting)
{
  // Modified diffuse lighting
  float d = 0.5*(1.0+dot(normal,lighting));
  return d*d;
}

subroutine( shadeModelType )
vec4 phongModel( vec3 col, vec3 norm )
{
  // Returned color set to ambient
  vec3 c = Material.ambient.xyz;
  
  // Modified diffuse lighting
  float d = Diffuse(norm,lightDir_frag);
   
  c +=  Material.diffuse.xyz * vec3(d,d,d);

  // Specular 
  vec3 r = reflect(lightDir_frag, norm);
  float l=0.5*(1.0+dot(r, eyeVec_frag));
  l=l*l;
  float s = pow(l, Material.shininess);
  c +=  Material.specular.xyz * vec3(s,s,s);
  
  return vec4(c,Material.ambient.a);
}


subroutine( shadeModelType )
vec4 phongVertexColorModel( vec3 col, vec3 norm )
{
  // Default : Ambient color
  vec3 final_color = col;
  
 // float lambertTerm = dot(norm,lightDir_frag);
 // if (lambertTerm > 0.0)
 // {
 //   // Diffuse color
	//final_color +=  col * vec3(lambertTerm,lambertTerm,lambertTerm);
	//// Specular color
	//vec3 R = reflect(-lightDir_frag, norm);
	//float specular = pow(max(dot(R, eyeVec_frag), 0.0), Material.shininess);
	//final_color +=  col * vec3(specular,specular,specular);
 // }
  return vec4(final_color,1.0);
}

subroutine( shadeModelType )
vec4 goochModel( vec3 col, vec3 norm )
{	
	vec3 warm = vec3(0.85,0.75,0.55);
	vec3 cold = vec3(0.45,0.55,0.8);

	// Reflected ray
	vec3 r= reflect(lightDir_frag, norm);

	float d=Diffuse(norm, lightDir_frag);

	float s=pow((1 + dot(r, eyeVec_frag)) / 2.0,20.0);

	vec3 diffuse= d*warm+(1-d)*cold ;
	vec3 specular= s*(warm+vec3(0.1,0.1,0.1))+(1-s)*(cold+vec3(0.1,0.1,0.1));
 
	return vec4(0.85*diffuse+0.15*specular, 1.0);
}

subroutine( shadeModelType ) // Tri planar mapping with colors
vec4 EricTriPlanar( vec3 col, vec3 norm )
{	
  vec3 g = abs(norm);
  g=pow(g,vec3(8.0));
  g /= (g[0] + g[1] + g[2]);

  vec3 color=vec3(0.5)+vec3(0.5,0.0,0.)*g.x+vec3(0.0,0.5,0.)*g.y+vec3(0.0,0.0,0.5)*g.z;
  return vec4(color, 1.0);
}

// Shade according to normals
// Shading range lies within (0.6,0.6,0.6) (1.0,1.0,1.0) range
subroutine( shadeModelType )
vec4 normalColorModel(vec3 col, vec3 norm)
{
  // Color computed directly from the normal
  vec4 c = vec4(0.2*(vec3(3.0,3.0,3.0)+2.0*norm),1.0);
	
  return c;
}

// Three dimensional grid
subroutine( shadeModelType )
vec4 gridModel( vec3 col, vec3 norm )
{
	float steps=10.0;
	vec3 colorLine = vec3(0.2,0.1,0.1);
	vec3 final_color;
	vec3 pos=vec3(abs(position_frag.x),abs(position_frag.y),abs(position_frag.z));
	
	if ((pos.x-float(int(pos.x/steps))*steps)<0.1) final_color = colorLine;
	else{
		if ((pos.y-float(int(pos.y/steps))*steps)<0.1) final_color = colorLine;
		else{
			if ((pos.z-float(int(pos.z/steps))*steps)<0.1) final_color = colorLine;
			else{
				final_color = Material.ambient.xyz;
		  
				float lambertTerm = dot(norm,lightDir_frag);
				if (lambertTerm > 0.0)
				{
					// Diffuse color
					final_color +=  Material.diffuse.xyz * vec3(lambertTerm,lambertTerm,lambertTerm);
					// Specular color
					vec3 R = reflect(-lightDir_frag, norm);
					float specular = pow(max(dot(R, eyeVec_frag), 0.0), Material.shininess);
					final_color +=  Material.specular.xyz * vec3(specular,specular,specular);
				}
			}
		}
	}
	return vec4(final_color,1.0);
}


subroutine( shadeModelType )
vec4 ignModel( vec3 col, vec3 norm )
{

	float z=position_frag.z;
	vec4 color_=vec4(col,1.0);
	vec4 final_color;

	const float heightRange = 100.0;
	const float offset = 0.0;
	const int nbStep = 4;
	const float rangeFullColored = 0.01;
	const float rangeGradientColored = 0.01;
	
	//how much will be multiplied the color in ambient area
	//const float coefBackLighting = 0.15;
	const float coefMultAmbient = 0.4; //0.53

	const vec4 colorLine    = vec4(205.0 / 255.0, 195.0 / 255.0, 195.0 / 255.0, 1.0);
	const vec4 colorNeutral = vec4(235.0 / 255.0, 235.0 / 255.0, 235.0 / 255.0, 1.0);
	const vec4 colorForest  = vec4(135.0 / 255.0, 165.0 / 255.0, 125.0 / 255.0, 1.0);

	// Get a float between 0.0 and 1.0 (0.0 = bottom of a 'slice', 1.0 = top)	
	float height = z;
	height = height - offset;
	
	height = height / heightRange;
	height = height * nbStep;
	height = fract(height);
	
	// Defining colors according to the position in 0.0 to 1.0
	float beginTr  = - height + rangeGradientColored;
	beginTr = clamp (beginTr, 0.0, rangeGradientColored);
	beginTr = beginTr / rangeGradientColored;

	float endTr = height - (1.0 - rangeGradientColored - rangeFullColored);
	endTr = clamp (endTr, 0.0, rangeGradientColored);
	endTr = endTr / rangeGradientColored;

	float ratioColorLine = max(beginTr, endTr);
	float ratioColorNeutral = 1.0 - ratioColorLine;

	// Diffuse color
	vec4 usedNeutralColor = (color_.b * colorForest) + ((1 - color_.b) * colorNeutral);
	vec4 colorDiffuse = (ratioColorLine * colorLine) + (ratioColorNeutral * usedNeutralColor);
	float dotLight = dot(norm, lightDir_frag);
	
	dotLight = clamp (dotLight, 0.0, 1.0);
	dotLight = abs(dotLight*0.90);

	//multAmbient goes from coefMultAmbient to 1.0
	float multAmbient = ((1.0 - coefMultAmbient) * dotLight) + coefMultAmbient;

	final_color = multAmbient * colorDiffuse;
	final_color.a = 1.0;
	
	return final_color;
}

subroutine( shadeModelType )
vec4 terrain3DModel( vec3 col, vec3 norm )
{
	vec3 worldNormal = normalize(norm);
	float Stretching = 0.3;
	
	// Generate planes blending coeffs
	vec3 absNormal   = abs(worldNormal);
	vec3 blendnormal = normalize(pow(absNormal, vec3(8.0)));
	float blendSqrX  = blendnormal.x * blendnormal.x;
	float blendSqrY  = blendnormal.y * blendnormal.y;
	float blendSqrZ  = blendnormal.z * blendnormal.z;
	
	// Generate planes uvs
	vec2 texCoordX = vec2(-position_frag.z, -position_frag.y);
	vec2 texCoordY = vec2(-position_frag.x, -position_frag.z);
	vec2 texCoordZ = vec2( position_frag.x, -position_frag.y);
	
	// Sample color maps
	vec4 color1 = vec4(0.0);
	color1 += texture2D(tex1, texCoordX.xy * Stretching) * blendSqrX;
	color1 += texture2D(tex1, texCoordY.xy * Stretching) * blendSqrY;
	color1 += texture2D(tex1, texCoordZ.xy * Stretching) * blendSqrZ;
	
	vec4 color2 = vec4(0.0);
	color2 += texture2D(tex2, texCoordX.xy * Stretching) * blendSqrX;
	color2 += texture2D(tex2, texCoordY.xy * Stretching) * blendSqrY;
	color2 += texture2D(tex2, texCoordZ.xy * Stretching) * blendSqrZ;
	
	vec4 color3 = vec4(0.0);
	color3 += texture2D(tex3, texCoordX.xy * Stretching) * blendSqrX;
	color3 += texture2D(tex3, texCoordY.xy * Stretching) * blendSqrY;
	color3 += texture2D(tex3, texCoordZ.xy * Stretching) * blendSqrZ;
	
	vec4 color4 = vec4(0.0);
	color4 += texture2D(tex4, texCoordX.xy * Stretching) * blendSqrX;
	color4 += texture2D(tex4, texCoordY.xy * Stretching) * blendSqrY;
	color4 += texture2D(tex4, texCoordZ.xy * Stretching) * blendSqrZ;
	
	vec4 mixcolor = vec4(0.0);
	mixcolor = mix(color1  , color2, col.x);
	mixcolor = mix(mixcolor, color3, col.y);
	mixcolor = mix(mixcolor, color4, col.z);		
			
	// Diffuse lighting 
	float d = Diffuse(worldNormal, lightDir_frag);

	vec4 final_color = d*mixcolor;
	final_color.a = 1.0;
	return final_color;	
}

// Texture mapping shader
subroutine( shadeModelType )
vec4 UVTextureModel( vec3 col, vec3 norm )
{
	vec3 worldNormal = normalize(norm);

	// Sample color maps
	vec4 color1 = texture2D(textmat, UV_frag);
	float alpha = color1.a;
	color1.a=1.0;
		
	// Diffuse lighting 
	float d = Diffuse(worldNormal, lightDir_frag);


	// Calculating and return the Final Color
	vec4 final_color = 0.4*color1+0.6*d*color1;
	//final_color = color1;

	final_color.a = alpha;
	return final_color;	
}

// ajout Eric Guerin (perlin noise 3D)
vec3 mod289(vec3 x) {
  return x - floor(x * (1.0 / 289.0)) * 289.0;
}

vec4 mod289(vec4 x) {
  return x - floor(x * (1.0 / 289.0)) * 289.0;
}

vec4 permute(vec4 x) {
     return mod289(((x*34.0)+1.0)*x);
}

vec4 taylorInvSqrt(vec4 r)
{
  return 1.79284291400159 - 0.85373472095314 * r;
}

float snoise(vec3 v)
  {
  const vec2 C = vec2(1.0/6.0, 1.0/3.0) ;
  const vec4 D = vec4(0.0, 0.5, 1.0, 2.0);

// First corner
  vec3 i = floor(v + dot(v, C.yyy) );
  vec3 x0 = v - i + dot(i, C.xxx) ;

// Other corners
  vec3 g = step(x0.yzx, x0.xyz);
  vec3 l = 1.0 - g;
  vec3 i1 = min( g.xyz, l.zxy );
  vec3 i2 = max( g.xyz, l.zxy );

  // x0 = x0 - 0.0 + 0.0 * C.xxx;
  // x1 = x0 - i1 + 1.0 * C.xxx;
  // x2 = x0 - i2 + 2.0 * C.xxx;
  // x3 = x0 - 1.0 + 3.0 * C.xxx;
  vec3 x1 = x0 - i1 + C.xxx;
  vec3 x2 = x0 - i2 + C.yyy; // 2.0*C.x = 1/3 = C.y
  vec3 x3 = x0 - D.yyy; // -1.0+3.0*C.x = -0.5 = -D.y

// Permutations
  i = mod289(i);
  vec4 p = permute( permute( permute(
             i.z + vec4(0.0, i1.z, i2.z, 1.0 ))
           + i.y + vec4(0.0, i1.y, i2.y, 1.0 ))
           + i.x + vec4(0.0, i1.x, i2.x, 1.0 ));

// Gradients: 7x7 points over a square, mapped onto an octahedron.
// The ring size 17*17 = 289 is close to a multiple of 49 (49*6 = 294)
  float n_ = 0.142857142857; // 1.0/7.0
  vec3 ns = n_ * D.wyz - D.xzx;

  vec4 j = p - 49.0 * floor(p * ns.z * ns.z); // mod(p,7*7)

  vec4 x_ = floor(j * ns.z);
  vec4 y_ = floor(j - 7.0 * x_ ); // mod(j,N)

  vec4 x = x_ *ns.x + ns.yyyy;
  vec4 y = y_ *ns.x + ns.yyyy;
  vec4 h = 1.0 - abs(x) - abs(y);

  vec4 b0 = vec4( x.xy, y.xy );
  vec4 b1 = vec4( x.zw, y.zw );

  //vec4 s0 = vec4(lessThan(b0,0.0))*2.0 - 1.0;
  //vec4 s1 = vec4(lessThan(b1,0.0))*2.0 - 1.0;
  vec4 s0 = floor(b0)*2.0 + 1.0;
  vec4 s1 = floor(b1)*2.0 + 1.0;
  vec4 sh = -step(h, vec4(0.0));

  vec4 a0 = b0.xzyw + s0.xzyw*sh.xxyy ;
  vec4 a1 = b1.xzyw + s1.xzyw*sh.zzww ;

  vec3 p0 = vec3(a0.xy,h.x);
  vec3 p1 = vec3(a0.zw,h.y);
  vec3 p2 = vec3(a1.xy,h.z);
  vec3 p3 = vec3(a1.zw,h.w);

//Normalise gradients
  vec4 norm = taylorInvSqrt(vec4(dot(p0,p0), dot(p1,p1), dot(p2, p2), dot(p3,p3)));
  p0 *= norm.x;
  p1 *= norm.y;
  p2 *= norm.z;
  p3 *= norm.w;

// Mix final noise value
  vec4 m = max(0.6 - vec4(dot(x0,x0), dot(x1,x1), dot(x2,x2), dot(x3,x3)), 0.0);
  m = m * m;
  return 42.0 * dot( m*m, vec4( dot(p0,x0), dot(p1,x1),
                                dot(p2,x2), dot(p3,x3) ) );
}

// VOIR COMMENT APPELER LE SHADER DANS APPVOID
subroutine( shadeModelType )
vec4 bumpNormal( vec3 col, vec3 norm )
{
  vec3 normalbump = vec3(0.4,0.4,0.4)*vec3(snoise(position_frag*10.0),snoise(position_frag*10.0),snoise(position_frag*10.0));
  
  vec3 shade1noise = vec3(snoise(position_frag*20.0),snoise(position_frag*20.0),snoise(position_frag*20.0));
  vec3 shade2noise = vec3(snoise(position_frag*50.0),snoise(position_frag*50.0),snoise(position_frag*50.0));

  // Default : Ambient color
 vec3  final_color = vec3(0.0,0.0,0.0);
  
  float lambertTerm = 0.5*(1.0+dot(norm,lightDir_frag));
  lambertTerm = lambertTerm*lambertTerm;
  //float lambertTerm = dot(normalize(norm+normalbump),lightDir_frag);
   // Interpolate
   vec3 t= vec3(1,1,1) + shade1noise * 0.05 + shade2noise * 0.025;
    // Diffuse color
	final_color +=  Material.ambient.xyz + (Material.diffuse.xyz* t-Material.ambient.xyz)  * vec3(lambertTerm,lambertTerm,lambertTerm);

	// Specular color
	vec3 R = reflect(lightDir_frag, normalize(norm+normalbump));
    float specular = pow(max(dot(R, eyeVec_frag), 0.0), Material.shininess);
	final_color +=  Material.specular.xyz * vec3(specular,specular,specular);
 
  return vec4(final_color,1.0);
}

subroutine( shadeModelType )
vec4 ThinSnow( vec3 col, vec3 norm )
{
  vec3 normalbump = 0.1*vec3(snoise(position_frag*70.0),snoise(position_frag*70.0),snoise(position_frag*70.0));

  // Default : Ambient color
  vec3 final_color = Material.specular.xyz;
  

  float lambertTerm = dot(norm+normalbump,lightDir_frag);
  if (lambertTerm > 0.0)
  {
    // Diffuse color
	final_color +=  Material.diffuse.xyz * vec3(lambertTerm,lambertTerm,lambertTerm);
	// Specular color
	vec3 R = reflect(-lightDir_frag, norm+normalbump);
    float specular = pow(max(dot(R, eyeVec_frag), 0.0), Material.shininess);
	final_color +=  Material.specular.xyz * vec3(specular,specular,specular);
  }
  
  float alpha;

  if (col.r>0.9999)
    alpha = 1.0;
  else if (col.r<0.001)
    alpha = 0.0;
  else if (snoise(position_frag*100.0)<(col.r*2.0-1.0) || snoise(position_frag*200.0)<(col.r*2.0-1.0))
  	alpha = 1.0;
  else
  	alpha = 0.3;
  return vec4(final_color,alpha);
}


// *******************************************************************************************************************************************************************************


struct lightSource
{
  vec3 position;
  vec3 diffuse;
  vec3 specular;
  float constantAttenuation, linearAttenuation, quadraticAttenuation;
  float spotCutoff, spotExponent;
  vec3 spotDirection;
};
const int numberOfLights = 4;
lightSource lights[numberOfLights];
lightSource light0 = lightSource(
  vec3(0.0,  0.0, -10000.0),
  vec3(0.9,  0.7,  0.8),
  vec3(0.1,  0.1,  0.1),
  0.0, 1.0, 0.0,
  180.0, 0.0,
  vec3(0.0, 0.0, 0.0)
);
lightSource light1 = lightSource(
    vec3(0.0, -10000.0,  10000.0),
    vec3(0.7,  0.7,  0.9),
    vec3(0.1,  0.1,  0.1),
    0.0, 1.0, 0.0,
    80.0, 10.0,
    vec3(0.0, 1.0, 0.0)
);
  
  lightSource light2 = lightSource(
    vec3(0.0,  10000.0,  10000.0),
    vec3(0.7,  0.7,  0.8),
    vec3(0.1,  0.1,  0.1),
    0.0, 1.0, 0.0,
    80.0, 10.0,
    vec3(0.0, 1.0, 0.0)
);

  lightSource light3 = lightSource(
    vec3(10000.0,  0.0,  10000.0),
    vec3(0.7,  0.7,  0.7),
    vec3(0.1,  0.1,  0.1),
    0.0, 1.0, 0.0,
    80.0, 10.0,
    vec3(0.0, 1.0, 0.0)
);

// Soft grey shading with cool warm tint
subroutine( shadeModelType )
vec4 GreyCoolWarm( vec3 col, vec3 norm )
{
  lights[0] = light0;
  lights[1] = light1;
  lights[2] = light2;
  lights[3] = light3;

  // Default : Ambient color
  vec3 final_color = Material.ambient.xyz;
  
  for (int index = 0; index < numberOfLights; index++) // for all light sources
  {
	vec3 positionToLightSource = vec3(lights[index].position - position_frag.xyz);
    vec3 lightDirection = normalize(positionToLightSource);

    // Smooth lighting coefficient
	float c = 0.5*(1.0+dot(norm,lightDirection));
    c*=c;
    // Diffuse color
	final_color +=  Material.diffuse.xyz * lights[index].diffuse * vec3(c,c,c);
  }
  return vec4(final_color,Material.diffuse.a);
}

// Soft grey shading with cool warm tint
subroutine( shadeModelType )
vec4 WireframeAspect( vec3 col, vec3 norm )
{
	lights[0] = light0;
	lights[1] = light1;
	lights[2] = light2;
	lights[3] = light3;

	// Default : Ambient color
	vec3 final_color = Material.ambient.xyz;
  
	for (int index = 0; index < numberOfLights; index++) // for all light sources
	{
		vec3 positionToLightSource = vec3(lights[index].position - position_frag.xyz);
		vec3 lightDirection = normalize(positionToLightSource);

		// Smooth lighting coefficient
		float c = 0.5*(1.0+dot(norm,lightDirection));
		c*=c;
		// Diffuse color
		final_color +=  Material.diffuse.xyz * lights[index].diffuse * vec3(c,c,c);
	}
	return vec4(final_color,Material.diffuse.a) + vec4(pow(1.0-ratio,3.0),0.0,0.0,1.0);
}

// *****************************************************************************************************************************************************
subroutine( shadeModelType )
vec4 Wireframe( vec3 col, vec3 norm )
{
    float d = min(dist[0],min(dist[1],dist[2]));
    float I = exp2(-1*d*d);
	return I*WIRE_COL + (1.0 - I)*FILL_COL;
}

void main()
{
	// Call shader model
	vec4 final_color = shadeModel( color_frag, normal_frag );
	
	// Display Wireframe
	if (wireframe==1) {
		float d = min(dist[0],min(dist[1],dist[2]));
		float I = exp2(-7.0*d*d);
		final_color = vec4( (I*0.85)*final_color.xyz +(1.0 - I)*final_color.xyz,1.0);
	}
	frag_color = final_color;
}


// *****************************************************************************************************************************************************
// Axel Terrain Shader
const vec3 AxelLight = normalize(vec3(1.0f, 1.0f, 5.0f));

vec3 UnpackNormal(vec3 packedNormal)
{
	vec3 unpacked = packedNormal.xyz * 2.0 - 1.0;
	// Normal are inverted on y-axis - may be Qt does something when we call
	// QGLWidget::convertToGLFormat(niceRockNormal);
	unpacked.y *= -1;
	return unpacked;
}

float Saturate(float v)
{
	return clamp(v, 0.0, 1.0);
}

vec3 GetTriplanarBlendingCoefficient(vec3 worldNormal)
{
	// Generate planes blending coeffs
	vec3 absNormal   = abs(worldNormal);
	vec3 blendnormal = normalize(pow(absNormal, vec3(8.0)));
	float blendSqrX  = blendnormal.x * blendnormal.x;
	float blendSqrY  = blendnormal.y * blendnormal.y;
	float blendSqrZ  = blendnormal.z * blendnormal.z;	
	return vec3(blendSqrX, blendSqrY, blendSqrZ);
}

vec3 GetTriplanarNormal(sampler2D tex, vec3 worldNormal, vec3 blendCoeff)
{
	// Generate planes uvs
	vec2 texCoordX = vec2(-position_frag.z, -position_frag.y);
	vec2 texCoordY = vec2(-position_frag.x, -position_frag.z);
	vec2 texCoordZ = vec2( position_frag.x, -position_frag.y);
	
	// Tangent space normal maps
	vec3 tnormalX = UnpackNormal(texture2D(tex, texCoordX * textureScaling).xyz);
	vec3 tnormalY = UnpackNormal(texture2D(tex, texCoordY * textureScaling).xyz);
	vec3 tnormalZ = UnpackNormal(texture2D(tex, texCoordZ * textureScaling).xyz);	
	
	// Swizzle world normals into tangent space and apply UDN blend.
	// These should get normalized, but it's very a minor visual
	// difference to skip it until after the blend.
	tnormalX = vec3(tnormalX.xy + worldNormal.zy, worldNormal.x);
	tnormalY = vec3(tnormalY.xy + worldNormal.xz, worldNormal.y);
	tnormalZ = vec3(tnormalZ.xy + worldNormal.xy, worldNormal.z);
	
	// Swizzle tangent normals to match world orientation and triblend
	return normalize(
		tnormalX.zyx * blendCoeff.x +
		tnormalY.xzy * blendCoeff.y +
		tnormalZ.xyz * blendCoeff.z
    );
}

vec4 GetTriplanarAlbedo(sampler2D tex, vec3 blendCoeff)
{
	// Generate planes uvs
	vec2 texCoordX = vec2(-position_frag.z, -position_frag.y);
	vec2 texCoordY = vec2(-position_frag.x, -position_frag.z);
	vec2 texCoordZ = vec2( position_frag.x, -position_frag.y);
	
	vec4 texCol = vec4(0.0);
	texCol += texture2D(tex, texCoordX.xy * textureScaling) * blendCoeff.x;
	texCol += texture2D(tex, texCoordY.xy * textureScaling) * blendCoeff.y;
	texCol += texture2D(tex, texCoordZ.xy * textureScaling) * blendCoeff.z;
	return texCol;
}

subroutine(shadeModelType)
vec4 AxelTerrain(vec3 col, vec3 norm)
{
	vec3 worldNormal = normalize(norm);
	
	// Trick : normalMappingStrength to 0.0f means we want to show shaded geology
	if (normalMappingStrength == 0.0f)
	{
		vec4 final_color = vec4(col, 1.0);
		final_color.a = 1.0;
		return final_color;
	}
	// Standard Terrain shading - with our without triplanar normal mapping
	else
	{
		vec3 blendCoeff = GetTriplanarBlendingCoefficient(worldNormal);
		
		// Albedo & Normal mix
		vec3 triplanarNormal1 = GetTriplanarNormal(tex6, worldNormal, blendCoeff);
		vec3 triplanarNormal2 = GetTriplanarNormal(tex8, worldNormal, blendCoeff);
	
		vec4 triplanarAlbedo1 = GetTriplanarAlbedo(tex5, blendCoeff);
		vec4 triplanarAlbedo2 = GetTriplanarAlbedo(tex7, blendCoeff);
		
		vec4 mixColor  	= mix(triplanarAlbedo1, triplanarAlbedo2, col.x);
		vec3 mixNormal 	= mix(triplanarNormal1, triplanarNormal2, col.x);
		mixNormal 		= mix(worldNormal, mixNormal, normalMappingStrength);
		
		// Lighting
		float ndotl = 0.0;
		if (applyNormalMapping == 1)
			ndotl = Saturate(dot(mixNormal, normalize(light_dir_axel)));
		else
			ndotl = Saturate(dot(worldNormal, normalize(light_dir_axel)));

		vec4 final_color = ndotl * mixColor;
		final_color.a = 1.0;
		return final_color;
	}
}
