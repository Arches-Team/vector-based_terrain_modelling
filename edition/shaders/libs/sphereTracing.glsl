#version 430 core

#define M_PI 3.14159

//If the SDF and PRIMITIVES where not defined, create a default object to compile the shader.
#ifndef PRIMITIVES
#define PRIMITIVES
#endif


//default sdf, a BIG RED CROSS
#ifndef SDF
#define UPLEFT vec3(-0.7071, 0, 0.7071)
#define UPRIGHT vec3(0.7071, 0, 0.7071)
#define SDF max(min(length(p-UPLEFT*(dot(p, UPLEFT)))-1, length(p-UPRIGHT*(dot(p, UPRIGHT)))-1), length(p)-10)
#define NOSDF
#endif

//default cost, the cost of the cross
#ifndef SDFCOST
#define SDFCOST 3
#endif


#ifndef GRAD
#define GRAD normalize(vec3(\
	(sdf(p + vec3(epsilon, 0, 0)) - sdf(p - vec3(epsilon, 0, 0))) / (2*epsilon),\
	(sdf(p + vec3(0, epsilon, 0)) - sdf(p - vec3(0, epsilon, 0))) / (2*epsilon),\
	(sdf(p + vec3(0, 0, epsilon)) - sdf(p - vec3(0, 0, epsilon))) / (2*epsilon)\
	))
#endif

#ifndef MATERIAL
#define MATERIAL material(vec3(1.0, 1.0, 1.0), 0.8, 4.);
#endif


uniform mat4 invMatrix;
uniform ivec2 iResolution;
uniform vec2 zfarnear;
uniform vec3 cameraPos;
uniform int maxSteps;
uniform int shading;
uniform vec3 sampleFieldNormal;
uniform vec3 sampleFieldPoint;
uniform float sampleFieldFreq;
uniform int maxCost;
uniform float kg;
uniform float kn;
uniform int aa;

uniform vec3 light_dir;
uniform vec3 light_color;
uniform vec3 ambiant_color;
uniform int nb_sample;
uniform float seed;
uniform int nb_bounce;
uniform ivec2 texture_size;



float hash13(vec3 p)
{
	p = fract(p* .1031);
	p += dot(p, p.zyx + 31.32);
	return fract((p.x + p.y) * p.z);
}



#define NORMAL 0
#define NBSTEPS 1
#define PLANE_SDF 2
#define PLANE_GRAD 4
#define COST 3
#define SEXY 5

const float epsilon = 0.001;

struct material{
	vec3 diffuse;//color of the object
	float kd; //diffuse coeficient (1 diffuse, 0 mirror)
	float ns; //specular coeficient (0diffuse, infinity : perfect mirror)
};


#ifdef VERTEX_SHADER

void main(void)
{
	vec4 vertices[4] = vec4[4](vec4(-1.0, -1.0, 1.0, 1.0),
                               vec4( 1.0, -1.0, 1.0, 1.0),
                               vec4(-1.0,  1.0, 1.0, 1.0),
                               vec4( 1.0,  1.0, 1.0, 1.0));
    vec4 pos = vertices[gl_VertexID];
    gl_Position = pos;

} 
#endif


#ifdef FRAGMENT_SHADER

out vec4 fragment;

uniform sampler2D readColor;

float m_seed = seed;

float random()
{
	m_seed = hash13(vec3( gl_FragCoord.xy, m_seed));
	return m_seed;
}

vec3 random_direction(vec3 n)
{
	float r1 = random();
	float r2 = random();

	float sqrt_r2 = sqrt(r2*(1-r2));
	float x = 2*cos(2*M_PI*r1)*sqrt_r2;
	float y = 2*sin(2*M_PI*r1)*sqrt_r2;
	float z = 1 - 2*r2;
	vec3 d = vec3(x, y, z);
	return dot(d, n) > 0 ? d : -d;
}

mat3 world_matrix(vec3 n)
{
	float signe= sign(n.z);
    float a= -1.0f / (signe + n.z);
    float d= n.x * n.y * a;
    vec3 t= vec3(1.0f + signe * n.x * n.x * a, signe * d, -signe * n.x);
    vec3 b= vec3(d, signe + n.y * n.y * a, -n.y);
	return mat3(t, b, n);
}

vec3 random_direction_brdf(vec3 incoming_d, vec3 n, float kd, float ns)
{
	float u1 = random();
	float u2 = random();
	float u3 = random();
	mat3 w = world_matrix(n);
	if(u1 <= kd)
    {
        // terme diffus
        // genere une direction cos theta / pi, cf GI compendium, eq 35
        float phi= float(2*M_PI) * u3;
        float cos_theta= sqrt(u2);
        float sin_theta= sqrt(1 - cos_theta*cos_theta);
            
        return w*vec3(cos(phi) * sin_theta, sin(phi) * sin_theta, cos_theta); //need to rotate to be along n
    }
    else
    {
        // terme reflechissant
        // genere une direction h
        // genere une direction cos^n theta / pi, cf GI compendium, eq 35+
        float phi= float(2*M_PI) * u3;
        float cos_theta= pow(u2, 1 / float(ns +1));
        float sin_theta= sqrt(1 - cos_theta*cos_theta);
            
        vec3 h= w*vec3(cos(phi) * sin_theta, sin(phi) * sin_theta, cos_theta);
        // genere une direction reflechie, l= reflect(o | h)
        return -incoming_d + 2 * dot(h, incoming_d) * h;    
    }
}

float pdf_brdf(vec3 d, vec3 incoming_d, vec3 n, float kd, float ns){
		float cos_theta= dot(n, d);
        if(cos_theta <= 0)
            return 0.;
        
        // pdf du terme diffus 
        float diffuse= cos_theta / M_PI;
        
        // pdf du terme reflechissant
        vec3 h= normalize(d - incoming_d);
        float cos_theta_h= dot(n, h);
        if(cos_theta_h <= 0)
            return 0;
        if(dot(d, h) <= 0)
            return 0;
        
        float s= (ns + 1.) / (2.*M_PI) * pow(cos_theta_h, ns) / (4. * dot(d, h));
        // le terme 1 / dot(l, h) est introduit par le changement de variable : direction h vers direction reflechie 
        
        // pdf de la mixture
        return kd*diffuse + (1 - kd)*s;
//return 0.;
}

PRIMITIVES


float sdf(vec3 p)
{	
	return SDF;
}

int sdfCost(vec3 p)
{
	return SDFCOST;
}


vec3 grad(vec3 p)
{
	return GRAD;
}

material getMaterial(vec3 p, vec3 n)
{
	return MATERIAL;
}



bool Raymarch(vec3 o, vec3 d, int STEPS, out float t, out int nb_steps, out int cost)
{
	t = 0;
	nb_steps = 0;
	cost = 0;
	float tmp_prev;
	float tmp = sdf(o);
	cost += sdfCost(o);
	for(nb_steps = 0; nb_steps < STEPS; nb_steps++)
	{
		tmp_prev = tmp;
		tmp = sdf(o);
		cost += sdfCost(o);
		if(tmp*tmp_prev <= 0) //changed signe
			return true;
		
		if(t > zfarnear[0]-zfarnear[1])
			return false;
		
		t+=abs(tmp)+epsilon;
		o = o + d*(abs(tmp)+epsilon);
	}
	return false;
}

vec4 coolWarm(float u)
{
  u = clamp(u, 0, 1);
  vec4 cool = vec4(97.0/256.0, 130.0/256.0, 234.0/256.0, 1.0);
  vec4 white = vec4(221.0/256.0, 221.0/256.0, 221.0/256.0, 1.0);
  vec4 warm = vec4(220.0/256.0, 94.0/256.0, 75.0/256.0, 1.0);

  if (u < 0.5)
  {
	u = 2*u;
    	return mix(cool, white, u);
  }
  else
  {
	u = 2*u-1;
    	return mix(white, warm, u);
  }
}

// Compute sky color 
// d  Ray direction
vec3 SkyShadeBlue(in vec3 d)
{
  	// light direction
	//vec3 lig = normalize(vec3( 0.3,0.5, 0.6));
	float sun = 0.5*(dot(light_dir,d)+1.0);
	vec3 color = vec3(0.35,0.45,0.65)*(0.75-0.25*d.z);
	color += vec3(0.65,0.6,0.55)*pow( sun, 18.0 );
	return color;
}

float brdf(in vec3 incoming_d, in vec3 bounce_d, in vec3 n, in float kd, in float ns){
	float diffuse = 1/M_PI;
	vec3 h = normalize(-incoming_d + bounce_d); // incoming is pointing toward the surface and bounce is pointing outward
	float specular = (ns+8)/(8*M_PI)*pow(dot(h, n), ns);

	if(dot(incoming_d, n) > 0) return 0.;
	if(dot(bounce_d, n) < 0) return 0.;
	if(dot(h, n) < 0) return 0.;
	return kd*diffuse + (1-kd)*specular;
}

void main()
{
	vec2 coord = 2*gl_FragCoord.xy/iResolution - 1;


	//vec4(coord.xy, 0, 1) : homogeneous starting point of the ray in the projected space
	//vec4(coord.xy, 1, 1) : homogeneous ending point of the ray in the projected space

	vec4 oh = invMatrix*vec4(coord.xy, 0, 1); //homogeneous starting point of the ray in the object space
	vec4 eh = invMatrix*vec4(coord.xy, 1, 1); //homogeneous ending point of the ray in the object space

	vec3 o= oh.xyz / oh.w; //starting point of the ray in the object space
	vec3 e= eh.xyz / eh.w; //ending point of the ray in the object space

	vec3 d = e - o;
	d = normalize(d); //ray

	//compute the pixel 3D size:
	vec2 coord_pixel = 2*(gl_FragCoord.xy+ivec2(1, 1))/iResolution - 1;
	vec4 oh_pixel = invMatrix*vec4(coord_pixel.xy, 0, 1);
	vec3 o_pixel = oh_pixel.xyz / oh_pixel.w;
	float pixelsize = length(o_pixel-o)/(sqrt(2)*length(cameraPos-o)); //sqrt(2) to have the radius of the circle in the pixel


	float t = 0;
	int nb_steps = 0;
	int cost = 0;

	if(shading == PLANE_SDF || shading == PLANE_GRAD)
	{
		

		t = dot(sampleFieldPoint - o, sampleFieldNormal)/dot(d, sampleFieldNormal);
		if(dot(d,sampleFieldNormal) == 0 || t <= 0)
		{
			fragment = vec4(0, 0.01, 0, 1.0);
			gl_FragDepth = 0.99;
		}
		else
		{
			o = o + t*d;
			gl_FragDepth = ((1.0f / length(cameraPos-o)) - (1.0f / zfarnear[1])) / ((1.0f / zfarnear[0]) - (1.0f / zfarnear[1]));

			if(shading==PLANE_SDF)
			{
				float s = sdf(o);
				float v = cos(abs(s*sampleFieldFreq))*0.5 + 0.5;
				fragment = vec4(s < 0 ? v : 0, 0, s >=0 ? v : 0, 1.0);
			}
			else
			{
				vec3 n = grad(o);
				fragment = vec4((n+vec3(1))*0.5, 1);
			}
		}
			
	}
	else
	{
		int bounce_max = shading == SEXY ? nb_bounce : 1;
		if(shading == SEXY && nb_sample == -1)
			bounce_max = 0;
		vec3 color_total = vec3(0., 0., 0.); // for sexy shader, the color of the all path
		float depth = 0;
		vec3 accumulated_brdf = vec3(1., 1., 1.); //for sexy shader, the accumulated attenuation of the different bounces
		for(int b = 0; b < bounce_max; b++)
		{

			if(Raymarch(o, d, maxSteps, t, nb_steps, cost))
			{
				//hit the surface
		
				o = o + t*d;

				vec3 n = vec3(0., 0., 0.);
				if(shading == NORMAL || shading == SEXY)
					n = normalize(grad(o));
					//anti aliasing on the normal
					/*float relativePixelsize = pixelsize*length(cameraPos-o);

					vec3 n = vec3(0);

					vec3 b = normalize(cross(d, vec3(1, 1, 1)));
					vec3 t = normalize(cross(d, b));

					for(int i = -aa; i <= aa; i++)
						for(int j = -aa; j <= aa; j++)
							if(i*i + j*j <= aa*aa)
								n += grad(o+relativePixelsize*(i*b+j*t)/max(1, aa));*/

					
					//n = normalize(n-max(0, dot(d, n))*d);
				if(shading == NORMAL)
					fragment = vec4(0.2 * (vec3(3.0) + 2.0 * n), 1.0);


				if(shading == SEXY){
					o += n*epsilon; // to avoid self intersection
					material local_material = getMaterial(o, n);

					vec3 diffuse_color = local_material.diffuse;
					float kd = local_material.kd;
					float ns = local_material.ns;


					float tmp;
					int tmp_nb_step;
					int tmp_cost;
					bool inShadow = Raymarch(o+epsilon*n, light_dir, 100, tmp, tmp_nb_step, tmp_cost);
					
					float pdf_sun = 1; //allways cast a ray in the same direction, the proba of this ray is 1
					float brdf_sun = brdf(d, light_dir, n, kd, ns)*max(0, dot(n, light_dir));

					if(!inShadow)
						color_total += diffuse_color*light_color*accumulated_brdf*(brdf_sun/pdf_sun);

					//generate a bounced ray
					vec3 old_d = d;
					d = random_direction(n);
					//d = random_direction_brdf(d, n, kd, ns); // not working


					float pdf_bounce = 1/(2*M_PI);
					//float pdf_bounce = pdf_brdf(d, old_d, n, kd, ns);
					vec3 brdf_bounce = diffuse_color*brdf(old_d, d, n, kd, ns)*max(0, dot(n, d));

					accumulated_brdf *= brdf_bounce/pdf_bounce;
				}
				
				if(shading == COST){
					int nb_aa = 0;
					for(int i = -aa; i <= aa; i++)
						for(int j = -aa; j <= aa; j++)
							if(i*i + j*j <= aa*aa)
								nb_aa++;
					cost += 6*nb_aa*sdfCost(o); //time 6 because one grad computation use 6 sdf computation.
					fragment = coolWarm(cost/float(maxCost));
				}
				
				if(shading == NBSTEPS)
					fragment = vec4(nb_steps/float(maxSteps), nb_steps/float(maxSteps), nb_steps/float(maxSteps), 1);

				if(b == 0){
					
					depth = ((1.0f / length(cameraPos-o)) - (1.0f / zfarnear[1])) / ((1.0f / zfarnear[0]) - (1.0f / zfarnear[1]));
					if(shading != SEXY)
						gl_FragDepth = depth;
				}
			}
			else
			{

				if(shading == SEXY)
				{
					if(b == 0)
						depth = 1.01; //To not overwrite the skybox
						//discard;

						//vec3 m_sky_color = SkyShadeBlue(d);
					color_total += ambiant_color*accumulated_brdf;
					break;
				}

				if(shading == NORMAL)
					discard;

				if(shading == NBSTEPS)
					fragment = vec4(nb_steps/float(maxSteps), nb_steps/float(maxSteps), nb_steps/float(maxSteps), 1);
				if(shading == COST)
					fragment = coolWarm(cost/float(maxCost));
				gl_FragDepth = 0.99;
				
			}
		}
		if(shading == SEXY)
		{

			vec4 oldColor = texture(readColor, vec2(gl_FragCoord.xy)/texture_size);

			vec4 new_color = (vec4(color_total, depth)+(oldColor*float(nb_sample)))/(1.+float(nb_sample));


			if(nb_sample == -1){
				fragment = vec4(oldColor.rgb, 1.);
				gl_FragDepth = oldColor.a;
			}
			else{

				fragment = new_color;
				gl_FragDepth = 0.;
				if(any(isnan(new_color)) || any(isinf(new_color)))
					fragment = oldColor;
				}
		}
			
	}
	//2*gl_FragCoord.xy/iResolution - 1
	//if(coord.x < 0.)
	
	//if(2*gl_FragCoord.x < iResolution.x)
	//	fragment = vec4(0, 1, 0, 1);
		

#ifdef NOSDF
fragment = vec4(1, 0, 0, 1);
#endif
}

#endif
