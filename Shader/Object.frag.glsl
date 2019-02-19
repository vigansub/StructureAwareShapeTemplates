#version 400 core

#define FRAG_COLOR 0
#define MAX_LIGHTS 10

layout(location = FRAG_COLOR) out vec4 FragColor;



uniform struct Light 
{
    vec3 position;
    vec3 direction;
    vec3 intensity;
    vec3 attenuation;
    vec4 cone;
    int type;
    vec3 color;
} 
lights[MAX_LIGHTS];

uniform struct Material
{
    vec3 Ka;
    vec3 Kd;
    vec3 Ks;
    float Ns;
    int hasTex;
} material;

uniform vec3 camPos;
uniform float shadowIntensity;
uniform int renderMode;
uniform float alpha;
uniform int applyShadow;
uniform int numLights;
uniform int isSelected;
uniform int colorMode;

uniform sampler2D tex; 
uniform sampler2D shadowMap[MAX_LIGHTS];

in vec4 VertPosition;
in vec4 VertColor;
in vec4 VertTexture;
in vec4 VertNormal;
in vec4 VertShadowCoord[MAX_LIGHTS];

uniform float ScaleFactor = 0.8;

const float C1 = 0.429043;
const float C2 = 0.511664;
const float C3 = 0.743125;
const float C4 = 0.886227;
const float C5 = 0.247708;

////Old Town Square
//const vec3 L00 = vec3( 0.871297, 0.875222, 0.864470);
//const vec3 L1m1 = vec3( 0.175058, 0.245335, 0.312891);
//const vec3 L10 = vec3( 0.034675, 0.036107, 0.037362);
//const vec3 L11 = vec3(-0.004629, -0.029448, -0.048028);
//const vec3 L2m2 = vec3(-0.120535, -0.121160, -0.117507);
//const vec3 L2m1 = vec3( 0.003242, 0.003624, 0.007511);
//const vec3 L20 = vec3(-0.028667, -0.024926, -0.020998);
//const vec3 L21 = vec3(-0.077539, -0.086325, -0.091591);
//const vec3 L22 = vec3(-0.161784, -0.191783, -0.219152); 

////Grace Cathedral
//const vec3 L00 = vec3(0.79, 0.44, 0.54);
//const vec3 L1m1 = vec3(0.39, 0.35, 0.60);
//const vec3 L10 = vec3(-0.34, -0.18, -0.27);
//const vec3 L11 = vec3(-0.29, -0.06, 0.01);
//const vec3 L2m2 = vec3(-0.11, -0.05, -0.12);
//const vec3 L2m1 = vec3(-0.26, -0.22, -0.47);
//const vec3 L20 = vec3(-0.16, -0.09, -0.15);
//const vec3 L21 = vec3(0.56, 0.21, 0.14);
//const vec3 L22 = vec3(0.21, -0.5, -0.30); 

//Eucalyptus Grove
const vec3 L00  = vec3(0.38, 0.43, 0.45);
const vec3 L1m1 = vec3(0.29, 0.36, 0.41);
const vec3 L10  = vec3(0.04, 0.03, 0.01);
const vec3 L11  = vec3(-0.1, -0.1, -0.09);
const vec3 L2m2 = vec3(-0.06, -0.06, 0.04);
const vec3 L2m1 = vec3(0.01, 0.01, -0.05);
const vec3 L20  = vec3(-0.09, -0.13, -0.15);
const vec3 L21  = vec3(-0.06, -0.05, -0.04);
const vec3 L22  = vec3(0.02, 0.0, -0.05); 

////St. Peter's Basilica
//const vec3 L00 = vec3(0.36, 0.26, 0.23);
//const vec3 L1m1 = vec3(0.18, 0.14, 0.13);
//const vec3 L10 = vec3(-0.02, -0.01, 0.0);
//const vec3 L11 = vec3(0.03, 0.02, 0.0);
//const vec3 L2m2 = vec3(0.02, 0.01, 0.0);
//const vec3 L2m1 = vec3(-0.05, -0.03, -0.01);
//const vec3 L20 = vec3(-0.09, -0.08, -0.07);
//const vec3 L21 = vec3(0.01, 0.0, 0.0);
//const vec3 L22 = vec3(-0.08, -0.03, 0.0); 

////Uffizi Gallery
//const vec3 L00 = vec3(0.32, 0.31, 0.35);
//const vec3 L1m1 = vec3(0.37, 0.37, 0.43);
//const vec3 L10 = vec3(0.0, 0.0, 0.0);
//const vec3 L11 = vec3(-0.01, -0.01, -0.01);
//const vec3 L2m2 = vec3(-0.02, -0.02, -0.03);
//const vec3 L2m1 = vec3(-0.01, -0.01, -0.01);
//const vec3 L20 = vec3(-0.28, -0.28, -0.32);
//const vec3 L21 = vec3(0.00, 0.00, 0.00);
//const vec3 L22 = vec3(-0.24, -0.24, -0.28); 

////Gallileo's Tomb
//const vec3 L00 = vec3(1.04, 0.76, 0.71);
//const vec3 L1m1 = vec3(0.44, 0.34, 0.34);
//const vec3 L10 = vec3(-0.22, -0.18, -0.17);
//const vec3 L11 = vec3(0.71, 0.54, 0.56);
//const vec3 L2m2 = vec3(0.64, 0.50, 0.52);
//const vec3 L2m1 = vec3(-0.12, -0.09, -0.08);
//const vec3 L20 = vec3(-0.37, -0.28, -0.29);
//const vec3 L21 = vec3(-0.17, -0.13, -0.13);
//const vec3 L22 = vec3(0.55, 0.42, 0.42); 

////Vine Street Kitchen
//const vec3 L00 = vec3(0.64, 0.67, 0.73);
//const vec3 L1m1 = vec3(0.28, 0.32, 0.33);
//const vec3 L10 = vec3(0.42, 0.60, 0.77);
//const vec3 L11 = vec3(-0.05, -0.04, -0.02);
//const vec3 L2m2 = vec3(-0.10, -0.08, -0.05);
//const vec3 L2m1 = vec3(0.25, 0.39, 0.53);
//const vec3 L20 = vec3(0.38, 0.54, 0.71);
//const vec3 L21 = vec3(0.06, 0.01, -0.02);
//const vec3 L22 = vec3(-0.03, -0.02, -0.03); 

////Breezeway
//const vec3 L00 = vec3(0.32, 0.36, 0.38);
//const vec3 L1m1 = vec3(0.37, 0.41, 0.45);
//const vec3 L10 = vec3(-0.01, -0.01, -0.01);
//const vec3 L11 = vec3(-0.10, -0.12, -0.12);
//const vec3 L2m2 = vec3(-0.13, -0.15, -0.17);
//const vec3 L2m1 = vec3(-0.01, -0.02, 0.02);
//const vec3 L20 = vec3(-0.07, -0.08, -0.09);
//const vec3 L21 = vec3(0.02, 0.03, 0.003);
//const vec3 L22 = vec3(-0.29, -0.32, -0.36); 

////Campus Sunset 
//const vec3 L00  = vec3(0.79, 0.94, 0.98);
//const vec3 L1m1 = vec3(0.44, 0.56, 0.70);
//const vec3 L10  = vec3(-0.10, -0.18, -0.27);
//const vec3 L11  = vec3(0.45, 0.38, 0.20);
//const vec3 L2m2 = vec3(0.18, 0.14, 0.05);
//const vec3 L2m1 = vec3(-0.14, -0.22, -0.31);
//const vec3 L20  = vec3(-0.39, -0.40, -0.36);
//const vec3 L21  = vec3(0.09, 0.07, 0.04);
//const vec3 L22  = vec3(0.67, 0.67, 0.52); 

////Funston Beach Sunset
//const vec3 L00 = vec3(0.68, 0.69, 0.70);
//const vec3 L1m1 = vec3(0.32, 0.37, 0.44);
//const vec3 L10 = vec3(-0.17, -0.17, -0.17);
//const vec3 L11 = vec3(-0.45, -0.42, -0.34);
//const vec3 L2m2 = vec3(-0.17, -0.17, -0.15);
//const vec3 L2m1 = vec3(-0.08, -0.09, -0.10);
//const vec3 L20 = vec3(-0.03, -0.02, -0.01);
//const vec3 L21 = vec3(0.16, 0.14, 0.10);
//const vec3 L22 = vec3(0.37, 0.31, 0.20); 

vec3 sphericalHarmonic()
{
    vec3 tnorm = VertPosition.xyz;

    vec3 diffuse = C1 * L22 * (tnorm.x * tnorm.x - tnorm.y * tnorm.y) + 
    C3 * L20 * tnorm.z * tnorm.z +
    C4 * L00 -
    C5 * L20 +
    2.0 * C1 * L2m2 * tnorm.x * tnorm.y +
    2.0 * C1 * L21 * tnorm.x * tnorm.z +
    2.0 * C1 * L2m1 * tnorm.y * tnorm.z +
    2.0 * C2 * L11 * tnorm.x +
    2.0 * C2 * L1m1 * tnorm.y +
    2.0 * C2 * L10 * tnorm.z;

    diffuse *= ScaleFactor; 

    return diffuse;
}

float lightAttenuation(vec3 P, vec3 L, vec3 attenuation, float intensity)
{
    float d = distance(P, L);
    float r = 1 / (attenuation.x + attenuation.y* d + attenuation.z * d * d) * intensity;
    return r;
}

float dualConeSpotlight(vec3 P, vec3 LP, vec3 LD, float outerCone, float innerCone)
{
    vec3 V = normalize(P - LP);
    float d = dot(V, LD);
    return smoothstep(outerCone, innerCone, d);
}

float VSM(vec4 smcoord, sampler2D sm)
{
	vec3 coords = smcoord.xyz / smcoord.w ;    

    if(smcoord.z < 1)
        return 1;

    float depth = coords.z;
    
    vec4 depthBlurrred = texture(sm, coords.xy);

    float depthSM = depthBlurrred.x;
	float sigma2  = depthBlurrred.y;

    float realDepth = texture(sm, coords.xy).x;

	sigma2 -= depthSM * depthSM;

	float bias = 0.0001;

	float dist = depth - depthSM;
	float P = sigma2 / ( sigma2 + dist * dist );
	float lit = max( P, ( depth - bias ) <= depthSM ? 1.0 : 0.0);
	lit = min(1.0, lit);   

    return mix(shadowIntensity, 1.0, lit);

    return lit;
}

vec3 applyLight(Light light, Material material, vec3 P, vec3 N, vec3 V)
{
    vec3 texColor = vec3(1, 1, 1);    
    vec3 ambient =  material.Ka * light.color;
    vec3 diffuse = vec3(0);
    vec3 specular = vec3(0);

    vec3 L = normalize(light.position - P);
       
    if(material.hasTex == 1)
    {
        texColor = texture(tex, vec2(VertTexture.s, 1-VertTexture.t)).rgb;
		//texColor = vec3(0.5, 0.0, 0.5);
    }    

    float d = max(dot(N, L), 0);        
    float s = 0;

    if(d >= 0)
    {
        vec3 H = normalize(L+V);
        s = pow(max(dot(N, H), 0), material.Ns);
        specular = material.Ks * light.color * s * 1;
    }

    vec3 selecColor = vec3(1);
    if(isSelected == 1)
    {
        selecColor = vec3(1, 0, 0);
    }
	
    diffuse  = material.Kd * selecColor * texColor.rgb * light.color * d;    

    float t = dualConeSpotlight(P, light.position, -light.direction, 21.0, 32.0);
    float a = lightAttenuation(P, light.position.xyz-P, light.attenuation, light.intensity.x);

    vec3 color = ambient + (diffuse + specular) * a;

    return color;
}

void main()
{
    vec4 color = vec4(0, 0, 0, 1);

    float shadow = 1.0;

    vec3 P = VertPosition.xyz;
    vec3 N = normalize(VertNormal.xyz);
    vec3 V = normalize(camPos.xyz-P);
   
    if(applyShadow == 0)    
    {
        for(int i=0; i<numLights; i++)
        {            
            shadow *= VSM(VertShadowCoord[i], shadowMap[i]); 
        }   
    }

    for(int i=0; i<numLights; i++)
    {
        color.rgb += applyLight(lights[i], material, P, N, V) * shadow;
    }   
    color.a = alpha;

	if(colorMode == 0)
	{
		FragColor = vec4(color.xyz * 1.0, 1);	
	}
	else
	{
		vec3 L = normalize(lights[0].position - P);
		float ndotl = max(0, dot(N, L));
		
		FragColor = vec4(VertColor.rgb * ndotl, 1);
		FragColor = vec4(vec3(0.5, 0.5, 0.5) * ndotl, 1);
		//FragColor = vec4(vec3(0.3, 0.1, 0.1) * ndotl, 1);
	}

	//FragColor = vec4(1, 0, 0, 1);
	//FragColor = vec4(VertTexture.st, 0, 1);
}
