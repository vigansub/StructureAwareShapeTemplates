#version 400 core

in vec4 VertPosition;
in vec4 VertNormal;
in vec4 VertColor;
in vec4 VertTexture;

uniform vec4 color;
uniform vec3 lightPos;
uniform float transparency;
uniform int colorMode;

void main()
{
   vec4 color = color;

   float ndotl = dot(normalize(VertNormal.xyz), normalize(lightPos - VertPosition.xyz));

   color = vec4(1, 1, 1, 1);

   if(colorMode == 0)
	 color = vec4(1.0f, 0.6f, 0.3f, 1.0f);
   if(colorMode == 1)
	 color = vec4(0.3f, 0.3f, 1.0f, 1.0f);
   if(colorMode == 2)
	 color = vec4(1.0f, 1.0f, 0.3f, 1.0f);
   if(colorMode == 3)
	 color = vec4(1.0f, 0.3f, 0.3f, 1.0f);
   if(colorMode == 3)
	 color = vec4(0.3f, 1.0f, 0.0f, 1.0f);
   
   color.xyz *= max(0.3, ndotl);
   
   gl_FragColor = vec4(color.xyz, transparency);	
}

/*
case 0:
    color = vec4(1.0f, 0.5f, 0.0f, 1.0f);
    break;
case 1:
    color = vec4(0.0f, 0.0f, 1.0f, 1.0f);
    break;
case 2:
    color = vec4(0.0f, 1.0f, 0.0f, 1.0f);
    break;
case 4:
    color = vec4(1.0f, 0.3f, 0.3f, 1.0f);
    break;
case 3:
    color = vec4(1.0f, 0.5f, 0.0f, 1.0f);
    break;
case 5:
    color = vec4(0.0f, 0.5f, 1.0f, 1.0f);
    break;
case 6:
    color = vec4(0.5f, 0.5f, 1.0f, 1.0f);
    break;
default:
    break;
	*/