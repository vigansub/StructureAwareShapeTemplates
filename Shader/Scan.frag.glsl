//Author: Sören Pirk
//Date: 22.01.2013

#version 400 core

in vec4 VertPosition;
in vec4 VertNormal;
in vec4 VertColor;
in vec4 VertTexture;
in vec4 MVPPosition;

uniform vec3 lightPos;

void main()
{
   vec4 color = VertColor;
   gl_FragData[0] = vec4(MVPPosition.xyz, 1);   
   gl_FragData[1] = vec4(VertNormal.xyz, 1);	
}
