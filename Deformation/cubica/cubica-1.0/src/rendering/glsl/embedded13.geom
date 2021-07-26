#version 120 
#extension GL_EXT_geometry_shader4 : enable

#define NORMAL_LENGTH 0.02

varying in vec3 Light[];
varying in vec3 Normal[];

varying out vec3 oNormal;
varying out vec3 oLight;

void main()
{
  gl_Position = gl_PositionIn[0] - NORMAL_LENGTH * vec4(Normal[0], 0.0);
  gl_FrontColor = gl_FrontColorIn[0];
  oNormal = Normal[0];
  oLight = Light[0];
  EmitVertex();

  gl_Position = gl_PositionIn[0] + NORMAL_LENGTH * vec4(Normal[0], 0.0);
  gl_FrontColor = gl_FrontColorIn[0];
  oNormal = Normal[0];
  oLight = Light[0];
  EmitVertex();

  EndPrimitive();
  
  /////////////////////////////////////////////////////////////
}

