
#define BASIS_LENGTH 0.05

varying in vec4 r1[];
varying in vec4 r2[];
varying in vec4 r3[];

varying in vec3 c1[];
varying in vec3 c2[];
varying in vec3 c3[];


void main()
{  
  gl_Position = gl_PositionIn[0];
  gl_FrontColor = vec4(c1[0], 1.0);
  EmitVertex();
  gl_Position = gl_PositionIn[0] + BASIS_LENGTH * r1[0];
  gl_FrontColor = vec4(c1[0], 1.0);
  EmitVertex();
  EndPrimitive();

  gl_Position = gl_PositionIn[0];
  gl_FrontColor = vec4(c2[0], 1.0);
  EmitVertex();
  gl_Position = gl_PositionIn[0] + BASIS_LENGTH * r2[0];
  gl_FrontColor = vec4(c2[0], 1.0);
  EmitVertex();
  EndPrimitive();

  gl_Position = gl_PositionIn[0];
  gl_FrontColor = vec4(c3[0], 1.0);
  EmitVertex();
  gl_Position = gl_PositionIn[0] + BASIS_LENGTH * r3[0];
  gl_FrontColor = vec4(c3[0], 1.0);
  EmitVertex();
  EndPrimitive();

  ///////////////////////////////////////////
}

