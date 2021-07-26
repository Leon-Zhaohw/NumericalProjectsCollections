
#define BASIS_LENGTH 0.02

varying in vec3 n1Rot[];
varying in vec3 n2Rot[];
varying in vec3 n3Rot[];

varying in vec3 n1Rest[];
varying in vec3 n2Rest[];
varying in vec3 n3Rest[];

void main()
{  
  gl_Position = gl_PositionIn[0];
  gl_FrontColor = vec4(0.5, 0.0, 0.0, 1.0);
  EmitVertex();
  gl_Position = gl_PositionIn[0] + BASIS_LENGTH * vec4(n1Rest[0], 0.0);
  gl_FrontColor = vec4(0.5, 0.0, 0.0, 1.0);
  EmitVertex();
  gl_Position = gl_PositionIn[0] + BASIS_LENGTH * vec4(n1Rot[0], 0.0);
  gl_FrontColor = vec4(1.0, 0.0, 0.0, 1.0);
  EmitVertex();
  EndPrimitive();

  gl_Position = gl_PositionIn[0];
  gl_FrontColor = vec4(0.0, 0.5, 0.0, 1.0);
  EmitVertex();
  gl_Position = gl_PositionIn[0] + BASIS_LENGTH * vec4(n2Rest[0], 0.0);
  gl_FrontColor = vec4(0.0, 0.5, 0.0, 1.0);
  EmitVertex();
  gl_Position = gl_PositionIn[0] + BASIS_LENGTH * vec4(n2Rot[0], 0.0);
  gl_FrontColor = vec4(0.0, 1.0, 0.0, 1.0);
  EmitVertex();
  EndPrimitive();

  gl_Position = gl_PositionIn[0];
  gl_FrontColor = vec4(0.0, 0.0, 0.5, 1.0);
  EmitVertex();
  gl_Position = gl_PositionIn[0] + BASIS_LENGTH * vec4(n3Rest[0], 0.0);
  gl_FrontColor = vec4(0.0, 0.0, 0.5, 1.0);
  EmitVertex();
  gl_Position = gl_PositionIn[0] + BASIS_LENGTH * vec4(n3Rot[0], 0.0);
  gl_FrontColor = vec4(0.0, 0.0, 1.0, 1.0);
  EmitVertex();
  EndPrimitive();// */

  ///////////////////////////////////////////
/*  gl_Position = gl_PositionIn[0];
  gl_FrontColor = vec4(1.0, 0.0, 0.0, 1.0);
  EmitVertex();
  gl_Position = gl_PositionIn[0] + BASIS_LENGTH * vec4(n1Rot[0], 0.0);
  gl_FrontColor = vec4(1.0, 0.0, 0.0, 1.0);
  EmitVertex();
  EndPrimitive();

  gl_Position = gl_PositionIn[0];
  gl_FrontColor = vec4(0.0, 1.0, 0.0, 1.0);
  EmitVertex();
  gl_Position = gl_PositionIn[0] + BASIS_LENGTH * vec4(n2Rot[0], 0.0);
  gl_FrontColor = vec4(0.0, 1.0, 0.0, 1.0);
  EmitVertex();
  EndPrimitive();

  gl_Position = gl_PositionIn[0];
  gl_FrontColor = vec4(0.0, 0.0, 1.0, 1.0);
  EmitVertex();
  gl_Position = gl_PositionIn[0] + BASIS_LENGTH * vec4(n3Rot[0], 0.0);
  gl_FrontColor = vec4(0.0, 0.0, 1.0, 1.0);
  EmitVertex();
  EndPrimitive();// */
}

