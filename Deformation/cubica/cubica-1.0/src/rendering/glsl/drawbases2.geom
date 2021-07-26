
varying in vec4 n1[];
varying in vec4 n2[];
varying in vec4 n3[];
varying in vec4 n4[];

void main()
{  
  ///////////////////////////////////////////
  gl_Position = gl_PositionIn[0];
  gl_FrontColor = vec4(0.25, 0.5, 0.25, 1.0);
  //gl_FrontColor = vec4(0.25, 0.0, 0.0, 1.0);
  EmitVertex();
  gl_Position = (gl_PositionIn[0] + n1[0]) / 2.0;
  gl_FrontColor = vec4(0.5, 0.5, 0.5, 1.0);
  //gl_FrontColor = vec4(1.0, 0.0, 0.0, 1.0);
  EmitVertex();
  EndPrimitive();

  gl_Position = gl_PositionIn[0];
  gl_FrontColor = vec4(0.5, 0.5, 0.25, 1.0);
  //gl_FrontColor = vec4(0.0, 0.25, 0.0, 1.0);
  EmitVertex();
  gl_Position = (gl_PositionIn[0] + n2[0]) / 2.0;
  gl_FrontColor = vec4(0.5, 0.5, 0.5, 1.0);
  //gl_FrontColor = vec4(0.0, 1.0, 0.0, 1.0);
  EmitVertex();
  EndPrimitive();

  gl_Position = gl_PositionIn[0];
  gl_FrontColor = vec4(0.5, 0.25, 0.5, 1.0);
  //gl_FrontColor = vec4(0.05, 0.05, 0.3, 1.0);
  EmitVertex();
  gl_Position = (gl_PositionIn[0] + n3[0]) / 2.0;
  gl_FrontColor = vec4(0.5, 0.5, 0.5, 1.0);
  //gl_FrontColor = vec4(0.1, 0.1, 1.0, 1.0);
  EmitVertex();
  EndPrimitive();

  gl_Position = gl_PositionIn[0];
  gl_FrontColor = vec4(0.25, 0.5, 0.5, 1.0);
  //gl_FrontColor = vec4(0.25, 0.25, 0.0, 1.0);
  EmitVertex();
  gl_Position = (gl_PositionIn[0] + n4[0]) / 2.0;
  gl_FrontColor = vec4(0.5, 0.5, 0.5, 1.0);
  //gl_FrontColor = vec4(1.0, 1.0, 0.0, 1.0);
  EmitVertex();
  EndPrimitive();
}

