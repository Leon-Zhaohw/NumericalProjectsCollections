
varying in vec4 Point0[];
varying in vec4 Point1[];
varying in vec4 Point2[];

void main()
{
  gl_Position = Point0[0];
  gl_FrontColor = vec4(0.4);
  EmitVertex();
  gl_Position = Point1[0];
  gl_FrontColor = vec4(0.4);
  EmitVertex();
  gl_Position = Point2[0];
  gl_FrontColor = vec4(0.4);
  EmitVertex();
  EndPrimitive();
}
