/////////////////////////////////////////////////
// Very basic OpenGL fragment shader
/////////////////////////////////////////////////

void main()
{
  
  //gl_FragColor = vec4(vec3(Depth), 1.0);
  gl_FragColor = gl_Color;
}
