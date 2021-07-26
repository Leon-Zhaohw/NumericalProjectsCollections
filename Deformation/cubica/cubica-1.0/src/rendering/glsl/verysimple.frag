/////////////////////////////////////////////////
// Very basic OpenGL fragment shader
/////////////////////////////////////////////////

varying float Depth;

void main()
{
  vec2 norm;
  norm.x = dFdx(Depth);
  norm.y = dFdy(Depth);
 
  //normalize(norm);
  //norm += vec2(1.0);
  //norm += vec2(1.0);
  norm = abs(norm);
  vec3 col = vec3(1.0 / (10.0 * length(norm)));
  
  //gl_FragColor = vec4(vec3(Depth), 1.0);
  gl_FragColor = vec4(col, 1.0);
}
