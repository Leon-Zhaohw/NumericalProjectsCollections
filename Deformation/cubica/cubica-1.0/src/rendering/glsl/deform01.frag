/////////////////////////////////////////////////
// Very basic OpenGL fragment shader
/////////////////////////////////////////////////

//uniform sampler2D Umatrix;
//uniform vec4 q[5];
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
  //vec3 col = vec3(1.0 / (10.0 * length(norm)));

  vec3 col = gl_Color.rgb;
  
  //gl_FragColor = vec4(vec3(Depth), 1.0);
  gl_FragColor = vec4(col, 1.0);
}
