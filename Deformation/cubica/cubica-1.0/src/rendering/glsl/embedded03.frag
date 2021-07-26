/////////////////////////////////////////////////
// OpenGL fragment shader for rendering
// embedded tet mesh with NO SHADOWS
/////////////////////////////////////////////////

varying vec3 Light;
varying vec3 Normal;

void main()
{ 
  //normalize(Light);
  vec3 N = normalize(Normal);
  vec3 L = normalize(Light);
  float NdotL = max(dot(N, L), 0.0);
  NdotL = (NdotL * 0.75) + 0.25;
 
  gl_FragColor = vec4(vec3(NdotL), 1.0);
  //gl_FragColor = vec4(vec3(1.0), 1.0);
  //gl_FragColor = vec4(abs(N), 1.0);
}
