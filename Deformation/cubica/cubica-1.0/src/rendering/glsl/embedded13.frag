/////////////////////////////////////////////////
// OpenGL fragment shader for rendering
// embedded tet mesh with NO SHADOWS
/////////////////////////////////////////////////

varying vec3 oLight;
varying vec3 oNormal;

void main()
{ 
  //normalize(Light);
  vec3 N = normalize(oNormal);
  vec3 L = normalize(oLight);
  float NdotL = max(dot(N, L), 0.0);
  NdotL = (NdotL + 0.25) * 0.8;
 
  //gl_FragColor = vec4(vec3(NdotL), 1.0);
  gl_FragColor = gl_Color;
  //gl_FragColor = vec4(abs(N), 1.0);
}
