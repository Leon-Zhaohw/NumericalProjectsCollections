/////////////////////////////////////////////////
// Very basic OpenGL fragment shader
/////////////////////////////////////////////////

varying vec3 Light;
varying vec3 Normal;

void main()
{ 
  //normalize(Light);
  vec3 norm = normalize(Normal);
  float NdotL = max(dot(norm, Light), 0.0);
  NdotL = (NdotL + 0.25) * 0.8;
 
  gl_FragColor = vec4(vec3(NdotL), 1.0);
  //gl_FragColor = vec4(abs(norm), 1.0);
}
