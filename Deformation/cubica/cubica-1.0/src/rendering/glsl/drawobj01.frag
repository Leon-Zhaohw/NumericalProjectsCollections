/////////////////////////////////////////////////
// OpenGL fragment shader for rendering
// basic OBJ with no shadows
/////////////////////////////////////////////////

// Vector from pixel position to light (should be ~unit length)
varying vec3 Light;
// Normal vector of pixel position (should be ~unit length)
varying vec3 Normal;
// Vector from 'eye' to pixel position (not unit length)
varying vec4 Pos;

void main()
{
  // Normalize these (somewhat redundant)
  vec3 N = normalize(Normal);
  vec3 L = normalize(Light);
  vec3 P = normalize(Pos.xyz);

  // Half-Vector between eye and light vectors (for specular calculation)
  vec3 halfVector = normalize(L + P);

  // cosine of the angle between the normal and the light vectors
  float NdotL = max(dot(N, L), 0.0);
  // cosine of the angle between the normal and the half-vector
  float NdotHV = max(dot(N, halfVector), 0.0);
  // Power factor
  float pf = 0.0;
  // Shininess, used for the specular exponent. 
  float shininess = max(gl_FrontMaterial.shininess, 10.0);

  // Specular doesn't show up if the surface is completely in shadow
  if (NdotL > 0.0)
    pf = pow(NdotHV, shininess);

  // Add all light and surface color contributions
  vec3 color = gl_LightSource[0].ambient.rgb * (1.0 - NdotL) * gl_FrontMaterial.ambient.rgb;
  color += gl_LightSource[0].diffuse.rgb * NdotL * gl_FrontMaterial.diffuse.rgb;
  color += gl_LightSource[0].specular.rgb * pf * gl_FrontMaterial.specular.rgb;
 
  // Output final color  
  gl_FragColor = vec4(color, 1.0);
}

