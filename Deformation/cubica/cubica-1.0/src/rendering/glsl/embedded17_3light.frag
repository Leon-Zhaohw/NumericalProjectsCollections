/////////////////////////////////////////////////
// OpenGL fragment shader for rendering
// deformed-embedded vertices and normals 
// from the texture into the 3d space.
// Used in: _meshRenderNoShadowShader
/////////////////////////////////////////////////

// This is added in GLSL_SCENE_VIEWER::init(), when loading the shader
//#define TOTAL_LIGHTS 5

// Vector from pixel position to light (should be ~unit length)
varying vec3 Light[TOTAL_LIGHTS];
// Normal vector of pixel position (should be ~unit length)
varying vec3 Normal;
// Vector from 'eye' to pixel position (not unit length)
varying vec4 Pos;

void main()
{
  // Normalize these (somewhat redundant)
  vec3 N = normalize(Normal);
  vec3 P = normalize(Pos.xyz);

  vec3 colorTotal;

  int i;
  for ( i = 0; i < TOTAL_LIGHTS; i++ )
  {
    vec3 L = normalize(Light[i]);
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
    vec3 color = gl_LightSource[i].ambient.rgb * (1.0 - NdotL) * gl_FrontMaterial.ambient.rgb;
    color += gl_LightSource[i].diffuse.rgb * NdotL * gl_FrontMaterial.diffuse.rgb;
    color += gl_LightSource[i].specular.rgb * pf * gl_FrontMaterial.specular.rgb;

    colorTotal += color;
  }
 
  // Output final color  
  gl_FragColor = vec4(colorTotal, 1.0);
}

