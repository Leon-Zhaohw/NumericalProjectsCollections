/////////////////////////////////////////////////
// OpenGL fragment shader for rendering
// deformed-embedded vertices and normals 
// from the texture into the 3d space.
// Used in: _meshRenderNoShadowShader
/////////////////////////////////////////////////

#define TOTAL_LIGHTS 5

// Vector from pixel position to light (should be ~unit length)
varying vec3 Light[TOTAL_LIGHTS];
// Normal vector of pixel position (should be ~unit length)
varying vec3 Normal;
// Vector from 'eye' to pixel position (not unit length)
varying vec4 Pos;

varying vec3 eyeVec;

void main()
{
  // Normalize these (somewhat redundant)
  vec3 N = normalize(Normal);

  vec3 colorTotal = gl_FrontLightModelProduct.sceneColor;

  int i;
  for ( i = 0; i < TOTAL_LIGHTS; i++ )
  {
    vec3 L = normalize(Light[i]);

    // cosine of the angle between the normal and the light vectors
    float NdotL = max(dot(N, L), 0.0);

    float lambertTerm = dot(N, L);

    if ( lambertTerm > 0.0 )
    {
      colorTotal += gl_LightSource[i].diffuse *
                    gl_FrontMaterial.diffuse *
                    lambertTerm;

      vec3 E = normalize(eyeVec);
      vec3 R = reflect(-L, N);

      float specular = pow(max(dot(R, E), 0.0), gl_FrontMaterial.shininess);

      colorTotal += gl_LightSource[i].specular *
                    gl_FrontMaterial.specular *
                    specular;
    }
  }
 
  // Output final color  
  gl_FragColor = vec4(colorTotal, 1.0);
}

