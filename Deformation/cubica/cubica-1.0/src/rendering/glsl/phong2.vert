/////////////////////////////////////////////////
// OpenGL vertex shader for rendering
// deformed-embedded vertices and normals 
// from the texture into the 3d space.
// Used in: _meshRenderWithShadowShader
//          _meshRenderNoShadowShader
/////////////////////////////////////////////////

#define TOTAL_LIGHTS 5

varying vec3 Light[TOTAL_LIGHTS];
varying vec3 Normal;
varying vec4 Pos;
varying vec3 eyeVec;

void main()
{
  gl_Position = ftransform();
  Pos = gl_Vertex;

  Normal = gl_Normal;

  vec4 ePos = gl_ModelViewMatrix * Pos;

  eyeVec = -ePos.xyz;

  int i;
  for ( i = 0; i < TOTAL_LIGHTS; i++ )
  {
    if (gl_LightSource[i].position.w == 1.0)
    {
      Light[i] = gl_LightSource[i].position.xyz - ePos.xyz;
    } else {
      Light[i] = -gl_LightSource[i].position.xyz;
    }
    if (dot(Light[i],Light[i]) > 0.0)
      Light[i] = normalize(Light[i]);
  }
  
  Normal = gl_NormalMatrix * Normal;
}

