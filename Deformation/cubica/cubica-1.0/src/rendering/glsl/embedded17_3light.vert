/////////////////////////////////////////////////
// OpenGL vertex shader for rendering
// deformed-embedded vertices and normals 
// from the texture into the 3d space.
// Used in: _meshRenderWithShadowShader
//          _meshRenderNoShadowShader
/////////////////////////////////////////////////

uniform sampler2D embeddedVertexTex;
uniform sampler2D embeddedNormalTex;

// This is added in GLSL_SCENE_VIEWER::init(), when loading the shader
//#define TOTAL_LIGHTS 5

varying vec3 Light[TOTAL_LIGHTS];
varying vec3 Normal;
varying vec4 Pos;

void main()
{
  Pos = texture2D(embeddedVertexTex, gl_Vertex.xy);
  if (Pos.w > 0.0)
    Pos = Pos / Pos.w;

  Normal = texture2D(embeddedNormalTex, gl_Vertex.xy).rgb;
  if (dot(Normal,Normal) > 0.0)
    Normal = normalize(Normal);

  vec4 ePos = gl_ModelViewMatrix * Pos;

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
  gl_Position = gl_ModelViewProjectionMatrix * Pos;
}

