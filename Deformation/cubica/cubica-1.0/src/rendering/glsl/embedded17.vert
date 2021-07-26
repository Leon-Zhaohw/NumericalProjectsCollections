/////////////////////////////////////////////////
// OpenGL vertex shader for rendering
// deformed-embedded vertices and normals 
// from the texture into the 3d space.
// Used in: _meshRenderWithShadowShader
//          _meshRenderNoShadowShader
/////////////////////////////////////////////////

uniform sampler2D embeddedVertexTex;
uniform sampler2D embeddedNormalTex;

varying vec3 Light;
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

  if (gl_LightSource[0].position.w == 1.0)
  {
    Light = gl_LightSource[0].position.xyz - ePos.xyz;
  } else {
    Light = -gl_LightSource[0].position.xyz;
  }
  if (dot(Light,Light) > 0.0)
    Light = normalize(Light);
  
  Normal = gl_NormalMatrix * Normal;
  gl_Position = gl_ModelViewProjectionMatrix * Pos;
}

