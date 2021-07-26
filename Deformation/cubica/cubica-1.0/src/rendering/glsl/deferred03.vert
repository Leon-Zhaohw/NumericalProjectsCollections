/////////////////////////////////////////////////
// OpenGL vertex shader for rendering
// deformed-embedded vertices and normals 
// from the texture into the 3d space.
// Used in: _meshRenderDeferredShader
/////////////////////////////////////////////////

uniform sampler2D embeddedVertexTex;
uniform sampler2D embeddedNormalTex;

varying vec2 Depth;
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
  Normal = gl_NormalMatrix * Normal;

  vec4 ePos = gl_ModelViewMatrix * Pos;

  Depth = (gl_ProjectionMatrix * ePos).zw;
  Pos.xyz = ePos.xyz;
  gl_Position = gl_ProjectionMatrix * ePos;
 
}

