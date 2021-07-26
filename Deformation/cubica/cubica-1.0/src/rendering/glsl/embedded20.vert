/////////////////////////////////////////////////
// OpenGL vertex shader for rendering
// deformed-embedded vertices and normals 
// from the texture into the 3d space.
// Used in: _meshRenderBasicShader
/////////////////////////////////////////////////

uniform sampler2D embeddedVertexTex;

void main()
{
  vec4 Pos = texture2D(embeddedVertexTex, gl_Vertex.xy);
  if (Pos.w > 0.0)
    Pos = Pos / Pos.w;

  Pos = gl_ModelViewProjectionMatrix * Pos;

  gl_FrontColor = gl_Color;
  gl_Position = Pos;
 
}

