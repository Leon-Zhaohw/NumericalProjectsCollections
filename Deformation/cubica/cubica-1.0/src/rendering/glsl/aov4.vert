/////////////////////////////////////////////////
// Vertex shader for ambient occlusion volumes.
/////////////////////////////////////////////////

uniform sampler2D embeddedVertexTex;
uniform float delta;

varying float vDelta;

void main()
{
  vec4 Pos = texture2D(embeddedVertexTex, gl_Vertex.xy);
  if (Pos.w > 0.0)
    Pos = Pos / Pos.w;

  //vDelta = gl_Vertex.z;
  vDelta = delta;
  //vDelta = 0.5;

  gl_Position = gl_ModelViewMatrix * Pos;
}

