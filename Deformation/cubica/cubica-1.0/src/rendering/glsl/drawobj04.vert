/////////////////////////////////////////////////
// Vertex shader for ambient occlusion volumes,
// for an obj.
/////////////////////////////////////////////////

uniform float delta;

varying float vDelta;

void main()
{
  vec4 Pos = gl_Vertex; 
  if (Pos.w > 0.0)
    Pos = Pos / Pos.w;

  //vDelta = gl_Vertex.z;
  vDelta = delta;
  //vDelta = 0.5;

  gl_Position = gl_ModelViewMatrix * Pos;
}

