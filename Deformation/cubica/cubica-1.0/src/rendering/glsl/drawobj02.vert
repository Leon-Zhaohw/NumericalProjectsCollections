/////////////////////////////////////////////////
// OpenGL vertex shader for rendering
// the obj into 3d space. 
/////////////////////////////////////////////////

varying vec2 Depth;

void main()
{
  vec4 Pos = gl_Vertex;
  if (Pos.w > 0.0)
    Pos = Pos / Pos.w;

  Pos = gl_ModelViewProjectionMatrix * Pos;

  Depth = Pos.zw;
  gl_Position = Pos;
 
}

