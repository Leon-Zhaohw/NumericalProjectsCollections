/////////////////////////////////////////////////
// OpenGL vertex shader for rendering
// the obj with current color 
/////////////////////////////////////////////////

void main()
{
  vec4 Pos = gl_Vertex;
  if (Pos.w > 0.0)
    Pos = Pos / Pos.w;

  Pos = gl_ModelViewProjectionMatrix * Pos;

  gl_FrontColor = gl_Color;
  gl_Position = Pos;
 
}

