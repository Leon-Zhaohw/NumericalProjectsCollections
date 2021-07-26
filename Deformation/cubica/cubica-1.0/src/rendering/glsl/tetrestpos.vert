/////////////////////////////////////////////////
// Very simple shader to take a vertex position
// and a texture coordinate, and write the position
// into the texture at the coord.
/////////////////////////////////////////////////

void main()
{
  vec4 screenPos = vec4(gl_MultiTexCoord0.st, 0.0, 1.0);

  gl_FrontColor = gl_Vertex;
  gl_Position = gl_ModelViewProjectionMatrix * screenPos;
}

