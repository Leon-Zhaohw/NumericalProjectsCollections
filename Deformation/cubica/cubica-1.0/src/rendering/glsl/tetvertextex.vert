/////////////////////////////////////////////////
// OpenGL vertex shader 
/////////////////////////////////////////////////

uniform sampler2D tetMeshTex;
uniform sampler2D tetNormalTex;

void main()
{
  // vertex position stores barycentric coords inside a tetrahedral
  vec2 coords = gl_Vertex.xy;

  // Grab positions of four tetMesh vertices
  vec4 pos = texture2D(tetMeshTex, coords);
  vec3 col = texture2D(tetNormalTex, coords).rgb;
  if (length(col) > 0.0)
    col = abs(normalize(col));
  else
    col = abs(col);

  // Color and output position
  //gl_FrontColor = vec4(col, 1.0);
  gl_FrontColor = vec4(0,0,0,1);
  gl_Position = gl_ModelViewProjectionMatrix * pos;

}

