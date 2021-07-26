/////////////////////////////////////////////////
// OpenGL vertex shader for rendering
// obj into deferred render targets. 
/////////////////////////////////////////////////

varying vec2 Depth;
varying vec3 Normal;
varying vec4 Pos;

void main()
{
  Pos = gl_Vertex; 
  if (Pos.w > 0.0)
    Pos = Pos / Pos.w;

  Normal = gl_Normal;
  if (dot(Normal,Normal) > 0.0)
    Normal = normalize(Normal);
  Normal = gl_NormalMatrix * Normal;

  vec4 ePos = gl_ModelViewMatrix * Pos;

  Depth = (gl_ProjectionMatrix * ePos).zw;
  Pos.xyz = ePos.xyz;
  gl_Position = gl_ProjectionMatrix * ePos;
 
}

