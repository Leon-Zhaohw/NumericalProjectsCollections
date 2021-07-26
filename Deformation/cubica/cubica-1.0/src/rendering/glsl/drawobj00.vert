/////////////////////////////////////////////////
// OpenGL vertex shader for rendering
// basic OBJ vertices and normals 
/////////////////////////////////////////////////

varying vec3 Light;
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

