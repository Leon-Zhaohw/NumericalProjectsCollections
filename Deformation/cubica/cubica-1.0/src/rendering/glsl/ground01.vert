//////
// Shader to draw the ground plane with shadow mapping
/////

varying vec3 Light;
varying vec3 Normal;
varying vec4 Pos;

void main()
{
  Normal = normalize(gl_NormalMatrix * normalize(gl_Normal));
  Pos = vec4(gl_Vertex.xyz, 1.0);
 
  vec4 ePos = gl_ModelViewMatrix * gl_Vertex;

  if (gl_LightSource[0].position.w == 1.0)
  {
    Light = gl_LightSource[0].position.xyz - ePos.xyz;
  } 
  else
  {
    Light = -gl_LightSource[0].position.xyz;
  }
  Light = normalize(Light);

//	gl_TexCoord[0] = gl_MultiTexCoord0;
	gl_Position = ftransform();
}


