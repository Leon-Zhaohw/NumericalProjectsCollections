//////
// Shader to draw the ground plane to a shadow map
/////

varying vec2 vDepth;
varying vec3 vPos;

void main()
{ 
  vec4 eyepos = gl_ModelViewMatrix * gl_Vertex;
  vDepth = (gl_ProjectionMatrix * eyepos).zw;
  vPos = eyepos.xyz;

	gl_Position = ftransform();
}



