//////
// Shader to draw the ground plane to a shadow map
/////

varying vec2 Depth;

void main()
{ 
  vec4 pos = gl_ModelViewProjectionMatrix * gl_Vertex;
  Depth = pos.zw;

	gl_Position = ftransform();
}



