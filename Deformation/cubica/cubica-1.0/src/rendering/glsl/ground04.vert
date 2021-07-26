
varying vec4 Origin;
varying vec3 Down;

void main()
{
  //vec4 Pos = gl_Vertex;
  //Pos.z += 7.0 * delta;

  Origin = gl_ModelViewMatrix * vec4(0.0,0.0,0.15,1.0);
  Down = gl_NormalMatrix * vec3(0.0, 0.0,-1.0);
  
  //gl_Position = gl_ModelViewProjectionMatrix * Pos;
  gl_Position = gl_Vertex + vec4(0.0,0.0,-0.5,0.0);
}

