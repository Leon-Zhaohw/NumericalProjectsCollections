/////////////////////////////////////////////////
// Fragment shader for calculating deformed normals and tangents
/////////////////////////////////////////////////

varying vec3 Normal;
varying vec3 Tangent;

void main()
{
  gl_FragData[0].rgb = Normal;
  gl_FragData[1].rgb = Tangent;
  //gl_FragData[1].rgb = vec3(0.0, 0.1, 0.0);
}

