/////////////////////////////////////////////////
// 
/////////////////////////////////////////////////

varying vec3 Rotation1;
varying vec3 Rotation2;
varying vec3 Rotation3;

void main()
{
  // Store our calculated n1 and n2 rest frame vectors in the two render targets.
  gl_FragData[0].rgb = Rotation1.xyz;
  gl_FragData[1].rgb = Rotation2.xyz;
  gl_FragData[2].rgb = Rotation3.xyz;
}
