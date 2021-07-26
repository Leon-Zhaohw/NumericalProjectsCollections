/////////////////////////////////////////////////
// 
/////////////////////////////////////////////////

varying vec3 n1Frame;
varying vec3 n2Frame;

void main()
{
  // Store our calculated n1 and n2 rest frame vectors in the two render targets.
  gl_FragData[0].rgb = n1Frame.xyz;
  gl_FragData[1].rgb = n2Frame.xyz;
}
