/////////////////////////////////////////////////
// OpenGL fragment shader for deferred rendering
// face normal, per-pixel linear depth, and world position
/////////////////////////////////////////////////

varying vec2 fDepth;
varying vec3 fPos;
varying vec3 fNormal;

#define APPROX_Z 0.00143

/*void main()
{ 
  // Bring linear depth values from into (0, 1) range
	float frag = (0.5 * Depth.x / Depth.y) + 0.5;

  float dx = dFdx(Depth.x);
  float dy = dFdy(Depth.x);
  vec3 norm = vec3(dx,dy, APPROX_Z);

  norm = normalize(norm);
  //norm = 0.5 * (norm + vec3(1.0));
  

  gl_FragData[0] = vec4(norm, frag);
  gl_FragData[1] = vec4(Pos, 1.0);
  //gl_FragColor = vec4(vec3(frag), 1.0);
  //gl_FragColor = vec4(0.6, 0.2, 0.4, 1.0);
}*/


void main()
{ 
  // Bring linear depth values from into (0, 1) range
	float frag = (0.5 * fDepth.x / fDepth.y) + 0.5;

  vec3 norm = fNormal;
  if (dot(norm,norm) > 0.0)
    norm = normalize(norm);

  //norm = 0.5 * (norm + vec3(1.0));

  gl_FragData[0] = vec4(norm, frag);
  gl_FragData[1] = vec4(fPos, 1.0);
}

