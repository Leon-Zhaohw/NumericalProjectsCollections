/////////////////////////////////////////////////
// OpenGL fragment shader for rendering
// embedded tet mesh INTO A SHADOWMAP
/////////////////////////////////////////////////

// Two input variables used to modify the pseudo polygon offset
//uniform float factor;
//uniform float units;

varying vec2 Depth;
//varying float Epsilon;
#define EXP_FACTOR_K 80.0

void main()
{ 
  // Bring linear depth values from into (0, 1) range
	float frag = (0.5 * Depth.x / Depth.y) + 0.5;

	// The moment calculation is used with variance shadowmapping
	vec2 Moments;
	Moments.x = frag;
	float dx = dFdx(Moments.x);
	float dy = dFdy(Moments.x);
	
	Moments.y = Moments.x * Moments.x + 0.25 * (dx*dx + dy*dy);	

	// Calculate pseudo polygon offset
	//float o = frag * (factor - 1.0) + 0.01 * units;

  //float depth = frag;// - Epsilon * o;
  
  // calculated fragment depth + offset
  gl_FragColor.r = exp2(EXP_FACTOR_K * frag);
  gl_FragColor.g = Moments.x;
  gl_FragColor.b = Moments.y;

  //gl_FragColor = vec4(vec3(frag), 1.0);
  //gl_FragColor = vec4(0.6, 0.2, 0.4, 1.0);
}

