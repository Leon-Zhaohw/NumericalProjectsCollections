/////////////////////////////////////////////////
// OpenGL fragment shader for rendering
// obj INTO A SHADOWMAP
/////////////////////////////////////////////////

varying vec2 Depth;
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

  gl_FragColor.r = exp2(EXP_FACTOR_K * frag);
  gl_FragColor.g = Moments.x;
  gl_FragColor.b = Moments.y;

}

