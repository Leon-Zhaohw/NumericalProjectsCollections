/////////////////////////////////////////////////
// Very basic OpenGL vertex shader
/////////////////////////////////////////////////

varying float Depth;

void main()
{
	// Transform from object space to eye space
	vec4 P = gl_ModelViewMatrix * gl_Vertex;

  // Do projection to get the linear eye depth
	vec2 EyePos = (gl_ProjectionMatrix * P).zw;
  Depth = EyePos.x / EyePos.y;
  Depth *= 100.0 * abs(EyePos.x);

	gl_Position = ftransform();
}
