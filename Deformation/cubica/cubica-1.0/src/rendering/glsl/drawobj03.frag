/////////////////////////////////////////////////
// OpenGL fragment shader for rendering
// obj into deferred render targets.
// Used in: _meshRenderDeferredShader
/////////////////////////////////////////////////

varying vec2 Depth;
varying vec3 Normal;
varying vec4 Pos;

void main()
{
  // Bring linear depth values from into (0, 1) range
	float frag = (0.5 * Depth.x / Depth.y) + 0.5;

  vec3 norm = Normal;
  if (dot(norm,norm) > 0.0)
    norm = normalize(norm);

  gl_FragData[0] = vec4(norm, frag);
  gl_FragData[1] = Pos;



}


