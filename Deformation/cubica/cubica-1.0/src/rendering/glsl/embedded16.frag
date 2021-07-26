/////////////////////////////////////////////////
// OpenGL fragment shader for rendering
// deformed-embedded vertices and normals into
// the embedded vertex textures.
// Expects additive blending to be on.

// Used in: _embeddedMeshToTextureShader
/////////////////////////////////////////////////


varying vec4 fPos;
varying vec3 fNormal;
varying float fDelta;

void main()
{
  gl_FragData[0] = fPos;
  gl_FragData[1] = vec4(fNormal, fDelta);

}

