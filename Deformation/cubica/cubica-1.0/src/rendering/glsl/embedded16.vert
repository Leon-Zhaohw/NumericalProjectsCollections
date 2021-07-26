/////////////////////////////////////////////////
// OpenGL vertex shader for rendering
// deformed-embedded vertices and normals into
// the embedded vertex textures

// Used in: _embeddedMeshToTextureShader
/////////////////////////////////////////////////

uniform sampler2D tetMeshTex;

varying vec3 vPos;
varying float vDelta;

void main()
{
  // gl_Vertex  stores barycentric coords inside a tetrahedral
  vec3 coords = gl_Vertex.xyz;

  // Grab positions of four tetMesh vertices
  vec3 v0 = texture2D(tetMeshTex, gl_MultiTexCoord0.st).rgb;
  vec3 v1 = texture2D(tetMeshTex, gl_MultiTexCoord0.pq).rgb;
  vec3 v2 = texture2D(tetMeshTex, gl_MultiTexCoord1.st).rgb;
  vec3 v3 = texture2D(tetMeshTex, gl_MultiTexCoord1.pq).rgb;

  // Construct the original position from the deformed verts
  float subtract = 1.0 - coords.x - coords.y - coords.z;
  vPos = subtract * v3 + v0 * coords.x + v1 * coords.y + v2 * coords.z;

  // Per vertex variable delta stored in third component.
  vDelta = gl_MultiTexCoord2.z;
 
  // Position is texture coordinate of output point. Expects viewport to be set to NDC.
  gl_Position = gl_ProjectionMatrix * vec4(gl_MultiTexCoord2.st, 0.0, 1.0);
  //gl_Position = gl_ProjectionMatrix * vec4(0.5,0.5,0.0,1.0);
}

