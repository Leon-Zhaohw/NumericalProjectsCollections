/////////////////////////////////////////////////
// Vertex shader for ambient occlusion volumes.
// Takes texture coords of three points, from the
// transformed vertex texture tetMeshTex, and passes
// the points (and other info) on.
/////////////////////////////////////////////////

uniform sampler2D tetMeshTex;
//uniform vec2 offset;

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
  vec3 oPos = subtract * v3 + v0 * coords.x + v1 * coords.y + v2 * coords.z;

  gl_Position = gl_ModelViewMatrix * vec4(oPos, 1.0);
}

