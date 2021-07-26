/////////////////////////////////////////////////
// OpenGL vertex shader for deferred rendering
// face normal, per-pixel linear depth, and world position
/////////////////////////////////////////////////

uniform sampler2D tetMeshTex;

varying vec2 Depth;
varying vec3 Pos;
//varying float Epsilon;

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
  vec3 objpos = subtract * v3 + v0 * coords.x + v1 * coords.y + v2 * coords.z;

  vec4 eyepos = gl_ModelViewMatrix * vec4(objpos, 1.0);
  
  Depth = (gl_ProjectionMatrix * eyepos).zw;
  Pos = eyepos.xyz;
  gl_Position = gl_ProjectionMatrix * eyepos;
}

