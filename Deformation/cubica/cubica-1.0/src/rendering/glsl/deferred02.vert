/////////////////////////////////////////////////
// OpenGL vertex shader for rendering
// embedded tet mesh WITH SHADOWS
/////////////////////////////////////////////////

uniform sampler2D tetMeshTex;
uniform sampler2D tetRotation1Tex;
uniform sampler2D tetRotation2Tex;
uniform sampler2D tetRotation3Tex;

varying vec2 fDepth;
varying vec3 fPos;
varying vec3 fNormal;
//varying vec2 Depth;
//varying vec3 Pos;

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
  vec3 Pos = subtract * v3 + v0 * coords.x + v1 * coords.y + v2 * coords.z;

  vec3 r0 = texture2D(tetRotation1Tex, gl_MultiTexCoord0.st).rgb;
  vec3 r1 = texture2D(tetRotation1Tex, gl_MultiTexCoord0.pq).rgb;
  vec3 r2 = texture2D(tetRotation1Tex, gl_MultiTexCoord1.st).rgb;
  vec3 r3 = texture2D(tetRotation1Tex, gl_MultiTexCoord1.pq).rgb;
  
  v0.x = dot(r0, gl_Normal);
  v1.x = dot(r1, gl_Normal);
  v2.x = dot(r2, gl_Normal);
  v3.x = dot(r3, gl_Normal);
  
  r0 = texture2D(tetRotation2Tex, gl_MultiTexCoord0.st).rgb;
  r1 = texture2D(tetRotation2Tex, gl_MultiTexCoord0.pq).rgb;
  r2 = texture2D(tetRotation2Tex, gl_MultiTexCoord1.st).rgb;
  r3 = texture2D(tetRotation2Tex, gl_MultiTexCoord1.pq).rgb;

  v0.y = dot(r0, gl_Normal);
  v1.y = dot(r1, gl_Normal);
  v2.y = dot(r2, gl_Normal);
  v3.y = dot(r3, gl_Normal);

  r0 = texture2D(tetRotation3Tex, gl_MultiTexCoord0.st).rgb;
  r1 = texture2D(tetRotation3Tex, gl_MultiTexCoord0.pq).rgb;
  r2 = texture2D(tetRotation3Tex, gl_MultiTexCoord1.st).rgb;
  r3 = texture2D(tetRotation3Tex, gl_MultiTexCoord1.pq).rgb;

  v0.z = dot(r0, gl_Normal);
  v1.z = dot(r1, gl_Normal);
  v2.z = dot(r2, gl_Normal);
  v3.z = dot(r3, gl_Normal);

  vec3 norm = subtract * v3 + v0 * coords.x + v1 * coords.y + v2 * coords.z;
  fNormal = normalize(gl_NormalMatrix * normalize(norm));

  vec4 ePos = gl_ModelViewMatrix * vec4(Pos, 1.0);

  fDepth = (gl_ProjectionMatrix * ePos).zw;
  fPos = ePos.xyz;

  gl_Position = gl_ProjectionMatrix * ePos;
}

