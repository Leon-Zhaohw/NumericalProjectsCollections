/////////////////////////////////////////////////
// OpenGL vertex shader for rendering
// embedded tet mesh INTO A SHADOWMAP
/////////////////////////////////////////////////

uniform sampler2D tetMeshTex;
//uniform sampler2D tetNormalTex;
//uniform sampler2D tetTangentTex;

varying vec2 Depth;
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

 /* v0 = texture2D(tetNormalTex, gl_MultiTexCoord0.st).rgb;
  v1 = texture2D(tetNormalTex, gl_MultiTexCoord0.pq).rgb;
  v2 = texture2D(tetNormalTex, gl_MultiTexCoord1.st).rgb;

  if (length(v0) > 0.0)
    v0 = normalize(v0);
  if (length(v1) > 0.0)
    v1 = normalize(v1);
  if (length(v2) > 0.0)  
    v2 = normalize(v2);

  vec3 tetnorm = coords.x * v0 + coords.y * v1 + coords.z * v2;
  if (length(tetnorm) > 0.0)
    tetnorm = normalize(tetnorm);
    
  v0 = texture2D(tetTangentTex, gl_MultiTexCoord0.st).rgb;
  v1 = texture2D(tetTangentTex, gl_MultiTexCoord0.pq).rgb;
  v2 = texture2D(tetTangentTex, gl_MultiTexCoord1.st).rgb;

  if (length(v0) > 0.0)
    v0 = normalize(v0);
  if (length(v1) > 0.0)
    v1 = normalize(v1);
  if (length(v2) > 0.0)  
    v2 = normalize(v2);

  vec3 tangent = v0 * coords.x + v1 * coords.y + v2 * coords.z;
  if (length(tangent) > 0.0)
    tangent = normalize(tangent);

  vec3 bitangent = cross(tetnorm, tangent);
  
  vec3 norm;
  norm.x = dot(tangent,   gl_Normal);
  norm.y = dot(bitangent, gl_Normal);
  norm.z = dot(tetnorm,   gl_Normal);
  norm = gl_NormalMatrix * normalize(norm);*/

  vec4 eyepos = gl_ModelViewMatrix * vec4(objpos, 1.0);

  /*vec3 V = normalize(eyepos.xyz / eyepos.w);
  float NdotV = dot(norm, V);
  Epsilon = 1.0 - max(0.0, NdotV);*/
  
  Depth = (gl_ProjectionMatrix * eyepos).zw;
  gl_Position = gl_ProjectionMatrix * eyepos;
}

