/////////////////////////////////////////////////
// OpenGL vertex shader for rendering
// embedded tet mesh with NO SHADOWS
/////////////////////////////////////////////////

uniform sampler2D tetMeshTex;
uniform sampler2D tetNormalTex;
uniform sampler2D tetTangentTex;

varying vec3 Light;
varying vec3 Normal;

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
  vec3 pos = subtract * v3 + v0 * coords.x + v1 * coords.y + v2 * coords.z;

  v0 = texture2D(tetNormalTex, gl_MultiTexCoord0.st).rgb;
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
  Normal = normalize(gl_NormalMatrix * norm);
  //Normal = normalize(norm);

  //norm = gl_NormalMatrix * norm;

  vec4 temppos = gl_ModelViewMatrix * vec4(pos, 1.0);

  if (gl_LightSource[0].position.w == 1.0)
  {
    Light = gl_LightSource[0].position.xyz - temppos.xyz;
  } else {
    Light = -gl_LightSource[0].position.xyz;
  }
  Light = normalize(Light);
  //float NdotL = max(dot(norm, Light), 0.0);
  //NdotL = (NdotL + 0.25) * 0.8;
    
  // Color and output position
//  gl_FrontColor = gl_Color;
  //gl_FrontColor.rgb = NdotL * vec3(1.0);
  //gl_FrontColor.a = 1.0;

  //vec3 diff = abs(gl_Normal - tangent);

  //gl_FrontColor = vec4(abs(gl_Normal), 1.0);
  //gl_FrontColor = vec4(diff, 1.0);
  //gl_FrontColor = vec4(abs(tangent), 1.0);

  gl_Position = gl_ModelViewProjectionMatrix * vec4(pos, 1.0);
}

