/////////////////////////////////////////////////
// OpenGL vertex shader for rendering
// embedded tet mesh with NO SHADOWS
/////////////////////////////////////////////////

uniform sampler2D tetMeshTex;
uniform sampler2D tetRotation1Tex;
uniform sampler2D tetRotation2Tex;
uniform sampler2D tetRotation3Tex;

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

  //vec3 norm = v0 * coords.x + v1 * coords.y + v2 * coords.z;
  vec3 norm = subtract * v3 + v0 * coords.x + v1 * coords.y + v2 * coords.z;

 // float error = length(gl_Normal - v0);

 // if (error < 0.01)
 //   Normal = vec3(1.0, 0.0, 0.0);
 // else
 //   Normal = vec3(0.0, 1.0, 0.0);
  //gl_FrontColor.rgb = abs(normalize(norm));
  //gl_FrontColor.a = 1.0;
  //Normal = normalize(gl_NormalMatrix * normalize(norm)); //normalize(gl_NormalMatrix * normalize(norm));
  Normal = normalize(gl_NormalMatrix * normalize(norm));
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

  //pos = pos - 0.01 * norm;

  gl_Position = gl_ModelViewProjectionMatrix * vec4(pos, 1.0);
}

