/////////////////////////////////////////////////
// Vertex shader for ambient occlusion volumes.
// Takes texture coords of three points, from the
// transformed vertex texture tetMeshTex, and passes
// the points (and other info) on.
/////////////////////////////////////////////////

uniform sampler2D tetMeshTex;
uniform vec2 offset;

// Triangle points
varying vec4 P0;
varying vec4 P1;
varying vec4 P2;

// Edge (0..2) and negative face (3) normals
varying vec4 M0;
varying vec4 M1;
varying vec4 M2;
varying vec4 M3;

void main()
{
  P0.xyz = texture2D(tetMeshTex, gl_Vertex.xy + offset).rgb;
  P1.xyz = texture2D(tetMeshTex, gl_Vertex.zw + offset).rgb;
  P2.xyz = texture2D(tetMeshTex, gl_MultiTexCoord0.xy + offset).rgb;

  vec3 n1 = P1.xyz - P0.xyz;
  vec3 n2 = P2.xyz - P0.xyz;
  M3 = vec4(normalize(cross(n2, n1)), 0.0);

  /*vec3 n1 = P1.xyz - P0.xyz;
  vec3 n2 = P2.xyz - P0.xyz;
  vec3 norm = normalize(cross(n1, n2));*/

  // Debug stuff
  //P0.xyz -= 0.005*M3;
  //P1.xyz -= 0.005*M3;
  //P2.xyz -= 0.005*M3;

 // vec3 norm = -M3.xyz;
 // P0.xyz -= 0.0025*norm;
 // P1.xyz -= 0.0025*norm;
 // P2.xyz -= 0.0025*norm;


  M0 = vec4(normalize(cross(n1, M3.xyz)), 0.0);
  M2 = vec4(normalize(cross(M3.xyz, n2)), 0.0);

  n1 = P2.xyz - P1.xyz;  
  M1 = vec4(normalize(cross(n1, M3.xyz)), 0.0);
  
  P0 = gl_ModelViewMatrix * vec4(P0.xyz, 1.0);
  P1 = gl_ModelViewMatrix * vec4(P1.xyz, 1.0);
  P2 = gl_ModelViewMatrix * vec4(P2.xyz, 1.0);

  M0 = gl_ModelViewMatrix * vec4(M0.xyz, 0.0);
  M1 = gl_ModelViewMatrix * vec4(M1.xyz, 0.0);
  M2 = gl_ModelViewMatrix * vec4(M2.xyz, 0.0);
  M3 = gl_ModelViewMatrix * vec4(M3.xyz, 0.0);

  /*M0 = gl_NormalMatrix * M0;
  M1 = gl_NormalMatrix * M1;
  M2 = gl_NormalMatrix * M2;
  M3 = gl_NormalMatrix * M3;*/

  gl_Position = vec4(0.0);
}

