

uniform sampler2D tetMeshTex;
uniform vec2 offset;

varying vec4 n1;
varying vec4 n2;
varying vec4 n3;
varying vec4 n4;

void main()
{
  vec2 coord = gl_Vertex.xy + offset;
  vec2 neighbors[4];
  neighbors[0] = gl_MultiTexCoord0.st;
  neighbors[1] = gl_MultiTexCoord0.pq;
  neighbors[2] = gl_MultiTexCoord1.st;
  neighbors[3] = gl_MultiTexCoord1.pq;
  
  // Grab position of 'this' vertex in tetMesh
  vec3 p = texture2D(tetMeshTex, coord).rgb;
  // Grab positions of four neighbor vertices
  n1.xyz = texture2D(tetMeshTex, neighbors[0]+offset).rgb;
  n2.xyz = texture2D(tetMeshTex, neighbors[1]+offset).rgb;
  n3.xyz = texture2D(tetMeshTex, neighbors[2]+offset).rgb;
  n4.xyz = texture2D(tetMeshTex, neighbors[3]+offset).rgb;

  n1 = gl_ModelViewProjectionMatrix * vec4(n1.xyz, 1.0);
  n2 = gl_ModelViewProjectionMatrix * vec4(n2.xyz, 1.0);
  n3 = gl_ModelViewProjectionMatrix * vec4(n3.xyz, 1.0);
  n4 = gl_ModelViewProjectionMatrix * vec4(n4.xyz, 1.0);

  gl_Position = gl_ModelViewProjectionMatrix * vec4(p, 1.0);
}

