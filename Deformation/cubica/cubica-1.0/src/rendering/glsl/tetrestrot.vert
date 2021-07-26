/////////////////////////////////////////////////
// Using 'this' vertex position, and four neighbor
// vertex positions, construct two vectors of an
// orthonormal frame
/////////////////////////////////////////////////

uniform sampler2D tetMeshTex;
uniform vec2 offset;

varying vec3 n1Frame;
varying vec3 n2Frame;

void main()
{
  vec2 neighbors[4];
  neighbors[0] = gl_MultiTexCoord0.st;
  neighbors[1] = gl_MultiTexCoord0.pq;
  neighbors[2] = gl_MultiTexCoord1.st;
  neighbors[3] = gl_MultiTexCoord1.pq;
  
  // Grab position of 'this' vertex in tetMesh
  vec3 p = texture2D(tetMeshTex, gl_Vertex.xy+offset).rgb;
  // Grab positions of four neighbor vertices
  vec3 v1 = normalize(texture2D(tetMeshTex, neighbors[0]+offset).rgb - p);
  vec3 v2 = normalize(texture2D(tetMeshTex, neighbors[1]+offset).rgb - p);
  vec3 v3 = normalize(texture2D(tetMeshTex, neighbors[2]+offset).rgb - p);
  vec3 v4 = normalize(texture2D(tetMeshTex, neighbors[3]+offset).rgb - p);

  n1Frame = normalize(v1 + v2 + v3);
  n2Frame = normalize(cross(n1Frame, v4));  

  vec4 pos = vec4(gl_Vertex.xy, 0.0, 1.0);
  gl_Position = gl_ModelViewProjectionMatrix * pos;
}

