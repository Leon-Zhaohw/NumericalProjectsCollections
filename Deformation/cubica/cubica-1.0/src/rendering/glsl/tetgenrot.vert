/////////////////////////////////////////////////
// Using 'this' vertex position, and four neighbor
// vertex positions, construct two vectors of an
// orthonormal frame
/////////////////////////////////////////////////

uniform sampler2D tetMeshTex;
uniform sampler2D n1RestFrame;
uniform sampler2D n2RestFrame;
uniform vec2 offset;

varying vec3 Rotation1;
varying vec3 Rotation2;
varying vec3 Rotation3;

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
  vec3 v1 = texture2D(tetMeshTex, neighbors[0]+offset).rgb - p;
  if (dot(v1,v1) > 0.0)    v1 = normalize(v1);
  vec3 v2 = texture2D(tetMeshTex, neighbors[1]+offset).rgb - p;
  if (dot(v2,v2) > 0.0)    v2 = normalize(v2);
  vec3 v3 = texture2D(tetMeshTex, neighbors[2]+offset).rgb - p;
  if (dot(v3,v3) > 0.0)    v3 = normalize(v3);
  vec3 v4 = texture2D(tetMeshTex, neighbors[3]+offset).rgb - p;
  if (dot(v4,v4) > 0.0)    v4 = normalize(v4);

  vec3 n1Rot = v1 + v2 + v3;
  if (dot(n1Rot,n1Rot) > 0.0)    n1Rot = normalize(n1Rot);
  vec3 n2Rot = cross(n1Rot, v4);
  if (dot(n2Rot,n2Rot) > 0.0)    n2Rot = normalize(n2Rot);
  vec3 n3Rot = cross(n1Rot, n2Rot);
  if (dot(n3Rot,n3Rot) > 0.0)    n3Rot = normalize(n3Rot);

  vec3 n1Rest = texture2D(n1RestFrame, coord).rgb;
  if (dot(n1Rest,n1Rest) > 0.0)    n1Rest = normalize(n1Rest);
  vec3 n2Rest = texture2D(n2RestFrame, coord).rgb;
  if (dot(n2Rest,n2Rest) > 0.0)    n2Rest = normalize(n2Rest);
  vec3 n3Rest = cross(n1Rest, n2Rest);
  if (dot(n3Rest,n3Rest) > 0.0)    n3Rest = normalize(n3Rest);

  // Construct rotation matrix R = Nrot * Nrest^T
  Rotation1.x = dot(n1Rot, n1Rest);
  Rotation1.y = dot(n1Rot, n2Rest);
  Rotation1.z = dot(n1Rot, n3Rest);

  Rotation2.x = dot(n2Rot, n1Rest);
  Rotation2.y = dot(n2Rot, n2Rest);
  Rotation2.z = dot(n2Rot, n3Rest);

  Rotation3.x = dot(n3Rot, n1Rest);
  Rotation3.y = dot(n3Rot, n2Rest);
  Rotation3.z = dot(n3Rot, n3Rest);// */

  /*Rotation1.x = dot(n1Rest, n1Rot);
  Rotation1.y = dot(n1Rest, n2Rot);
  Rotation1.z = dot(n1Rest, n3Rot);

  Rotation2.x = dot(n2Rest, n1Rot);
  Rotation2.y = dot(n2Rest, n2Rot);
  Rotation2.z = dot(n2Rest, n3Rot);

  Rotation3.x = dot(n3Rest, n1Rot);
  Rotation3.y = dot(n3Rest, n2Rot);
  Rotation3.z = dot(n3Rest, n3Rot);// */

  /*Rotation1 = n1Rot;
  Rotation2 = n2Rot;
  Rotation3 = n3Rot;// */
 
  vec4 pos = vec4(gl_Vertex.xy, 0.0, 1.0);
  gl_Position = gl_ModelViewProjectionMatrix * pos;
}

