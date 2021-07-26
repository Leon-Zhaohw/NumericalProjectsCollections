

uniform sampler2D tetMeshTex;
uniform sampler2D n1RestFrame;
uniform sampler2D n2RestFrame;
uniform vec2 offset;

varying vec3 n1Rot;
varying vec3 n2Rot;
varying vec3 n3Rot;

varying vec3 n1Rest;
varying vec3 n2Rest;
varying vec3 n3Rest;

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
  vec3 v1 = normalize(texture2D(tetMeshTex, neighbors[0]+offset).rgb - p);
  vec3 v2 = normalize(texture2D(tetMeshTex, neighbors[1]+offset).rgb - p);
  vec3 v3 = normalize(texture2D(tetMeshTex, neighbors[2]+offset).rgb - p);
  vec3 v4 = normalize(texture2D(tetMeshTex, neighbors[3]+offset).rgb - p);

  vec3 n1 = texture2D(n1RestFrame, coord).rgb;
  vec3 n2 = texture2D(n2RestFrame, coord).rgb;
  vec3 n3 = normalize(cross(n1, n2));
  
  n1Rest = (gl_ModelViewProjectionMatrix * vec4(n1, 0.0)).xyz;
  n2Rest = (gl_ModelViewProjectionMatrix * vec4(n2, 0.0)).xyz;
  n3Rest = (gl_ModelViewProjectionMatrix * vec4(n3, 0.0)).xyz;

  n1 = normalize(v1 + v2 + v3);
  n2 = normalize(cross(n1, v4));
  n3 = normalize(cross(n1, n2));

  n1Rot = (gl_ModelViewProjectionMatrix * vec4(n1, 0.0)).xyz;
  n2Rot = (gl_ModelViewProjectionMatrix * vec4(n2, 0.0)).xyz;
  n3Rot = (gl_ModelViewProjectionMatrix * vec4(n3, 0.0)).xyz;

  gl_Position = gl_ModelViewProjectionMatrix * vec4(p, 1.0);
}

