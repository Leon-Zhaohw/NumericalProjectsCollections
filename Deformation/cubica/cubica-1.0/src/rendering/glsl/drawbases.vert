

uniform sampler2D tetMeshTex;
uniform sampler2D Rotation1;
uniform sampler2D Rotation2;
uniform sampler2D Rotation3;
uniform vec2 offset;

varying vec4 r1;
varying vec4 r2;
varying vec4 r3;

varying vec3 c1;
varying vec3 c2;
varying vec3 c3;

void main()
{
  vec2 coord = gl_Vertex.xy + offset;
  /*vec2 neighbors[4];
  neighbors[0] = gl_MultiTexCoord0.st;
  neighbors[1] = gl_MultiTexCoord0.pq;
  neighbors[2] = gl_MultiTexCoord1.st;
  neighbors[3] = gl_MultiTexCoord1.pq;*/
  
  // Grab position of 'this' vertex in tetMesh
  vec3 p = texture2D(tetMeshTex, coord).rgb;
  
  // Grab rotation matrix
  c1 = normalize(texture2D(Rotation1, coord).rgb);
  c2 = normalize(texture2D(Rotation2, coord).rgb);
  c3 = normalize(texture2D(Rotation3, coord).rgb);

  r1 = gl_ModelViewMatrix * vec4(c1, 0.0);
  r2 = gl_ModelViewMatrix * vec4(c2, 0.0);
  r3 = gl_ModelViewMatrix * vec4(c3, 0.0);
    
  gl_Position = gl_ModelViewProjectionMatrix * vec4(p, 1.0);
}

