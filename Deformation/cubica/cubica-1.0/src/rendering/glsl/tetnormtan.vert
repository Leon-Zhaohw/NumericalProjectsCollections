/////////////////////////////////////////////////
// OpenGL vertex shader 
/////////////////////////////////////////////////

// #define TARGET 0   // inserted within the program

uniform sampler2D tetMeshTex;
uniform vec2 offset;

varying vec3 Normal;
varying vec3 Tangent;

void main()
{
  vec2 coord[3];
  coord[0] = gl_Vertex.st;
  coord[1] = gl_Vertex.pq;
  coord[2] = gl_MultiTexCoord0.st;
  
  // Grab positions of three tetMesh face vertices
  vec3 v0 = texture2D(tetMeshTex, coord[0]+offset).rgb;
  vec3 v1 = texture2D(tetMeshTex, coord[1]+offset).rgb;
  vec3 v2 = texture2D(tetMeshTex, coord[2]+offset).rgb;

  Tangent = v1 - v0;
  Tangent = normalize(Tangent);
  vec3 bitan = v2 - v0;
  
  Normal = cross(Tangent, bitan);
  Normal = normalize(Normal);

  vec4 pos = vec4(coord[TARGET], 0.0, 1.0);
 
  // Color and output position
  gl_Position = gl_ModelViewProjectionMatrix * pos;

}

