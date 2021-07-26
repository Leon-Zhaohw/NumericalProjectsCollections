/////////////////////////////////////////////////
// OpenGL geometry shader for deferred rendering
// face normal, per-pixel linear depth, and 
// per-pixel world position
/////////////////////////////////////////////////
// input: GL_TRIANGLES
// output: GL_TRIANGLE_STRIP

varying in vec2 Depth[];
varying in vec3 Pos[];

varying out vec2 fDepth;
varying out vec3 fPos;
varying out vec3 fNormal;

void main()
{
  vec3 n1 = Pos[1] - Pos[0];
  vec3 n2 = Pos[2] - Pos[0];
  vec3 norm = normalize(cross(n1,n2));

  int x;
  for (x = 0; x < 3; x++)
  {
    gl_Position = gl_PositionIn[x];
    fDepth = Depth[x];
    fPos = Pos[x];
    fNormal = norm;
    EmitVertex();
  }
  EndPrimitive();
}
