/////////////////////////////////////////////////
// OpenGL geometry shader for deferred rendering
// face normal, per-pixel linear depth, and 
// per-pixel world position
/////////////////////////////////////////////////
// input: GL_TRIANGLES
// output: GL_TRIANGLE_STRIP

varying in vec2 vDepth[];
varying in vec3 vPos[];

varying out vec2 Depth;
varying out vec4 Pos;
varying out vec3 Normal;

void main()
{
  vec3 n1 = vPos[1] - vPos[0];
  vec3 n2 = vPos[2] - vPos[0];
  vec3 norm = normalize(cross(n1,n2));

  int x;
  for (x = 0; x < 3; x++)
  {
    gl_Position = gl_PositionIn[x];
    Depth = vDepth[x];
    Pos = vec4(vPos[x], 1.0);
    Normal = norm;
    EmitVertex();
  }
  EndPrimitive();
}
