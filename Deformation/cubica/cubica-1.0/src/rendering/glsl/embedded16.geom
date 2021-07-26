/////////////////////////////////////////////////
// OpenGL geometry shader for rendering
// deformed-embedded vertices and normals into
// the embedded vertex textures
// Input: TRIANGLES
// Output: POINTS
// VerticesOut: 3

// Used in: _embeddedMeshToTextureShader
/////////////////////////////////////////////////

#define _INVERSE_PI_ 0.318309886183791

// Embedded vertex positions
varying in vec3 vPos[];
// Embedded variable deltas
varying in float vDelta[];

// Position and normal values to pass to fragment shader
varying out vec4 fPos;
varying out vec3 fNormal;
varying out float fDelta;

void main()
{
  // Grab the three triangle edges
  vec3 a = vPos[1] - vPos[0];
  if (dot(a,a) > 0.0)
    a = normalize(a);
  vec3 b = vPos[2] - vPos[0];
  if (dot(b,b) > 0.0)
    b = normalize(b);
  vec3 c = vPos[2] - vPos[1];
  if (dot(c,c) > 0.0)
    c = normalize(c);

  // Keep the largest delta
  float maxDelta = 0.0;

  maxDelta = max(maxDelta, vDelta[0]);
  maxDelta = max(maxDelta, vDelta[1]);
  maxDelta = max(maxDelta, vDelta[2]);

  // Weight of this triangles' contribution to the vertex normal
  float weight;

  // Calculate the triangles' face normal
  vec3 norm = cross(a,b);
  if (dot(norm,norm) > 0.0)
    norm = normalize(norm);

  // Pass vertex position and weighted normal (output as points)
  gl_Position = gl_PositionIn[0];
  fPos = vec4(vPos[0], 1.0);
  weight = acos(clamp(dot( a, b), -1.0, 1.0)) * _INVERSE_PI_;
  fNormal = norm * weight;
  fDelta = maxDelta;
  EmitVertex();
  EndPrimitive();

  gl_Position = gl_PositionIn[1];
  fPos = vec4(vPos[1], 1.0);
  weight = acos(clamp(dot( c,-a), -1.0, 1.0)) * _INVERSE_PI_;
  fNormal = norm * weight;
  fDelta = maxDelta;
  EmitVertex();
  EndPrimitive();

  gl_Position = gl_PositionIn[2];
  fPos = vec4(vPos[2], 1.0);
  weight = acos(clamp(dot(-b,-c), -1.0, 1.0)) * _INVERSE_PI_;
  fNormal = norm * weight;
  fDelta = maxDelta;
  EmitVertex();
  EndPrimitive();

}

