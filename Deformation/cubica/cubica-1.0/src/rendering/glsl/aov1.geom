/////////////////////////////////////////////////
// Geometry shader for ambient occlusion volumes.
// Takes the points of a triangle, and extrudes them
// into the volume.
/////////////////////////////////////////////////

// Uniform storing the maximum obscurance distance
//uniform float delta;
varying out float oDelta;

// Input variables, storing triangle positions, and edge/face normals
varying in vec4 P0[];
varying in vec4 P1[];
varying in vec4 P2[];

varying in vec4 M0[];
varying in vec4 M1[];
varying in vec4 M2[];
varying in vec4 M3[];

// Output variables to the fragment shader
varying out vec3 p0;
varying out vec3 p1;
varying out vec3 p2;

varying out vec3 m0;
varying out vec3 m1;
varying out vec3 m2;
varying out vec3 m3;

/*vec3 M0;
vec3 M1;
vec3 M2;
vec3 M3;*/

vec4 V0;
vec4 V1;
vec4 V2;

void copyVarying()
{
  p0 = P0[0].xyz;  p1 = P1[0].xyz;  p2 = P2[0].xyz;
  //p0 = V0.xyz;  p1 = V1.xyz;  p2 = V2.xyz;
  m0 = M0[0].xyz;  m1 = M1[0].xyz;  m2 = M2[0].xyz;   m3 = M3[0].xyz;
  //m0 = M0;      m1 = M1;      m2 = M2;       m3 = M3;
}

void main()
{
  float delta = 0.1;
  float delta2 = delta * 2.0;

 /* vec3 n1 = P1[0].xyz - P0[0].xyz;
  vec3 n2 = P2[0].xyz - P0[0].xyz;
  M3 = normalize(cross(n2, n1));
  M0 = normalize(cross(n1, M3));
  M2 = normalize(cross(M3, n2));

  n1 = P2[0].xyz - P1[0].xyz;  
  M1 = normalize(cross(n1, M3));*/
  
  // Construct the vertices in the plane, clamping the extension vectors to 2*delta
  V0 = -delta * ((M0[0]+M2[0]) / (1.0 + dot(M0[0],M2[0])));
  if (length(V0) > delta2)  {
    V0 = delta2 * normalize(V0);
  }
  V0 = gl_ProjectionMatrix * (V0 + P0[0]);
  
  V1 = -delta * ((M1[0]+M0[0]) / (1.0 + dot(M1[0],M0[0])));
  if (length(V1) > delta2)  {
    V1 = delta2 * normalize(V1);
  }
  V1 = gl_ProjectionMatrix * (V1 + P1[0]);
  
  V2 = -delta * ((M2[0]+M1[0]) / (1.0 + dot(M2[0],M1[0])));
  if (length(V2) > delta2)  {
    V2 = delta2 * normalize(V2);
  }
  V2 = gl_ProjectionMatrix * (V2 + P2[0]);

  /*V0 = P0[0];
  V1 = P1[0];
  V2 = P2[0];// */

  vec4 dM3 = -delta * gl_ProjectionMatrix * M3[0];

  oDelta = delta;
 
  // Plane B0..B2
  gl_Position = V0 + dM3;
  copyVarying();
//  gl_FrontColor = vec4(1.0, 0.3, 0.3, 1.0);
  EmitVertex();  
  gl_Position = V0;
  copyVarying();
//  gl_FrontColor = vec4(0.5, 0.3, 0.3, 1.0);
  EmitVertex();
  gl_Position = V1 + dM3;
  copyVarying();
//  gl_FrontColor = vec4(0.3, 1.0, 0.3, 1.0);
  EmitVertex();
  gl_Position = V1;
  copyVarying();
//  gl_FrontColor = vec4(0.3, 0.5, 0.3, 1.0);
  EmitVertex();
  gl_Position = V2 + dM3;
  copyVarying();
//  gl_FrontColor = vec4(0.3, 0.3, 1.0, 1.0);
  EmitVertex();
  gl_Position = V2;
  copyVarying();
//  gl_FrontColor = vec4(0.3, 0.3, 0.5, 1.0);
  EmitVertex();
  gl_Position = V0 + dM3;
  copyVarying();
//  gl_FrontColor = vec4(1.0, 0.3, 0.3, 1.0);
  EmitVertex();
  gl_Position = V0;
  copyVarying();
//  gl_FrontColor = vec4(0.5, 0.3, 0.3, 1.0);
  EmitVertex();
  EndPrimitive();

  // Plane B3 (above the triangle)
  gl_Position = V0 + dM3;
  copyVarying();
//  gl_FrontColor = vec4(1.0, 0.3, 0.3, 1.0);
  EmitVertex();
  gl_Position = V1 + dM3;
  copyVarying();
//  gl_FrontColor = vec4(0.3, 1.0, 0.3, 1.0);
  EmitVertex();
  gl_Position = V2 + dM3;
  copyVarying();
//  gl_FrontColor = vec4(0.3, 0.3, 1.0, 1.0);
  EmitVertex();
  EndPrimitive();

  // Plane B4 (containing the triangle)
  gl_Position = V0;
  copyVarying();
//  gl_FrontColor = vec4(0.5, 0.3, 0.3, 1.0);
  EmitVertex();
  gl_Position = V2;
  copyVarying();
//  gl_FrontColor = vec4(0.3, 0.3, 0.5, 1.0);
  EmitVertex();
  gl_Position = V1;
  copyVarying();
//  gl_FrontColor = vec4(0.3, 0.5, 0.3, 1.0);
  EmitVertex();
  EndPrimitive();
}
