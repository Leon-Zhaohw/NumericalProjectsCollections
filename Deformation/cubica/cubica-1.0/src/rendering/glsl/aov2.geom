/////////////////////////////////////////////////
// Geometry shader for ambient occlusion volumes.
// Takes the points of a triangle, and extrudes them
// into the volume.
/////////////////////////////////////////////////

// Uniform storing the maximum obscurance distance
//uniform float delta;

//varying in vec3 oPos[];
varying in float vDelta[];

// Output variables to the fragment shader
varying out vec3 p0;
varying out vec3 p1;
varying out vec3 p2;

varying out vec3 m0;
varying out vec3 m1;
varying out vec3 m2;
varying out vec3 m3;

varying out float oDelta;

vec3 P0;
vec3 P1;
vec3 P2;

vec4 M0;
vec4 M1;
vec4 M2;
vec4 M3;

void copyVarying()
{
  p0 = P0.xyz;  p1 = P1.xyz;  p2 = P2.xyz;
  m0 = M0.xyz;  m1 = M1.xyz;  m2 = M2.xyz;   m3 = M3.xyz;
  //m0 = M0;      m1 = M1;      m2 = M2;       m3 = M3;
}

void main()
{
  oDelta = 0.0;
  oDelta = max(oDelta, vDelta[0]);
  oDelta = max(oDelta, vDelta[1]);
  oDelta = max(oDelta, vDelta[2]);
 
  float delta2 = oDelta * 2.0;

  // Dividing by .w is unnecessary right now
  P0 = gl_PositionIn[0].xyz;// / gl_PositionIn[0].w;
  P1 = gl_PositionIn[1].xyz;// / gl_PositionIn[1].w;
  P2 = gl_PositionIn[2].xyz;// / gl_PositionIn[2].w;

  vec3 n1 = normalize(P1 - P0);
  vec3 n2 = normalize(P2 - P0);
  M3 = vec4(normalize(cross(n2, n1)), 0.0);
  M0 = vec4(normalize(cross(n1, M3.xyz)), 0.0);
  M2 = vec4(normalize(cross(M3.xyz, n2)), 0.0);

  n1 = P2 - P1;  
  M1 = vec4(normalize(cross(n1, M3.xyz)), 0.0);
   
  // Construct the vertices in the plane, clamping the extension vectors to 2*delta
  vec4 V0 = -oDelta * ((M0+M2) / (1.0 + dot(M0,M2)));
  if (length(V0) > delta2)  {
    V0 = delta2 * normalize(V0);
  }
  V0 = gl_ProjectionMatrix * (V0 + gl_PositionIn[0]);
  
  vec4 V1 = -oDelta * ((M1+M0) / (1.0 + dot(M1,M0)));
  if (length(V1) > delta2)  {
    V1 = delta2 * normalize(V1);
  }
  V1 = gl_ProjectionMatrix * (V1 + gl_PositionIn[1]);
  
  vec4 V2 = -oDelta * ((M2+M1) / (1.0 + dot(M2,M1)));
  if (length(V2) > delta2)  {
    V2 = delta2 * normalize(V2);
  }
  V2 = gl_ProjectionMatrix * (V2 + gl_PositionIn[2]);

  /*vec4 aV0 = gl_ProjectionMatrix * gl_PositionIn[0];
  vec4 aV1 = gl_ProjectionMatrix * gl_PositionIn[1];
  vec4 aV2 = gl_ProjectionMatrix * gl_PositionIn[2];// */

  vec4 dM3 = -oDelta * gl_ProjectionMatrix * M3;
 
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


/////////////////////////////////////////////
// Debugging stuff

/*  gl_Position = gl_ProjectionMatrix * gl_PositionIn[0];
  copyVarying();
  EmitVertex();  
  gl_Position = gl_ProjectionMatrix * gl_PositionIn[1];
  copyVarying();
  EmitVertex(); 
  gl_Position = gl_ProjectionMatrix * gl_PositionIn[2];
  copyVarying();
  EmitVertex();
  gl_Position = gl_ProjectionMatrix * gl_PositionIn[0];
  copyVarying();
  EmitVertex();  
  EndPrimitive();

  gl_Position = gl_ProjectionMatrix * 0.5 * (gl_PositionIn[0] + gl_PositionIn[1]);
  copyVarying();
  EmitVertex(); 
  gl_Position = gl_ProjectionMatrix * ((-oDelta * 0.3* M3) + (0.5 * (gl_PositionIn[0] + gl_PositionIn[1])) - (oDelta * (M0)));
  copyVarying();
  EmitVertex();
  EndPrimitive();


  gl_Position = gl_ProjectionMatrix * 0.5 * (gl_PositionIn[1] + gl_PositionIn[2]);
  copyVarying();
  EmitVertex(); 
  gl_Position = gl_ProjectionMatrix * ((-oDelta * 0.3* M3) + (0.5 * (gl_PositionIn[1] + gl_PositionIn[2])) - (oDelta * (M1)));
  copyVarying();
  EmitVertex();
EndPrimitive();

  gl_Position = gl_ProjectionMatrix * 0.5 * (gl_PositionIn[2] + gl_PositionIn[0]);
  copyVarying();
  EmitVertex(); 
  gl_Position = gl_ProjectionMatrix * ((-oDelta * 0.3* M3) + (0.5 * (gl_PositionIn[2] + gl_PositionIn[0])) - (oDelta * (M2)));
  copyVarying();
  EmitVertex();
EndPrimitive();


  gl_Position = gl_ProjectionMatrix * gl_PositionIn[0];
  copyVarying();
  EmitVertex(); 
  gl_Position = gl_ProjectionMatrix * (V0 + gl_PositionIn[0] + (-oDelta * M3));
  copyVarying();
  EmitVertex();  
  EndPrimitive();

  gl_Position = gl_ProjectionMatrix * gl_PositionIn[1];
  copyVarying();
  EmitVertex();
  gl_Position = gl_ProjectionMatrix * (V1 + gl_PositionIn[1] + (-oDelta * M3));
  copyVarying();
  EmitVertex();
  EndPrimitive();
  
  gl_Position = gl_ProjectionMatrix * gl_PositionIn[2];
  copyVarying();
  EmitVertex();
  gl_Position = gl_ProjectionMatrix * (V2 + gl_PositionIn[2] + (-oDelta * M3));
  copyVarying();
  EmitVertex();
  EndPrimitive();*/


/////////////////////////////////////////////

