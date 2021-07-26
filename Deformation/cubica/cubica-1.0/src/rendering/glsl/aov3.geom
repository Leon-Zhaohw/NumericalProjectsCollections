/////////////////////////////////////////////////
// Geometry shader for ambient occlusion volumes.
// Takes the points of a triangle, and extrudes them
// into the volume.
/////////////////////////////////////////////////

// Uniform storing the maximum obscurance distance
uniform float delta;

vec3 P0;
vec3 P1;
vec3 P2;

vec4 M0;
vec4 M1;
vec4 M2;
vec4 M3;

void main()
{
  float delta2 = delta * 2.0;

   // Dividing by .w is unnecessary right now
  P0 = gl_PositionIn[0].xyz / gl_PositionIn[0].w;
  P1 = gl_PositionIn[1].xyz / gl_PositionIn[1].w;
  P2 = gl_PositionIn[2].xyz / gl_PositionIn[2].w;

  vec3 n1 = normalize(P1 - P0);
  vec3 n2 = normalize(P2 - P0);
  M3 = vec4(normalize(cross(n2, n1)), 0.0);
  M0 = vec4(normalize(cross(n1, M3.xyz)), 0.0);
  M2 = vec4(normalize(cross(M3.xyz, n2)), 0.0);

  n1 = P2 - P1;  
  M1 = vec4(normalize(cross(n1, M3.xyz)), 0.0);
   
  // Construct the vertices in the plane, clamping the extension vectors to 2*delta
  vec4 V0 = -delta * ((M0+M2) / (1.0 + dot(M0,M2)));
  if (length(V0) > delta2)  {
    V0 = delta2 * normalize(V0);
  }
  V0 = gl_ProjectionMatrix * (V0 + gl_PositionIn[0]);
  
  vec4 V1 = -delta * ((M1+M0) / (1.0 + dot(M1,M0)));
  if (length(V1) > delta2)  {
    V1 = delta2 * normalize(V1);
  }
  V1 = gl_ProjectionMatrix * (V1 + gl_PositionIn[1]);
  
  vec4 V2 = -delta * ((M2+M1) / (1.0 + dot(M2,M1)));
  if (length(V2) > delta2)  {
    V2 = delta2 * normalize(V2);
  }
  V2 = gl_ProjectionMatrix * (V2 + gl_PositionIn[2]);

  vec4 midP = gl_ProjectionMatrix * ((gl_PositionIn[0] + gl_PositionIn[1] + gl_PositionIn[2]) / 3.0);

  vec4 dM3 = -delta * gl_ProjectionMatrix * M3;
  vec4 dM0 = 0.5 * delta * gl_ProjectionMatrix * M0;
  vec4 dM1 = 0.5 * delta * gl_ProjectionMatrix * M1;
  vec4 dM2 = 0.5 * delta * gl_ProjectionMatrix * M2;

  midP = midP + (0.25 * dM3);

  /*gl_Position = gl_ProjectionMatrix * gl_PositionIn[0];
  gl_FrontColor = vec4(1.0);
  EmitVertex();*/
  gl_Position = V0;
  gl_FrontColor = vec4(1.0);
  EmitVertex();  
  gl_Position = V1;
  gl_FrontColor = vec4(1.0);
  EmitVertex();
  gl_Position = V1 + dM3;
  gl_FrontColor = vec4(1.0);
  EmitVertex();
  gl_Position = V0 + dM3;
  gl_FrontColor = vec4(1.0);
  EmitVertex();
  gl_Position = V0;
  gl_FrontColor = vec4(1.0);
  EmitVertex();
  EndPrimitive();


  /*gl_Position = gl_ProjectionMatrix * gl_PositionIn[1];
  gl_FrontColor = vec4(1.0);
  EmitVertex();*/
  gl_Position = V1;
  gl_FrontColor = vec4(1.0);
  EmitVertex();  
  gl_Position = V2;
  gl_FrontColor = vec4(1.0);
  EmitVertex();
  gl_Position = V2 + dM3;
  gl_FrontColor = vec4(1.0);
  EmitVertex();
  gl_Position = V1 + dM3;
  gl_FrontColor = vec4(1.0);
  EmitVertex();
  gl_Position = V1;
  gl_FrontColor = vec4(1.0);
  EmitVertex();
  EndPrimitive();


  /*gl_Position = gl_ProjectionMatrix * gl_PositionIn[2];
  gl_FrontColor = vec4(1.0);
  EmitVertex();*/
  gl_Position = V2;
  gl_FrontColor = vec4(1.0);
  EmitVertex();  
  gl_Position = V0;
  gl_FrontColor = vec4(1.0);
  EmitVertex();
  gl_Position = V0 + dM3;
  gl_FrontColor = vec4(1.0);
  EmitVertex();
  gl_Position = V2 + dM3;
  gl_FrontColor = vec4(1.0);
  EmitVertex();
  gl_Position = V2;
  gl_FrontColor = vec4(1.0);
  EmitVertex();
  EndPrimitive();

  /*gl_Position = midP;
  gl_FrontColor = vec4(0.3, 0.9, 0.3, 1.0);
  EmitVertex();
  gl_Position = midP + dM0;
  gl_FrontColor = vec4(0.3, 0.9, 0.3, 1.0);
  EmitVertex();
  EndPrimitive();

  gl_Position = midP;
  gl_FrontColor = vec4(0.9, 0.3, 0.3, 1.0);
  EmitVertex();
  gl_Position = midP + dM1;
  gl_FrontColor = vec4(0.9, 0.3, 0.3, 1.0);
  EmitVertex();
  EndPrimitive();

  gl_Position = midP;
  gl_FrontColor = vec4(0.3, 0.3, 0.9, 1.0);
  EmitVertex();
  gl_Position = midP + dM2;
  gl_FrontColor = vec4(0.3, 0.3, 0.9, 1.0);
  EmitVertex();
  EndPrimitive();

  gl_Position = midP;
  gl_FrontColor = vec4(0.8, 0.8, 0.3, 1.0);
  EmitVertex();
  gl_Position = midP + (0.5*dM3);
  gl_FrontColor = vec4(0.8, 0.8, 0.3, 1.0);
  EmitVertex();
  EndPrimitive();*/
}




