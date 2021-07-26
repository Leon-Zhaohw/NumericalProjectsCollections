/////////////////////////////////////////////////
// Fragment shader used by ambiend occlusion volumes
/////////////////////////////////////////////////

#define INVERSE_2_PI  0.159154943091895
#define DOT_EPSILON  -0.001
#define EPSILON       0.00001

uniform sampler2D normalTex;
uniform sampler2D positionTex;
uniform vec2 screenSize;
//uniform float delta;

varying float oDelta;

varying vec3 p0;
varying vec3 p1;
varying vec3 p2;

varying vec3 m0;
varying vec3 m1;
varying vec3 m2;
varying vec3 m3;

//vec3 color;

float AOp3(in vec3 N, in vec3 P0, in vec3 P1, in vec3 P2)
{
  float aop = 0.0;
  float PidotPj;
  vec3 PixPj;
  
  PidotPj = clamp(dot(P0,P1),-1.0, 1.0);
  PixPj = cross(P1,P0);
  if (dot(PixPj,PixPj) > 0.0)
    PixPj = normalize(PixPj);
  aop += acos(PidotPj) * dot(N, PixPj);
  
  PidotPj = clamp(dot(P1,P2),-1.0, 1.0);
  PixPj = cross(P2,P1);
  if (dot(PixPj,PixPj) > 0.0)
    PixPj = normalize(PixPj);
  aop += acos(PidotPj) * dot(N, PixPj);

  PidotPj = clamp(dot(P2,P0),-1.0, 1.0);
  PixPj = cross(P0,P2);
  if (dot(PixPj,PixPj) > 0.0)
    PixPj = normalize(PixPj);
  aop += acos(PidotPj) * dot(N, PixPj);
  
  aop *= INVERSE_2_PI;

  return aop;
}


float AOp4(in vec3 N, in vec3 P0, in vec3 P1, in vec3 P2, in vec3 P3)
{
  float aop = 0.0;
  float PidotPj;
  vec3 PixPj;

  PidotPj = clamp(dot(P0,P1),-1.0, 1.0);
  PixPj = cross(P1,P0);
  if (dot(PixPj,PixPj) > 0.0)
    PixPj = normalize(PixPj);
  aop += acos(PidotPj) * dot(N, PixPj);

  PidotPj = clamp(dot(P1,P2),-1.0, 1.0);
  PixPj = cross(P2,P1);
  if (dot(PixPj,PixPj) > 0.0)
    PixPj = normalize(PixPj);
  aop += acos(PidotPj) * dot(N, PixPj);

  PidotPj = clamp(dot(P2,P3),-1.0, 1.0);
  PixPj = cross(P3,P2);
  if (dot(PixPj,PixPj) > 0.0)
    PixPj = normalize(PixPj);
  aop += acos(PidotPj) * dot(N, PixPj);

  PidotPj = clamp(dot(P3,P0),-1.0, 1.0);
  PixPj = cross(P0,P3);
  if (dot(PixPj,PixPj) > 0.0)
    PixPj = normalize(PixPj);
  aop += acos(PidotPj) * dot(N, PixPj);
  
  aop *= INVERSE_2_PI;

  return aop;
}


void main()
{
  float delta = oDelta;
  vec2 uv = (gl_FragCoord.xy) / screenSize;
  float invD = 1.0 / delta;
  
  // .rgb is face normal, .a is linear depth
  vec4 Nd = texture2D(normalTex, uv);
  // If alpha is 1.0, then assume this fragment has no prevous position/normal
  if (Nd.a == 1.0)
    discard;

  vec3 N = Nd.rgb;
  vec3 X = texture2D(positionTex, uv).rgb;

  //N.z = N.z * 5.0;
  //N = normalize(N);

  //color = vec3(1.0,1.0,1.0);

  vec3 P[4];
  P[0] = p0 - X;
  P[1] = p1 - X;
  P[2] = p2 - X;
  P[3] = vec3(0.0);
  
  // Clip to positive half space of tangent plane at this pixel
  float dots[3];
  float pLen;
  int pdots = 3;
  dots[0] = dot(P[0], N);
  pLen = dot(P[0],P[0]);
  if (pLen > 0.0)
    dots[0] *= inversesqrt(pLen);
  if (dots[0] < DOT_EPSILON)
  {
    pdots -= 1;
  }

  dots[1] = dot(P[1], N);
  pLen = dot(P[1],P[1]);
  if (pLen > 0.0)
    dots[1] *= inversesqrt(pLen);
  if (dots[1] < DOT_EPSILON)
  {
    pdots -= 1;
  }

  dots[2] = dot(P[2], N);
  pLen = dot(P[2],P[2]);
  if (pLen > 0.0)
    dots[2] *= inversesqrt(pLen);
  if (dots[2] < DOT_EPSILON)
  {
    pdots -= 1;
  }

  // pdots be the number of positive vertices
  if (pdots == 0) 
    discard;

  // Compute g(x) for this fragment
  float g = 1.0;      // init this to average alpha of the polygon if its available
  g *= max(0.0, min(1.0, (1.0 - dot(P[0], normalize(m0))*invD)));
  g *= max(0.0, min(1.0, (1.0 - dot(P[1], normalize(m1))*invD)));
  g *= max(0.0, min(1.0, (1.0 - dot(P[2], normalize(m2))*invD)));
  g *= max(0.0, min(1.0, (1.0 - dot(P[0], normalize(m3))*invD)));

  if (g == 0.0)
    discard;

  bool isTriangle = true;
  if (pdots == 1) 
  {    
    // Clip into triangle
    vec3 edge;
    float t;
    if (dots[0] >= DOT_EPSILON)
    {
     // color = vec3(1.0,0.0,0.0);
      edge = P[0] - P[1];
      t = (-dot(P[1], N)) / (EPSILON + dot(edge, N));
      P[1] += (t * edge);
      edge = P[0] - P[2];
      t = (-dot(P[2],N)) / (EPSILON + dot(edge, N));
      P[2] += (t * edge);
    } 
    else if (dots[1] >= DOT_EPSILON)
    {
      //color = vec3(1.0,1.0,0.0);

      edge = P[1] - P[0];
      t = (-dot(P[0],N)) / (EPSILON + dot(edge, N));
      P[0] += (t * edge);
      edge = P[1] - P[2];
      t = (-dot(P[2],N)) / (EPSILON + dot(edge, N));
      P[2] += (t * edge);
    } 
    else // if (dots[2] >= 0.0)
    {
      //color = vec3(1.0,0.0,1.0);
      edge = P[2] - P[0];
      t = (-dot(P[0],N)) / (EPSILON + dot(edge, N));
      P[0] += (t * edge);
      edge = P[2] - P[1];
      t = (-dot(P[1],N)) / (EPSILON + dot(edge, N));
      P[1] += (t * edge);
    }
  }

  if (pdots == 2) 
  {
    isTriangle = false;    
    // Clip into quad
    vec3 edge;
    float t;
    if (dots[0] < DOT_EPSILON)
    {
      //color = vec3(0.0,1.0,1.0);
      edge = P[2] - P[0];
      t = (-dot(P[0],N)) / (EPSILON + dot(edge, N));
      P[3] = P[0] + (t * edge);
      edge = P[1] - P[0];
      t = (-dot(P[0],N)) / (EPSILON + dot(edge, N));
      P[0] += (t * edge);
    }
    else if (dots[1] < DOT_EPSILON)
    {
      //color = vec3(0.0,0.0,1.0);
      P[3] = P[0];
      edge = P[0] - P[1];
      t = (-dot(P[1],N)) / (EPSILON + dot(edge, N));
      P[0] = P[1] + (t * edge);
      edge = P[2] - P[1];
      t = (-dot(P[1],N)) / (EPSILON + dot(edge, N));
      P[1] += (t * edge);
    } 
    else 
    {
      //color = vec3(0.5,0.5,0.5);
      edge = P[0] - P[2];
      t = (-dot(P[2],N)) / (EPSILON + dot(edge, N));
      P[3] = P[2] + (t * edge);
      edge = P[1] - P[2];
      t = (-dot(P[2],N)) / (EPSILON + dot(edge, N));
      P[2] += (t * edge);
    }
  }
 
  if (dot(P[0],P[0]) > 0.0)
    P[0] = normalize(P[0]);
  if (dot(P[1],P[1]) > 0.0)
    P[1] = normalize(P[1]);
  if (dot(P[2],P[2]) > 0.0)
    P[2] = normalize(P[2]);
  if (dot(P[3],P[3]) > 0.0)
    P[3] = normalize(P[3]);
  
  float AOp = 0.0;
  if (isTriangle)
  {
    AOp = AOp3(N, P[0], P[1], P[2]);
  } else {
    AOp = AOp4(N, P[0], P[1], P[2], P[3]);
  }

  float accessibility = clamp(AOp, 0.0, 1.0);
  
  gl_FragData[0] = vec4(N, accessibility*g);
  //vec3 foo = vec3(g);
  //gl_FragData[0] = vec4(foo, 1.0);
}


