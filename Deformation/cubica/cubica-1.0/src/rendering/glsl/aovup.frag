
uniform sampler2D aovLowTex;
uniform sampler2D defHiTex;
uniform vec2 screenSize;

#define EPSILON 0.01
#define NORMAL_WEIGHT_EXP 32.0

// aovLowTex Sample
vec4 lo;
// defHiTex Sample
vec4 hi;

// Weights
float wt;         // total
float wn;         // normal
float wz;         // depth
float wi;

void main()
{
  vec2 uv = gl_FragCoord.xy / screenSize;
  vec2 offset = vec2(-1.0) / screenSize;  
  float color = 0.0;
  wt = EPSILON;
  vec3 n;
  float nLen;
 
  hi = texture2D(defHiTex, uv);
  nLen = length(hi.xyz);
  if (nLen > 0.0)
    hi.xyz /= nLen;

  // Pixel 0,0
  lo = texture2D(aovLowTex, uv);
  /*n.xy = lo.gb;
  nLen = min(length(n.xy), 1.0);
  n.z = sqrt(1.0 - nLen);
  n = normalize(n);

  wn = clamp(pow(dot(hi.xyz, n.xyz), NORMAL_WEIGHT_EXP), EPSILON, 1.0);
  wz = 1.0 / (EPSILON + abs(lo.a - hi.a));
  //wi = wn * wz * 0.424889;
  wi = wn * wz * 0.25;
  wt = wt + wi;
  color += (lo.r * wi);
  //color = color + wi; // */

  // Pixel 1,0
  /*lo = texture2D(aovLowTex, uv + vec2(offset.s, 0.0));
  n.xy = lo.gb;
  nLen = min(length(n.xy), 1.0);
  n.z = sqrt(1.0 - nLen);
  n = normalize(n);

  wn = clamp(pow(dot(hi.xyz, n.xyz), NORMAL_WEIGHT_EXP), EPSILON, 1.0);
  wz = 1.0 / (EPSILON + abs(lo.a - hi.a));
  //wi = wn * wz * 0.212445;
  wi = wn * wz * 0.25;
  wt = wt + wi;
  color += (lo.r * wi); 
  
  //color = color + wi;
  //color = color + max(wi, 0.0); // */

  // Pixel 1,1
  /*lo = texture2D(aovLowTex, uv + vec2(offset.st));
  n.xy = lo.gb;
  nLen = min(length(n.xy), 1.0);
  n.z = sqrt(1.0 - nLen);
  n = normalize(n);
  
  wn = clamp(pow(dot(hi.xyz, n.xyz), NORMAL_WEIGHT_EXP), EPSILON, 1.0);
  wz = 1.0 / (EPSILON + abs(lo.a - hi.a));
  //wi = wn * wz * 0.150221;
  wi = wn * wz * 0.25;
  wt += wi;
  color += (lo.r * wi);
  //color += wi; // */

  // Pixel 0,1
  /*lo = texture2D(aovLowTex, uv + vec2(0.0, offset.t));
  n.xy = lo.gb;
  nLen = min(length(n.xy), 1.0);
  n.z = sqrt(1.0 - nLen);
  n = normalize(n);

  wn = clamp(pow(dot(hi.xyz, n.xyz), NORMAL_WEIGHT_EXP), EPSILON, 1.0);
  wz = 1.0 / (EPSILON + abs(lo.a - hi.a));
  //wi = wn * wz * 0.212445;
  wi = wn * wz * 0.25;
  wt += wi;
  color += (lo.r * wi);
  //color += wi; // */

  //color = color / wt;
  //gl_FragData[0] = vec4(color);
  gl_FragData[0] = vec4(lo.r);
}

