
uniform sampler2D defLowTex;
uniform sampler2D defHiTex;
uniform sampler2D aovTex;
uniform vec2 screenSize;

//#define _USE_FILTER_ 1

#define EPSILON 0.01
#define WEIGHT_NORM_CALC  clamp(pow(dot(hi.xyz, n.xyz), 32.0), EPSILON, 1.0)
#define WEIGHT_Z_CALC     clamp(pow(EPSILON / (EPSILON + abs(lo.a - hi.a)), 8.0), 0.0, 1.0)
//#define WEIGHT_Z_CALC_OLD (1.0 / (EPSILON + abs(lo.a - hi.a)))
#define NUM_FILTER_SAMPLES 12
#ifndef SAMPLE_RADIUS
#define SAMPLE_RADIUS 4.0
#endif

// aovLowTex Sample
vec4 lo;
// defHiTex Sample
vec4 hi;

// Weights
float wt;         // total
float wn;         // normal
float wz;         // depth
float wi;

vec2 vTaps[12];

void main()
{
  vTaps[0] = vec2(-0.326212,-0.40581)/screenSize;
  vTaps[1] = vec2(-0.840144,-0.07358)/screenSize;
  vTaps[2] = vec2(-0.695914, 0.457137)/screenSize;
  vTaps[3] = vec2(-0.203345, 0.620716)/screenSize;
  vTaps[4] = vec2( 0.96234, -0.194983)/screenSize;
  vTaps[5] = vec2( 0.473434,-0.480026)/screenSize;
  vTaps[6] = vec2( 0.519456, 0.767022)/screenSize;
  vTaps[7] = vec2( 0.185461,-0.893124)/screenSize;
  vTaps[8] = vec2( 0.507431, 0.064425)/screenSize;
  vTaps[9] = vec2( 0.89642,  0.412458)/screenSize;
  vTaps[10] = vec2(-0.32194, -0.932615)/screenSize;
  vTaps[11] = vec2(-0.791559,-0.59771)/screenSize;

  vec2 uv = gl_FragCoord.xy / screenSize;
  //vec2 offset = vec2(2.0) / screenSize;
  float sample;
  //vec3 color = vec3(0.0);
  float color = 0.0;
  wn = wz = 1.0;
  wt = EPSILON;
  vec3 n;
  float nLen;
  
 
  hi = texture2D(defHiTex, uv);
  nLen = length(hi.xyz);
  if (nLen > 0.0)
    hi.xyz /= nLen;

  // Center pixel
  sample = texture2D(aovTex, uv).a;

#ifdef _USE_FILTER_
  lo = texture2D(defLowTex, uv);
  n = lo.xyz;
  nLen = length(n);
  if (nLen > 0.0)
    n = normalize(n);

  wn = WEIGHT_NORM_CALC;
  wz = WEIGHT_Z_CALC;
#endif

  wi = wn * wz;
  wt = wt + wi;
  color += (sample * wi); 

#ifdef _USE_FILTER_
  int i;
  for (i = 0; i < NUM_FILTER_SAMPLES; i++)
  {
    sample = texture2D(aovTex, uv + SAMPLE_RADIUS*vTaps[i]).a;
    lo = texture2D(defLowTex, uv + SAMPLE_RADIUS*vTaps[i]);
    n = lo.xyz;
    nLen = length(n);
    if (nLen > 0.0)
      n = normalize(n);

    wn = WEIGHT_NORM_CALC;
    wz = WEIGHT_Z_CALC;
    wi = wn * wz * 0.35;
    wt = wt + wi;
    color += (sample * wi);
  }
#endif
  
  color = color / wt;

  //lo = texture2D(aovTex, uv);
  //color = texture2D(aovTex, uv).a;
  //color = log(1.0 + 10.0 * length(hi.xyz - lo.xyz));

  //gl_FragData[0] = vec4(vec3(color), 1.0);
  gl_FragData[0] = vec4(vec3(sample), 1.0);
  //gl_FragData[0] = vec4(vec3(texture2D(aovTex, uv).a), 1.0);

}

