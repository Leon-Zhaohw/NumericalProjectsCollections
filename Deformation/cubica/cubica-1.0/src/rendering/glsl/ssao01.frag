/////////////////////////////////////////////////
// Fragment shader to do screen space ambient occlusion.
// Expects texture to be rgb=normal, a=linear-depth
/////////////////////////////////////////////////
// Originally got this from:
// http://www.gamerendering.com/category/lighting/ssao-lighting/
/////////////////////////////////////////////////

uniform sampler2D texture;
#define NUM_SAMPLES 16

const float totStrength = 1.38;
const float strength = 0.02;
const float offset = 18.0;
const float falloff = 0.00000002;
const float rad = 0.02;
const float invSamples = 1.0/16.0;

vec3 pSphere[NUM_SAMPLES];

void main()
{
  pSphere[0]  = vec3(0.53812504, 0.18565957, -0.43192);
  pSphere[1]  = vec3(0.13790712, 0.24864247, 0.44301823);
  pSphere[2]  = vec3(0.33715037, 0.56794053, -0.005789503);
  pSphere[3]  = vec3(-0.6999805, -0.04511441, -0.0019965635);
  pSphere[4]  = vec3(0.06896307, -0.15983082, -0.85477847);
  pSphere[5]  = vec3(0.056099437, 0.006954967, -0.1843352);
  pSphere[6]  = vec3(-0.014653638, 0.14027752, 0.0762037);
  pSphere[7]  = vec3(0.010019933, -0.1924225, -0.034443386);
  pSphere[8]  = vec3(-0.35775623, -0.5301969, -0.43581226);
  pSphere[9]  = vec3(-0.3169221, 0.106360726, 0.015860917);
  pSphere[10] = vec3(0.010350345, -0.58698344, 0.0046293875);
  pSphere[11] = vec3(-0.08972908, -0.49408212, 0.3287904);
  pSphere[12] = vec3(0.7119986, -0.0154690035, -0.09183723);
  pSphere[13] = vec3(-0.053382345, 0.059675813, -0.5411899);
  pSphere[14] = vec3(0.035267662, -0.063188605, 0.54602677);
  pSphere[15] = vec3(-0.47761092, 0.2847911, -0.0271716);
		
  vec2 uv = gl_TexCoord[0].st;

  vec4 currentPixelSample = texture2D(texture,uv);
 
  float currentPixelDepth = currentPixelSample.a;
 
  // current fragment coords in screen space
  vec3 epos = vec3(uv.xy,currentPixelDepth);
  // get the normal of current fragment
  vec3 norm = currentPixelSample.xyz;
 
  float bl = 0.0;
  // adjust for the depth ( not sure if this is good..)
  float radD = rad/currentPixelDepth;
 
  vec3 ray, se, occNorm;
  float occluderDepth, depthDifference;
  float occlusion;
 
  for(int i=0; i<NUM_SAMPLES; ++i)
  {
    // get a vector (randomized inside of a sphere with radius 1.0) from a texture and reflect it
    ray = radD*pSphere[i];
   
    // if the ray is outside the hemisphere then change direction
    se = epos + sign(dot(ray,norm) )*ray;
 
    // get the depth of the occluder fragment
    vec4 occluderFragment = texture2D(texture,se.xy);
 
    // get the normal of the occluder fragment
    occNorm = occluderFragment.xyz;
 
    // if depthDifference is negative = occluder is behind current fragment
    depthDifference = currentPixelDepth-occluderFragment.a;
 
    // calculate the difference between the normals as a weight
 
    // the falloff equation, starts at falloff and is kind of 1/x^2 falling

    occlusion = step(falloff, depthDifference);
    occlusion *= 1.0-dot(occNorm,norm);
    occlusion *= 1.0-smoothstep(falloff,strength,depthDifference);

    bl += occlusion;
  }
 
  // output the result
  float ao = 1.0-totStrength*bl*invSamples;
  
  gl_FragColor = vec4(vec3(1.0), ao);
}
