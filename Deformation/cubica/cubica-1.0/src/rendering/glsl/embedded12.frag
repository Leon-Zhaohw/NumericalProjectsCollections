/////////////////////////////////////////////////
// OpenGL fragment shader for rendering
// embedded tet mesh WITH SHADOWS
/////////////////////////////////////////////////

uniform sampler2D shadowTex;
uniform float minVariance;
uniform float reduceBleed;

varying vec3 Light;
varying vec3 Normal;
varying vec3 Pos;

//#define PCF_NUM_SAMPLES 16
#define PCF_NUM_SAMPLES 4
#define EXP_FACTOR_K  80.0

vec2 PoissonDisk[PCF_NUM_SAMPLES];

//float Rotation;

// Sample the given shadow map
float LookUp(in vec4 ShadowCoord, in sampler2D Shadowmap, in vec2 Offset, in float zReceiver)
{
  vec4 Coord = ShadowCoord + vec4(Offset.x, Offset.y, 0.0, 0.0) * 0.5;
  float Depth = texture2DProj(Shadowmap, Coord).r;
  float Shade;
 // if (Depth < min(1.0, zReceiver))
 //   Depth = 0.0;
 // else
 //   Depth = 1.0;
 
 // float edge1 = min(1.0, zReceiver) - 0.05;
 // float edge2 = edge1 + 0.1;
 // Depth = smoothstep(edge1, edge2,Depth);
 // Depth = step(min(1.0, zReceiver), Depth);

  //Shade = exp2(EXP_FACTOR_K*Depth) * exp2(-EXP_FACTOR_K*zReceiver);
  //Depth = log2(Depth) / EXP_FACTOR_K;
  
  Shade = Depth * exp2(-EXP_FACTOR_K*zReceiver);

  if (Shade > 1.0)
    Shade = 1.0;

  return Shade;
}

// Shade at the given coordinate using variance shadow mapping
float ShadeVar(in vec4 ShadowCoord, in sampler2D Shadowmap, in float zReceiver)
{
  vec3 Depth = texture2DProj(Shadowmap, ShadowCoord).rgb;
	float Shade;

  float p;
  if (Depth.g < min(1.0, zReceiver))
    p = 0.0;
	else
	  p = 1.0;
		
	float Variance = Depth.b - (Depth.g*Depth.g);
	Variance = max(Variance, minVariance);
	
	float d = zReceiver - Depth.g;
	float p_max = Variance / (Variance + d*d);
	
	p_max = max(p, p_max);

  Shade = clamp((p_max - reduceBleed) / (1.0 - reduceBleed), 0.0, 1.0);

  return Shade;
}

// Shade at the given coordinate using exponential shadow mapping (with 2x2 PCF correction)
float ShadeExp(in vec4 ShadowCoord, in sampler2D Shadowmap, in float zReceiver, in float filterRadiusUV)
{
  vec3 Depth = texture2DProj(Shadowmap, ShadowCoord).rgb;
	float Shade = Depth.r * exp2(-EXP_FACTOR_K*zReceiver);

  if (Shade > 1.02)
  {
    Shade = 1.0;
  	for (int i = 0; i < PCF_NUM_SAMPLES; ++i)
	  {
		  vec2 offset = PoissonDisk[i] * filterRadiusUV;
  		Shade += LookUp(ShadowCoord, Shadowmap, offset, zReceiver);
    }
    Shade = Shade / (1.0 + float(PCF_NUM_SAMPLES));
  }

  return Shade;
}

// Shade at the given coordinate using exponential shadow mapping, fading to variance
// shadow mapping for corrections
float ShadeExpVar(in vec4 ShadowCoord, in sampler2D Shadowmap, in float zReceiver)
{
  vec3 Depth = texture2DProj(Shadowmap, ShadowCoord).rgb;
	float Shade = Depth.r * exp2(-EXP_FACTOR_K*zReceiver);

  float p;
	if (Depth.g < min(1.0, zReceiver))
    p = 0.0;
	else
	  p = 1.0;
		
	float Variance = Depth.b - (Depth.g*Depth.g);
	Variance = max(Variance, minVariance);
	
	float d = zReceiver - Depth.g;
	float p_max = Variance / (Variance + d*d);
	
	p_max = max(p, p_max);

  // Variance shadows suffer from light leaks, and so over estimate darkness
  // so reduce it. So lets use ESM to darken those areas, and reduce the darkening in general.
  Shade = min(Shade, clamp((p_max - reduceBleed) / (1.0 - reduceBleed), 0.0, 1.0));
  
  return Shade;
}

void main()
{ 
	/*PoissonDisk[0] = vec2(-0.94201624, -0.39906216);
	PoissonDisk[1] = vec2(0.94558609, -0.76890725);
	PoissonDisk[2] = vec2(-0.094184101, -0.92938870);
	PoissonDisk[3] = vec2(0.34495938, 0.29387760);
	PoissonDisk[4] = vec2(-0.91588581, 0.45771432);
	PoissonDisk[5] = vec2(-0.81544232, -0.87912464);
	PoissonDisk[6] = vec2(-0.38277543, 0.27676845);
	PoissonDisk[7] = vec2(0.97484398, 0.75648379);
	PoissonDisk[8] = vec2(0.44323325, -0.97511554);
	PoissonDisk[9] = vec2(0.53742981, -0.47373420);
	PoissonDisk[10] = vec2(-0.26496911, -0.41893023);
	PoissonDisk[11] = vec2(0.79197514, 0.19090188);
	PoissonDisk[12] = vec2(-0.24188840, 0.99706507);
	PoissonDisk[13] = vec2(-0.81409955, 0.91437590);
	PoissonDisk[14] = vec2(0.19984126, 0.78641367);
	PoissonDisk[15] = vec2(0.14383161, -0.14100790);*/
  
  PoissonDisk[0] = vec2( 0.25, 0.75);
	PoissonDisk[1] = vec2(-0.75, 0.25);
	PoissonDisk[2] = vec2(-0.25,-0.75);
	PoissonDisk[3] = vec2( 0.75,-0.25);

	vec4 ShadowCoord;
	float Shadow;
	vec3 Mask;

  vec3 L = normalize(Light);
  vec3 N = normalize(Normal);
  float NdotL = max(dot(N, L), 0.0);
  
  ShadowCoord = gl_TextureMatrix[3] * vec4(Pos.xyz, 1.0);
  
  //Shadow = ShadeExpVar(ShadowCoord, shadowTex, (ShadowCoord.z / ShadowCoord.w));
  Shadow = ShadeVar(ShadowCoord, shadowTex, (ShadowCoord.z / ShadowCoord.w));
  //Shadow = ShadeExp(ShadowCoord, shadowTex, (ShadowCoord.z / ShadowCoord.w), 0.01);
  //Shadow = LookUp(ShadowCoord, shadowTex, vec2(0.0), (ShadowCoord.z / ShadowCoord.w));

  //Shadow = ShadowCoord.z / ShadowCoord.w;
  //Shadow = Shadow - 0.01;
  
  //float Depth = texture2DProj(shadowTex, ShadowCoord).r;
  //if (Depth < min(1.0, Shadow))
  //  Depth = 0.0;
  //else
  //  Depth = 1.0;
 
 // Mask = min(Shadow, NdotL);
  /*if (Shadow < 0.0)
    Mask = vec3(1.0, 0.0, 0.0);
  else*/
    Mask = vec3(Shadow * NdotL);
       
  Mask = (Mask + vec3(0.1)) * 0.91;
 
  gl_FragColor = vec4(Mask, 1.0);
  //gl_FragColor = vec4(1.0, Shadow, Shadow, 1.0);

}

