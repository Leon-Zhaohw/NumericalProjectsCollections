/////////////////////////////////////////////////
// OpenGL fragment shader for rendering
// embedded tet mesh WITH SHADOWS
// Used in: _meshRenderWithShadowShader
/////////////////////////////////////////////////

uniform sampler2D shadowTex;
uniform sampler2D aovTex;
uniform float minVariance;
uniform float reduceBleed;

varying vec3 Light;
varying vec3 Normal;
varying vec4 Pos;

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

void main()
{ 
	vec4 ShadowCoord;
	float Shadow = 1.0;
	float Mask;

  vec3 L = normalize(Light);
  vec3 N = normalize(Normal);
  vec3 P = normalize(Pos.xyz);

  // Half-Vector between eye and light vectors (for specular calculation)
  vec3 halfVector = normalize(L + P);

  // cosine of the angle between the normal and the light vectors
  float NdotL = max(dot(N, L), 0.0);
  // cosine of the angle between the normal and the half-vector
  float NdotHV = max(dot(N, halfVector), 0.0);
  // Power factor
  float pf = 0.0;
  // Shininess, used for the specular exponent. 
  float shininess = max(gl_FrontMaterial.shininess, 10.0);

  // Specular doesn't show up if the surface is completely in shadow
  if (NdotL > 0.0)
    pf = pow(NdotHV, shininess);
  
  ShadowCoord = gl_TextureMatrix[4] * Pos;
  vec2 shadowOutside = ShadowCoord.xy / ShadowCoord.w;
  if (abs(shadowOutside.x) < 1.0 && abs(shadowOutside.y) < 1.01)
  //if (length(shadowOutside) <= 1.0)
  {
    Shadow = ShadeVar(ShadowCoord, shadowTex, (ShadowCoord.z / ShadowCoord.w));
  }

  // Mask is the total shading on the pixel, due to shadowing and surface orientation
  Mask = Shadow * NdotL;
  
  // Add all light and surface color contributions
  vec3 color = gl_LightSource[0].ambient.rgb * (1.0 - Mask) * gl_FrontMaterial.ambient.rgb;
  color += gl_LightSource[0].diffuse.rgb * Mask * gl_FrontMaterial.diffuse.rgb;
  color += gl_LightSource[0].specular.rgb * pf * gl_FrontMaterial.specular.rgb;
 
  gl_FragColor = vec4(color, 1.0);

}


