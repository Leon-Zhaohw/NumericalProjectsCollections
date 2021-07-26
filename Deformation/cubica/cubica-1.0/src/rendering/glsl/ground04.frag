uniform sampler2D normalTex;
uniform sampler2D positionTex;
uniform vec2 screenSize;
uniform float delta;

varying vec4 Origin;
varying vec3 Down;

void main()
{
  vec2 uv = (gl_FragCoord.xy) / screenSize;
   
  // .rgb is face normal, .a is linear depth
  vec4 Nd = texture2D(normalTex, uv);
  // If alpha is 1.0, then assume this fragment has no prevous position/normal
  if (Nd.a == 1.0)
    discard;

  vec3 N = Nd.rgb;
  vec3 X = texture2D(positionTex, uv).rgb;
  vec3 P0 = X - Origin.xyz;
  
  //float h = 2.0 * length(X.y - Origin.y);
  float h = 1.0 - dot(P0, -Down);
 
  float shade = 0.4 * (1.0 + dot(Down, N));
 
  //vec3 color;
  //color.rgb = vec3(shade * clamp(h, 0.0, 0.0));
  
  gl_FragColor = vec4(vec3(shade), 0.0);
  //gl_FragColor = vec4(vec3(0.0), shade);

}
