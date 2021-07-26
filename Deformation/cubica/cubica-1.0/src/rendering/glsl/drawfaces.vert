
#define NORMAL_OFFSET 0.02
#define SCALE_FACTOR 0.2

uniform sampler2D tetMeshTex;
uniform vec2 offset;

varying vec4 Point0;
varying vec4 Point1;
varying vec4 Point2;

void main()
{
  Point0.xyz = texture2D(tetMeshTex, gl_Vertex.xy + offset).rgb;
  Point1.xyz = texture2D(tetMeshTex, gl_Vertex.zw + offset).rgb;
  Point2.xyz = texture2D(tetMeshTex, gl_MultiTexCoord0.xy + offset).rgb;

/*  vec3 n1 = Point1.xyz - Point0.xyz;
  vec3 n2 = Point2.xyz - Point0.xyz;
  vec3 norm = normalize(cross(n1, n2));

  vec3 midPt = (Point0.xyz + Point1.xyz + Point2.xyz) / 3.0;

  vec3 v = Point0.xyz - midPt;
  Point0.xyz += SCALE_FACTOR * v;

  v = Point1.xyz - midPt;
  Point1.xyz += SCALE_FACTOR * v;

  v = Point2.xyz - midPt;
  Point2.xyz += SCALE_FACTOR * v;

  Point0.xyz -= NORMAL_OFFSET * norm;
  Point1.xyz -= NORMAL_OFFSET * norm;
  Point2.xyz -= NORMAL_OFFSET * norm;*/ 

  Point0 = gl_ModelViewProjectionMatrix * vec4(Point0.xyz, 1.0);
  Point1 = gl_ModelViewProjectionMatrix * vec4(Point1.xyz, 1.0);
  Point2 = gl_ModelViewProjectionMatrix * vec4(Point2.xyz, 1.0);

  gl_Position = vec4(0.0);
}

