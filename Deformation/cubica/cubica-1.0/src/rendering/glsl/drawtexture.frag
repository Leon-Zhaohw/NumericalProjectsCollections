uniform sampler2D tex;

void main()
{
  vec3 texcolor = texture2D(tex, gl_TexCoord[0].st).rgb;
  gl_FragColor = vec4(texcolor, 1.0);
  //gl_FragColor = vec4(1.0);
}
