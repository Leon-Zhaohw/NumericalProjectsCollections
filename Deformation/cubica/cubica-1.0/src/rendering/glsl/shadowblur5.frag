/////////////////////////////////////////////////
// Separable gaussian blur fragment shader (of 5x5 blur)
/////////////////////////////////////////////////

uniform sampler2D texture;
uniform float texDim;

void main()
{
	vec2 offset = vec2((1.0 / texDim), 0.0);	
	vec4 sum = vec4(0.0);
	
	sum += texture2D(texture, gl_TexCoord[0].st - 2.0 * offset);
  sum += 3.71 * texture2D(texture, gl_TexCoord[0].st - offset);
	sum += 5.86 * texture2D(texture, gl_TexCoord[0].st);
	sum += 3.71 * texture2D(texture, gl_TexCoord[0].st + offset);	
	sum += texture2D(texture, gl_TexCoord[0].st + 2.0 * offset);
	
	sum = sum / 15.28;
	
	gl_FragColor = sum;
}
