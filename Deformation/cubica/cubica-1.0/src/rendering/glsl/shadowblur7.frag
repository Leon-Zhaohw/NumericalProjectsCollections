/////////////////////////////////////////////////
// Separable gaussian blur fragment shader (of 7x7 blur)
/////////////////////////////////////////////////

uniform sampler2D texture;
uniform float texDim;

void main()
{
	vec2 offset = vec2((1.0 / texDim), 0.0);	
	vec4 sum = vec4(0.0);
	
	sum += texture2D(texture, gl_TexCoord[0].st - 3.0 * offset);
	sum += 10.0 * texture2D(texture, gl_TexCoord[0].st - 2.0 * offset);
	sum += 40.0 * texture2D(texture, gl_TexCoord[0].st - offset);
	sum += 64.0 * texture2D(texture, gl_TexCoord[0].st);
	sum += 40.0 * texture2D(texture, gl_TexCoord[0].st + offset);
	sum += 10.0 * texture2D(texture, gl_TexCoord[0].st + 2.0 * offset);
	sum += texture2D(texture, gl_TexCoord[0].st + 3.0 * offset);
	
	sum = sum / 166.0;
	
	gl_FragColor = sum;
}
