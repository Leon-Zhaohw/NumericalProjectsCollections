//////
// Shader to draw the ground plane with current color
/////

void main()
{ 
  gl_FrontColor = gl_Color;
	gl_Position = ftransform();
}



