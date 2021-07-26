//////
// Shader to draw the ground plane 
/////

void main()
{ 
  
	gl_Position = gl_ModelViewMatrix * gl_Vertex;
}



