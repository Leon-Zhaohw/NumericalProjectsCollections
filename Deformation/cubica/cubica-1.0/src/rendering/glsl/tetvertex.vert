
/////////////////////////////////////////////////
// OpenGL vertex shader which does subspace tet 
// mesh deformation given a surfaceU matrix 
// (packed into texture Umatrix), and vector q,
// where q is Q_SIZE*4 floats, and Umatrix has
// texel width Q_STEP.
/////////////////////////////////////////////////

//#define Q_SIZE 5
//#define Q_STEP 90.0 

uniform sampler2D Umatrix;
uniform vec4 q[Q_SIZE];

void main()
{
	// Transform from object space to eye space
	//vec4 P = gl_ModelViewMatrix * gl_Vertex;
  
  int i;
  float offsetMultiplier = 0.0;
  vec2 offsetVector = vec2((1.0 / Q_STEP), 0.0);
  //vec4 Ux, Uy, Uz;
  vec4 uVal;
  vec2 uBasisLoc = gl_MultiTexCoord0.st;
  float unconstrained = gl_MultiTexCoord1.p;
  vec4 screenPos = vec4(gl_MultiTexCoord1.st, 0.0, 1.0);
  vec3 xVal = vec3(0.0);

  // Want to read Q_SIZE*3 texels, Q_SIZE for each _x coord.
  for (i = 0; i < Q_SIZE; i++)
  {
    uVal = texture2D(Umatrix, uBasisLoc + (offsetMultiplier * offsetVector));
    offsetMultiplier += 1.0;

    xVal.x += dot(uVal, q[i]) * unconstrained;
    //xVal.x += dot(uVal, vec4(0.01));
  }
  for (i = 0; i < Q_SIZE; i++)
  {
    uVal = texture2D(Umatrix, uBasisLoc + (offsetMultiplier * offsetVector));
    offsetMultiplier += 1.0;

    xVal.y += dot(uVal, q[i]) * unconstrained;
    //xVal.y += dot(uVal, vec4(0.01));
  }
  for (i = 0; i < Q_SIZE; i++)
  {
    uVal = texture2D(Umatrix, uBasisLoc + (offsetMultiplier * offsetVector));
    offsetMultiplier += 1.0;

    xVal.z += dot(uVal, q[i]) * unconstrained;
    //xVal.z += dot(uVal, vec4(0.01));
  }

  vec4 vertPos = gl_Vertex;
  //vec4 C = vec4(1.0, 1.0, 1.0, 1.0);
  //vec4 C = texture2D(Umatrix, uBasisLoc + (0.0 * offsetVector));

  vertPos.xyz += xVal;
  //vec4 vertPos = gl_ModelViewMatrix * (gl_Vertex + vec4(xVal, 0.0));
  //vertPos = gl_ModelViewMatrix * vertPos;

  //gl_FrontColor = gl_ProjectionMatrix * vertPos;
  gl_FrontColor = vertPos;

  //gl_FrontColor = vec4(1.0, 0.0, 0.0, 1.0);
  
  gl_Position = gl_ModelViewProjectionMatrix * screenPos;

  //gl_Position = gl_ModelViewProjectionMatrix * vec4(1.0, 1.0, 0.0, 1.0);
}

