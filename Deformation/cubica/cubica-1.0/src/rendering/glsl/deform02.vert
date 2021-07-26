/////////////////////////////////////////////////
// OpenGL vertex shader which does subspace tet 
// mesh deformation given a surfaceU matrix 
// (packed into texture Umatrix), and vector q,
// where q is Q_SIZE*4 floats, and Umatrix has
// texel width Q_STEP.
/////////////////////////////////////////////////

//#define Q_SIZE 5
//#define Q_SIZE_F 5.0
//#define Q_STEP 90.0 

uniform sampler2D Umatrix;
uniform vec4 q[Q_SIZE];
varying float Depth;

void main()
{
	// Transform from object space to eye space
	//vec4 P = gl_ModelViewMatrix * gl_Vertex;

  // Note: Currently, object space _is_ world space, so we just
  // add x to it then transform to eye space. However, if
  // more transforms are used client side, you'll need to
  // pass in the inverse modelview matrix to go from
  // object -> eye -> world -> eye (or something like that).
  
  int i;
  float offsetMultiplier = 0.0;
  vec2 offsetVector = vec2((1.0 / Q_STEP), 0.0);
  //vec4 Ux, Uy, Uz;
  vec4 uVal;
  vec2 uBasisLoc = gl_MultiTexCoord0.st;
  float unconstrained = gl_MultiTexCoord1.p;
  vec3 xVal = vec3(0.0);

/*  for (i = 0; i < Q_SIZE; i++)
  {
    Ux = texture2D(Umatrix, uBasisLoc + (offsetMultiplier * offsetVector));
    Uy = texture2D(Umatrix, uBasisLoc + ((Q_SIZE_F+offsetMultiplier) * offsetVector));
    Uz = texture2D(Umatrix, uBasisLoc + ((Q_SIZE_F+Q_SIZE_F+offsetMultiplier) * offsetVector));
    offsetMultiplier += 1.0; // */

 /*   xVal.x += dot(Ux, q[i]);
    xVal.y += dot(Uy, q[i]);
    xVal.z += dot(Uz, q[i]);// */
 /*   xVal.x += dot(Ux, vec4(0.01));
    xVal.y += dot(Uy, vec4(0.01));
    xVal.z += dot(Uz, vec4(0.01)); // */
 // }

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

  vec4 P = gl_Vertex;
  vec3 N = abs(gl_Normal);
  vec4 C = vec4(N, 1.0);
  //vec4 C = texture2D(Umatrix, uBasisLoc + (0.0 * offsetVector));

  P += vec4(xVal, 0.0);
  //vec4 P = gl_ModelViewMatrix * (gl_Vertex + vec4(xVal, 0.0));
  P = gl_ModelViewMatrix * P;

  // Do projection to get the linear eye depth
	vec2 EyePos = (gl_ProjectionMatrix * P).zw;
  Depth = EyePos.x / EyePos.y;
  Depth *= 100.0 * abs(EyePos.x);

  gl_FrontColor = C;
  
  gl_Position = gl_ProjectionMatrix * P;
}

