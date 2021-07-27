//////////////////////////////////////////////////////////////////////
// This file is part of Closest Point Turbulence.
// 
// Closest Point Turbulence is free software: you can redistribute it 
// and/or modify it under the terms of the GNU General Public License 
// as published by the Free Software Foundation, either version 3 of 
// the License, or (at your option) any later version.
// 
// Closest Point Turbulence is distributed in the hope that it will 
// be useful, but WITHOUT ANY WARRANTY; without even the implied 
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Closest Point Turbulence. 
// If not, see <http://www.gnu.org/licenses/>.
// 
// Copyright 2013 Theodore Kim and Nils Thuerey
// 
//////////////////////////////////////////////////////////////////////

#include <cmath>

#include "glvu.h"
#include <VEC3.h>
#include "FIELD_3D.h"
#include "BOX.h"
#include "TIMER.h"

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/gl.h> // OpenGL itself.
#include <GL/glu.h> // GLU support library.
#include <GL/glut.h> // GLUT support library.
#endif

#include <iostream>
#include <IWAVE_3D.h>

using namespace std;

// OGL preview
GLVU glvu;

// by default have the simulation run when the OGL window opens
bool animate = true;

// track which frame we're on
int simulationFrame = 0;
int endFrame = 301;

// object where the simulation heavy lifting happens
IWAVE_3D iWave;

// amount to up-res the existing liquid sim
int upresFactor = 2;

// number of iWave substeps per coarse simulation step
int substeps = 10;

// file paths
string inputPath("./data/houdini_input/");
string outputPath("./data/houdini_output/");
string texturePath(outputPath + "textures");
string surfacePath(outputPath + "surfaces");

// start on frame 4, since the paddle actually doesn't
// appear before then
int firstFrame = 4;
int houdiniFrame = firstFrame;

// paddle parameters
VEC3F initialCenter(0.00749427, 0.260715, -0.00250574);
float velocities[] = {3.02,3.02,3.02,3.02,3.02,2.6,2.125,2.09,2.01,1.955,1.915,1.875,1.855,1.83,1.81,1.79,1.775,1.765,1.75,1.74,1.73,1.725,1.715,1.715,1.71,1.705,1.705,1.705,1.705,1.705,1.71,1.715,1.72,1.73,1.735,1.74,1.75,1.76,1.77,1.78,1.79,1.8,1.815,1.83,1.84,1.85,1.865,1.88,1.89,1.905,1.92,1.935,1.945,1.96,1.97,1.985,1.995,2.005,2.015,2.025,2.035,2.05,2.06,2.065,2.075,2.085,2.095,2.105,2.11,2.12,2.13,2.135,2.145,2.15,2.16,2.17,2.175,2.185,2.19,2.2,2.205,2.215,2.22,2.23,2.235,2.245,2.25,2.255,2.265,2.27,2.275,2.285,2.29,2.3,2.305,2.31,2.315,2.325,2.33,2.335,2.34,2.35,2.355,2.36,2.365,2.37,2.38,2.385,2.39,2.395,2.4,2.405,2.415,2.42,2.425,2.43,2.435,2.445,2.45,2.455,2.46,2.465,2.47,2.475,2.48,2.485,2.49,2.495,2.505,2.51,2.515,2.52,2.525,2.5325,2.5375,2.5425,2.5475,2.5525,2.5575,2.5625,2.5675,2.5725,2.5775,2.5825,2.59,2.595,2.6,2.605,2.61,2.615,2.62,2.625,2.63,2.635,2.64,2.645,2.65,2.6525,2.6575,2.6625,2.6675,2.6725,2.675,2.68,2.685,2.69,2.695,2.7,2.7025,2.7075,2.7125,2.7175,2.72,2.725,2.73,2.735,2.7375,2.7425,2.7475,2.7525,2.755,2.76,2.7625,2.765,2.77,2.7725,2.7775,2.7825,2.785,2.79,2.7925,2.7975,2.8,2.805,2.81,2.8125,2.8175,2.8225,2.83,2.835,2.84,2.845,2.8525,2.8575,2.8625,2.8675,2.875,2.88,2.885,2.89,2.895,2.9,2.905,2.91,2.9125,2.9175,2.92,2.9225,2.925,2.9275,2.93,2.93,2.9325,2.9325,2.9325,2.9325,2.93,2.9275,2.925,2.9225,2.9175,2.915,2.91,2.9025,2.8975,2.89,2.8825,2.8725,2.8625,2.8525,2.8425};
float xScale = 0.26;
float yScale = 1;
float zScale = 0.51;
VEC3F translation(0,0,0.67);

///////////////////////////////////////////////////////////////////////
// GL and GLUT callbacks
///////////////////////////////////////////////////////////////////////
void glutDisplay()
{
  glvu.BeginFrame();
    glEnable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    iWave.surfaceMesh().drawSolidTexturePreview(4, 0.5);

    // draw the paddle
    BOX paddle(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5);

    glColor4f(1,0,0,1);

    // draw the paddle itself
    glPushMatrix();
      float theta = -43.55;
      theta += (houdiniFrame - 3) * velocities[houdiniFrame];
      glRotatef(theta, 0,1,0);
      glTranslatef(translation[0], translation[1], translation[2]);
      glTranslatef(initialCenter[0], initialCenter[1], initialCenter[2]);

      glTranslatef(0,-0.5, 0);
      glScalef(xScale, yScale, zScale);
      paddle.draw();
    glPopMatrix();
  glvu.EndFrame();
}

///////////////////////////////////////////////////////////////////////
// step the simulation
///////////////////////////////////////////////////////////////////////
void stepEverything()
{
  TIMER functionTimer(__FUNCTION__);
  static int stepsSeen = 1;

  if (simulationFrame != 0 && simulationFrame % substeps == 0)
  {
    // write out the current simulation state
    TIMER fileWriteTimer("Dumping texture");
    iWave.writeHeight(texturePath, houdiniFrame);
    iWave.writeSurfaceMesh(surfacePath, houdiniFrame);
    houdiniFrame++;
  
    // only output the timings after each substep sequence
    TIMER::printTimingsPerFrame(stepsSeen);

    // load up the next Houdini frame
    iWave.loadNextHoudini12(houdiniFrame, inputPath, upresFactor);

    // seed turbulence according to curvature
    iWave.setSourceToHoudiniCurvature(true);

    // generate OGL textures for preview display
    iWave.refreshSurfaceTextures();
  }

  // if we're done, just bail
  if (houdiniFrame == endFrame)
    exit(0);

  // timestep the iWave simulation
  iWave.stepHoudiniFrozenCore();

  stepsSeen++;
  simulationFrame++;
}

///////////////////////////////////////////////////////////////////////
// animate and display new result
///////////////////////////////////////////////////////////////////////
void glutIdle()
{
  if (animate)
    stepEverything();

  glutPostRedisplay();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutKeyboard(unsigned char key, int x, int y)
{
  switch (key)
  {
    case 'a':
      animate = !animate;
      break;
    case ' ':
      stepEverything();
      break;
    case 'q':
      exit(0);
      break;
    default:
      break;
  }

  glvu.Keyboard(key,x,y);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutMouseClick(int button, int state, int x, int y)
{
  glvu.Mouse(button,state,x,y);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutMouseMotion(int x, int y)
{
  glvu.Motion(x,y);
}

//////////////////////////////////////////////////////////////////////////////
// open the GLVU window
//////////////////////////////////////////////////////////////////////////////
int glvuWindow()
{
  char title[] = "iWave 3D Viewer";
  glvu.Init(title,
            GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH,
            0, 0, 800, 800);

  glvuVec3f ModelMin(-10,-10,-10), ModelMax(10,10,10), 
        Eye(1.86418, 3.54351, 2.1454),  LookAtCntr(1.45312, 2.76776, 1.6666),  Up(-0.438467, 0.628726, -0.642223);

  float Yfov = 45;
  float Aspect = 1;
  float Near = 0.001f;
  float Far = 10.0f;
  glvu.SetAllCams(ModelMin, ModelMax, Eye, LookAtCntr, Up, Yfov, Aspect, Near, Far);

  glvuVec3f worldCenter;
  worldCenter[0] = iWave.surfaceMesh().vertexMean()[0];
  worldCenter[1] = iWave.surfaceMesh().vertexMean()[1];
  worldCenter[2] = iWave.surfaceMesh().vertexMean()[2];

  glvu.SetWorldCenter(worldCenter);

  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  GLfloat lightZeroPosition[] = {10.0, 4.0, 10.0, 1.0};
  GLfloat lightZeroColor[] = {0.8, 1.0, 1.0, 1.0};
  glLightfv(GL_LIGHT0, GL_POSITION, lightZeroPosition);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, lightZeroColor);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_COLOR_MATERIAL);
  glShadeModel(GL_SMOOTH);

  glClearColor(0.1,0.1,0.1,0);
  glutDisplayFunc(&glutDisplay);
  glutIdleFunc(&glutIdle);
  glutKeyboardFunc(&glutKeyboard);
  glutMouseFunc(&glutMouseClick);
  glutMotionFunc(&glutMouseMotion);

  glutMainLoop();

  // Control flow will never reach here
  return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  if (argc < 2)
  {
    cout << " Usage: " << argv[0] << " <upres factor>" << endl;
    return 0;
  }

  // hardcode the number of substeps -- you can change this to get
  // different effects
  substeps = 10;
  cout << " Using " << substeps << " substeps " << endl;

  // hard code the damping factor to 0.2 -- you can change this
  // to get different effects
  double alpha = 0.2;

  // get the upres factor
  upresFactor = atoi(argv[1]);
  if (upresFactor < 2)
  {
    cout << upresFactor << " is not a valid upres factor! " << endl;
    return 0;
  }
  cout << " Using upres factor of " << upresFactor << endl;

  // create the output paths
  char buffer[256];
  sprintf(buffer, "%i", upresFactor);
  string upresString(buffer);
  sprintf(buffer, "%i", substeps);
  string substepsString(buffer);
  cout << " Using texture path: " << texturePath.c_str() << endl;
  cout << " Using surface path: " << surfacePath.c_str() << endl;
  string mkdir;
  mkdir = string("mkdir ") + outputPath;
  system(mkdir.c_str());
  mkdir = string("mkdir ") + texturePath;
  system(mkdir.c_str());
  mkdir = string("mkdir ") + surfacePath;
  system(mkdir.c_str());

  // create the actual iWave simulation
  int simulationRes = 100;
  int filterWidth = 15;
  iWave = IWAVE_3D(simulationRes, simulationRes, simulationRes, filterWidth);

  // set the damping alpha
  cout << " Using damping alpha: " << alpha << endl;
  iWave.alpha() = alpha;

  // load up a Houdini frame
  iWave.loadHoudini12(houdiniFrame, inputPath, upresFactor);

  // generate some OGL textures for display
  iWave.refreshSurfaceTextures();

  glutInit(&argc, argv);
  glvuWindow();
  return 1;
}
