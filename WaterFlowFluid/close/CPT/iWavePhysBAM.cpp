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
//////////////////////////////////////////////////////////////////////

#include <cmath>

#include <glvu.h>
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

int simulationFrame = 0;
int physBAMFrame = 2;

int endFrame = 201;
bool animate = true;

int simulationRes = 101;
int filterWidth = 15;

double target = 4.0;
double width = 0.4;

GLVU glvu;
IWAVE_3D iWave;

int upresFactor = 4;
int substeps = 5;
string path("./data/PhysBAM_input/");

string outputPath("./data/PhysBAM_output/");
string texturePath(outputPath + "textures");
string surfacePath(outputPath + "surfaces");

int skipSourcing = 4;

///////////////////////////////////////////////////////////////////////
// draw coordinate axes
///////////////////////////////////////////////////////////////////////
void drawAxes()
{
  // draw coordinate axes
  glPushMatrix();
  glTranslatef(-0.1f, -0.1f, -0.1f);
  glLineWidth(3.0f);
  glBegin(GL_LINES);
    // x axis is red
    glColor4f(10.0f, 0.0f, 0.0f, 1.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glColor4f(10.0f, 0.0f, 0.0f, 0.0f);
    glVertex3f(1.0f, 0.0f, 0.0f);

    // y axis is green 
    glColor4f(0.0f, 10.0f, 0.0f, 1.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glColor4f(0.0f, 10.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 1.0f, 0.0f);
    
    // z axis is blue
    glColor4f(0.0f, 0.0f, 10.0f, 1.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glColor4f(0.0f, 0.0f, 10.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 1.0f);
  glEnd();
  glLineWidth(1.0f);
  glPopMatrix();
}

///////////////////////////////////////////////////////////////////////
// GL and GLUT callbacks
///////////////////////////////////////////////////////////////////////
void glutDisplay()
{
  glvu.BeginFrame();
    glEnable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    iWave.surfaceMesh().drawSolidTexturePreview(4, 0.5);
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
    TIMER fileWriteTimer("Dumping texture");
    iWave.writeHeight(texturePath, physBAMFrame);
    iWave.writeSurfaceMesh(surfacePath, physBAMFrame);
    physBAMFrame++;
    TIMER::printTimingsPerFrame(stepsSeen);

    iWave.loadNextPouringFrozenCore(physBAMFrame, path, upresFactor);
    iWave.setSourceToPouringCurvature(true);
 
    iWave.refreshSurfaceTextures();
  }

  if (physBAMFrame == endFrame)
    exit(0);

  iWave.stepPouringFrozenCore();

  stepsSeen++;
  simulationFrame++;
}

///////////////////////////////////////////////////////////////////////
// animate and display new result
///////////////////////////////////////////////////////////////////////
void glutIdle()
{
  if (animate)
  {
    stepEverything();
  }
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
        Eye(0.722933, 1.71411, 0.239635),  LookAtCntr(0.33598, 0.802006, 0.104214),  Up(-0.480711, 0.324866, -0.814481);

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
  if (argc != 2)
  {
    cout << " Usage: " << argv[0] << " <upres factor>" << endl;
    return 0;
  }

  // get the upres factor
  upresFactor = atoi(argv[1]);
  if (upresFactor < 2)
  {
    cout << upresFactor << " is not a valid upres factor! " << endl;
    return 0;
  }
  cout << " Using upres factor of " << upresFactor << endl;

  double alpha = 0.2;
  cout << " Using texture path: " << texturePath.c_str() << endl;
  cout << " Using surface path: " << surfacePath.c_str() << endl;

  int simulationRes = 100;
  int filterWidth = 15;

  string mkdir;
  mkdir = string("mkdir ") + outputPath;
  system(mkdir.c_str());
  mkdir = string("mkdir ") + texturePath;
  system(mkdir.c_str());
  mkdir = string("mkdir ") + surfacePath;
  system(mkdir.c_str());

  // path to dump data to
  iWave = IWAVE_3D(simulationRes, simulationRes, simulationRes, filterWidth);

  iWave.alpha() = alpha;
  iWave.loadPhysBAM(physBAMFrame, path, upresFactor);
  iWave.setSourceToPouringCurvature(true);
  
  iWave.refreshSurfaceTextures();

  VEC3F center = iWave.height().center();
  VEC3F lengths = iWave.height().lengths();
  cout << " Simulation center: " << center << endl;
  cout << " Simulation lengths: " << lengths << endl;

  TIMER::printTimings();

  glutInit(&argc, argv);
  glvuWindow();
  return 1;
}
