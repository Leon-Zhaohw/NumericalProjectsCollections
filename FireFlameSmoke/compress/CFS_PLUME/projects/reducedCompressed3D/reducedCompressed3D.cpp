/*
 This file is part of SSFR (Zephyr).
 
 Zephyr is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 Zephyr is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Zephyr.  If not, see <http://www.gnu.org/licenses/>.
 
 Copyright 2013 Theodore Kim
 */
#include "EIGEN.h"
#include <cmath>
#include "QUICKTIME_MOVIE.h"

#include <glvu.h>
#include <VEC3.h>
#include <iostream>
#if __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/freeglut.h>
#include <GL/glu.h>
#endif
#include "SUBSPACE_FLUID_3D_COMPRESSED_EIGEN.h"
#include "MATRIX.h"
#include "COMPRESSION.h"
#include "COMPRESSION_DATA.h"
#include "MATRIX_COMPRESSION_DATA.h"
#include "SIMPLE_PARSER.h"

using namespace std;

GLVU glvu;

// is the mouse pressed?
bool mouseClicked = false;

// what's the z coordinate of the last mouse click?
float clickZ;

// animate it this frame?
bool animate = true;

// currently capturing frames for a movie?
bool captureMovie = true;

// fluid being simulated
SUBSPACE_FLUID_3D_COMPRESSED_EIGEN* fluid = NULL;

// Quicktime movie to capture to
QUICKTIME_MOVIE movie;

void runOnce();
void runEverytime();
// helper function to write unique filenames
bool fileExists(const string& filename);
// helper function to write to quicktime movie
void writeToQuicktime();

// resolutions
int xRes = 0;
int yRes = 0;
int zRes = 0;

vector<VECTOR> snapshots;

// user configuration
string reducedPath;
string compressionPath;
string snapshotPath;
string previewReducedMovie;
int simulationSnapshots;

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

//////////////////////////////////////////////////////////////////////////////
// Translate screen space to world space
//
// Adapted from the NeHe page:
// http://nehe.gamedev.net/data/articles/article.asp?article=13
//////////////////////////////////////////////////////////////////////////////
VEC3F unproject(float x, float y, float z)
{
  GLint viewport[4];
	GLdouble modelview[16];
	GLdouble projection[16];

	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);
	glGetIntegerv(GL_VIEWPORT, viewport);

  double worldX, worldY, worldZ;
	gluUnProject(x, viewport[3] - y, z,
               modelview, projection, viewport, 
               &worldX, &worldY, &worldZ);

  return VEC3F(worldX, worldY, worldZ);
}

///////////////////////////////////////////////////////////////////////
// GL and GLUT callbacks
///////////////////////////////////////////////////////////////////////
void glutDisplay()
{
  glvu.BeginFrame();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);

    glPushMatrix();
      glTranslatef(0.5, 0.5, 0.5);
      /*
      /////////////////////////////////////////////////////////////////
      // take difference here
      
      auto density = fluid->density();
      auto ground_density = ground->density();
      auto diff = density - ground_density;
      diff.draw();
      diff.drawBoundingBox();
      */
      fluid->density().draw();
      fluid->density().drawBoundingBox();
      /////////////////////////////////////////////////////////////////
    glPopMatrix();
   

  glvu.EndFrame();
  // if we're recording a movie, capture a frame
  if (captureMovie) {
    movie.addFrameGL();
  }
}

///////////////////////////////////////////////////////////////////////
// animate and display new result
///////////////////////////////////////////////////////////////////////
void glutIdle()
{
  runEverytime();

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
    case 'q':
      exit(0);
      break;
    case 'm':
        // if we were already capturing a movie
        if (captureMovie)
        {
         // write out the movie
         movie.writeMovie("movie.mov");

        // reset the movie object
        movie = QUICKTIME_MOVIE();

       // stop capturing frames
       captureMovie = false;
        }
        else
        {
         cout << " Starting to capture movie. " << endl;
        captureMovie = true;
        }
       break; 
    default:
      break;
  }
  glvu.Keyboard(key,x,y);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutSpecial(int key, int x, int y)
{
  switch (key)
  {
    case GLUT_KEY_LEFT:
      break;
    case GLUT_KEY_RIGHT:
      break;
    case GLUT_KEY_UP:
      break;
    case GLUT_KEY_DOWN:
      break;
  }
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
  char title[] = "3D Viewer";
  glvu.Init(title,
            GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH,
            0, 0, 800, 800);

  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  GLfloat lightZeroPosition[] = {10.0, 4.0, 10.0, 1.0};
  GLfloat lightZeroColor[] = {1.0, 1.0, 1.0, 1.0};
  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, 1);
  glLightfv(GL_LIGHT0, GL_POSITION, lightZeroPosition);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, lightZeroColor);

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); 
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_DEPTH_TEST);
  glShadeModel(GL_SMOOTH);
  glClearColor(0,0,0,0);

  glutDisplayFunc(&glutDisplay);
  glutIdleFunc(&glutIdle);
  glutKeyboardFunc(&glutKeyboard);
  glutSpecialFunc(&glutSpecial);
  glutMouseFunc(&glutMouseClick);
  glutMotionFunc(&glutMouseMotion);

  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  // ADJ: modify the Eye vector to change the scene's perspective
  glvuVec3f ModelMin(-10,-10,-10), ModelMax(10,10,10), 
        Eye(1.5,0.5,2), LookAtCntr(0.5,0.5,0.5),  Up(0,1,0);

  float Yfov = 45;
  float Aspect = 1;
  float Near = 0.001f;
  float Far = 10.0f;
  glvu.SetAllCams(ModelMin, ModelMax, Eye, LookAtCntr, Up, Yfov, Aspect, Near, Far);

  glvuVec3f center(0.5, 0.5, 0.5);
  glvu.SetWorldCenter(center);

  glutMainLoop();

  // Control flow will never reach here
  return EXIT_SUCCESS;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
  // read in the cfg file
  if (argc != 2)
  {
    cout << " Usage: " << argv[0] << " *.cfg" << endl;
    return 0;
  }
  SIMPLE_PARSER parser(argv[1]);

  int xRes = parser.getInt("xRes", 48);
  int yRes = parser.getInt("yRes", 64);
  int zRes = parser.getInt("zRes", 48);
  reducedPath = parser.getString("reduced path", "./data/reduced.dummy/");
  compressionPath = parser.getString("compression path", "./data/reduced.dummy/");
  printf(" Using compression path: %s\n", compressionPath.c_str());
  snapshotPath = parser.getString("snapshot path", "./data/dummy/");
  previewReducedMovie = parser.getString("preview movie", "./movies/movie.mov");
  cout << "Reduced movie is written to: " << previewReducedMovie << endl;
  simulationSnapshots = parser.getInt("simulation snapshots", 20);
  Real vorticity = parser.getFloat("vorticity", 0);

  cout << " Using vorticity: " << vorticity << endl;

  unsigned int boundaries[6];
  boundaries[0] = parser.getInt("front", 1);
  boundaries[1] = parser.getInt("back", 1);
  boundaries[2] = parser.getInt("left", 1);
  boundaries[3] = parser.getInt("right", 1);
  boundaries[4] = parser.getInt("top", 0);
  boundaries[5] = parser.getInt("bottom", 0);

  string names[] = {"front", "back", "left", "right", "top", "bottom"};
  for (int x = 0; x < 6; x++)
  {
    cout << " Boundary on " << names[x].c_str() << "\tis set to " << flush;
    if (boundaries[x] == 0)
      cout << "Neumann " << endl;
    else
      cout << "Dirichlet " << endl;
  }

  double discardThreshold = parser.getFloat("discard threshold", 1e-9);
  cout << " Using discard threshold: " << discardThreshold << endl;

  bool usingIOP = parser.getBool("iop", 0);
  cout << "Using IOP: " << usingIOP << endl;

  bool fastPow = parser.getBool("fastPow", false);
  cout << " Fast pow: " << fastPow << endl;

  bool debug = parser.getBool("debug", 0);
  cout << "Debug: " << debug << endl;


	fluid = new SUBSPACE_FLUID_3D_COMPRESSED_EIGEN(xRes, yRes, zRes, reducedPath, &boundaries[0], usingIOP);
  fluid->setCompressionPath(compressionPath);
  fluid->loadReducedIOP(string(""));
  // For debugging, use loadReducedIOPAll here and stepMovingDebug in runEverytime.
  // fluid->loadReducedIOPAll(string(""));
  
  puts("finished loadReducedIOP");

  fluid->initCompressionData();
  fluid->fullRankPath() = snapshotPath;
  fluid->vorticityEps() = vorticity;


  // set the FIELD_3D static
  FIELD_3D::usingFastPow() = fastPow;
  
  // ground = new FLUID_3D_MIC(xRes, yRes, zRes, 0);
 
  TIMER::printTimings();
 
  glutInit(&argc, argv);
  glvuWindow();

  return 1;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

void runOnce()
{
  fluid->initCompressionData();
}

void runEverytime()
{
  if (animate)
  {
    static int step = 0;
    if (step == 0) { 
      fluid->writeCompressedErrorMatrixDims(simulationSnapshots);
      // fluid->loadSubspaceComparisonMatrix();
    }

    cout << " Simulation step " << 1 + step << endl;
    fluid->addSmokeColumn();
    fluid->stepPlume();

    // write to pbrt
    char buffer[256];
    string pbrtPath = compressionPath + string("pbrt/");
    sprintf(buffer, "%splume.compressed.%04i.pbrt", pbrtPath.c_str(), step);
    FIELD_3D::exportPbrt(fluid->density(), buffer);

    fluid->appendCompressedSubspaceVectors();
    // fluid->compareSubspace(step);

    // Use the Debug routine with loadReducedIOPAll above

    // if (step % 10 == 0)
    {
      VECTOR::printVertical = false;
      TIMER::printTimingsPerFrame(step);
    }
    
    // check if we're done
    if (step == simulationSnapshots - 1) 
    {
      string path = compressionPath + string("compression.relative.L2.subspace.error.vector");
      EIGEN::write(path.c_str(), fluid->l2Error());
      path = compressionPath + string("compression.relative.L2.ground.error.vector");
      // convert from vector<Real> to VECTOR
      VECTOR groundRelativeError(fluid->velocityErrorRelative());
      groundRelativeError.write(path);

      path = compressionPath + string("compression.density.relative.L2.ground.error.vector");
      VECTOR densityRelativeError(fluid->densityErrorRelative());
      densityRelativeError.write(path);

      TIMER::printTimings();
      // if we were already capturing a movie
      if (captureMovie) {
        // writeToQuicktime();
        movie.writeMovie(previewReducedMovie.c_str());

        // stop capturing frames

        captureMovie = false;
      }
      exit(0);
    }
   
    step++;
  }
}

//////////////////////////////////////////////////////////////////////
// Helper functions
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

// check if a file exists
//////////////////////////////////////////////////////////////////////
bool fileExists(const string& filename)
{
  FILE* file;
  file = fopen(filename.c_str(), "rb");
  
  if (file == NULL)
    return false;

  fclose(file);
  return true;
}

//////////////////////////////////////////////////////////////////////
// write to unique quicktime movie file
//////////////////////////////////////////////////////////////////////
void writeToQuicktime()
{
  // write out the movie
  int i = 0;
  char buffer[256];
  sprintf(buffer, "movieObstacleSubspaceCompressed%i.mov", i);
  string movieString(buffer);

  while (fileExists(movieString)) {
    i++;
    sprintf(buffer, "movieObstacleSubspaceCompressed%i.mov", i);
    movieString = string(buffer);
  }
  movie.writeMovie(movieString.c_str());
  // reset the movie object
  movie = QUICKTIME_MOVIE();
}
