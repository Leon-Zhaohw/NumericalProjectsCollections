/*
This file is part of Cubica.
 
Cubica is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Cubica is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Cubica.  If not, see <http://www.gnu.org/licenses/>.
*/
//------------------------------------------------------------------------------
// GL interface elements are from:
//------------------------------------------------------------------------------
// GLVU : Copyright 1997 - 2002 
//        The University of North Carolina at Chapel Hill
//------------------------------------------------------------------------------
// Permission to use, copy, modify, distribute and sell this software and its 
// documentation for any purpose is hereby granted without fee, provided that 
// the above copyright notice appear in all copies and that both that copyright 
// notice and this permission notice appear in supporting documentation. 
// Binaries may be compiled with this software without any royalties or 
// restrictions. 
//
// The University of North Carolina at Chapel Hill makes no representations 
// about the suitability of this software for any purpose. It is provided 
// "as is" without express or implied warranty.

//#include <GL/glew.h>
#include <glvu.hpp>
#include <snapshot.hpp>
#include <iostream>
#include <TET_MESH.h>
#include <SUBSPACE_TET_MESH.h>
#include <SPARSE_MATRIX.h>
#include <STVK.h>
#include <MOONEY_RIVLIN.h>
#include <ARRUDA_BOYCE.h>
#include <NEO_HOOKEAN.h>
#include <INVERTIBLE.h>
#include <SIMPLE_PARSER.h>
#include <SUBSPACE_INTEGRATOR.h>
//#include <GLSL_SHADER.h>
#include <SKELETON.h>

using namespace std;

//////////////////////////////////////////////////////////////////////////////
// Meshes and integrators
//////////////////////////////////////////////////////////////////////////////
OBJ *embeddedMesh;
SKELETON* skeleton;
TET_MESH* tetMesh;
string tetMeshName("");
string dataPath("");
string posePath("");

//////////////////////////////////////////////////////////////////////////////
// Various function declarations
//////////////////////////////////////////////////////////////////////////////
void generateSurfaceBasis();
void updateSurfaceMesh();
void initShaders(const char *vertexFilename, const char *fragmentFilename);
void setupLighting();

//////////////////////////////////////////////////////////////////////////////
// GLOBALS
//////////////////////////////////////////////////////////////////////////////
GLVU glvu;
static int windowStartX = 0;
static int windowStartY = 0;
static int windowWidth = 700;
static int windowHeight = 700;
static bool animate = false;
static int currentBone = 0;

//////////////////////////////////////////////////////////////////////////////
// Shader stuff for rendering
//////////////////////////////////////////////////////////////////////////////
GLenum shaderProgram;
GLenum vertexShader;
GLenum fragmentShader;

//////////////////////////////////////////////////////////////////////////////
// Material display constants
//////////////////////////////////////////////////////////////////////////////
GLfloat mat_shininess[] = { 120.0 };

//GLfloat mat_diffuse[] = { 0.54, 0.0, 0.0, 1.0 };
GLfloat mat_diffuse[] = { 0.5, 0.2, 0.3, 1.0 };
//GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat mat_specular[] = { 0.0, 0.0, 0.0, 0.0 };

GLfloat cursor_diffuse[] = { 0.75, 0.75, 0.75, 1.0 };

struct timeval myTime;
uint32_t oldTime;
uint32_t newTime;
uint32_t frameCount;

//////////////////////////////////////////////////////////////////////////////
// For controlling external forces.
// If a point selected for external forces on the mesh is
// farther away than this distance from the actual clicked
// point in world space - no forces are applied
//////////////////////////////////////////////////////////////////////////////
Real maxForceRadius = -1.0;

//////////////////////////////////////////////////////////////////////////////
// USER-PROVIDED DRAWING ROUTINE
//////////////////////////////////////////////////////////////////////////////
void userDisplayFunc()
{
  // get camera settings
  Camera* camera = glvu.GetCurrentCam();
  Vec3f Eye, Lookat, Up;
  camera->GetLookAtParams(&Eye, &Lookat, &Up);

  // repack camera settings into arrays
  float eye[] = {Eye.x, Eye.y, Eye.z};
  float look[3];
  look[0] = Lookat.x - eye[0];
  look[1] = Lookat.y - eye[1];
  look[2] = Lookat.z - eye[2];
  float magnitude = 1.0f / sqrt(look[0] * look[0] + look[1] * look[1] + look[2] * look[2]);
  look[0] *= magnitude;
  look[1] *= magnitude;
  look[2] *= magnitude;

  glvu.BeginFrame();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

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

    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);

    glColor4f(1,1,1,1);
    //tetMesh->drawSurfaceFaces();

    skeleton->drawTetWeights();
    //skeleton->drawSkinning();

    glDisable(GL_DEPTH_TEST);
      glLineWidth(5.0);
      glColor4f(0,10,0,1);
      skeleton->drawOdeBones();

      glColor4f(10,0,0,1);
      tetMesh->drawConstrainedNodes();
    glEnable(GL_DEPTH_TEST);

  glvu.EndFrame();
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

//////////////////////////////////////////////////////////////////////////////
// USER-PROVIDED KEYBOARD HANDLING ROUTINE
//////////////////////////////////////////////////////////////////////////////
void userKeyboardFunc(unsigned char Key, int x, int y)
{
  switch(Key)
  {
    case 'a':
      animate = !animate;
      break;

    case '1':
      currentBone++;
      if (currentBone >= skeleton->totalBones() - 1)
        currentBone = 0;

      cout << " Viewing bone " << currentBone << endl;
      break;

    case 'q':
    case 'Q':
      exit(0);
      break;

    case 'C':
      {
        string constrainedName = tetMeshName + string(".constrained");
        skeleton->constrainOdeBoneTets(tetMesh, constrainedName);
      }
      break;
    // make bone lengths shorted
    case 'S':
      {
        vector<Real>& lengths = skeleton->boneLengths();
        for (unsigned int x = 0; x < lengths.size(); x++)
          lengths[x] *= 0.75;
      }
      break;
    case 'b':
      skeleton->buildOdeSkinning(tetMesh);
      break;
  };

  glutPostRedisplay();
  if (Key != '=')
    glvu.Keyboard(Key,x,y);
}

//////////////////////////////////////////////////////////////////////////////
// USER-PROVIDED IDLE ROUTINE
//////////////////////////////////////////////////////////////////////////////
void userIdleFunc()
{
  static int steps = 0;
  char buffer[256];

  if (animate)
  {
    sprintf(buffer, "%04i", steps);
    string skeletonFile = posePath + string("ode.motion.") + string(buffer) + string(".skeleton");
    skeleton->loadOdeFrame(skeletonFile.c_str());
    skeleton->updateOdeSkinning(false);

    steps++;
  }

  glutPostRedisplay();
}

//////////////////////////////////////////////////////////////////////////////
// USER-PROVIDED MOUSE HANDLING ROUTINE
//////////////////////////////////////////////////////////////////////////////
void userMouseFunc(int button, int state, int x, int y)
{
  // pass through to default handler
  glvu.Mouse(button,state,x,y);
}

//////////////////////////////////////////////////////////////////////////////
// USER-PROVIDED MOUSE MOTION ROUTINE
//////////////////////////////////////////////////////////////////////////////
void userMotionFunc(int x, int y)
{
  // pass through to default handler
  glvu.Motion(x,y);
}

//////////////////////////////////////////////////////////////////////////////
// open the GLVU window
//////////////////////////////////////////////////////////////////////////////
int glvuWindow( Vec3f bboxCenter, Real eyeDistanceScale = 1.0 )
{
  glvu.Init("GLVU Window",
            GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH,
            windowStartX, windowStartY, windowWidth, windowHeight);

  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  GLfloat lightZeroPosition[] = {10.0, 4.0, 10.0, 1.0};
  GLfloat lightZeroColor[] = {0.8, 1.0, 0.8, 1.0};
  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, 1);
  glLightfv(GL_LIGHT0, GL_POSITION, lightZeroPosition);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, lightZeroColor);

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); 
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_COLOR_MATERIAL);
  //glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_LINE_SMOOTH);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glShadeModel(GL_SMOOTH);
  //glShadeModel(GL_FLAT);
  glClearColor(0,0,0,0);

  glutDisplayFunc(userDisplayFunc);
  glutMouseFunc(userMouseFunc);
  glutMotionFunc(userMotionFunc);
  glutKeyboardFunc(userKeyboardFunc);
  glutIdleFunc(userIdleFunc);

  Vec3f ModelMin(-10,-10,-10), ModelMax(10,10,10), 
        Eye(0.5,0.5,3), LookAtCntr(0.5,0.5,0.5), Up(0,1,0);
  float Yfov = 45;
  float Aspect = 1;
  float Near = 0.001f;
  float Far = 10.0f;
  glvu.SetAllCams(ModelMin, ModelMax, Eye,LookAtCntr,Up, Yfov,Aspect, Near,Far);

  // reset center
  Vec3f center(0.5f, 0.5f, 0.5f);
  glvu.SetWorldCenter(center);
  
  glutMainLoop();

  // Control flow will never reach here
  return EXIT_SUCCESS;
}

//////////////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  if (argc < 2)
  {
    cout << " USAGE: " << argv[0] << " *.cfg <run headless>" << endl;
    cout << " The last argument is optional. If it is set to anything," << endl;
    cout << " it will prevent a GL window from popping up at the end, " << endl;
    cout << " so you will not get the inspect the results. " << endl;
    return 0;
  }

  // run the program headless?
  bool headless = (argc >= 3) ? true : false;

  // input parameters
  string outputPath("");
  string triangleMeshName = string("bunny_watertight.obj");

  // read in different parameters
  string configName(argv[1]);
  SIMPLE_PARSER configFile(configName); 
  outputPath = configFile.getString("output path", outputPath);
  dataPath   = configFile.getString("data path", dataPath);
  posePath   = configFile.getString("pose path", posePath);
  triangleMeshName = configFile.getString("triangle mesh name", triangleMeshName);

  cout << " triangle mesh name: " << triangleMeshName.c_str() << endl;

  tetMeshName = outputPath + triangleMeshName + string(".tetmesh");
  string constrainedName = tetMeshName + string(".constrained");
  FILE* file = fopen(constrainedName.c_str(), "rb");
  if (file != NULL)
  {
    tetMeshName = constrainedName;
    fclose(file);
  }

  cout << " Reading tet mesh " << tetMeshName.c_str() << endl;
  tetMesh = new TET_MESH(tetMeshName.c_str(), NULL, 0, false);

  // load in the skeleton
  string skeletonFile = posePath + string("ode.motion.0000.skeleton");
  skeleton = new SKELETON(skeletonFile, true);

  skeleton->tetMesh() = tetMesh;
  skeleton->buildOdeSkinning(tetMesh);

  Real moveSpeed = 1.0;
  Real eyeDistanceScale = 1.0;
  glvu.SetMoveSpeed( moveSpeed * glvu.GetMoveSpeed() );

  if (tetMesh->constrained() == false)
  {
    string constrainedName = tetMeshName + string(".constrained");
    skeleton->constrainOdeBoneTets(tetMesh, constrainedName);
  }

  if (!headless)
  {
    glutInit(&argc, argv);
    Vec3f boxCenter(0.5, 0.5, 0.5);
    glvuWindow( boxCenter, eyeDistanceScale );
  }
}
