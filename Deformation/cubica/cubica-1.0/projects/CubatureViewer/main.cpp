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

using namespace std;
SUBSPACE_TET_MESH* tetMesh;
SUBSPACE_INTEGRATOR* integrator;

//////////////////////////////////////////////////////////////////////////////
// GLOBALS
//////////////////////////////////////////////////////////////////////////////
GLVU glvu;
int windowStartX = 0;
int windowStartY = 0;
int windowWidth = 700;
int windowHeight = 700;
bool mouseClicked = false;
float clickZ;
bool animate = false;
bool regressionTesting = false;
string outputPath;

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

    glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
    tetMesh->drawSurfaceFaces();
    glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
    glPointSize(10.0f);
    tetMesh->drawConstrainedNodes();

    if (mouseClicked)
      integrator->drawForceVector();
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
// Dump several frames to disk after a poke
//////////////////////////////////////////////////////////////////////////////
void pokeTest()
{
  integrator->poke();
  for (int x = 0; x < 25; x++)
  {
    string filename("./matlab/bunny_full.m");
    string varname("full");
    tetMesh->verticesToMatlab(filename.c_str(), varname.c_str(), x);
    for (int y = 0; y < 5; y++)
      integrator->stepInvertibleImplicit(true, false);
  }
}

//////////////////////////////////////////////////////////////////////////////
// Dump several frames to disk over a gradual pull
//////////////////////////////////////////////////////////////////////////////
void pullTest()
{
  for (int x = 0; x < 300; x++)
  {
    string filename("./matlab/bunny_full_pull.m");
    string binaryname("./matlab/bunny_full_pull.binary");
    string varname("full");
    //tetMesh->verticesToMatlab(filename.c_str(), varname.c_str(), x);
    tetMesh->verticesToBinary(binaryname.c_str(), x);
    for (int y = 0; y < 1; y++)
    {
      integrator->stepInvertibleQuasistatic(true);
      integrator->pull();
    }
  }
  glutPostRedisplay();
}

//////////////////////////////////////////////////////////////////////////////
// USER-PROVIDED KEYBOARD HANDLING ROUTINE
//////////////////////////////////////////////////////////////////////////////
void userKeyboardFunc(unsigned char Key, int x, int y)
{
  static int steps = 0;
  
  switch(Key)
  {
    case 'a':
      animate = !animate;
      break;
    case 'p':
      cout << __FILE__ << " " << __LINE__ << " POKE: " << endl;
      integrator->poke();
      break;

    case 'q':
    case 'Q':
      exit(0);
      break;

    case 'd':
      //tetMesh->verticesToMatlab("bunny_full.m", "full");
      pullTest();
      break;

    case ' ':
      integrator->printTimingBreakdown();
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

  if (animate)
  {
    if (regressionTesting)
    {
      integrator->changeGravity(VEC3F(0, 0, 1), 1);
      integrator->addGravity();
    }
    bool debugPrint = false;
    integrator->stepInvertibleImplicit(true, debugPrint);
    tetMesh->updateSurfaceMesh();

    if (steps == 100 && regressionTesting)
    {
      string filename = outputPath + string("regression.vector");
      cout << " Writing out regression result " << filename.c_str() << endl;
      integrator->velocity().write(filename.c_str());
      exit(0);
    }

    //integrator->stepExplicit();
    //integrator->stepInvertibleSemiImplicit(true);
    //integrator->stepInvertibleQuasistatic(true);
    //integrator->pull();

    //if (steps == 20) animate = !animate;
    steps++;
  }
  glutPostRedisplay();
}

//////////////////////////////////////////////////////////////////////////////
// USER-PROVIDED MOUSE HANDLING ROUTINE
//////////////////////////////////////////////////////////////////////////////
void userMouseFunc(int button, int state, int x, int y)
{
  int Modifiers = glutGetModifiers();
  if (button == GLUT_LEFT_BUTTON && 
      state == GLUT_DOWN &&
      Modifiers & GLUT_ACTIVE_SHIFT)
  {
    // retrieve and store the depth of this click
    GLint viewport[4];
	  glGetIntegerv(GL_VIEWPORT, viewport);
    glReadPixels(x, viewport[3] - y, 
                 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &clickZ);

    // get the world space coordinate
    VEC3F point = unproject(x,y,clickZ);

    // hand the coordinate to the tet mesh
    integrator->click(point);
    mouseClicked = true;
    return;
  }
  if (button == GLUT_LEFT_BUTTON && 
      state == GLUT_UP)
  {
    mouseClicked = false;
    integrator->unclick();
    return;
  }

  // pass through to default handler
  glvu.Mouse(button,state,x,y);
}

//////////////////////////////////////////////////////////////////////////////
// USER-PROVIDED MOUSE MOTION ROUTINE
//////////////////////////////////////////////////////////////////////////////
void userMotionFunc(int x, int y)
{
  if (mouseClicked)
  {
    VEC3F point = unproject(x,y,clickZ);
    integrator->drag(point);
    return;
  }

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
 
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);
  glShadeModel(GL_SMOOTH);
  glClearColor(0,0,0,0);

  glutDisplayFunc(userDisplayFunc);
  glutMouseFunc(userMouseFunc);
  glutMotionFunc(userMotionFunc);
  glutKeyboardFunc(userKeyboardFunc);
  glutIdleFunc(userIdleFunc);

  Vec3f ModelMin(-10,-10,-10), ModelMax(10,10,10), 
        Eye(0.5,0.5,3), LookAtCntr(0.5,0.5,0.5), Up(0,1,0);

  LookAtCntr = bboxCenter;

  Vec3f EyeDirection = Eye - LookAtCntr;
  EyeDirection *= eyeDistanceScale;
  Eye = LookAtCntr + EyeDirection;

  float Yfov = 45;
  float Aspect = 1;
  float Near = 0.001f;
  float Far = 10.0f;
  glvu.SetAllCams(ModelMin, ModelMax, Eye,LookAtCntr,Up, Yfov,Aspect, Near,Far);

  // reset center
  //Vec3f center(0.5f, 0.5f, 0.5f);
  glvu.SetWorldCenter(bboxCenter);
  
  glutMainLoop();

  // Control flow will never reach here
  return EXIT_SUCCESS;
}

//////////////////////////////////////////////////////////////////////////////
// Read a material file and return a pointer to it
//////////////////////////////////////////////////////////////////////////////
MATERIAL* readMaterial(SIMPLE_PARSER& parser)
{
  return SIMPLE_PARSER::READ_MATERIAL( parser );
}

//////////////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  if (argc <= 1)
  {
    cout << " USAGE: " << argv[0] << " *.cfg " << endl;
    return 0;
  }
  // input parameters
  string triangleMeshPath;
  string triangleMeshName;
  string tetMeshName;

  bool inexactNewton = false;
  if ( argc > 2 )
  {
    inexactNewton = ( atoi( argv[2] ) != 0 );
  }

  // read in different parameters
  string configName(argv[1]);
  SIMPLE_PARSER configFile(configName); 
  triangleMeshPath = configFile.getString("triangle mesh path", triangleMeshPath);
  triangleMeshName = configFile.getString("triangle mesh name", triangleMeshName);
  outputPath       = configFile.getString("output path", outputPath);

  if (configFile.defined("tet mesh name"))
    tetMeshName = configFile.getString("tet mesh name", tetMeshName);
  else
    tetMeshName = outputPath + triangleMeshName + string(".tetmesh");

  // read in how many materials there are
  int totalMaterials = 0;
  totalMaterials = configFile.getInt("total materials", totalMaterials);
  if (totalMaterials == 0)
  {
    cout << " NO MATERIALS SPECIFIED!!!!" << endl;
    exit(1);
  }
  
  // read in the actual materials
  MATERIAL** materials = new MATERIAL*[totalMaterials];
  for (int x = 0; x < totalMaterials; x++)
  {
    // read in the config file name for the material
    char buffer[256];
    sprintf(buffer, "material %i", x);
    string materialString(buffer);
    string materialFilename;
    materialFilename = configFile.getString(materialString.c_str(), materialFilename);

    // open the config file
    SIMPLE_PARSER materialFile(materialFilename);
    
    // get the material
    MATERIAL* material = readMaterial(materialFile);
    
    // set the invertible wrapper if necessary
    bool invertible = false;
    invertible = configFile.getBool("invertible", invertible);
    if (invertible)
    {
      cout << " Setting material to invertible" << endl;
      material = new INVERTIBLE(material);
    }

    materials[x] = material;
  }

  // output precision
  cout << " Precision: " << sizeof(Real) * 8 << " bit" << endl;

  Real meshScale           = configFile.getFloat("mesh scale", 1.0);
  
  cout << "Got mesh scale " << meshScale << endl;

  cout << "==================================================" << endl;
  cout << " Reading mesh file " << tetMeshName.c_str() << endl;
  cout << "==================================================" << endl;
  tetMesh = new SUBSPACE_TET_MESH(tetMeshName.c_str(), materials, totalMaterials,
                                  false, NULL, NULL, meshScale);

  Real bbox[6];
  tetMesh->boundingBox( bbox );

  cout << "Bounding box " << endl;
  cout << "\t [ " << bbox[0] << "\t" << bbox[1] << endl;
  cout << "\t [ " << bbox[2] << "\t" << bbox[3] << endl;
  cout << "\t [ " << bbox[4] << "\t" << bbox[5] << endl;
  cout << endl << endl;

  Vec3f boxCenter( ( bbox[0] + bbox[1] ) * 0.5,
                   ( bbox[2] + bbox[3] ) * 0.5,
                   ( bbox[4] + bbox[5] ) * 0.5 );
 
  Real timestep = 1.0f / 60.0f;
  timestep = configFile.getFloat("timestep", timestep);
  cout << endl;
  cout << " Using timestep: " << timestep << endl;
  //integrator = new SUBSPACE_INTEGRATOR(tetMesh, timestep, 0.01, 0.01);

  //Real *whatIsGoingOn = new Real[30];

  Real rayleighAlpha = 0.003;
  Real rayleighBeta = 0.003;

  rayleighAlpha = configFile.getFloat("rayleigh alpha", rayleighAlpha);
  rayleighBeta = configFile.getFloat("rayleigh beta", rayleighBeta);

  Real deformationTolerance = 0.005;
  integrator = new SUBSPACE_INTEGRATOR(tetMesh, timestep, rayleighAlpha,
                                       rayleighBeta, 9.8,
                                       // Inexact newton solves with
                                       // a 10% deformation tolerance
                                       inexactNewton, deformationTolerance);

  //integrator->maxNewtonSteps() = 10;
  integrator->solverEps() = 0.01;
  integrator->maxNewtonSteps() = 10;
  cout << " Max Newton steps: " << integrator->maxNewtonSteps() << endl;

  // factor to increase the mouse-input force by
  double forceMultiplier = 0.1;
  forceMultiplier = configFile.getFloat("force multiplier", forceMultiplier);
  integrator->forceMultiplier() = forceMultiplier;
  cout << " Using force multiplier: " << forceMultiplier << endl;

  Real moveSpeed = 1.0;
  Real eyeDistanceScale = 1.0;

  moveSpeed = configFile.getFloat("navigation speed", moveSpeed);
  eyeDistanceScale = configFile.getFloat("view distance", eyeDistanceScale);

  // are we doing a regression test?
  regressionTesting = configFile.getBool("regression testing", regressionTesting);
  if (regressionTesting)
    animate = true;

  glvu.SetMoveSpeed( moveSpeed * glvu.GetMoveSpeed() );

  glutInit(&argc, argv);
  glvuWindow( boxCenter, eyeDistanceScale );
}
