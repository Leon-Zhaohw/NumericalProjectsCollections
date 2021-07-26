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

#include <iostream>
#include <cstring>
#include <SIMPLE_PARSER.h>
#include "OBJ.h"
#include "BOX.h"
#include "SPHERE.h"
#include "CYLINDER.h"
#include "ISO_STUFFER.h"
#include "SKELETON.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////////
// GLOBALS
//////////////////////////////////////////////////////////////////////////////
int res = 32;
ISO_STUFFER* isoStuffer;
bool animate = false;
int DIV = 1048576;
char *divisor = "M";
int WIDTH = 4;
int startingFree;
float zSlice = 0.5f;

// In Windows, pop up a preview window
#ifdef USING_GLVU
#include <glvu.hpp>
#include <snapshot.hpp>

OBJ objFile;
GLVU glvu;
SURFACE* compound;
Real meshScale;
VEC3F meshCenter;

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

    glColor4f(1.0f, 1.0f, 1.0f, 0.5f);
    glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
    isoStuffer->drawFinalTets();

    glColor4f(10.0f, 0.0f, 0.0f, 1.0f);
    glPointSize(5.0f);
    isoStuffer->drawConstrainedNodes();
    glColor4f(1.0f, 1.0f, 1.0f, 1.0f);

    //glColor4f(1,1,1,1);
    //objFile.draw();

  glvu.EndFrame();
}

//////////////////////////////////////////////////////////////////////////////
// USER-PROVIDED KEYBOARD HANDLING ROUTINE
//////////////////////////////////////////////////////////////////////////////
void userKeyboardFunc(unsigned char Key, int x, int y)
{
  static int steps = 0;
  
  switch(Key)
  {
    case 'q':
    case 'Q':
      exit(0);
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
  glutPostRedisplay();
}

//////////////////////////////////////////////////////////////////////////////
// open the GLVU window
//////////////////////////////////////////////////////////////////////////////
int glvuWindow()
{
  glvu.Init("GLVU Window",
            GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH,
            0,0,700,700);

  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  GLfloat lightZeroPosition[] = {10.0, 4.0, 10.0, 1.0};
  GLfloat lightZeroColor[] = {0.8, 1.0, 0.8, 1.0};
  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, 1);
  glLightfv(GL_LIGHT0, GL_POSITION, lightZeroPosition);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, lightZeroColor);
 
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_COLOR_MATERIAL);
  //glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);
  glShadeModel(GL_SMOOTH);
  glClearColor(0,0,0,0);

  glutDisplayFunc(userDisplayFunc);
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
#endif

//////////////////////////////////////////////////////////////////////////////
// rescale the OBJ mesh
//////////////////////////////////////////////////////////////////////////////
void rescaleObj()
{
  // get the bounding box
  VEC3F mins;
  VEC3F maxs;
  objFile.BoundingBox(mins, maxs);
  meshCenter = (mins + maxs) * 0.5;

  // get the scaling
  meshScale = maxs[0] - mins[0];
  for (int x = 0; x < 3; x++)
    if (maxs[x] - mins[x] > meshScale)
      meshScale = maxs[x] - mins[x];

  meshScale = 0.95 / meshScale;

  // scale all the vertices
  vector<VEC3F>& vertices = objFile.vertices;
  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    vertices[x] -= meshCenter;
    vertices[x] *= meshScale;
    vertices[x] += VEC3F(0.5, 0.5, 0.5);
  }
}

//////////////////////////////////////////////////////////////////////////////
// Rescale the model and the isostuffing results back to the original size
//////////////////////////////////////////////////////////////////////////////
void rescaleResults()
{
  // rescale model
  Real invScale = 1.0 / meshScale;

  vector<VEC3F>& vertices = objFile.vertices;
  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    vertices[x] -= VEC3F(0.5, 0.5, 0.5);
    vertices[x] *= invScale;
    vertices[x] += meshCenter;
  }

  vector<VEC3F*> unconstrainedNodes = isoStuffer->unconstrainedNodes();
  vector<VEC3F*> constrainedNodes = isoStuffer->constrainedNodes();

  for (unsigned int x = 0; x < unconstrainedNodes.size(); x++)
  {
    VEC3F& vertex = *unconstrainedNodes[x];
    vertex -= VEC3F(0.5, 0.5, 0.5);
    vertex *= invScale;
    vertex += meshCenter;
  }
  for (unsigned int x = 0; x < constrainedNodes.size(); x++)
  {
    VEC3F& vertex = *constrainedNodes[x];
    vertex -= VEC3F(0.5, 0.5, 0.5);
    vertex *= invScale;
    vertex += meshCenter;
  }
}

//////////////////////////////////////////////////////////////////////////////
// generate test from implicit functions
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  if (argc < 2)
  {
    cout << " Usage: " << argv[0] << " <*.cfg> <run headless?>" << endl;
    cout << " The last argument is optional. If it is set to anything," << endl;
    cout << " it will prevent a GL window from popping up at the end, " << endl;
    cout << " so you will not get the inspect the results. " << endl;
    return 0;
  }

  // run the program headless?
  bool headless = (argc >= 3) ? true : false;

  // set default parameters
  int bccRes = 64;
  string triangleMeshPath = string("../examples/");
  string triangleMeshName = string("bunny_watertight.obj");

  // read in different parameters
  string configName("");
  if (argc > 1)
    configName = string(argv[1]);
  SIMPLE_PARSER configFile(configName); 

  string outputPath = string("./meshes/tet/ragdoll/");
  string dataPath = string("./data/ragdoll/");
  triangleMeshPath = configFile.getString("triangle mesh path", triangleMeshPath);
  triangleMeshName = configFile.getString("triangle mesh name", triangleMeshName);

  bccRes      = configFile.getInt("bccRes", bccRes);
  outputPath  = configFile.getString("output path", outputPath);
  dataPath    = configFile.getString("data path", dataPath);

  // start the isostuffer
  cout << "==================================================" << endl;
  cout << "RUNNING ISOSURFACE STUFFING" << endl;
  cout << "==================================================" << endl;
  cout << "BCC grid resolution: " << bccRes << "^3" << endl;

  // see if the mesh is an OBJ or PLY file
  bool plyFile = false;
  const char* data = triangleMeshName.data();
  data += triangleMeshName.length() - 3;
  if (strcmp(data, "ply") == 0)
    plyFile = true;

  // load the triangle mesh
  string triangleMeshFullpath = triangleMeshPath + triangleMeshName;
  if (plyFile)
    objFile.LoadPly(triangleMeshFullpath);
  else
    objFile.Load(triangleMeshFullpath);
  objFile.ComputeVertexNormals();
  objFile.setBCCRes(bccRes);

  // scale the obj down to the unit cube
  rescaleObj();

  // create acceleration structures
  cout << " Creating acceleration structures ... " << endl;
  objFile.createAccelGrid();
  objFile.createDistanceGrid(bccRes);
  cout << " done. " << endl;

  // generate the mesh
  isoStuffer = new ISO_STUFFER(bccRes, bccRes, bccRes);
  isoStuffer->useExistingCaches() = false;
  isoStuffer->generateLimitedTets(objFile);
  isoStuffer->generateInsideTets(objFile);
  isoStuffer->constrainNone();

  // scale the obj and the tets back up
  rescaleResults();

  // make the output dir
  string mkdir("mkdir ");
  mkdir += outputPath;
  cout << "mkdir command: " << mkdir.c_str() << endl;
  system(mkdir.c_str());
  string finalMeshFile = outputPath + triangleMeshName + string(".tetmesh");
  isoStuffer->writeFile(finalMeshFile.c_str());

#ifdef USING_GLVU
  if (!headless)
  {
    glutInit(&argc, argv);
    glvuWindow();
  }
#endif
  return 0;
}
