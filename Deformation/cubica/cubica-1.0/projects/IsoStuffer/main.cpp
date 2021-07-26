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
//#include "obj.h"
#include "OBJ.h"
#include "BOX.h"
#include "SPHERE.h"
#include "CYLINDER.h"
#include "COMPOUND.h"
#include "ISO_STUFFER.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////////
// GLOBALS
//////////////////////////////////////////////////////////////////////////////
//int res = 64;
int res = 32;
//int res = 16;
ISO_STUFFER* isoStuffer;
OBJ objFile;
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

GLVU glvu;

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
    //objFile.draw();
    glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
    //isoStuffer->drawTets();
    isoStuffer->drawFinalTets();
    //isoStuffer->drawFinalTetsZSlice(zSlice);
    //isoStuffer->drawSurfaceTetsZSlice(zSlice);
    //isoStuffer->drawInsideTetsZSlice(zSlice);
    //isoStuffer->drawOutsideTets();

    glColor4f(10.0f, 0.0f, 0.0f, 1.0f);
    glPointSize(5.0f);
    isoStuffer->drawConstrainedNodes();
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
    /*
    case '!':
      isoStuffer->generateCase1Only();
      break;
    case '@':
      isoStuffer->generateCase4Only();
      break;
    case '#':
      isoStuffer->generateCase5Only();
      break;
    case '$':
      isoStuffer->generateCase6Only();
      break;
    case '%':
      isoStuffer->generateCase7Only();
      break;
    case '^':
      isoStuffer->generateCase9Only();
      break;
    case '&':
      isoStuffer->generateCase10Only();
      break;
      */
    case '*':
      isoStuffer->generateFinalTets();
      break;
    case ',':
      zSlice -= 1.0f / 256.0f;
      break;
    case '.':
      zSlice += 1.0f / 256.0f;
      break;
    case 'a':
      animate = !animate;
      break;
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
  glEnable(GL_CULL_FACE);
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

//////////////////////////////////////////////////////////////////////////////
// generate test from an OBJ file 
//////////////////////////////////////////////////////////////////////////////
int tetsFromObj(int argc, char* argv[])
{
  objFile.Load("../examples/bunny_watertight.obj");
  //objFile.LoadPly("../examples/happy_vrip_res4.ply");
  //objFile.LoadPly("../examples/dragon_vrip_res4.ply");
  objFile.ComputeVertexNormals();
  objFile.normalize(isoStuffer->xRes());
  objFile.setBCCRes(isoStuffer->xRes());

  // center the mesh about (0.5, 0.5, 0.5)
  VEC3 minVecd;
  VEC3 maxVecd;
  objFile.BoundingBox(minVecd, maxVecd);
  VEC3F minVec(minVecd);
  VEC3F maxVec(maxVecd);
  VEC3F halfVec(0.5f, 0.5f, 0.5f);
  VEC3F center;
  center[0] = (maxVec[0] - minVec[0]) * halfVec[0] + minVec[0];
  center[1] = (maxVec[1] - minVec[1]) * halfVec[1] + minVec[1];
  center[2] = (maxVec[2] - minVec[2]) * halfVec[2] + minVec[2];
  objFile.translate(halfVec - center);

  // create acceleration structures
  objFile.createAccelGrid();
  objFile.createDistanceGrid(isoStuffer->xRes());

  // generate the final mesh
  isoStuffer->useExistingCaches() = false;
  isoStuffer->generateLimitedTets(objFile);
  isoStuffer->generateInsideTets(objFile);
  isoStuffer->constrainMinZ();

  // save the final mesh
  char buffer[256];
  sprintf(buffer, "%i", res);
  string filename = string("bunny.") + string(buffer) + string(".binary.mesh");
  isoStuffer->writeFile(filename.c_str());
  isoStuffer->writeFileLegacy("bunny.32.legacy.mesh");

  glvuWindow();

  return 0;
}

//////////////////////////////////////////////////////////////////////////////
// generate test from implicit functions
//////////////////////////////////////////////////////////////////////////////
int tetsFromImplicit(int argc, char* argv[])
{
  int invRes = 1.0f / res;
  BOX box(0.95, 0.05,
          //0.46, 0.45, // good 100 setting
          0.48, 0.45, // good 32 setting
          //0.55, 0.45, // good 8 setting
          //5.125*invRes + 0.05, 5.125*invRes -0.01,
          0.95, 0.05);
  SPHERE sphere(0.1f);
  CYLINDER cylinder;

  COMPOUND compound;

  compound.addSurface(new SPHERE(0.25f));
  compound.addSurface(new CYLINDER());

  cout << " Generating all tets ... " << endl;
  isoStuffer->useExistingCaches() = false;
  isoStuffer->generateAllTets();
  cout << " Generating inside tets ... " << endl;
  isoStuffer->generateInsideTets(box);

  cout << " Generating constrained nodes ... " << endl;
  float center[] = {0.5, 0.5};
  float caps[] = {0.45, 0.45333};
  CYLINDER constraint(center, 0.1, caps);
  isoStuffer->generateConstrainedNodes(constraint);
  cout << " Outputting final tets file ... " << endl;
  isoStuffer->writeFile("finalTets.mesh");

  glvuWindow();
  
  return 0;
}

//////////////////////////////////////////////////////////////////////////////
// generate test from implicit functions
//////////////////////////////////////////////////////////////////////////////
int tetsBalloon(int argc, char* argv[])
{
  float invRes = 1.0f / res;
  COMPOUND compound;

  float bottom = invRes;
  float top = 4.0f * invRes;

  compound.addSurface(new BOX(0.95, 0.05,
                             top, bottom,
                             0.95, 0.05));
  float thickness = invRes;
  compound.subtractSurface(new BOX(0.95 - thickness, 0.05 + thickness,
                                  top  - thickness, bottom + thickness,
                                  0.95 - thickness, 0.05 + thickness));

  cout << " Generating all tets ... " << endl;
  isoStuffer->useExistingCaches() = false;
  isoStuffer->generateAllTets();
  cout << " Generating inside tets ... " << endl;
  isoStuffer->generateInsideTets(compound);

  cout << " Generating constrained nodes ... " << endl;
  float center[] = {0.5, 0.5};
  float caps[] = {bottom, bottom + invRes};
  CYLINDER constraint(center, 0.1, caps);
  isoStuffer->generateConstrainedNodes(constraint);
  cout << " Outputting final tets file ... " << endl;
  isoStuffer->writeFile("finalTets.mesh");

  glvuWindow();
  
  return 0;
}

//////////////////////////////////////////////////////////////////////////////
// generate test from implicit functions
//////////////////////////////////////////////////////////////////////////////
int tetsBalloonL(int argc, char* argv[])
{
  float invRes = 1.0f / res;
  COMPOUND compound;

  float bottom = invRes;
  float top = 4.0f * invRes;
  float thickness = invRes;
  float jitter = 0.25 * invRes;

  cout << " outside top : " << top << " bottom: " << bottom << endl;

  compound.addSurface(new BOX((1 - 2.0f * invRes), 0.5,
                             top, bottom,
                             (1 - 2.0f * invRes), invRes));

  compound.addSurface(new BOX((1 - 2.0f * invRes), invRes,
                             top, bottom,
                             (1 - 2.0f * invRes), 0.5));

  compound.subtractSurface(new BOX((1 - 2.0f * invRes) - thickness, 0.5 + thickness,
                                  top  - thickness, bottom + thickness,
                                  (1 - 2.0f * invRes) - thickness, invRes + thickness));
  compound.subtractSurface(new BOX((1 - 2.0f * invRes) - thickness, invRes + thickness,
                                  top  - thickness, bottom + thickness,
                                  (1 - 2.0f * invRes) - thickness, 0.5 + thickness));

  cout << " Generating limited tets ... " << endl;
  //isoStuffer->useExistingCaches() = false;
  isoStuffer->generateAllTets();
  //isoStuffer->generateLimitedTets(compound);
  cout << " Generating inside tets ... " << endl;
  isoStuffer->generateInsideTets(compound);

  cout << " Generating constrained nodes ... " << endl;
  float center[] = {0.75, 0.75};
  float caps[] = {bottom, bottom + invRes};
  CYLINDER constraint(center, 0.1, caps);
  isoStuffer->generateConstrainedNodes(constraint);
  cout << " Outputting final tets file ... " << endl;
  isoStuffer->writeFile("finalTets.mesh");

  glvuWindow();
  
  return 0;
}

//////////////////////////////////////////////////////////////////////////////
// generate test from implicit functions
//////////////////////////////////////////////////////////////////////////////
int tetsBalloonH(int argc, char* argv[])
{
  float invRes = 1.0f / res;
  COMPOUND compound;
  float oneThird = (int)(res / 3) / (float)res;

  float bottom = invRes;
  float top = 4.0f * invRes;
  float thickness = invRes;
  float jitter = 0.25 * invRes;

  float ribHeight = 1.0f;
  float ribWidth = 2.0f;

  cout << " outside top : " << top << " bottom: " << bottom << endl;

  // right stick
  compound.addSurface(new BOX((1 - 3.0f * invRes), 2 * oneThird,
                             top, bottom,
                             (1 - 2.0f * invRes), invRes));
  // left stick
  compound.addSurface(new BOX(oneThird, invRes,
                             top, bottom,
                             (1 - 2.0f * invRes), invRes));

  // middle stick
  compound.addSurface(new BOX((1 - 3.0f * invRes), invRes,
                             top, bottom,
                             2 * oneThird, oneThird));

  // middle top rib
  compound.addSurface(new BOX((1 - 3.0f * invRes), invRes,
                             top + invRes * ribHeight, bottom,
                             2 * oneThird, 2 * oneThird - invRes * ribWidth));

  // middle bottom rib
  compound.addSurface(new BOX((1 - 3.0f * invRes), invRes,
                             top + invRes * ribHeight, bottom,
                             oneThird + invRes * ribWidth, oneThird));

  // right stick right rib
  compound.addSurface(new BOX((1 - 3.0f * invRes), (1 - 3.0f * invRes) - invRes * ribWidth,
                             top + invRes * ribHeight, bottom,
                             (1 - 2.0f * invRes), invRes));

  // right stick top rib
  compound.addSurface(new BOX((1 - 3.0f * invRes), 2 * oneThird,
                             top + invRes * ribHeight, bottom,
                             (1 - 2.0f * invRes), (1 - 2.0f * invRes) - invRes * ribWidth));

  // right stick bottom rib
  compound.addSurface(new BOX((1 - 3.0f * invRes), 2 * oneThird,
                             top + invRes * ribHeight, bottom,
                             invRes + invRes * ribWidth, invRes));

  // right stick left rib
  compound.addSurface(new BOX(2 * oneThird + invRes * ribWidth, 2 * oneThird, 
                             top + invRes * ribHeight, bottom,
                             (1 - 2.0f * invRes), invRes));

  // left stick right rib
  compound.addSurface(new BOX(oneThird, oneThird - invRes * ribWidth,
                             top + invRes * ribHeight, bottom,
                             (1 - 2.0f * invRes), invRes));
  
  // left stick left rib
  compound.addSurface(new BOX(invRes + invRes * ribWidth, invRes,
                             top + invRes * ribHeight, bottom,
                             (1 - 2.0f * invRes), invRes));
  // left stick top rib
  compound.addSurface(new BOX(oneThird, invRes,
                             top + invRes * ribHeight, bottom,
                             (1 - 2.0f * invRes), (1 - 2.0f * invRes) - invRes * ribWidth));
  // left stick bottom rib
  compound.addSurface(new BOX(oneThird, invRes,
                             top + invRes * ribHeight, bottom,
                             invRes + invRes * ribWidth, invRes));

  // carve out the inside
  // right stick
  compound.subtractSurface(new BOX((1 - 3.0f * invRes) - invRes, 2 * oneThird + invRes,
                             top - invRes, bottom + invRes,
                             (1 - 2.0f * invRes) - invRes, invRes + invRes));
  // left stick
  compound.subtractSurface(new BOX(oneThird - invRes, invRes + invRes,
                                   top - invRes, bottom + invRes,
                                   (1 - 2.0f * invRes) - invRes, invRes + invRes));
  // middle stick
  compound.subtractSurface(new BOX((1 - 3.0f * invRes) - invRes, invRes + invRes,
                                   top - invRes, bottom + invRes,
                                   2 * oneThird - invRes, oneThird + invRes));

  cout << " Generating limited tets ... " << endl;
  isoStuffer->useExistingCaches() = false;
  isoStuffer->generateAllTets();
  //isoStuffer->generateLimitedTets(compound);
  cout << " Generating inside tets ... " << endl;
  isoStuffer->generateInsideTets(compound);

  cout << " Generating constrained nodes ... " << endl;
  BOX constraint((1 - 3.0f * invRes), (1 - 4.0f * invRes),
                             bottom + 1.5f * invRes, 0,
                             (1 - 2.0f * invRes), invRes);
      
  isoStuffer->generateConstrainedNodes(constraint);
  cout << " Outputting final tets file ... " << endl;
  isoStuffer->writeFile("balloon.h.binary.mesh");
  glvuWindow();
  
  return 0;
}

//////////////////////////////////////////////////////////////////////////////
// generate a flagpole
//////////////////////////////////////////////////////////////////////////////
int tetsFlagpole(int argc, char* argv[])
{
  float invRes = 1.0f / res;
  COMPOUND compound;

  float center[] = {0.5, 0.5};
  float caps[] = {invRes, 0.8};
  compound.addSurface(new CYLINDER(center, invRes, caps));

  float sphereCenter[] = {0.5, 0.8, 0.5};
  compound.addSurface(new SPHERE(invRes * 2, sphereCenter));

  cout << " Generating all tets ... " << endl;
  isoStuffer->useExistingCaches() = false;
  isoStuffer->generateAllTets();
  cout << " Generating inside tets ... " << endl;
  isoStuffer->generateInsideTets(compound);

  cout << " Generating constrained nodes ... " << endl;
  caps[0] = invRes;
  caps[1] = 1.25 * invRes;
  CYLINDER constraint(center, 0.1, caps);
  isoStuffer->generateConstrainedNodes(constraint);
  cout << " Outputting final tets file ... " << endl;

  isoStuffer->centerMesh();
  isoStuffer->scaleMesh(4.0f);
  isoStuffer->writeFile("flagpole.mesh");

  glvuWindow();
  
  return 0;
}

//////////////////////////////////////////////////////////////////////////////
// generate an accordian mesh
//////////////////////////////////////////////////////////////////////////////
int tetsAccordian(int argc, char* argv[])
{
  float invRes = 1.0f / res;
  COMPOUND compound;

  float center[] = {0.5, 0.5};
  float caps[] = {invRes, 0.5f};
  compound.addSurface(new CYLINDER(center, 0.2f, caps));

  caps[0] = 2 * invRes;
  caps[1] = 0.5f - invRes;
  compound.subtractSurface(new CYLINDER(center, 0.2f - 2 * invRes, caps));

  for (int x = 2; x < res / 2; x+=3)
  {
    caps[0] = x * invRes;
    caps[1] = (x+1) * invRes;
    compound.addSurface(new CYLINDER(center, 0.3f - invRes, caps));
  }

  cout << " Generating all tets ... " << endl;
  isoStuffer->useExistingCaches() = false;
  isoStuffer->generateAllTets();
  cout << " Generating inside tets ... " << endl;
  isoStuffer->generateInsideTets(compound);

  cout << " Generating constrained nodes ... " << endl;
  caps[0] = invRes;
  caps[1] = 1.25 * invRes;
  CYLINDER constraint(center, 0.5, caps);
  isoStuffer->generateConstrainedNodes(constraint);
  cout << " Outputting final tets file ... " << endl;

  isoStuffer->centerMesh();
  isoStuffer->scaleMesh(4.0f);
  isoStuffer->writeFile("accordian.64.mesh");

  glvuWindow();
  
  return 0;
}

//////////////////////////////////////////////////////////////////////////////
// generate an accordian mesh
//////////////////////////////////////////////////////////////////////////////
int tetsSquareAccordian(int argc, char* argv[])
{
  float invRes = 1.0f / res;
  COMPOUND compound;

  float center[] = {0.5, 0.5};
  float caps[] = {invRes, 0.5f};
  compound.addSurface(new BOX(1.0f - 10 * invRes, 10 * invRes,
                             0.5, invRes,
                             1.0f - 10 * invRes, 10 * invRes));

  caps[0] = 2 * invRes;
  caps[1] = 0.5f - invRes;
  compound.subtractSurface(new BOX(1.0f - 11 * invRes, 11 * invRes,
                             0.5 - invRes, invRes,
                             1.0f - 11 * invRes, 11 * invRes));

  for (int x = 2; x < res / 2; x+=3)
  {
    caps[0] = x * invRes;
    caps[1] = (x+1) * invRes;
    compound.addSurface(new BOX(1.0f - 8 * invRes, 8 * invRes,
                               (x+1) * invRes, x * invRes,
                               1.0f - 8 * invRes, 8 * invRes));
  }

  cout << " Generating all tets ... " << endl;
  isoStuffer->useExistingCaches() = false;
  isoStuffer->generateAllTets();
  cout << " Generating inside tets ... " << endl;
  isoStuffer->generateInsideTets(compound);

  cout << " Generating constrained nodes ... " << endl;
  caps[0] = 0;
  caps[1] = 1.25 * invRes;
  CYLINDER constraint(center, 0.5, caps);
  isoStuffer->generateConstrainedNodes(constraint);
  cout << " Outputting final tets file ... " << endl;

  isoStuffer->centerMesh();
  isoStuffer->scaleMesh(4.0f);
  isoStuffer->writeFile("accordian.square.32.mesh");
  isoStuffer->writeFileLegacy("accordian.square.32.legacy.mesh");

  glvuWindow();
  
  return 0;
}

//////////////////////////////////////////////////////////////////////////////
// generate tets from obj and implicit function
//////////////////////////////////////////////////////////////////////////////
int tetsEar(int argc, char* argv[])
{
  //objFile.Load("face.obj");
  objFile.Load("hand.obj");
  objFile.ComputeVertexNormals();
  objFile.normalize(isoStuffer->xRes());
  objFile.setBCCRes(isoStuffer->xRes());

  VEC3 minVecd;
  VEC3 maxVecd;
  objFile.BoundingBox(minVecd, maxVecd);
  VEC3F minVec(minVecd);
  VEC3F maxVec(maxVecd);
  VEC3F halfVec(0.5f, 0.5f, 0.5f);
  VEC3F center;
  center[0] = (maxVec[0] - minVec[0]) * halfVec[0] + minVec[0];
  center[1] = (maxVec[1] - minVec[1]) * halfVec[1] + minVec[1];
  center[2] = (maxVec[2] - minVec[2]) * halfVec[2] + minVec[2];
  cout << " computed center: " << center << endl;
  cout << " translation: " << halfVec - center << endl;
  objFile.translate(halfVec - center);
  cout << "Building accelGrid ...";
  objFile.createAccelGrid();
  cout << " done." << endl;
  cout << "Building distance grid ...";
  objFile.createDistanceGrid(isoStuffer->xRes());
  cout << " done." << endl;

  float invRes = 1.0f / res;
  COMPOUND compound;

  float bottom = invRes;
  float top = 4.0f * invRes;
  compound.addSurface(&objFile);

  cout << " Generating all tets ... " << endl;
  isoStuffer->useExistingCaches() = false;
  isoStuffer->generateAllTets();
  cout << " Generating inside tets ... " << endl;
  isoStuffer->generateInsideTets(compound);

  BOX constraint(0,0,0,0);
  isoStuffer->generateConstrainedNodes(constraint);
  isoStuffer->writeFile("hand.mesh");
  isoStuffer->writeOBJ("hand.tets.obj");

  glvuWindow();
  
  return 0;
}

#endif

//////////////////////////////////////////////////////////////////////////////
// generate a cube
//////////////////////////////////////////////////////////////////////////////
int tetsCube(int argc, char* argv[])
{
  // set default parameters
  int bccRes = 32;
  string triangleMeshPath = string("../examples/");
  string triangleMeshName = string("bunny_watertight.obj");
  string outputPath = string("../tests/bunny/");
  string constrainedAxis = string("negative z");
  bool safeMeshing = true;

  // read in different parameters
  string configName("");
  if (argc > 1)
    configName = string(argv[1]);
  SIMPLE_PARSER configFile(configName); 
  bccRes           = configFile.getInt("bccRes", bccRes);
  triangleMeshPath = configFile.getString("triangle mesh path", triangleMeshPath);
  triangleMeshName = configFile.getString("triangle mesh name", triangleMeshName);
  outputPath       = configFile.getString("output path", outputPath);
  constrainedAxis  = configFile.getString("constrained axis", outputPath);
  safeMeshing      = configFile.getBool("safe meshing", safeMeshing);

  // try to create the output directory
  string mkdir("mkdir ");
  mkdir = mkdir + outputPath;
  system(mkdir.c_str());

  // start the isostuffer
  cout << "==================================================" << endl;
  cout << "RUNNING ISOSURFACE STUFFING" << endl;
  cout << "==================================================" << endl;
  cout << "BCC grid resolution: " << bccRes << "^3" << endl;

  // don't load the mesh, just create a cube
  float invRes = 4.0f / bccRes;
  BOX box(1.0 - invRes, invRes, 1.0 - invRes, invRes, 1.0 - invRes, invRes);

   // generate the mesh
  isoStuffer = new ISO_STUFFER(bccRes, bccRes, bccRes);
  isoStuffer->useExistingCaches() = false;
  if (safeMeshing)
    isoStuffer->generateAllTets();
  else
    isoStuffer->generateLimitedTets(box);
  isoStuffer->generateInsideTets(box);
 
  if (constrainedAxis.compare(string("positive and negative z")) == 0)
  {
    cout << " Pinning maximum and minimum nodes along Z axis." << endl;
    isoStuffer->constrainMinMaxZ();
  }
  else if (constrainedAxis.compare(string("negative z")) == 0)
  {
    cout << " Pinning minimum nodes along Z axis." << endl;
    isoStuffer->constrainMinZ();
  }
  else if (constrainedAxis.compare(string("positive z")) == 0)
  {
    cout << " Pinning maximum nodes along Z axis." << endl;
    isoStuffer->constrainMaxZ();
  }
  else if (constrainedAxis.compare(string("positive and negative y")) == 0)
  {
    cout << " Pinning maximum and minimum nodes along Y axis." << endl;
    isoStuffer->constrainMinMaxY();
  }
  else if (constrainedAxis.compare(string("negative y")) == 0)
  {
    cout << " Pinning minimum nodes along Y axis." << endl;
    isoStuffer->constrainMinY();
  }
  else if (constrainedAxis.compare(string("positive y")) == 0)
  {
    cout << " Pinning maximum nodes along Y axis." << endl;
    isoStuffer->constrainMaxY();
  }
  else if (constrainedAxis.compare(string("positive and negative x")) == 0)
  {
    cout << " Pinning maximum and minimum nodes along X axis." << endl;
    isoStuffer->constrainMinMaxX();
  }
  else if (constrainedAxis.compare(string("negative x")) == 0)
  {
    cout << " Pinning minimum nodes along X axis." << endl;
    isoStuffer->constrainMinX();
  }
  else if (constrainedAxis.compare(string("positive x")) == 0)
  {
    cout << " Pinning maximum nodes along X axis." << endl;
    isoStuffer->constrainMaxX();
  }
  else if (constrainedAxis.compare(string("none")) == 0)
  {
    cout << " Generating unconstrained mesh." << endl;
    isoStuffer->constrainNone();
  }
  else
  {
    cout << "**** No constraint axis was specified! ****" << endl;
  }

  // save the final mesh
  string finalMeshFile = outputPath + triangleMeshName + string(".tetmesh");
  isoStuffer->writeFile(finalMeshFile.c_str());
  
  // stomp any previous vertex and face cache that may have been created
  string faceFile = finalMeshFile + string(".surfacefaces");
  string vertexFile = finalMeshFile + string(".surfacevertices");
  faceFile = string("rm ") + faceFile;
  vertexFile = string("rm ") + vertexFile;
  system(faceFile.c_str());
  system(vertexFile.c_str());

#if _WIN32
  glvuWindow();
#endif
  
  return 0;
}

//////////////////////////////////////////////////////////////////////////////
// Do a conjoined normalization
//////////////////////////////////////////////////////////////////////////////
void normalize(OBJ* first, OBJ* second, int res)
{
  // do a conjoined normalization
  vector<VEC3>& firstVertices = first->vertices;
  vector<VEC3>& secondVertices = second->vertices;

  // first get the center of mass
  VEC3F centerOfMass;
  for (unsigned int x = 0; x < firstVertices.size(); x++)
    centerOfMass += firstVertices[x];
  for (unsigned int x = 0; x < secondVertices.size(); x++)
    centerOfMass += secondVertices[x];
  centerOfMass *= 1.0 / (firstVertices.size() + secondVertices.size());

  // translate everything to the center of mass
  for (unsigned int x = 0; x < firstVertices.size(); x++)
    firstVertices[x] -= centerOfMass;
  for (unsigned int x = 0; x < secondVertices.size(); x++)
    secondVertices[x] -= centerOfMass;

  // find the maximum magnitude
  double maxVal = 0.0f;
  for (unsigned int x = 0; x < firstVertices.size(); x++)
  {
    maxVal = (fabs(firstVertices[x][0]) > maxVal) ? fabs(firstVertices[x][0]) : maxVal;
    maxVal = (fabs(firstVertices[x][1]) > maxVal) ? fabs(firstVertices[x][1]) : maxVal;
    maxVal = (fabs(firstVertices[x][2]) > maxVal) ? fabs(firstVertices[x][2]) : maxVal;
  }
  for (unsigned int x = 0; x < secondVertices.size(); x++)
  {
    maxVal = (fabs(secondVertices[x][0]) > maxVal) ? fabs(secondVertices[x][0]) : maxVal;
    maxVal = (fabs(secondVertices[x][1]) > maxVal) ? fabs(secondVertices[x][1]) : maxVal;
    maxVal = (fabs(secondVertices[x][2]) > maxVal) ? fabs(secondVertices[x][2]) : maxVal;
  }

  // scale everything
  //double scale = 0.5 - 1.25 / res;
  //double scale = 0.5 - 2.0 / res;
  double scale = 0.5 - 4.0 / res;
  for (unsigned int x = 0; x < firstVertices.size(); x++)
    firstVertices[x] *= scale / maxVal;
  for (unsigned int x = 0; x < secondVertices.size(); x++)
    secondVertices[x] *= scale / maxVal;

  // translate everything to 0.5, 0.5, 0.5
  VEC3F half(0.5, 0.5, 0.5);
  for (unsigned int x = 0; x < firstVertices.size(); x++)
    firstVertices[x] += half;
  for (unsigned int x = 0; x < secondVertices.size(); x++)
    secondVertices[x] += half;
}

//////////////////////////////////////////////////////////////////////////////
// mesh the head
//////////////////////////////////////////////////////////////////////////////
void meshHead(int res)
{
  //int bccRes = 64;
  int bccRes = res;

  // load both of the files
  OBJ* headFile = new OBJ();
  headFile->Load("head.obj");
  headFile->ComputeVertexNormals();

  OBJ* skeletonFile = new OBJ();
  skeletonFile->Load("internal.obj");
  skeletonFile->ComputeVertexNormals();

  // do a conjoined normalization
  normalize(headFile, skeletonFile, bccRes);

  // initialize everything else now
  headFile->setBCCRes(bccRes);
  headFile->createAccelGrid();
  //headFile->createDistanceGrid(bccRes);

  skeletonFile->setBCCRes(bccRes);
  skeletonFile->createAccelGrid();
  //skeletonFile->createDistanceGrid(bccRes);
 
  COMPOUND compound;
  compound.addSurface(headFile);
  compound.subtractSurface(skeletonFile);

  compound.bccRes() = bccRes;
  compound.meshingHead() = true;

  isoStuffer = new ISO_STUFFER(bccRes, bccRes, bccRes);

  //isoStuffer->generateAllTets();
  isoStuffer->generateLimitedTets(compound);

  isoStuffer->generateInsideTets(compound);
  isoStuffer->constrainNone();
  isoStuffer->writeFile("head.tetmesh");

  isoStuffer->printTimingBreakdown();
}

//////////////////////////////////////////////////////////////////////////////
// generate test from implicit functions
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  /*
  int res = atoi(argv[1]);
  if (argc != 2)
  {
    cout << " Usage: IsoStuffer <res>" << endl;
    exit(1);
  }
  meshHead(res);
  exit(0);
  */

  /*
  tetsCube(argc, argv);
  exit(0);
  */
 
  // set default parameters
  int bccRes = 32;
  string triangleMeshPath = string("../examples/");
  string triangleMeshName = string("bunny_watertight.obj");
  string outputPath = string("../tests/bunny/");
  string constrainedAxis = string("negative z");
  bool safeMeshing = true;

  // read in different parameters
  string configName("");
  if (argc > 1)
    configName = string(argv[1]);
  SIMPLE_PARSER configFile(configName); 
  bccRes           = configFile.getInt("bccRes", bccRes);
  triangleMeshPath = configFile.getString("triangle mesh path", triangleMeshPath);
  triangleMeshName = configFile.getString("triangle mesh name", triangleMeshName);
  outputPath       = configFile.getString("output path", outputPath);
  constrainedAxis  = configFile.getString("constrained axis", outputPath);
  safeMeshing      = configFile.getBool("safe meshing", safeMeshing);
  bool normalizeMesh = configFile.getBool("normalize mesh", true);

  // try to create the output directory
  string mkdir("mkdir ");
  mkdir = mkdir + outputPath;
  system(mkdir.c_str());

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
  if (normalizeMesh)
    objFile.normalize(bccRes);
  objFile.setBCCRes(bccRes);

  // center the mesh about (0.5, 0.5, 0.5)
  VEC3 minVecd;
  VEC3 maxVecd;
  objFile.BoundingBox(minVecd, maxVecd);
  VEC3F minVec(minVecd);
  VEC3F maxVec(maxVecd);
  VEC3F halfVec(0.5f, 0.5f, 0.5f);
  VEC3F center;
  center[0] = (maxVec[0] - minVec[0]) * halfVec[0] + minVec[0];
  center[1] = (maxVec[1] - minVec[1]) * halfVec[1] + minVec[1];
  center[2] = (maxVec[2] - minVec[2]) * halfVec[2] + minVec[2];
  objFile.translate(halfVec - center);
  
  // create acceleration structures
  cout << " Creating acceleration structures ... " << endl;
  objFile.createAccelGrid();
  objFile.createDistanceGrid(bccRes);
  cout << " done. " << endl;

  // generate the mesh
  isoStuffer = new ISO_STUFFER(bccRes, bccRes, bccRes);
  isoStuffer->useExistingCaches() = false;
  if (safeMeshing)
    isoStuffer->generateAllTets();
  else
  {
    isoStuffer->generateLimitedTets(objFile);
  }
  isoStuffer->generateInsideTets(objFile);

  if (constrainedAxis.compare(string("positive and negative z")) == 0)
  {
    cout << " Pinning maximum and minimum nodes along Z axis." << endl;
    isoStuffer->constrainMinMaxZ();
  }
  else if (constrainedAxis.compare(string("negative z")) == 0)
  {
    cout << " Pinning minimum nodes along Z axis." << endl;
    isoStuffer->constrainMinZ();
  }
  else if (constrainedAxis.compare(string("positive z")) == 0)
  {
    cout << " Pinning maximum nodes along Z axis." << endl;
    isoStuffer->constrainMaxZ();
  }
  else if (constrainedAxis.compare(string("positive and negative y")) == 0)
  {
    cout << " Pinning maximum and minimum nodes along Y axis." << endl;
    isoStuffer->constrainMinMaxY();
  }
  else if (constrainedAxis.compare(string("negative y")) == 0)
  {
    cout << " Pinning minimum nodes along Y axis." << endl;
    isoStuffer->constrainMinY();
  }
  else if (constrainedAxis.compare(string("positive y")) == 0)
  {
    cout << " Pinning maximum nodes along Y axis." << endl;
    isoStuffer->constrainMaxY();
  }
  else if (constrainedAxis.compare(string("positive and negative x")) == 0)
  {
    cout << " Pinning maximum and minimum nodes along X axis." << endl;
    isoStuffer->constrainMinMaxX();
  }
  else if (constrainedAxis.compare(string("negative x")) == 0)
  {
    cout << " Pinning minimum nodes along X axis." << endl;
    isoStuffer->constrainMinX();
  }
  else if (constrainedAxis.compare(string("positive x")) == 0)
  {
    cout << " Pinning maximum nodes along X axis." << endl;
    isoStuffer->constrainMaxX();
  }
  else if (constrainedAxis.compare(string("none")) == 0)
  {
    cout << " Generating unconstrained mesh." << endl;
    isoStuffer->constrainNone();
  }
  else
  {
    cout << "**** No constraint axis was specified! ****" << endl;
  }

  // save the final mesh
  mkdir = string("mkdir ");
  mkdir += outputPath;
  cout << "mkdir command: " << mkdir.c_str() << endl;
  system(mkdir.c_str());
  string finalMeshFile = outputPath + triangleMeshName + string(".tetmesh");
  isoStuffer->writeFile(finalMeshFile.c_str());
  //finalMeshFile = outputPath + triangleMeshName + string(".legacy");
  //isoStuffer->writeFileLegacy(finalMeshFile.c_str());

  string unconstrainedMeshFile = outputPath + triangleMeshName
                               + string(".tetmesh.unconstrained");
  isoStuffer->writeFile(unconstrainedMeshFile.c_str());
  
  // stomp any previous vertex and face cache that may have been created
  string faceFile = finalMeshFile + string(".surfacefaces");
  string vertexFile = finalMeshFile + string(".surfacevertices");
  faceFile = string("rm ") + faceFile;
  vertexFile = string("rm ") + vertexFile;
  system(faceFile.c_str());
  system(vertexFile.c_str());

  bool regressionTesting = false;
  regressionTesting = configFile.getBool("regression testing", regressionTesting);

  // if regression testing, go headless
  if (regressionTesting) return 0;

#ifdef USING_GLVU
  glutInit(&argc, argv);
  glvuWindow();
#endif
  return 0;
}
