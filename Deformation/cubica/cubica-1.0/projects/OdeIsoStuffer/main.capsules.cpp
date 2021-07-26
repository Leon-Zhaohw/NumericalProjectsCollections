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

    glPushMatrix();
      glTranslatef(1,1,1);
      compound->draw();
    glPopMatrix();
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
// build the compound rag doll
//////////////////////////////////////////////////////////////////////////////
void buildCompoundRagdoll(SKELETON* skeleton, COMPOUND* compound)
{
  vector<BONE*>& bones = skeleton->bones();

  for (unsigned int x = 0; x < bones.size(); x++)
  {
    Real length = skeleton->boneLength(x);
    Real radius = skeleton->boneRadius(x);
    float center[] = {0,0,0};
    float caps[] = {-length * 0.5, length * 0.5};

    CYLINDER* cylinder = new CYLINDER(center, radius, caps);

    // make it point along z axis
    QUATERNION& rotation = cylinder->cylinderRotation();
    rotation = QUATERNION(MATRIX3::rotation(VEC3F(1,0,0), M_PI / 2)) * rotation;

    // add the ragdoll rotation
    rotation = bones[x]->quaternion() * rotation;

    // add the ragdoll translation
    VEC3F& translation = cylinder->cylinderTranslation();
    translation = bones[x]->translation();

    compound->addSurface(cylinder);

    VEC3F beginVertex(0,0, length * 0.5);
    VEC3F endVertex(0,0, -length * 0.5);
    MATRIX3 sphereRotation = bones[x]->quaternion().toExplicitMatrix3x3();
    beginVertex = sphereRotation * beginVertex + translation;
    endVertex = sphereRotation * endVertex + translation;

    // add the end cap spheres
    float leftCenter[] = {beginVertex[0], beginVertex[1], beginVertex[2]};
    SPHERE* leftCap = new SPHERE(radius, leftCenter);

    float rightCenter[] = {endVertex[0], endVertex[1], endVertex[2]};
    SPHERE* rightCap = new SPHERE(radius, rightCenter);

    compound->addSurface(leftCap);
    compound->addSurface(rightCap);
  }

  // get the bounding box
  VEC3F mins;
  VEC3F maxs;
  compound->boundingBox(mins, maxs);
  meshCenter = (mins + maxs) * 0.5;

  // get the scaling
  meshScale = maxs[0] - mins[0];
  for (int x = 0; x < 3; x++)
    if (maxs[x] - mins[x] > meshScale)
      meshScale = maxs[x] - mins[x];

  meshScale = 0.95 / meshScale;

  // scale after the fact
  vector<SURFACE*>& surfaces = compound->surfaces();
  for (unsigned int x = 0; x < surfaces.size(); x++)
  {
    if (surfaces[x]->type().compare("SPHERE") == 0)
    {
      SPHERE* sphere = (SPHERE*)surfaces[x];
      sphere->radius() *= meshScale;

      sphere->center() -= meshCenter;
      sphere->center() *= meshScale;
      sphere->center() += VEC3F(0.5, 0.5, 0.5);
    } 
    if (surfaces[x]->type().compare("CYLINDER") == 0)
    {
      CYLINDER* cylinder = (CYLINDER*)surfaces[x];
      Real radius = cylinder->radius();
      radius *= meshScale;
      cylinder->setRadius(radius);

      cylinder->caps()[0] *= meshScale;
      cylinder->caps()[1] *= meshScale;

      cylinder->cylinderTranslation() -= meshCenter;
      cylinder->cylinderTranslation() *= meshScale;
      cylinder->cylinderTranslation() += VEC3F(0.5, 0.5, 0.5);
    } 
  }
}

//////////////////////////////////////////////////////////////////////////////
// Rescale the model and the isostuffing results back to the original size
//////////////////////////////////////////////////////////////////////////////
void rescaleResults()
{
  // rescale model
  Real invScale = 1.0 / meshScale;

  vector<SURFACE*>& surfaces = ((COMPOUND*)compound)->surfaces();
  for (unsigned int x = 0; x < surfaces.size(); x++)
  {
    if (surfaces[x]->type().compare("SPHERE") == 0)
    {
      SPHERE* sphere = (SPHERE*)surfaces[x];
      sphere->radius() *= invScale;

      sphere->center() -= VEC3F(0.5, 0.5, 0.5);
      sphere->center() *= invScale;
      sphere->center() += meshCenter;
    } 
    if (surfaces[x]->type().compare("CYLINDER") == 0)
    {
      CYLINDER* cylinder = (CYLINDER*)surfaces[x];
      Real radius = cylinder->radius();
      radius *= invScale;
      cylinder->setRadius(radius);

      cylinder->caps()[0] *= invScale;
      cylinder->caps()[1] *= invScale;

      cylinder->cylinderTranslation() -= VEC3F(0.5, 0.5, 0.5);
      cylinder->cylinderTranslation() *= invScale;
      cylinder->cylinderTranslation() += meshCenter;
    } 
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
  if (argc != 2)
  {
    cout << " Usage: " << argv[0] << " <*.cfg>" << endl;
    return 0;
  }

  // set default parameters
  int bccRes = 64;

  // read in different parameters
  string configName("");
  if (argc > 1)
    configName = string(argv[1]);
  SIMPLE_PARSER configFile(configName); 

  string outputPath = string("./meshes/tet/ragdoll/");
  string dataPath = string("./data/ragdoll/");

  bccRes      = configFile.getInt("bccRes", bccRes);
  outputPath  = configFile.getString("output path", outputPath);
  dataPath    = configFile.getString("data path", dataPath);

  // start the isostuffer
  cout << "==================================================" << endl;
  cout << "RUNNING ISOSURFACE STUFFING" << endl;
  cout << "==================================================" << endl;
  cout << "BCC grid resolution: " << bccRes << "^3" << endl;

  // read in the skeleton
  string skeletonFile = dataPath + string("ode.motion.0000.skeleton");
  SKELETON* skeleton = new SKELETON(skeletonFile, true);

  // build the compound and send it to the stuffer
  compound = new COMPOUND();

  buildCompoundRagdoll(skeleton, (COMPOUND*)compound);

  // generate the mesh
  isoStuffer = new ISO_STUFFER(bccRes, bccRes, bccRes);
  isoStuffer->useExistingCaches() = false;
  isoStuffer->generateLimitedTets(*compound);
  isoStuffer->generateInsideTets(*compound);

  cout << " Generating unconstrained mesh." << endl;
  isoStuffer->constrainNone();

  // rescale everything back up
  rescaleResults();

  // make the output dir
  string mkdir("mkdir ");
  mkdir += outputPath;
  cout << "mkdir command: " << mkdir.c_str() << endl;
  system(mkdir.c_str());

  string finalMeshFile = outputPath + string("ode.tetmesh");
  isoStuffer->writeFile(finalMeshFile.c_str());

#ifdef USING_GLVU
  glutInit(&argc, argv);
  glvuWindow();
#endif
  return 0;
}
