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
// "as is" without jxpress or implied warranty.

#include <glvu.hpp>
#include <snapshot.hpp>
#include <iostream>
#include <TET_MESH.h>
#include <SPARSE_MATRIX.h>
#include <STVK.h>
#include <MOONEY_RIVLIN.h>
#include <ARRUDA_BOYCE.h>
#include <NEO_HOOKEAN.h>
#include <INVERTIBLE.h>
#include <SIMPLE_PARSER.h>
#include <FULLSPACE_INTEGRATOR.h>
#include <SKELETON.h>
#include "SUBSPACE_MULTIBODY_INTEGRATOR.h"
#include "PARTITIONED_SKINNED_SUBSPACE_TET_MESH.h"
#include <BD_TREE.h>

using namespace std;
vector<BD_TREE*> bdTrees;
vector<vector<pair<SURFACE*, int> > > collisionPairs;
PARTITIONED_SKINNED_SUBSPACE_TET_MESH* partitionedMesh;
SUBSPACE_MULTIBODY_INTEGRATOR* integrator;
SKELETON* skeleton = NULL;
OBJ* objMesh = NULL;
vector<BOX> boxes;
vector<OBJ> boxObjs;

#define USING_YNEG 0
#define USING_YPOS 0
#define USING_ZNEG 0
#define USING_ZPOS 0
#define USING_PACHINKO 0
#define USING_TEETH 1

//////////////////////////////////////////////////////////////////////////////
// GLOBALS
//////////////////////////////////////////////////////////////////////////////
GLVU glvu;
int windowStartX = 0;
int windowStartY = 0;
//int windowWidth = 700;
//int windowHeight = 700;
int windowWidth = 800;
int windowHeight = 600;
//int windowWidth = 200;
//int windowHeight = 200;
bool mouseClicked = false;
bool animate = true;
//bool animate = false;
bool singleStep = false;
float clickZ;
VEC3F* closest;
int drawPartition = 1;
Real springConst = 2.0;
bool whichDraw = false;

vector<VEC3F*> windVertices;
int bccRes = 32;
Real noiseMultiplier = 0.1;
int maxNewtonSteps = 1;
string previewVideo("preview.avi");
string renderPath("./renders/tmp/");
string odeDataPath("./");
string dataPath("./");
string posePath("./");
string outputPath("");
int steps = 0;
int substeps = 0;
int singleDomain = 2;

int startingFrame = 0;
int totalFrames = 1000;
Real gravityMagnitude = 9.81;
Real collisionStiffness = 2000;
Real collisionDamping = 0.1;
//Real collisionStiffness = 0;
//Real collisionDamping = 0;

bool regressionTesting = false;

////////////////////////////////////////////////////////////////
// Print integer to a zero-padded string
//////////////////////////////////////////////////////////////////
static std::string itoPaddedString(int frame)
{
  char buffer[256];
  sprintf(buffer, "%i", frame);

  std::string number = std::string(buffer);
  if (frame < 10) number = std::string("0") + number;
  if (frame < 100) number = std::string("0") + number;
  if (frame < 1000) number = std::string("0") + number;

  return number;
}

//////////////////////////////////////////////////////////////////////////////
// Dump a frame with this number to this directory
//////////////////////////////////////////////////////////////////////////////
void screenshot(string renderPath, int frame)
{
  FILE *fp;

  string number = itoPaddedString(frame);
  char FileName[256];
  sprintf(FileName,"%srenderGL.%s.ppm", renderPath.c_str(), number.c_str());
  cout << " Dumping frame " << FileName << endl;

  GLint OldReadBuffer;
  glGetIntegerv(GL_READ_BUFFER,&OldReadBuffer);
  glReadBuffer(GL_FRONT);

  GLint OldPackAlignment;
  glGetIntegerv(GL_PACK_ALIGNMENT,&OldPackAlignment); 
  glPixelStorei(GL_PACK_ALIGNMENT,1);

  int WW = glutGet((GLenum)GLUT_WINDOW_WIDTH);
  int WH = glutGet((GLenum)GLUT_WINDOW_HEIGHT);
  int NumPixels = WW*WH;
  GLubyte* Pixels = new GLubyte[NumPixels*3];
  if (Pixels==NULL) { printf("UNABLE TO ALLOC PIXEL READ ARRAY!\n"); return; }
  glReadPixels(0,0,WW,WH,GL_RGB,GL_UNSIGNED_BYTE,Pixels);

  fp = fopen(FileName, "wb");
  fprintf(fp, "P6\n%d %d\n255\n", WW, WH);
  fwrite(Pixels,1,NumPixels*3,fp);
  fclose(fp);
  delete[] Pixels;

  glPixelStorei(GL_PACK_ALIGNMENT,OldPackAlignment);
  glReadBuffer((GLenum)OldReadBuffer);
}

//////////////////////////////////////////////////////////////////////////////
// Read a material file and return a pointer to it
//////////////////////////////////////////////////////////////////////////////
MATERIAL* readMaterial(SIMPLE_PARSER& parser)
{
  // set material
  MATERIAL* material = NULL;
  string materialType;
  materialType = parser.getString("material type", materialType);

  if (materialType.compare("stvk") == 0)
  {
    double lambda = 10.0;
    double mu = 50.0;
    lambda = parser.getFloat("stvk lambda", lambda);
    mu = parser.getFloat("stvk mu", mu);
    material = new STVK(lambda, mu);
    cout << "==================================================" << endl;
    cout << " Material is St-VK, lambda = " << lambda << " mu = " << mu << endl;
  }
  else if (materialType.compare("mooney-rivlin") == 0)
  {
    double mu01 = 100.0;
    double mu10 = 500.0;
    double k = 100000.0;
    mu01 = parser.getFloat("mooney-rivlin mu01", mu01);
    mu10 = parser.getFloat("mooney-rivlin mu10", mu10);
    k = parser.getFloat("mooney-rivlin k", k);
    material = new MOONEY_RIVLIN(mu01, mu10, k);
    cout << "==================================================" << endl;
    cout << " Material is Mooney-Rivlin, mu01 = " << mu01 << " mu10 = " << mu10 << " k = " << k << endl;
  }
  else if (materialType.compare("arruda-boyce") == 0)
  {
    double nkTheta = 5000.0;
    double N  = 5.0;
    double k = 1000.0;
    nkTheta = parser.getFloat("arruda-boyce nktheta", nkTheta);
    N = parser.getFloat("arruda-boyce n", N);
    k = parser.getFloat("arruda-boyce k", k);
    material = new ARRUDA_BOYCE(nkTheta, N, k);
    cout << "==================================================" << endl;
    cout << " Material is Arruda-Boyce, nkTheta = " << nkTheta << " N = " << N << " k = " << k << endl;
  }
  else if (materialType.compare("neo-hookean") == 0)
  {
    double mu = 50.0;
    double lambda = 10.0;
    mu = parser.getFloat("neo-hookean mu", mu);
    lambda = parser.getFloat("neo-hookean lambda", lambda);
    material = new NEO_HOOKEAN(mu, lambda);
    cout << "==================================================" << endl;
    cout << " Material is Neo-Hookean, mu = " << mu << " lambda = " << lambda << endl;
  }
  else
  {
    cout << " *** Material type undefined! *** " << endl;
    exit(1);
  }

  return material;
}

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

    /*
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
    */

    glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
    //partitionedMesh->drawBlendedMesh();
    //partitionedMesh->drawEmbedding();

    glPushMatrix();
      //glTranslatef(-0.5, 0,0);
      partitionedMesh->drawSurfaceFaces();

      //partitionedMesh->drawAllTets();
    glPopMatrix();
 
    glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
    glPointSize(10.0f);

    //objMesh->draw();
    //partitionedMesh->drawBlendedMesh();

    //glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
    //for (unsigned int x = 0; x < boxes.size(); x++)
    //  boxes[x].draw();

    //glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
    glColor4f(0.5, 0.5, 0.5, 1.0);
    for (unsigned int x = 0; x < boxObjs.size(); x++)
      boxObjs[x].draw();

    // always draw the axes on top
    glDisable(GL_DEPTH_TEST);
      glClear(GL_DEPTH);
      //skeleton->drawOdeBones();

      glColor4f(10.0f, 10.0f, 10.0f, 1.0f);
      //partitionedMesh->drawSprings();
      //partitionedMesh->drawRigidFrames();
      glColor4f(10.0f, 0.0f, 0.0f, 1.0f);
      //partitionedMesh->drawConstrainedNodes();
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
//////////////////////////////////////////////////////////////////////////////
void projectPosition()
{
  for (int x = 0; x < partitionedMesh->partitions(); x++)
  {
    SUBSPACE_TET_MESH* subMesh = (SUBSPACE_TET_MESH*)partitionedMesh->mesh(x);
    SUBSPACE_INTEGRATOR* subIntegrator = integrator->integrator(x);
    subMesh->recoverX();
    VECTOR position = subMesh->TET_MESH::x();

    VECTOR groundTruth = position;

    VECTOR projected = subMesh->U() ^ position;
    subMesh->q() = projected;
    subIntegrator->position() = projected;
    subMesh->SUBSPACE_TET_MESH::updateFullMesh();

    //VECTOR diff = groundTruth - (subMesh->U() * projected);
    //cout << " partition " << x << " diff norm: " << diff.norm2() << endl;
    //cout << " projected: " << projected << endl;

    //cout << " U: " << subMesh->U() << endl;
    //exit(0);
  }
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void zeroPositions()
{
  for (int x = 0; x < partitionedMesh->partitions(); x++)
  {
    SUBSPACE_TET_MESH* subMesh = (SUBSPACE_TET_MESH*)partitionedMesh->mesh(x);
    subMesh->q() *= 0;
    subMesh->updateFullMesh();
  }
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
    case 's':
      singleStep = true;
      animate = true;
      break;

    case '1':
      //partitionedMesh->loadMocapFrame(112);
      //loadFullSimulationFrame(111);
      partitionedMesh->loadSkeletonFrame(odeDataPath, 111);
      //projectPosition();
      break;
    case '2':
      singleDomain++;
      if (singleDomain >= partitionedMesh->partitions())
        singleDomain = 0;

      cout << " Selected domain " << singleDomain << endl;
      break;
    case '0':
      projectPosition();
      break;

    case '9':
      zeroPositions();
      break;
    
    case ' ':
      integrator->printTimingBreakdown();
      break;

    case 'v':
      {
        Camera* camera = glvu.GetCurrentCam();
        Vec3f eye;
        Vec3f lookat;
        Vec3f up;
        camera->GetLookAtParams(&eye, &lookat, &up);
        cout << " Eye(" << eye[0] << ", " << eye[1] << ", " << eye[2] << "), " ;
        cout << " LookAtCntr(" << lookat[0] << ", " << lookat[1] << ", " << lookat[2] << "), " ;
        cout << " Up(" << up[0] << ", " << up[1] << ", " << up[2] << ");" << endl;
      }
      break;
    case 'q':
    case 'Q':
      integrator->printTimingBreakdown();
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
  if (animate)
  {
    // if it's the first time, resize the collision vector
    if (collisionPairs.size() == 0)
      collisionPairs.resize(partitionedMesh->partitions());

    int totalPairs = 0;
    for (int x = 0; x < partitionedMesh->partitions(); x++)
    {
      vector<pair<SURFACE*, int> >& pairs = collisionPairs[x];
      pairs.clear();

      for (unsigned int y = 0; y < boxes.size(); y++)
      {
        bdTrees[x]->intersect(&boxes[y], pairs);
        totalPairs += pairs.size();
      }
    }
    cout << " Total collisions: " << totalPairs << endl;

    int frame = startingFrame + (steps % totalFrames);

    // load up the bone transforms and move the constrained nodes
    partitionedMesh->loadOdeSkeletonFrame(odeDataPath, frame);

#if USING_YNEG
    VEC3F down(0,-1,0);
#endif
#if USING_YPOS
    VEC3F down(0,1,0);
#endif
#if USING_ZNEG
    VEC3F down(0,0,-1);
#endif
#if USING_ZPOS
    VEC3F down(0,0,1);
#endif
#if USING_PACHINKO
    VEC3F down(0,-1,0);
#endif
#if USING_TEETH
    VEC3F down(0,-1,0);
#endif

    //Real magnitude = integrator->forceMultiplier();
    Real magnitude = gravityMagnitude;
    integrator->changeGravity(down, magnitude);
    integrator->addGravity();

    //integrator->stepReducedSkinnedDynamic(true);
    //integrator->stepReducedSkinnedDynamicOMP(true);
    //integrator->stepReducedSkinnedQuasistatic(true);
    //integrator->stepReducedSkinnedQuasistaticOMP(true);

    //integrator->stepReducedSkinnedQuasistaticWithCollisions(collisionPairs, true);
    //integrator->stepReducedSkinnedQuasistaticWithCollisionsOMP(collisionPairs, true);

    integrator->stepReducedSkinnedDynamicWithCollisionsOMP(collisionPairs, true);

    //partitionedMesh->updateSubspaceMeshes();

    /*
    partitionedMesh->addRigidsToSurface();
    objMesh->LoadTetSurfaceMesh(partitionedMesh->blendedSurfaceMesh());
    partitionedMesh->subtractRigidsFromSurface();
    */

    // always point camera at the middle of the torso
    Camera* camera = glvu.GetCurrentCam();
    Vec3f eye;
    Vec3f lookat;
    Vec3f up;
    camera->GetLookAtParams(&eye, &lookat, &up);
    Vec3f direction = eye - lookat;
    direction *= 4;
    VEC3F center = skeleton->bone(0)->translation();
    lookat[0] = center[0];
    lookat[1] = center[1];
    lookat[2] = center[2];
    eye = lookat + direction;
    camera->LookAt(eye, lookat, up);

    steps++;

    if (steps % 100 == 0)
      integrator->printTimingBreakdown();

    if (steps == 200 && regressionTesting)
    {
      // dump out state for regression analysis
      string regressionFile = outputPath + string("regression.vector");
      partitionedMesh->writeState(regressionFile.c_str());

      exit(0);
    }

    if (steps == 2000)
    {
      integrator->printTimingBreakdown();
      exit(0);
    }

    if (singleStep) 
    {
      singleStep = false;
      animate = false;
    }
  }
  partitionedMesh->updateSubspaceMeshes();
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
    //integrator->click(point);
    integrator->click(point, singleDomain);
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
int glvuWindow()
{
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  glvu.Init("GLVU Window",
            GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH,
            windowStartX, windowStartY, windowWidth, windowHeight);

  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  GLfloat lightZeroPosition[] = {10.0, 4.0, 10.0, 1.0};
  GLfloat lightZeroColor[] = {0.8, 1.0, 1.0, 1.0};
  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, 1);
  glLightfv(GL_LIGHT0, GL_POSITION, lightZeroPosition);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, lightZeroColor);
 
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);
  glShadeModel(GL_SMOOTH);
  //glClearColor(0,0,0,0);
  glClearColor(1,1,1,1);

  glutDisplayFunc(userDisplayFunc);
  glutMouseFunc(userMouseFunc);
  glutMotionFunc(userMotionFunc);
  glutKeyboardFunc(userKeyboardFunc);
  glutIdleFunc(userIdleFunc);

  Vec3f ModelMin(-10,-10,-10), ModelMax(10,10,10), 
#if USING_ZPOS
        Eye(2.30811, -4.96579, -0.240963),  LookAtCntr(2.06507, -4.00842, -0.0848854),  Up(-0.185364, 0.112098, -0.976255);
#endif
#if USING_YNEG
        Eye(3.54132, 0.550256, -2.96742),  LookAtCntr(2.87068, 0.249502, -2.28935),  Up(-0.28289, 0.948729, 0.141011);
#endif
#if USING_PACHINKO
        Eye(0.0297215, -0.221213, -3.88846),  LookAtCntr(-0.00337817, -0.221785, -2.88901),  Up(-0.0577071, 0.998334, -0.00134061);
#endif
#if USING_TEETH
        Eye(0.0297215, -0.221213, -3.88846),  LookAtCntr(-0.00337817, -0.221785, -2.88901),  Up(-0.0577071, 0.998334, -0.00134061);
#endif
#if USING_ZNEG
        Eye(2.30811, -4.96579, -0.240963),  LookAtCntr(2.06507, -4.00842, -0.0848854),  Up(-0.185364, 0.112098, -0.976255);
        Up *= -1;
#endif
  float Yfov = 45;
  float Aspect = 1;
  float Near = 0.001f;
  float Far = 10.0f;
  glvu.SetAllCams(ModelMin, ModelMax, Eye,LookAtCntr,Up, Yfov,Aspect, Near,Far);

  // reset center
  //Vec3f center(0.5f, 0.5f, 0.5f);
  Vec3f center(0.5f, 0.5f, 0.0f);
  glvu.SetWorldCenter(center);
  
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  glutMainLoop();

  // Control flow will never reach here
  return EXIT_SUCCESS;
}

//////////////////////////////////////////////////////////////////////////////
// Create the stair obstacles
//////////////////////////////////////////////////////////////////////////////
void createStairs()
{
  Real shrink = 0.0;

  // original
#if USING_ZPOS
  // set up the boxes in the scene
  float leftSpacing = 1.25;
  float downSpacing = 1.25;
  for (int x = 0; x < 30; x++)
  {
    VEC3F size(1 - shrink,10 - shrink,0.25 - shrink);
    VEC3F position(-0.5 + x * leftSpacing ,1.9, 2 + x * downSpacing);
    boxes.push_back(BOX(position, size));
  }
#endif
#if USING_YNEG
  // set up the boxes in the scene
  float leftSpacing = 1.25;
  float downSpacing = 1.25;
  for (int x = 0; x < 30; x++)
  {
    VEC3F size(10 - shrink,0.25 - shrink,1 - shrink);
    VEC3F position(0.5,-1.2 - x * downSpacing, 1 - x * leftSpacing);
    boxes.push_back(BOX(position, size));
  }
#endif
#if USING_YPOS
  // set up the boxes in the scene
  float leftSpacing = 1.35;
  float downSpacing = 1.25;
  for (int x = 0; x < 30; x++)
  {
    VEC3F size(10 - shrink,0.25 - shrink,1 - shrink);
    VEC3F position(0.5, 2 + x * downSpacing, 1 - x * leftSpacing);
    boxes.push_back(BOX(position, size));
  }
#endif
#if USING_ZNEG
  // set up the boxes in the scene
  float leftSpacing = 1.25;
  float downSpacing = 1.25;
  for (int x = 0; x < 30; x++)
  {
    VEC3F size(1 - shrink,10 - shrink,0.25 - shrink);
    VEC3F position(-0.5 + x * leftSpacing ,1.9, -1 - x * downSpacing);
    boxes.push_back(BOX(position, size));
  }
#endif
#if USING_PACHINKO
  //float downSpacing = 1.25;
  float downSpacing = 1.65;
  float leftSpacing = downSpacing * 0.5;

  for (int y = 0; y < 10; y++)
    for (int x = 0; x < 10; x++)
    {
      VEC3F size(0.125,0.125,5);
      // (left/right) (up/down) (near/far)

      float jitter = (y % 2) ? 0.5 * downSpacing : 0;

      //VEC3F position(10 - y * leftSpacing,-1.2 - x * downSpacing - jitter, 1);
      VEC3F position(10 - (y + 6) * leftSpacing,-1.2 - x * downSpacing - jitter, 1);

      boxes.push_back(BOX(position, size));
    }
#endif
#if USING_TEETH
  float downSpacing = 1.65;
  float leftSpacing = downSpacing * 0.75;

  for (int y = 0; y < 2; y++)
    for (int x = 0; x < 10; x++)
    {
      VEC3F size(0.55,0.125,5);
      // (left/right) (up/down) (near/far)

      float jitter = (y % 2) ? 0.5 * downSpacing : 0;

      VEC3F position(0.5 -y * leftSpacing,-1.2 - x * downSpacing - jitter, 1);

      boxes.push_back(BOX(position, size));
    }

  {
  // left wall 
  VEC3F size(0.1, 20, 5);
  VEC3F position(0.85 ,-10, 1);
  boxes.push_back(BOX(position, size));

  }

  {
  // right wall 
  VEC3F size(0.1, 20, 5);
  VEC3F position(-1.1 ,-10, 1);
  boxes.push_back(BOX(position, size));

  }

  //bottom wall
  { 
  VEC3F size(1.5, 1, 5);
  VEC3F position((-1 + 0.75) * 0.5 ,-20, 1);
  boxes.push_back(BOX(position, size));

  } 
#endif
  for (unsigned int x = 0; x < boxes.size(); x++)
  {
    boxes[x].collisionStiffness() = collisionStiffness;
    boxes[x].collisionDamping() = collisionDamping;
  }

  // pipe the stairs to OBJs
  boxObjs.resize(boxes.size());
  for (unsigned int x = 0; x < boxes.size(); x++)
  {
    //boxObjs.push_back(OBJ());
    //OBJ& obj = boxObjs.back();
    OBJ& obj = boxObjs[x];

    VEC3F mins = boxes[x].boxMins();
    VEC3F maxs = boxes[x].boxMaxs();

    Real xMinus = mins[0];
    Real yMinus = mins[1];
    Real zMinus = mins[2];
    Real xPlus = maxs[0];
    Real yPlus = maxs[1];
    Real zPlus = maxs[2];

    vector<TRIANGLE*> surfaceTriangles;

    VEC3F v000(xMinus, yMinus, zMinus); 
    VEC3F v100(xPlus, yMinus, zMinus); 
    VEC3F v010(xMinus, yPlus, zMinus); 
    VEC3F v110(xPlus, yPlus, zMinus); 
    VEC3F v001(xMinus, yMinus, zPlus); 
    VEC3F v101(xPlus, yMinus, zPlus); 
    VEC3F v011(xMinus, yPlus, zPlus); 
    VEC3F v111(xPlus, yPlus, zPlus); 

    surfaceTriangles.push_back(new TRIANGLE(&(v010), &(v110), &(v100)));
    surfaceTriangles.push_back(new TRIANGLE(&(v100), &(v000), &(v010)));

    surfaceTriangles.push_back(new TRIANGLE(&(v001), &(v101), &(v111)));
    surfaceTriangles.push_back(new TRIANGLE(&(v111), &(v011), &(v001)));

    surfaceTriangles.push_back(new TRIANGLE(&(v000), &(v100), &(v101)));
    surfaceTriangles.push_back(new TRIANGLE(&(v101), &(v001), &(v000)));

    surfaceTriangles.push_back(new TRIANGLE(&(v011), &(v111), &(v110)));
    surfaceTriangles.push_back(new TRIANGLE(&(v110), &(v010), &(v011)));

    surfaceTriangles.push_back(new TRIANGLE(&(v001), &(v011), &(v010)));
    surfaceTriangles.push_back(new TRIANGLE(&(v010), &(v000), &(v001)));

    surfaceTriangles.push_back(new TRIANGLE(&(v100), &(v110), &(v111)));
    surfaceTriangles.push_back(new TRIANGLE(&(v111), &(v101), &(v100)));

    obj.LoadTetSurfaceMesh(surfaceTriangles);

    for (unsigned int y = 0; y < surfaceTriangles.size(); y++)
      delete surfaceTriangles[y];
  }
}

//////////////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  if (argc != 2)
  {
    cout << " USAGE: " << argv[0] << " *.cfg" << endl;
    return 0;
  }
  // input parameters
  string triangleMeshPath("");
  string triangleMeshName("");
  string tetMeshName("");

  // read in different parameters
  string configName(argv[1]);
  SIMPLE_PARSER configFile(configName); 
  triangleMeshPath = configFile.getString("triangle mesh path", triangleMeshPath);
  triangleMeshName = configFile.getString("triangle mesh name", triangleMeshName);
  outputPath       = configFile.getString("output path", outputPath);
  renderPath       = configFile.getString("render path", renderPath);
  dataPath         = configFile.getString("data path", "./");
  posePath         = configFile.getString("pose path", "./");
  previewVideo     = configFile.getString("preview video", previewVideo);

  // make the render directory
  string mkdir = string("mkdir ") + renderPath;
  system(mkdir.c_str());

  if (configFile.defined("tet mesh name"))
    tetMeshName = configFile.getString("tet mesh name", tetMeshName);
  else
    tetMeshName = outputPath + triangleMeshName + string(".tetmesh");
  tetMeshName = tetMeshName + string(".constrained");

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
  bool invertible = false;
  invertible = configFile.getBool("invertible", invertible);
  for (int x = 0; x < totalMaterials; x++)
  {
    // read in the config file name for the material
    char buffer[256];
    sprintf(buffer, "material %i", x);
    string materialString(buffer);
    string materialFilename("");
    materialFilename = configFile.getString(materialString.c_str(), materialFilename);

    // open the config file
    SIMPLE_PARSER materialFile(materialFilename);
    
    // get the material
    MATERIAL* material = readMaterial(materialFile);
    
    // set the invertible wrapper if necessary
    if (invertible)
    {
      material = new INVERTIBLE(material);
      cout << " Setting material to invertible" << endl;
      float inversionEpsilon = configFile.getFloat("inversion epsilon", ((INVERTIBLE*)material)->epsilon());
      cout << " Inversion epsilon: " << inversionEpsilon << endl;
    }

    materials[x] = material;
  }

  // load in the skinning weights
  string skinningName = configFile.getString("attachment name", "");
  string skeletonName = configFile.getString("skeleton name", "");
  string mocapName = configFile.getString("mocap name", "");

  string skeletonFile = posePath + string("ode.motion.0000.skeleton");
  skeleton = new SKELETON(skeletonFile, true);

  int partitions = skeleton->totalBones();
  cout << " Using " << partitions << " partitions " << endl;

  // output precision
  cout << " Precision: " << sizeof(Real) * 8 << " bit" << endl;
  
  cout << "==================================================" << endl;
  cout << " Reading mesh file " << tetMeshName.c_str() << endl;
  cout << "==================================================" << endl;

  //Real timestep = 1.0f / 120.0f;
  Real timestep = 1.0f / 240.0f;
  timestep = configFile.getFloat("timestep", timestep);

  cout << " Using timestep: " << timestep << endl;

  springConst = 2.0;
  springConst = configFile.getFloat("spring constant", springConst);
  cout << " Using interface spring constant: " << springConst << endl;

  Real dampingConst = 0.01;
  dampingConst = configFile.getFloat("damping constant", dampingConst);
  cout << " Using interface damping constant: " << dampingConst << endl;

  Real rayleighAlpha = 0.01;
  Real rayleighBeta = 0.00001;

  rayleighAlpha = configFile.getFloat("rayleigh alpha", rayleighAlpha);
  rayleighBeta = configFile.getFloat("rayleigh beta", rayleighBeta);
  
  Real interfaceDampingConst = 0.0;
  interfaceDampingConst = configFile.getFloat("interface damping constant", interfaceDampingConst);
  cout << " Using interface damping constant: " << interfaceDampingConst << endl;
  
  partitionedMesh = new PARTITIONED_SKINNED_SUBSPACE_TET_MESH(skeleton, tetMeshName.c_str(), materials, totalMaterials, partitions, springConst, true, true, tetMeshName);
  integrator = new SUBSPACE_MULTIBODY_INTEGRATOR(partitionedMesh, timestep, springConst, dampingConst, rayleighAlpha, rayleighBeta);
  partitionedMesh->updateFullMeshes();

  //partitionedMesh->printConnectivity();
  //exit(0);

  /*
  // try piping the surface mesh to an OBJ
  objMesh = new OBJ();
  partitionedMesh->addRigidsToSurface();
  objMesh->LoadTetSurfaceMesh(partitionedMesh->blendedSurfaceMesh());
  partitionedMesh->subtractRigidsFromSurface();
  */

  // precompute sandwiches that use the ODE rest poses
  integrator->precomputeOdeSandwiches(dataPath);
  //integrator->precomputeOdeSandwiches("./data/ode.armadillo.poses/");

  VECTOR::printVertical = false;
  for (int x = 0; x < partitions; x++)
  {
    integrator->integrator(x)->position() = ((SUBSPACE_TET_MESH*)partitionedMesh->mesh(x))->q();
    integrator->integrator(x)->positionOld() = ((SUBSPACE_TET_MESH*)partitionedMesh->mesh(x))->q();
  }

  // factor to increase the mouse-input force by
  double forceMultiplier = 0.1;
  forceMultiplier = configFile.getFloat("force multiplier", forceMultiplier);
  for (int x = 0; x < partitions; x++)
  {
    integrator->forceMultiplier() = forceMultiplier;
    integrator->integrator(x)->forceMultiplier() = forceMultiplier;
  }

  // set the Newton parameters
  //int newtonDigits = 3;
  int newtonDigits = 2;
  newtonDigits = configFile.getInt("Newton digits of precision", newtonDigits);
  Real precision = pow(10.0, (Real)(-newtonDigits));
  integrator->solverEps() = precision;

  //integrator->maxNewtonSteps() = 10;
  integrator->maxNewtonSteps() = 3;

  // get the gravity magnitude
  gravityMagnitude = configFile.getFloat("gravity magnitude", gravityMagnitude);
  cout << " gravity magnitude: " << gravityMagnitude << endl;

  // set the gravity direction
  VEC3F down(0,-1,0);
  integrator->changeGravity(down, gravityMagnitude);

  // get the simulation frame numbers
  startingFrame = configFile.getInt("starting snapshot", startingFrame);
  totalFrames = configFile.getInt("snapshots", totalFrames);

  collisionStiffness = configFile.getFloat("collision stiffness", collisionStiffness);
  cout << " Using collision stiffness " << collisionStiffness << endl;

  // are we regression testing?
  regressionTesting = configFile.getBool("regression testing", regressionTesting);

#if USING_PACHINKO
  totalFrames = 1500;
  odeDataPath = string("./data/armadillo.ragdoll.pachinko/");
#endif
#if USING_TEETH
  totalFrames = 2000;
  //collisionStiffness = 100;
  odeDataPath = string("./data/armadillo.ragdoll.teeth/");
#endif
#if USING_YNEG
  totalFrames = 1000;
  odeDataPath = string("./data/armadillo.ragdoll.yNeg/");
#endif

  // build the BD-Trees
  int maxDepth = 12;
  for (int x = 0; x < partitions; x++)
  {
    SUBSPACE_TET_MESH* tetMesh = (SUBSPACE_TET_MESH*)(partitionedMesh->mesh(x));
    BD_TREE* newTree = new BD_TREE(tetMesh, maxDepth);
    newTree->unconstrained() = partitionedMesh->unconstrained(x);
    bdTrees.push_back(newTree);
  }

  // build the stairs
  createStairs();

  // gather interface statistics
  int totalInterfaces = 0;
  int totalClones = 0;
  for (int x = 0; x < partitions; x++)
    for (int y = x; y < partitions; y++)
    {
      int clones = partitionedMesh->totalClones(x,y);
      if (clones > 0)
        totalInterfaces++;
      totalClones += clones;
    }

  int totalCubaturePoints = 0;
  for (int x = 0; x < partitions; x++)
    totalCubaturePoints += ((SUBSPACE_TET_MESH*)partitionedMesh->mesh(x))->totalKeyTets();

  cout << " Total interfaces: " << totalInterfaces << endl;
  cout << " Total clones: " << totalClones << endl;
  cout << " Total key tets: " << totalCubaturePoints << endl;

  cout << "====================" << endl;
  cout << " Max Newton steps: " << integrator->maxNewtonSteps() << endl;
  cout << " Newton digits: " << newtonDigits << endl;
  cout << " Timestep: " << timestep << endl;
  cout << " Total ODE frames: " << totalFrames << endl;
  cout << "====================" << endl;

  glutInit(&argc, argv);
  glvuWindow();
}
