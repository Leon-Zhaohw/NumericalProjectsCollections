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
#include <SPARSE_MATRIX.h>
#include <STVK.h>
#include <MOONEY_RIVLIN.h>
#include <ARRUDA_BOYCE.h>
#include <NEO_HOOKEAN.h>
#include <INVERTIBLE.h>
#include <SIMPLE_PARSER.h>
#include <FULLSPACE_INTEGRATOR.h>
#include <SKELETON.h>
#include <BOX.h>
#include <SPHERE_TREE.h>

using namespace std;

//////////////////////////////////////////////////////////////////////////////
// Meshes and integrators
//////////////////////////////////////////////////////////////////////////////
OBJ* objMesh = NULL;
SKELETON* skeleton= NULL;
TET_MESH* tetMesh = NULL;
FULLSPACE_INTEGRATOR* integrator = NULL;
string tetMeshName("");
string dataPath("./");
string posePath("./");
string renderPath("./");

//////////////////////////////////////////////////////////////////////////////
// Various function declarations
//////////////////////////////////////////////////////////////////////////////
void generateSurfaceBasis();
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
//static bool animate = false;
static bool animate = true;
static int currentBone = 0;
static int totalFrames = 0;
static Real gravityMagnitude = 9.81;

static int frame = 0;

static int startingFrame = 0;
static int endingFrame = -1;

//Real collisionStiffness = 10000;
Real collisionStiffness = 1000;
Real collisionDamping = 0.1;

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

GLfloat mat_diffuse[] = { 0.5, 0.2, 0.3, 1.0 };
GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };

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

    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);

    glColor4f(1,1,1,1);
    tetMesh->drawSurfaceFaces();

    glDisable(GL_DEPTH_TEST);
      glLineWidth(1.0);
      skeleton->drawOdeBones();

      glColor4f(1,0,0,1);
      tetMesh->drawConstrainedNodes();
    glEnable(GL_DEPTH_TEST);

    glColor4f(1,1,1,1);

    /*
    for (unsigned int x = 0; x < boxes.size(); x++)
      boxes[x].draw();
      */

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
  static bool firstStep = true;

  if (animate)
  {
    // update the frame
    char buffer[256];
    sprintf(buffer, "%04i", frame);
    string skeletonFile = posePath + string("ode.motion.") + string(buffer) + string(".skeleton");
    skeleton->loadOdeFrame(skeletonFile.c_str());
    skeleton->updateOdeSkinning(false);

    // set position to zero so the solve doesn't have to struggle its way
    // back from the previous solution
    if (integrator->position().size() > 0)
    {
      integrator->position() *= 0;
      tetMesh->x() = integrator->position();
    }

    integrator->stepSkinnedQuasistatic();
    tetMesh->updateFullMesh();

    // dump out simulation data
    string integratorFile = dataPath + string("full.integrator.");
    integratorFile += itoPaddedString(frame);
    integratorFile += string(".state");
    integrator->writeState(integratorFile);
    tetMesh->writeDeformedMesh(frame, dataPath);

    //if (frame >= totalFrames)
    if (frame >= endingFrame)
      exit(0);
    /*
    if (frame >= totalFrames)
    {
      integrator->printTimingBreakdown();

      animate = !animate;

      // convert everything to jpg
      cout << " Converting to JPG  ... "; flush(cout);
      string convert("mogrify -flip -format jpg -quality 100 ");
      convert += renderPath + string("*.ppm");
      system(convert.c_str());
      cout << " done." << endl;

      cout << " Stomping PPMs ... "; flush(cout);
      string stomp("rm ");
      stomp += renderPath + string("*.ppm");
      system(stomp.c_str());
      cout << " done." << endl;

      cout << " Calling FFMpeg" << endl;

      string ffmpeg("ffmpeg -f image2 -i ");
      ffmpeg += renderPath + string("renderGL.%04d.jpg -r 100 -qcomp 0 -y -b 4000k ");
      ffmpeg += renderPath + string("groundtruth.avi");
      system(ffmpeg.c_str());

      cout << " Stomping JPGs ... "; flush(cout);
      stomp = string("rm ");
      stomp += renderPath + string("*.jpg");
      system(stomp.c_str());
      cout << " done." << endl;
      exit(1);
    }
    */

    frame++;
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
        //Eye(0.5,0.5,3), LookAtCntr(0.5,0.5,0.5), Up(0,1,0);
        //Eye(0.5,0.5,5), LookAtCntr(0.5,0.5,0.5), Up(0,1,0);
        //Eye(0.335209, 1.50887, 4.01166),  LookAtCntr(0.343554, 1.42032, 3.01562),  Up(-0.0176097, 0.995903, -0.0886895);
        //Eye(0.446743, 0.623206, 2.37347),  LookAtCntr(0.472317, 0.557065, 1.37599),  Up(0.0379395, 0.997156, -0.0651464);
        Eye(0.0933924, -0.393629, -2.00392),  LookAtCntr(0.044833, -0.263856, -1.01357),  Up(-0.00896117, 0.991427, -0.130353);
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
  if (argc <= 1)
  {
    cout << " USAGE: " << argv[0] << " *.cfg <starting frame number> <ending frame number>" << endl;
    return 0;
  }
  // input parameters
  string triangleMeshPath("");
  string triangleMeshName("");
  string outputPath("");

  // read in different parameters
  string configName(argv[1]);
  SIMPLE_PARSER configFile(configName); 
  triangleMeshPath = configFile.getString("triangle mesh path", triangleMeshPath);
  triangleMeshName = configFile.getString("triangle mesh name", triangleMeshName);
  outputPath       = configFile.getString("output path", outputPath);
  renderPath       = configFile.getString("render path", renderPath);
  dataPath         = configFile.getString("data path", dataPath);
  posePath         = configFile.getString("pose path", posePath);

  // read in how many materials there are
  int totalMaterials = 0;
  totalMaterials = configFile.getInt("total materials", totalMaterials);
  if (totalMaterials == 0)
  {
    cout << " NO MATERIALS SPECIFIED!!!!" << endl;
    return false;
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
    MATERIAL* material = SIMPLE_PARSER::READ_MATERIAL(materialFile);
    
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

  // read in the tet mesh
  if (configFile.defined("tet mesh name"))
    tetMeshName = configFile.getString("tet mesh name", tetMeshName);
  else
    tetMeshName = outputPath + triangleMeshName + string(".tetmesh");

  // see if the constrained version exists
  cout << " tet mesh name: " << tetMeshName.c_str() << endl;
  string constrainedName = tetMeshName + string(".constrained");
  FILE* file = fopen(constrainedName.c_str(), "rb");
  if (file != NULL)
    tetMeshName = constrainedName;
  else
  {
    cout << " Skeleton has not been embedded! Run OdeTetViewer and press 'C'!" << endl;
    exit(0);
  }
  fclose(file);

  tetMesh = new TET_MESH(tetMeshName.c_str(), materials, totalMaterials, true);

  // output precision
  cout << " Precision: " << sizeof(Real) * 8 << " bit" << endl;

  // Get the mesh to embed
  int bccRes = 32;
  bccRes = configFile.getInt("bccRes", bccRes);

  // load up the skeleton
  string skeletonFile = posePath + string("ode.motion.0000.skeleton");
  skeleton = new SKELETON(skeletonFile, true);
  skeleton->tetMesh() = tetMesh;
  skeleton->buildOdeSkinning(tetMesh);
 
  // create the integrator
  Real rayleighAlpha = 0.01;
  Real rayleighBeta = 0.00001;

  Real timestep = 1.0 / 120.0;
  timestep = configFile.getFloat("timestep", timestep);

  cout << " Using timestep: " << timestep << endl;
  integrator = new FULLSPACE_INTEGRATOR(tetMesh, timestep,
                                        rayleighAlpha, rayleighBeta);

  // old settings
  integrator->PCGDigits() = 8;
  integrator->maxNewtonSteps() = 75;
  //integrator->PCGDigits() = 6;
  //integrator->maxNewtonSteps() = 10;
  integrator->solverEps() = 1e-2;
  
  integrator->solver().useMINRES();
  integrator->solver().useJacobi();

  // factor to increase the mouse-input force by
  double forceMultiplier = 0.1;
  forceMultiplier = configFile.getFloat("force multiplier", forceMultiplier);
  integrator->forceMultiplier() = forceMultiplier;

  gravityMagnitude = configFile.getFloat("gravity magnitude", gravityMagnitude);

  Real moveSpeed = 1.0;
  Real eyeDistanceScale = 1.0;

  moveSpeed = configFile.getFloat("navigation speed", moveSpeed);
  eyeDistanceScale = configFile.getFloat("view distance", eyeDistanceScale);

  maxForceRadius = configFile.getFloat("max force radius", maxForceRadius);

  cout << "Using maxForceRadius = " << maxForceRadius << endl;

  // get the total number of frames to simulate
  totalFrames = configFile.getInt("snapshots", totalFrames);
  cout << " Simulating " << totalFrames << " frames " << endl;
  startingFrame = 0;
  endingFrame = totalFrames;

  // add stairs to the scene
  collisionStiffness = configFile.getFloat("collision stiffness", collisionStiffness);
  collisionDamping = configFile.getFloat("collision damping", collisionDamping);

  cout << " Using collision stiffness: " << collisionStiffness << endl;
  cout << " Using collision damping: " << collisionDamping << endl;
  //createStairs();

  // create the render and data directories
  string mkdirRender = string("mkdir ") + renderPath;
  system(mkdirRender.c_str());
  string mkdirData = string("mkdir ") + dataPath;
  system(mkdirData.c_str());

  // if there are user specified frame numbers, use them
  if (argc == 4)
  {
    startingFrame = atoi(argv[2]);
    endingFrame = atoi(argv[3]);

    frame = startingFrame;
    totalFrames = endingFrame - startingFrame;
  }

  glvu.SetMoveSpeed( moveSpeed * glvu.GetMoveSpeed() );

  glutInit(&argc, argv);

  Vec3f boxCenter(0.5, 0.5, 0.5);
  glvuWindow( boxCenter, eyeDistanceScale );
}
