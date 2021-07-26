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
// GLVU : CopyrightPartition 1997 - 2002 
//        The University of North Carolina at Chapel Hill
//------------------------------------------------------------------------------
// Permission to use, copy, modify, distribute and sell this software and its 
// documentation for any purpose is hereby granted without fee, provided that 
// the above copyrightPartition notice appear in all copies and that both that copyrightPartition 
// notice and this permission notice appear in supporting documentation. 
// Binaries may be compiled with this software without any royalties or 
// restrictions. 
//
// The University of North Carolina at Chapel Hill makes no representations 
// about the suitability of this software for any purpose. It is provided 
// "as is" without express or implied warranty.

#include <glvu.hpp>
#include <TET_MESH.h>
#include <PARTITIONED_TET_MESH.h>
#include <SIMPLE_PARSER.h>
#include <OBJ.h>
#include <SKELETON.h>

OBJ* objMesh = NULL;
SKELETON* skeleton= NULL;
TET_MESH* tetMesh = NULL;
PARTITIONED_TET_MESH* partitionedMesh = NULL;

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
int partitionToDraw = -1;
int totalPartitions = 16;
VEC3F* clickedVertex = NULL;

int leftPartition = 0;
int rightPartition = 1;
bool drawLimitedExploded = false;

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
    /*
    if (drawLimitedExploded)
      partitionedMesh->drawExploded(leftPartition, rightPartition);
    else
    */
    partitionedMesh->drawExploded();

    // draw the clicked node
    glPointSize(5.0f);
    if (clickedVertex != NULL)
    {
      glDisable(GL_DEPTH_TEST);
      glColor4f(100.0f, 100.0f, 100.0f, 100.0f);
        vector<int>& membership = tetMesh->tetMembership(clickedVertex);
        TET& tet = tetMesh->tets()[membership[0]];
        tet.drawTriangles();

        // get the one ring
        vector<VEC3F*> oneRing;
        tetMesh->oneRing(membership[0], oneRing);
        glBegin(GL_POINTS);
          glColor4f(100.0f, 0.0f, 0.0f, 100.0f);
          for (int x = 0; x < oneRing.size(); x++)
            glVertex3dv(*oneRing[x]);
        glEnd();
      glEnable(GL_DEPTH_TEST);
    }

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
  static int steps = 0;
  
  switch(Key)
  {
    case 'e':
      drawLimitedExploded = !drawLimitedExploded;
      break;

    case '1':
      leftPartition++;
      if (leftPartition == rightPartition)
        leftPartition++;
      leftPartition = leftPartition % totalPartitions;
      rightPartition = rightPartition % totalPartitions;
      cout << " Left partition: " << leftPartition << endl;
      cout << " Right partition: " << rightPartition << endl;
      break;

    case '2':
      rightPartition++;
      if (leftPartition == rightPartition)
        rightPartition++;
      leftPartition = leftPartition % totalPartitions;
      rightPartition = rightPartition % totalPartitions;
      cout << " Left partition: " << leftPartition << endl;
      cout << " Right partition: " << rightPartition << endl;
      break;

    case 'q':
    case 'Q':
      exit(0);
      break;

    case '=':
      partitionToDraw++;
      if (partitionToDraw == totalPartitions)
        partitionToDraw = -1;
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
    clickedVertex = tetMesh->closestSurfaceNode(point);
    mouseClicked = true;
    return;
  }
  if (button == GLUT_LEFT_BUTTON && 
      state == GLUT_UP)
  {
    mouseClicked = false;
    //clickedVertex = NULL;
    //return;
  }

  // pass through to default handler
  glvu.Mouse(button,state,x,y);
}

//////////////////////////////////////////////////////////////////////////////
// open the GLVU window
//////////////////////////////////////////////////////////////////////////////
int glvuWindow()
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
  if (argc <= 1)
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
  string tetMeshName("");
  string triangleMeshName("");
  string triangleMeshPath("");
  string posePath("");

  // read in different parameters
  string configName(argv[1]);
  SIMPLE_PARSER configFile(configName); 
  triangleMeshName = configFile.getString("triangle mesh name", triangleMeshName);
  outputPath       = configFile.getString("output path", outputPath);
  triangleMeshPath = configFile.getString("triangle mesh path", triangleMeshPath);
  posePath         = configFile.getString("pose path", posePath);

  if (configFile.defined("tet mesh name"))
    tetMeshName = configFile.getString("tet mesh name", tetMeshName);
  else
    tetMeshName = outputPath + triangleMeshName + string(".tetmesh");

  // see if the constrained version exists
  string constrainedName = tetMeshName + string(".constrained");
  FILE* probe = fopen(constrainedName.c_str(), "rb");
  if (probe != NULL)
    tetMeshName = constrainedName;
  else
  {
    cout << " Skeleton has not been embedded! Run PinocchioTetViewer and press 'c'!" << endl;
    exit(0);
  }
  fclose(probe);

  // output precision
  cout << " Precision: " << sizeof(Real) * 8 << " bit" << endl;
  
  cout << "==================================================" << endl;
  cout << " Reading mesh file " << tetMeshName.c_str() << endl;
  cout << "==================================================" << endl;
  tetMesh = new TET_MESH(tetMeshName.c_str(), NULL, 0);

  // stomp any previous partitioning
  string rm("rm ");
  rm = rm + tetMeshName + ".partition*";
  cout << " Removing any previous partitioning: " << rm.c_str() << endl;
  system(rm.c_str());
  rm = string("rm ") + tetMeshName + ".dimacs*";
  cout << " Removing any graph coloring: " << rm.c_str() << endl;
  system(rm.c_str());

  // load up the obj and the skeleton
  Real meshScale           = configFile.getFloat("mesh scale", 1.0);
  string originalMeshFile = triangleMeshPath + triangleMeshName;
  objMesh = new OBJ();
  objMesh->Load(originalMeshFile.c_str(), meshScale);
  objMesh->ComputeVertexNormals();

  // load in the skeleton
  string skeletonFile = posePath + string("ode.motion.0000.skeleton");
  skeleton = new SKELETON(skeletonFile, true);
  skeleton->tetMesh() = tetMesh;
  skeleton->buildOdeSkinning(tetMesh);

  cout << "==================================================" << endl;
  cout << " Performing ODE partitioning " << endl;
  cout << "==================================================" << endl;
 
  // do the partitioning
  skeleton->odeSkinningPartition();

  vector<TET>& tets = tetMesh->tets();
  int maxPartition = 0;
  for (unsigned int x = 0; x < tets.size(); x++)
    if (tets[x].partition > maxPartition)
      maxPartition = tets[x].partition;
  totalPartitions = maxPartition + 1;

  // write out final partitioned mesh
  tetMesh->writePartitions();
  partitionedMesh = new PARTITIONED_TET_MESH(tetMeshName.c_str(), NULL, 0, totalPartitions);
 
  vector<int> totalValence;
  totalValence.resize(totalPartitions);
  for (int x = 0; x < totalPartitions; x++)
    for (int y = 0; y <= x; y++)
      if (partitionedMesh->neighbors(x,y))
      {
        totalValence[x]++;
        totalValence[y]++;
      }

  for (int x = 0; x < totalPartitions; x++)
    cout << " Partition " << x << " neighbors: " << totalValence[x] << endl;

  cout << " Drawing ... " << endl;
  flush(cout);
  if (!headless)
  {
    glutInit(&argc, argv);
    glvuWindow();
  }
  delete tetMesh;
  delete partitionedMesh;
}
