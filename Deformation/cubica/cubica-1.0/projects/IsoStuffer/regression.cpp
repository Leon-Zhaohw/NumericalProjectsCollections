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
#include <SIMPLE_PARSER.h>
#include "obj.h"
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
// generate test from implicit functions
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
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
  objFile.createAccelGrid();
  objFile.createDistanceGrid(bccRes);

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
  
  return 0;
}
