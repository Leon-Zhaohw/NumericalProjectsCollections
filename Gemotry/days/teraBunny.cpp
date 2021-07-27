// OpenGL Graphics includes
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/freeglut.h>
#include <GL/glu.h>
#endif
#include <float.h>

// Includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstdio>

#include <iostream>
#include "TIMER.h"
#include "QUATERNION.h"
#include "POLYNOMIAL_4D.h"
#include "TRIANGLE_MESH.h"
#include "SIMPLE_PARSER.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////////
// Meshes and integrators
//////////////////////////////////////////////////////////////////////////////
TRIANGLE_MESH* triangleMesh = NULL;
Real marchingSurface = 0;
string outputPrefix("./temp/temp");

float minValue;
float transformValue = 0;
float maxValue;
int steps;
int maxIterations = 1;
int res = 97;

Real isosurface = 0;
double checkpointFrequency = 60 * 60;

//////////////////////////////////////////////////////////////////////////////
// Shape information
//////////////////////////////////////////////////////////////////////////////
int xRes, yRes, zRes;
VEC3F center, lengths;
POLYNOMIAL_4D top, bottom;
double expScaling;
FIELD_3D fractal;
FIELD_3D distanceField;
FIELD_3D curvatureField;
Real quaternionSlice;

///////////////////////////////////////////////////////////////////////
// read in optimized shape
///////////////////////////////////////////////////////////////////////
void read(const string& filename)
{
  TIMER functionTimer(__FUNCTION__);
  FILE* file;
  file = fopen(filename.c_str(), "rb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " OPTIMIZE_3D read failed! " << endl;
    cout << " Could not open file " << filename.c_str() << endl;
    exit(0);
  }
  cout << " Reading file " << filename.c_str() << " ... " << flush;

  // read dimensions
  fread((void*)&xRes, sizeof(int), 1, file);
  fread((void*)&yRes, sizeof(int), 1, file);
  fread((void*)&zRes, sizeof(int), 1, file);
  cout << " Read in res: " << xRes << " " << yRes << " " << zRes << endl;

  center.read(file);
  lengths.read(file);

  top.read(file);
  bottom.read(file);
  cout << " Top powers sum to:    " << top.powerSum() << endl;
  cout << " Bottom powers sum to: " << bottom.powerSum() << endl;
  cout << " Scaled Top powers sum to:    " << top.scaledPowerSum() << endl;
  cout << " Scaled Bottom powers sum to: " << bottom.scaledPowerSum() << endl;

  fread((void*)&expScaling, sizeof(double), 1, file);

  int size;
  fread((void*)&size, sizeof(int), 1, file);

  fractal.read(file);
  distanceField.read(file);
  curvatureField.read(file);

  quaternionSlice = 0;

  fclose(file);
  cout << " done. " << endl;
}

///////////////////////////////////////////////////////////////////////
// read in a slim version of the optimized shape
///////////////////////////////////////////////////////////////////////
void readSlimfile(const string& filename)
{
  TIMER functionTimer(__FUNCTION__);
  FILE* file;
  file = fopen(filename.c_str(), "rb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " OPTIMIZE_3D slim read failed! " << endl;
    cout << " Could not open file " << filename.c_str() << endl;
    exit(0);
  }
  cout << " Reading slim file " << filename.c_str() << " ... " << flush;

  // read dimensions
  fread((void*)&xRes, sizeof(int), 1, file);
  fread((void*)&yRes, sizeof(int), 1, file);
  fread((void*)&zRes, sizeof(int), 1, file);
  cout << " Read in res: " << xRes << " " << yRes << " " << zRes << endl;

  center.read(file);
  lengths.read(file);
  cout << " Read in center: " << center << endl;
  cout << " Read in lengths: " << lengths << endl;

  top.read(file);
  bottom.read(file);

  fread((void*)&expScaling, sizeof(double), 1, file);
  quaternionSlice = 0;

  fclose(file);
  cout << " done. " << endl;
  
  cout << " Top powers sum to:    " << top.powerSum() << endl;
  cout << " Bottom powers sum to: " << bottom.powerSum() << endl;
  cout << " Scaled Top powers sum to:    " << top.scaledPowerSum() << endl;
  cout << " Scaled Bottom powers sum to: " << bottom.scaledPowerSum() << endl;
  cout << " Exponential scaling: " << expScaling << endl;
  cout << " Top mean root:    " << top.meanRootPosition() << endl;
  cout << " Bottom mean root: " << bottom.meanRootPosition() << endl;
}

///////////////////////////////////////////////////////////////////////
// write out a text version of the optimized shape
///////////////////////////////////////////////////////////////////////
void writeTextfile(const string& filename)
{
  TIMER functionTimer(__FUNCTION__);
  FILE* file;
  file = fopen(filename.c_str(), "w");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " OPTIMIZE_3D slim write failed! " << endl;
    cout << " Could not open file " << filename.c_str() << endl;
    exit(0);
  }
  cout << " Writing text file " << filename.c_str() << " ... " << flush;

  // write dimensions
  fprintf(file, "res = %i %i %i\n", xRes, yRes, zRes);

  VEC3F& fractalCenter = fractal.center();
  VEC3F& fractalLengths = fractal.lengths();
  fprintf(file, "center = %.20e %.20e %.20e\n", fractalCenter[0], fractalCenter[1], fractalCenter[2]);
  fprintf(file, "lengths = %.20e %.20e %.20e\n", fractalLengths[0], fractalLengths[1], fractalLengths[2]);
  top.writeTextfile(file);
  bottom.writeTextfile(file);
  fprintf(file, "exp scaling = %.20e\n", expScaling);
  fclose(file);
  cout << " done. " << endl;
}

///////////////////////////////////////////////////////////////////////
// read in a text version of the optimized shape
///////////////////////////////////////////////////////////////////////
void readTextfile(const string& filename)
{
  TIMER functionTimer(__FUNCTION__);
  FILE* file;
  file = fopen(filename.c_str(), "r");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " OPTIMIZE_3D slim read failed! " << endl;
    cout << " Could not open file " << filename.c_str() << endl;
    exit(0);
  }
  cout << " Reading text file " << filename.c_str() << " ... " << flush;

  // read dimensions
  fscanf(file, "res = %i %i %i\n", &xRes, &yRes, &zRes);
  cout << " Read in res: " << xRes << " " << yRes << " " << zRes << endl;

  fscanf(file, "center = %lf %lf %lf\n", &center[0], &center[1], &center[2]);
  fscanf(file, "lengths = %lf %lf %lf\n", &lengths[0], &lengths[1], &lengths[2]);
  cout << " Read in center: " << center << endl;
  cout << " Read in lengths: " << lengths << endl;

  top.readTextfile(file);
  bottom.readTextfile(file);

  fscanf(file, "exp scaling = %lf\n", &expScaling);
  quaternionSlice = 0;
  fclose(file);
  cout << " done. " << endl;
  
  cout << " Top powers sum to:    " << top.powerSum() << endl;
  cout << " Bottom powers sum to: " << bottom.powerSum() << endl;
  cout << " Scaled Top powers sum to:    " << top.scaledPowerSum() << endl;
  cout << " Scaled Bottom powers sum to: " << bottom.scaledPowerSum() << endl;
  cout << " Exponential scaling: " << expScaling << endl;
  cout << " Top mean root:    " << top.meanRootPosition() << endl;
  cout << " Bottom mean root: " << bottom.meanRootPosition() << endl;
}

///////////////////////////////////////////////////////////////////////
// write out a slimmer version of the optimized shape
///////////////////////////////////////////////////////////////////////
void writeSlimfile(const string& filename)
{
  TIMER functionTimer(__FUNCTION__);
  FILE* file;
  file = fopen(filename.c_str(), "wb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " OPTIMIZE_3D slim write failed! " << endl;
    cout << " Could not open file " << filename.c_str() << endl;
    exit(0);
  }
  cout << " Writing slim file " << filename.c_str() << " ... " << flush;

  // write dimensions
  fwrite((void*)&xRes, sizeof(int), 1, file);
  fwrite((void*)&yRes, sizeof(int), 1, file);
  fwrite((void*)&zRes, sizeof(int), 1, file);

  VEC3F& fractalCenter = fractal.center();
  VEC3F& fractalLengths = fractal.lengths();
  fractalCenter.write(file);
  fractalLengths.write(file);

  top.write(file);
  bottom.write(file);
  fwrite((void*)&expScaling, sizeof(double), 1, file);
  quaternionSlice = 0;

  fclose(file);
  cout << " done. " << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void computeFractalHuge(int res, const std::string objFilename)
{
  TIMER functionTimer(__FUNCTION__);
 
  if (triangleMesh) delete triangleMesh;
  string cachePath;
  char buffer[256];
  sprintf(buffer, ".res.%i.iterations.%i.cache", res, maxIterations);
  string cacheFilename = outputPrefix + string(buffer);

  triangleMesh = new TRIANGLE_MESH(fractal.center(),
                                   fractal.lengths(),
                                   VEC3I(res,res,res),
                                   top, 
                                   bottom, 
                                   expScaling, 
                                   maxIterations, 
                                   quaternionSlice, 
                                   isosurface,
                                   fractal.rotation(), 
                                   cacheFilename,
                                   objFilename,
                                   checkpointFrequency);
  TIMER::printTimings();
}

//////////////////////////////////////////////////////////////////////////////
// For the tooth: 
//    0 to 1.0 along (0,0,1,0) gives interesting results
//
//    0 to -1.0 along (0,0,1,0) also gives interesting results,
//      not the same as 0 to 1.0
//
//    0 to 2.0 along (0,1,0,0) makes it disappear eventually, though
//      it doesn't seem as interesting 
//
//    0 to -0.75 along (0,1,0,0) also makes it disappear, and it gets interesting
//
//    0 to 1.0 along (1,0,0,0) makes it disappear
//
// BUNNY:
//    0 to 2.0 along (1,0,0,0) makes it disappear
//    0 to 2.0 along (0,1,0,0) makes it disappear, even more quickly than x
//    0 to 2.0 along (0,0,1,0) makes it disappear, also very quickly
//////////////////////////////////////////////////////////////////////////////
void marchTranslations(const VEC3F& translation)
{
  float delta = transformValue;
  VEC3F      vDelta = translation;

  QUATERNION qDelta(vDelta[0], vDelta[1], vDelta[2], 0);
  cout << " Using the translation delta: " << delta <<  " quaternion version: " << qDelta << endl;

  //const VEC3F centerOriginal = optimize3D.fractal().center();
  const VEC3F centerOriginal = fractal.center();

  // back up the original root positions
  const POLYNOMIAL_4D topOriginal = top;
  const POLYNOMIAL_4D bottomOriginal = bottom;

  for (unsigned int i = 0; i < topOriginal.roots().size(); i++)
    top.rootsMutable()[i] = topOriginal.roots()[i] + qDelta;
  for (unsigned int i = 0; i < bottomOriginal.roots().size(); i++)
    bottom.rootsMutable()[i] = bottomOriginal.roots()[i] + qDelta;

  //optimize3D.fractal().center() = centerOriginal + vDelta;
  fractal.center() = centerOriginal + vDelta;

  cout << " New delta: " << qDelta << endl;

  // build the final output filename
  char buffer[256];
  sprintf(buffer, ".res.%i.iterations.%i.obj", res, maxIterations);
  string outputFilename = outputPrefix + string(buffer);

  // compute the fractal
  computeFractalHuge(res, outputFilename);
}

//////////////////////////////////////////////////////////////////////////////
// Read in the optimized shape from the original QUIJIBO format
//////////////////////////////////////////////////////////////////////////////
void readQuijiboFile(const string& filename, SIMPLE_PARSER& parser, VEC3F& rootTranslation)
{
  read(filename.c_str());
  writeSlimfile(filename + string(".slim"));
  writeTextfile(filename + string(".txt"));
  
  VEC3F& fractalCenter = fractal.center();
  VEC3F& fractalLengths = fractal.lengths();
  QUATERNION& rotation = fractal.rotation();

  VEC3 centerNew  = parser.getVector3("field center", fractalCenter);
  VEC3 lengthsNew = parser.getVector3("field lengths", fractalLengths);
  QUATERNION rotationNew = parser.getQuaternion("rotation", rotation);
  checkpointFrequency = parser.getInt("checkpoint frequency", checkpointFrequency);

  cout << " Original center: " << fractalCenter << endl;
  cout << " Original lengths: " << fractalLengths << endl;
  cout << " Original rotation: " << rotation << endl;
  cout << " New center: " << centerNew << endl;
  cout << " New lengths: " << lengthsNew << endl;
  cout << " New rotation: " << rotationNew << endl;
  cout << " Checkpoint frequency (in seconds): " << checkpointFrequency << endl;

  bool recenter = parser.getBool("recenter", false);

  // need to translate the center by the root translation
  if (recenter)
  {
    cout << " Recentering " << centerNew << " to ";
    centerNew -= rootTranslation;
    cout << centerNew << endl;
  }
  else
    cout << " Not recentering" << endl;

  fractalCenter = centerNew;
  fractalLengths = lengthsNew;
  rotation = rotationNew;

  cout << " Allocating initial fields ... " << flush;
  fractal = FIELD_3D(100,100,100, fractalCenter, fractalLengths);
  fractal.rotation() = rotationNew;
  cout << " done. " << endl;
}

//////////////////////////////////////////////////////////////////////////////
// Read in a flat text version of the optimized shape
//////////////////////////////////////////////////////////////////////////////
void readTextfile(const string& filename, SIMPLE_PARSER& parser, VEC3F& rootTranslation)
{
  readTextfile(filename + string(".txt"));
  
  QUATERNION rotation(1,0,0,0);
  VEC3 centerNew  = parser.getVector3("field center", center);
  VEC3 lengthsNew = parser.getVector3("field lengths", lengths);
  QUATERNION rotationNew = parser.getQuaternion("rotation", rotation);
  checkpointFrequency = parser.getInt("checkpoint frequency", checkpointFrequency);

  cout << " New center: " << centerNew << endl;
  cout << " New lengths: " << lengthsNew << endl;
  cout << " New rotation: " << rotationNew << endl;
  cout << " Checkpoint frequency (in seconds): " << checkpointFrequency << endl;

  bool recenter = parser.getBool("recenter", false);

  // need to translate the center by the root translation
  if (recenter)
  {
    cout << " Recentering " << centerNew << " to ";
    centerNew -= rootTranslation;
    cout << centerNew << endl;
  }
  else
    cout << " Not recentering" << endl;
  rotation = rotationNew;

  cout << " Allocating initial fields ... " << flush;
  fractal = FIELD_3D(100,100,100, center,lengths);
  fractal.rotation() = rotationNew;
  cout << " done. " << endl;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  if (argc != 2)
  {
    cout << " USAGE: " << argv[0] << " <config file>" << endl;
    return 0;
  }

  string cfg(argv[1]);
  SIMPLE_PARSER parser(cfg);

  string filename       = parser.getString("optimization file", string(""), true);
  VEC3F rootTranslation = parser.getVector3("root translation", VEC3(0,0,0));
  outputPrefix     = parser.getString("output prefix", string("./temp/temp"));
  res              = parser.getInt("field res", 50);
  maxIterations    = parser.getInt("max iterations", 3);

  cout << " Using input file:       " << filename.c_str() << endl;
  cout << " Using root translation: " << rootTranslation << endl;
  cout << " Using output prefix:    " << outputPrefix.c_str() << endl;
  cout << " Using resolution:       " << res << endl;
  cout << " Using max iterations:   " << maxIterations << endl;

  //readQuijiboFile(filename, parser, rootTranslation);
  readTextfile(filename, parser, rootTranslation);

  cout << " Rendering root Z translations " << endl;
  marchTranslations(rootTranslation);

  return 0;
}
