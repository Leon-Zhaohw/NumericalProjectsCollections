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

//////////////////////////////////////////////////////////////////////////////
// GLOBALS
//////////////////////////////////////////////////////////////////////////////
int bccRes = 32;
Real noiseMultiplier = 0.1;
int maxNewtonSteps = 1;
string previewVideo("preview.avi");
string renderPath("./renders/tmp/");
string dataPath("./");
string posePath("./");

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
  string outputPath("");
  string tetMeshName("");

  // read in different parameters
  string configName(argv[1]);
  SIMPLE_PARSER configFile(configName); 
  triangleMeshPath = configFile.getString("triangle mesh path", triangleMeshPath);
  triangleMeshName = configFile.getString("triangle mesh name", triangleMeshName);
  outputPath       = configFile.getString("output path", outputPath);
  renderPath       = configFile.getString("render path", renderPath);
  dataPath         = configFile.getString("data path", dataPath);
  posePath         = configFile.getString("pose path", posePath);
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

  Real springConst = 2.0;
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

  // precompute sandwiches that use the ODE rest poses
  integrator->precomputeOdeSandwiches(dataPath);
}
