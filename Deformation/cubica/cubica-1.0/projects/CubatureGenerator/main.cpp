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
#include <TET_MESH.h>
#include <SUBSPACE_TET_MESH.h>
#include <CUBATURE_GENERATOR.h>
#include <SPARSE_MATRIX.h>
#include <STVK.h>
#include <MOONEY_RIVLIN.h>
#include <ARRUDA_BOYCE.h>
#include <INVERTIBLE.h>
#include <NEO_HOOKEAN.h>
#include <SIMPLE_PARSER.h>

using namespace std;
SUBSPACE_TET_MESH* tetMesh;

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

  string configName(argv[1]);
  SIMPLE_PARSER configFile(configName);

  // input parameters
  string triangleMeshPath;
  string triangleMeshName;
  string outputPath;
  string tetMeshName;

  // read in different parameters
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
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;

  // set training parameters
  int rank                 = configFile.getInt("rank", 20);
  int trainingSamples      = configFile.getInt("training samples", 100);
  double trainingMagnitude = configFile.getFloat("training magnitude", 1.0);
  int maxKeyTets           = configFile.getInt("max key tets", 1000);
  int candidatesPerTry     = configFile.getInt("candidates per try", 100);
  double errorTolerance    = configFile.getFloat("error tolerance", 0.01);
  int randSeed             = configFile.getInt("randomization seed", 123456);

  Real meshScale           = configFile.getFloat("mesh scale", 1.0);

  cout << "Got mesh scale " << meshScale << endl;

  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  Real *somethingIsVeryWrongHere = new Real[30];

  cout << "==================================================" << endl;
  cout << " Reading mesh file " << tetMeshName.c_str() << endl;
  cout << "==================================================" << endl;
  tetMesh = new SUBSPACE_TET_MESH(tetMeshName.c_str(), materials, totalMaterials,
                                  false, NULL, NULL, meshScale);
  CUBATURE_GENERATOR generator(tetMesh);
  generator.maxKeyTets() = maxKeyTets;
  generator.candidatesPerTry() = candidatesPerTry;
  generator.errorTolerance() = errorTolerance;
  generator.randSeed(randSeed);

  // check that a basis was loaded
  if (!tetMesh->basisLoaded())
  {
    cout << "==================================================" << endl;
    cout << " Generating basis of rank = " << rank << endl;
    cout << "==================================================" << endl;
    SPARSE_MATRIX& stiffness = tetMesh->generateSparseStiffnessMatrix();

    string outputFile = outputPath + "stiffnessMatrix.m";
    stiffness.writeToMatlab(outputFile, "stiffness");

    generator.generateLMABasis(rank);
    //generator.generateLMABasisArpack(rank);
    //exit(0);
  }

  // check that a training set was not already loaded
  //if (generator.trainingSetSize() == 0)
  {
    cout << "==================================================" << endl;
    cout << " Generating training samples" << endl;
    cout << "==================================================" << endl;
    generator.generateTrainingSamples(trainingSamples, trainingMagnitude);
  }
  cout << "==================================================" << endl;
  cout << " Training cubature" << endl;
  cout << "==================================================" << endl;
  generator.generateCubature();
  generator.writeCubature();

  delete tetMesh;
  delete materials;
}
