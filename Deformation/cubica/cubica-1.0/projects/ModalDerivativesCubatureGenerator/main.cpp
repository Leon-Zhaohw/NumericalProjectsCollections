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
// Perform PCA on the modal derivatives and LMA basis
//////////////////////////////////////////////////////////////////////////////
void modalDerivativesPCA(MATRIX& modalDerivatives, VECTOR& modalEigenvalues)
{
  // get the U basis
  MATRIX& LMA = tetMesh->U();
  VECTOR& eigenvalues = tetMesh->eigenvalues();

  // create a vector of all the combined eigenvalues
  VECTOR weights(eigenvalues.size() + modalEigenvalues.size());
  for (int x = 0; x < eigenvalues.size(); x++)
    weights[x] = eigenvalues[x];
  for (int x = 0; x < modalEigenvalues.size(); x++)
    weights[eigenvalues.size() + x] = modalEigenvalues[x];

  // create a matrix of all the combined modes
  MATRIX UBasis(LMA.rows(), LMA.cols() + modalDerivatives.cols());
  for (int x = 0; x < LMA.rows(); x++)
  {
   for (int y = 0; y < LMA.cols(); y++)
     UBasis(x,y) = LMA(x,y);
   for (int y = 0; y < modalDerivatives.cols(); y++)
     UBasis(x, LMA.cols() + y) = modalDerivatives(x,y);
  }

  // reweight according to eigenvalues so that PCA will see it
  UBasis.normalizeColumns();
  for (int j = 0; j < UBasis.cols(); j++)
    for (int i = 0; i < UBasis.rows(); i++)
      // weight the LMA modes and the modal derivatives differently
      if (j < eigenvalues.size())
        UBasis(i,j) *= weights[0] / weights[j];
      else
        UBasis(i,j) *= weights[0] * weights[0] / weights[j];

  VECTOR values;
  MATRIX components;
  UBasis.PCA(components, values);

  int significant = 0;
  vector<int> significantIndices;
  for (int y = 0; y < values.size(); y++)
    if (values(y) > 1e-4 && !isnan(values(y)))
    {
      significant++;
      significantIndices.push_back(y);
    }
  cout << " The mesh has " << significant << " of " << values.size() << " significant components" << endl;
  if (significant > LMA.cols())
    significant = LMA.cols();
  cout << " Using " << significant << " of these components" << endl;
  cout << " PCA values: " << values << endl;

  // build the new basis
  MATRIX newBasis(UBasis.rows(), significant);
  for (int j = 0; j < UBasis.rows(); j++)
    for (int i = 0; i < significant; i++)
      newBasis(j,i) = components(j,significantIndices[i]);

  // build the new principal values
  VECTOR newValues(significant);
  for (int i = 0; i < significant; i++)
    newValues[i] = values[significantIndices[i]];
 
  // square and invert them for training set generation 
  for (int i = 0; i < significant; i++)
    newValues[i] = 1.0 / (newValues[i] * newValues[i]);

  tetMesh->U() = newBasis;
  tetMesh->eigenvalues() = newValues;

  cout << " new eigenvalues: " << newValues << endl;
}

//////////////////////////////////////////////////////////////////////////////
// This version uses the (slow) Slepc solver. A faster version is available
// in ModalDerivatives
//////////////////////////////////////////////////////////////////////////////
void computeModalDerivatives(MATRIX& modalDerivatives, VECTOR& packedEigenvalues)
{
  string tetMeshName = tetMesh->filename();

  MATRIX& U = tetMesh->U();
  VECTOR& eigenvalues = tetMesh->eigenvalues();

  // generate stiffness just once
  cout << " Generating stiffness matrix ... "; flush(cout);
  SPARSE_MATRIX& temp = tetMesh->generateSparseStiffnessMatrix(true);
  SPARSE_PETSC_MATRIX stiffness(temp.rows(), temp.cols());
  stiffness.setSparsity(temp);
  stiffness.equals(temp);
  stiffness.useICC();
  stiffness.usePCG();
  stiffness.eps() = 1e-7;
  stiffness.maxIterations() = 1000;
  cout << "done." << endl;

  // precompute how many of derivatives we're going to compute
  int totalDerivatives = 0;
  int totalPossible = 0;
  Real discard = 0.01;
  //Real discard = 0.001;
  for (int x = 0; x < U.cols(); x++)
    for (int y = x; y < U.cols(); y++)
    {
      Real weight = eigenvalues[0] * eigenvalues[0] / (eigenvalues[x] * eigenvalues[y]);
      if (weight > discard)
        totalDerivatives++;
      totalPossible++;
    }
  cout << "Discarding all derivatives less than " << discard << endl;
  cout << "Keeping " << totalDerivatives << " of " << totalPossible << " possible derivatives " << endl;

  // compute modal derivatives
  modalDerivatives.resizeAndWipe(U.rows(), totalDerivatives);
  vector<Real> modalEigenvalues;
  int totalCols = U.cols();
  int current = 0;
  cout << " Generating modal derivatives ... " << endl;
  for (int x = 0; x < totalCols; x++)
    for (int y = x; y < totalCols; y++)
    {
      Real weight = eigenvalues[0] * eigenvalues[0] / (eigenvalues[x] * eigenvalues[y]);
      if (weight < discard) continue;

      VECTOR firstVector = U.getColumn(x);
      VECTOR secondVector = U.getColumn(y);

      VECTOR hessianProduct = tetMesh->hessianProduct(firstVector, secondVector);
      hessianProduct *= -1.0;

      VECTOR solution = hessianProduct;
  
      // solve via PCG
      stiffness.solveCG(solution, hessianProduct);

      // store the mode and its eigenvalue
      modalDerivatives.setColumn(solution, current);
      modalEigenvalues.push_back(eigenvalues[x] * eigenvalues[y]);

      cout << current << " of " << totalDerivatives << " modal derivatives computed (" << 100 * ((Real)current / totalDerivatives) << "%)" << endl;
      current++;
    }

  // pack the modes into a matrix for output
  string derivativeName = tetMeshName;
  derivativeName = derivativeName + (".modal.derivatives.matrix");
  cout << " Writing modal derivatives to " << derivativeName.c_str() << endl;
  modalDerivatives.write(derivativeName.c_str());

  // pack the eigenvalues into a vector for output
  packedEigenvalues = VECTOR(modalEigenvalues);
  string eigenName = tetMeshName;
  eigenName = eigenName + (".modal.eigenvalues.vector");
  cout << " Writing modal eigenvalues to " << eigenName.c_str() << endl;
  packedEigenvalues.write(eigenName.c_str());
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
// 
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  if (argc <= 1)
  {
    cout << " USAGE: " << argv[0] << " *.cfg " << endl;
    return 0;
  }

  // input parameters
  string triangleMeshPath;
  string triangleMeshName;
  string outputPath;
  string tetMeshName;

  // read in different parameters
  string configName(argv[1]);
  SIMPLE_PARSER configFile(configName); 
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

  // set training parameters
  int rank                 = configFile.getInt("rank", 20);
  int trainingSamples      = configFile.getInt("training samples", 100);
  double trainingMagnitude = configFile.getFloat("training magnitude", 1.0);
  int maxKeyTets           = configFile.getInt("max key tets", 1000);
  int candidatesPerTry     = configFile.getInt("candidates per try", 100);
  double errorTolerance    = configFile.getFloat("error tolerance", 0.01);
  int randSeed             = configFile.getInt("randomization seed", 123456);

  cout << "==================================================" << endl;
  cout << " Reading mesh file " << tetMeshName.c_str() << endl;
  cout << "==================================================" << endl;
  tetMesh = new SUBSPACE_TET_MESH(tetMeshName.c_str(), materials, totalMaterials);
  CUBATURE_GENERATOR generator(tetMesh);
  generator.maxKeyTets() = maxKeyTets;
  generator.candidatesPerTry() = candidatesPerTry;
  generator.errorTolerance() = errorTolerance;
  generator.randSeed(randSeed);

  MATRIX modalDerivatives;
  VECTOR modalEigenvalues;

  // check that a basis was loaded
  if (!tetMesh->basisLoaded())
  {
    cout << "==================================================" << endl;
    cout << " Generating basis of rank = " << rank << endl;
    cout << "==================================================" << endl;
    SPARSE_MATRIX& stiffness = tetMesh->generateSparseStiffnessMatrix();

    string outputFile = outputPath + "stiffnessMatrix.m";
    stiffness.writeToMatlab(outputFile, "stiffness");

    generator.generateLMABasis(rank, false);
    computeModalDerivatives(modalDerivatives, modalEigenvalues);
  }
  else
  {
    // load the modal derivatives
    string derivativeName = tetMeshName;
    derivativeName = derivativeName + (".modal.derivatives.matrix");
    cout << " Looking for modal derivatives in " << derivativeName.c_str() << " ... "; flush(cout);
    bool fileFound = modalDerivatives.read(derivativeName.c_str());
    if (fileFound)
    {
      cout << " found." << endl;
      derivativeName = tetMeshName;
      derivativeName = derivativeName + (".modal.eigenvalues.vector");
      cout << " Looking for modal eigenvalues in " << derivativeName.c_str() << " ... "; flush(cout);
      fileFound = modalEigenvalues.read(derivativeName.c_str());
      if (fileFound)
        cout << "found." << endl;
      else
      {
        cout << "not found." << endl;
        cout << " No modal derivative eigenvalues found! Rerun ModalDerivatives! " << endl;
        exit(0);
      }
    }
    else
    {
      cout << " No modal derivatives found! Rerun ModalDerivatives! " << endl;
      exit(0);
    }
  }

  // see if the modal derivatives have been PCA'd yet
  string pcaCheck = tetMeshName + (".pca.peformed");
  FILE* check;
  check = fopen(pcaCheck.c_str(), "r");
  if (check == NULL)
  {
    // perform PCA on the modal derivatives
    modalDerivativesPCA(modalDerivatives, modalEigenvalues);
  
    // make sure PCA doesn't get applied twice by accident
    fclose(check);
    check = fopen(pcaCheck.c_str(), "wb");
    int checked = 1;
    fwrite((void*)&checked, sizeof(int), 1, check);
    fclose(check);

    // stomp the old basis
    string valuename = tetMeshName + string(".eigenvalues.vector");
    string vectorname = tetMeshName + string(".eigenvectors.matrix");
    cout << " Writing eigenvalues to file " << valuename.c_str() << endl;
    cout << " Writing eigenvectors to file " << vectorname.c_str() << endl;
    tetMesh->eigenvalues().write(valuename.c_str());
    tetMesh->U().write(vectorname.c_str());
  }

  // check that a training set was not already loaded
  cout << "==================================================" << endl;
  cout << " Generating training samples" << endl;
  cout << "==================================================" << endl;
  generator.generateTrainingSamples(trainingSamples, trainingMagnitude);

  cout << "==================================================" << endl;
  cout << " Training cubature" << endl;
  cout << "==================================================" << endl;
  generator.generateCubature();
  generator.writeCubature();

  delete tetMesh;
  delete materials;
}
