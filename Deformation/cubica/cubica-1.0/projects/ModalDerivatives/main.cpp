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

#include <umfpack.h>

using namespace std;
TET_MESH* tetMesh;

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void computeModalDerivatives(string tetMeshName)
{
  // load LMA modes
  MATRIX U;
  string LMAvectors = tetMeshName + string(".eigenvectors.matrix");
  FILE* file = fopen(LMAvectors.c_str(), "rb");
  if (file == NULL)
  {
    cout << " No eigenvalues found! Run ArpackLMA first! " << endl;
    exit(1);
  }
  fclose(file);
  U.read(LMAvectors.c_str());

  // load LMA eigenvalues
  VECTOR eigenvalues;
  string LMAvalues = tetMeshName + string(".eigenvalues.vector");
  file = fopen(LMAvalues.c_str(), "rb");
  if (file == NULL)
  {
    cout << " No eigenvalues found! Run ArpackLMA first! " << endl;
    exit(1);
  }
  fclose(file);
  eigenvalues.read(LMAvalues.c_str());

  // generate stiffness just once
  cout << " Generating stiffness matrix ... "; flush(cout);
  SPARSE_MATRIX& stiffness = tetMesh->generateSparseStiffnessMatrix(true);
  cout << "done." << endl;

  // populate Umfpack version of stiffness matrix
  cout << " Passing matrix to Umfpack ..."; flush(cout);
  vector<int> rows;
  vector<int> cols;
  vector<Real> vals;

  stiffness.entries(rows, cols, vals);
  int nz = vals.size();
  int n = stiffness.cols();

  int* Ap = new int[n + 1];
  int* Ai = new int[nz];
  double* Ax = new double[nz];

  int currentRow = -1;
  int entryAp = 0;
  for (int x = 0; x < nz; x++)
  {
    Ai[x] = cols[x];
    Ax[x] = vals[x];
    if (rows[x] != currentRow)
    {
      Ap[entryAp] = x;
      currentRow = rows[x];
      entryAp++;
    }
  }
  Ap[n] = nz;
  cout << " done." << endl;
  
  // clearing stiffness matrix
  stiffness.clear();

  // factorize stiffness matrix
  cout << " Factoring matrix with Umfpack ..."; flush(cout);
  void* symbolic;
  void* numeric;
  umfpack_di_symbolic(n, n, Ap, Ai, Ax, &symbolic, NULL, NULL);
  umfpack_di_numeric(Ap, Ai, Ax, symbolic, &numeric, NULL, NULL);
  umfpack_di_free_symbolic(&symbolic);
  cout << " done." << endl;

  // precompute how many of derivatives we're going to compute
  int totalDerivatives = 0;
  int totalPossible = 0;
  Real discard = 0.01;
  //Real discard = 0.001;
  //Real discard = 0.0001;
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
  
  // if it's computing more than a thousand derivatives, clamp it
  if (totalDerivatives > 1000)
  {
    cout << " Clamping total derivatives to 1000 to save time " << endl;
    totalDerivatives = 1000;
  }

  // compute modal derivatives
  MATRIX modalDerivatives(U.rows(), totalDerivatives);
  vector<Real> modalEigenvalues;
  int totalCols = U.cols();
  int current = 0;
  cout << " Generating modal derivatives ... " << endl;
  for (int x = 0; x < totalCols; x++)
  {
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
      umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, solution.data(), hessianProduct.data(), numeric, NULL, NULL);

      // store the mode and its eigenvalue
      modalDerivatives.setColumn(solution, current);
      modalEigenvalues.push_back(eigenvalues[x] * eigenvalues[y]);

      cout << current << " of " << totalDerivatives << " modal derivatives computed (" << 100 * ((Real)current / totalDerivatives) << "%)" << endl;
      current++;

      if (current == totalDerivatives) break;
    }
    if (current == totalDerivatives) break;
  }

  // pack the modes into a matrix for output
  string derivativeName = tetMeshName;
  derivativeName = derivativeName + (".modal.derivatives.matrix");
  cout << " Writing modal derivatives to " << derivativeName.c_str() << endl;
  modalDerivatives.write(derivativeName.c_str());

  // pack the eigenvalues into a vector for output
  VECTOR packedEigenvalues(modalEigenvalues);
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
  //runTest();
  //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  //exit(0);

  if (argc <= 1)
  {
    cout << " USAGE: " << argv[0] << " *.cfg " << endl;
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

  if (configFile.defined("tet mesh name"))
    tetMeshName = configFile.getString("tet mesh name", tetMeshName);
  else
    tetMeshName = outputPath + triangleMeshName + string(".tetmesh");

  bool skinned = false;
  skinned = configFile.getBool("skinned", false);
  if (skinned)
  {
    tetMeshName = tetMeshName + string(".constrained");
    cout << " Model is skinned! Using tet mesh filename " << tetMeshName.c_str() << endl;
  }
  else
    cout << " Model is not skinned! " << endl;

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
    string materialFilename("");
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

  cout << "==================================================" << endl;
  cout << " Reading mesh file " << tetMeshName.c_str() << endl;
  cout << "==================================================" << endl;
  tetMesh = new TET_MESH(tetMeshName.c_str(), materials, totalMaterials, false);

  // compute the modal derivatives
  computeModalDerivatives(tetMeshName);

  // remove the PCA flag
  string pcaCheck = tetMeshName + (".pca.peformed");
  string rm("rm ");
  rm = rm + pcaCheck;
  system(rm.c_str());

  delete tetMesh;
  delete materials;
}
