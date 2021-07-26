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

#include <arpackf.h>
#include <umfpack.h>

#if 0
extern "C" {
#include <taucs.h>
}
#endif

using namespace std;
SUBSPACE_TET_MESH* tetMesh;

#ifndef mod
#define mod(x,n) ((x) % (n))
#endif

//////////////////////////////////////////////////////////////////////
// Generate a basis for the mesh using linear modal analysis (LMA)
// with Arpack and Umfpack
//////////////////////////////////////////////////////////////////////
int generateLMABasisArpackUmfpack(int totalModes)
{
  int unconstrainedDiscard = 10;

  // if it's unconstrained, solve for 6 extra modes
  if (!tetMesh->constrained())
    totalModes += unconstrainedDiscard;

  int N = tetMesh->stiffnessMatrix().rows();

  SPARSE_MATRIX& B = tetMesh->massMatrix();
  SPARSE_MATRIX& K = tetMesh->stiffnessMatrix();

  // Getting a Bus Error on this on OSX, 12/2/10
  /*
  {
    cout << "Saving sparse mass and stiffness matrices" << endl;

    string meshFilename = tetMesh->filename();
    string massFileName = meshFilename + string(".mass.bcsm");
    string stiffnessFileName = meshFilename + string(".stiffness.bcsm");

    B.writeToBinary( massFileName );
    K.writeToBinary( stiffnessFileName );
  }
  */

  vector<int> rows;
  vector<int> cols;
  vector<Real> vals;

  K.entries(rows, cols, vals);
  int nz = vals.size();
  //int m = K.rows();
  int n = K.cols();

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
  
  // clearing stiffness matrix
  K.clear();

  void* symbolic;
  void* numeric;
  umfpack_di_symbolic(n, n, Ap, Ai, Ax, &symbolic, NULL, NULL);
  umfpack_di_numeric(Ap, Ai, Ax, symbolic, &numeric, NULL, NULL);
  umfpack_di_free_symbolic(&symbolic);

  int problemSize = N;
  int totalEigenvalues = totalModes;

  cout << " problem size: " << N << endl;
  cout << " total eigenvalues: " << totalModes  << endl;
  
  ARint ido = 0;
  char which[] = "LM";
  double tol = 1e-8;
  double* residual = new double[problemSize];
  ARint ncv = 4 * totalEigenvalues;
  if (ncv > problemSize) ncv = problemSize;
  ARint ldv = problemSize;
  double* v = new double[ldv*ncv];
  ARint* iparam = new ARint[11];
  for (int x = 0; x < 11; x++)
    iparam[x] = 0;
  iparam[0] = 1;
  iparam[2] = 3 * problemSize;
  iparam[3] = 1;
  iparam[4] = totalEigenvalues;

  char bmat[2] = "G";
  iparam[6] = 3;

  ARint* ipntr = new ARint[11];
  for (int x = 0; x < 11; x++)
    ipntr[x] = 0;

  double* workd = new double[3*problemSize];
  double* workl = new double[ncv*(ncv+8)];
  ARint lworkl = ncv*(ncv+8);
  ARint info = 1;
  ARint rvec = 1;  // Changed from above
  ARint* select = new ARint[ncv];
  double* d = new double[2*ncv];
  double sigma = 0;
  ARint ierr;

  for (int x = 0; x < 3 * problemSize; x++)
    workd[x] = 0;
  for (int x = 0; x < ncv*(ncv+8); x++)
    workl[x] = 0;
  for (int x = 0; x < 2*ncv; x++)
    d[x] = 0;
  for (int x = 0; x < problemSize; x++)
    residual[x] = x+1;
  for (int x = 0; x < ldv*ncv; x++)
    v[x] = 0;
  for (int x = 0; x < ncv; x++)
    select[x] = 0;

  VECTOR xVec(problemSize);
  VECTOR result(problemSize);
  do {
    F77NAME(dsaupd)(&ido, bmat, &problemSize, &which[0], 
      &totalEigenvalues, &tol, residual, 
	    &ncv, v, &ldv, iparam, ipntr, workd, workl,
	    &lworkl, &info);

    double* work1 = workd + ipntr[0] - 1;
    double* work2 = workd + ipntr[1] - 1;
    double* work3 = workd + ipntr[2] - 1;

    if (ido == -1)
    {
      for (int x = 0; x < problemSize; x++)
        xVec[x] = work1[x];

      result = B * xVec;
      umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, xVec.data(), result.data(), numeric, NULL, NULL);
      for (int x = 0; x < problemSize; x++)
        work2[x] = xVec[x];
    }

    if (ido == 1)
    { 
      for (int x = 0; x < problemSize; x++)
        xVec[x] = work3[x];

      umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, result.data(), xVec.data(), numeric, NULL, NULL);
      for (int x = 0; x < problemSize; x++)
        work2[x] = result[x];
    }

    if (ido == 2)
    { 
      for (int x = 0; x < problemSize; x++)
        xVec[x] = work1[x];

      result = B * xVec;
      for (int x = 0; x < problemSize; x++)
        work2[x] = result[x];
    }
  } while (ido != 99);

  if (info < 0) {
         cout << "Error with dsaupd, info = " << info << "\n";
         cout << "Check documentation in dsaupd\n\n";
  } else {
    char All[] = "All";
    F77NAME(dseupd)(&rvec, &All[0], select, d, v, &ldv, &sigma, bmat,
	    &problemSize, which, &totalEigenvalues, &tol, residual, &ncv, v, &ldv,
	    iparam, ipntr, workd, workl, &lworkl, &ierr);

    if (ierr != 0) {
      cout << "Error with dseupd, info = " << ierr << "\n";
      cout << "Check the documentation of dseupd.\n\n";
    } else if (info == 1) {
      cout << "Maximum number of iterations reached.\n\n";
    } else if (info == 3) {
      cout << "No shifts could be applied during implicit\n";
      cout << "Arnoldi update, try increasing NCV.\n\n";
    }

    // get all the eigenvalues, see if any need to be discarded
    VECTOR allEigs(totalEigenvalues);
    for (int i = 0; i < totalEigenvalues; i++)
      allEigs[i] = d[i];
    int start = 0;
    for (int i = 0; i < totalEigenvalues; i++)
      if (fabs(allEigs[i]) < 1e-1)
        start = i + 1;
      else
        i = totalEigenvalues + 1;

    // retrieve and store the eigensystem
    MATRIX& UBasis = tetMesh->U();
    VECTOR& eigenvalues = tetMesh->eigenvalues();

    // decide how to resize based on if the mesh is constrained
    if (tetMesh->constrained())
    {
      eigenvalues.resizeAndWipe(totalModes);
      UBasis.resizeAndWipe(N, totalModes);
    }
    else
    {
      eigenvalues.resizeAndWipe(totalModes - start);
      UBasis.resizeAndWipe(N, totalModes - start);
    }

    // bail if we didn't pad the unconstraineds enough
    if (start > unconstrainedDiscard)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " The mesh has more than " << unconstrainedDiscard << " unconstrained modes!!! Boost the number of bogus modes to solve for and rerun!" << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      exit(0);
    }

    for (int i = 0; i < totalEigenvalues - start; i++)
      eigenvalues[i] = d[i + start];

    for (int i = 0; i < totalEigenvalues - start; i++) 
      for (int j = 0; j < problemSize; j++) 
        UBasis(j,i) = v[(i + start)*problemSize+j];

    // make each column magnitude 1 
    //UBasis.orthogonalize();

    string meshFilename = tetMesh->filename();
    string valuename = meshFilename + string(".eigenvalues.vector");
    string vectorname = meshFilename + string(".eigenvectors.matrix");
    cout << " Writing eigenvalues to file " << valuename.c_str() << endl;
    cout << " Writing eigenvectors to file " << vectorname.c_str() << endl;
    eigenvalues.write(valuename.c_str());
    UBasis.write(vectorname.c_str());
 
    cout << " Found eigenvalues: " << endl;
    if (eigenvalues.size() < 10)
      cout << eigenvalues << endl;
    else
    {
      cout << " (Only printing first 10 ...)" << endl;
      cout << "["<< endl;
      for (int x = 0; x < 10; x++)
        cout << eigenvalues[x] << endl;
      cout << "]" << endl;
    }

    if (!tetMesh->constrained())
    {
      cout << " Discarded these rigid modes (these should all be near zero): " << endl;
      for (int i = 0; i < start; i++)
        cout << d[i] << endl;
    }
  }

  umfpack_di_free_numeric(&numeric);
  delete[] Ap;
  delete[] Ai;
  delete[] Ax;

  delete residual;
  delete v;
  delete iparam;
  delete ipntr;
  delete workd;
  delete workl;
  delete select;
  delete d;

  // see how many valid eigenvalues were produced
  Real threshold = 1e-8;
  VECTOR& eigenvalues = tetMesh->eigenvalues();
  int totalValid = 0;
  for (int x = 0; x < eigenvalues.size(); x++)
    if (eigenvalues[x] > threshold)
      totalValid++;

  return totalValid;
}

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

  bool skinned = false;
  skinned = configFile.getBool("skinned", false);
  if (skinned)
  {
    tetMeshName = tetMeshName + string(".constrained");
    cout << " Model is skinned! Using tet mesh filename " << tetMeshName.c_str() << endl;
  }
  else
  {
    cout << " Model is not skinned! " << endl;
  }

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

  bool useVolumeScale      = configFile.getBool("volume scale", false);
  bool useDensity          = configFile.getBool("use density", false);

  cout << "==================================================" << endl;
  cout << " Stomping any previous eigenanalysis" << endl;
  cout << "==================================================" << endl;
  string valuename = string("rm ") + tetMeshName + string(".eigenvalues.vector");
  string vectorname = string("rm ") + tetMeshName + string(".eigenvectors.matrix");
  string surfacename = string("rm ") + tetMeshName + string(".surfaceU.matrix");
  string trainingset = string("rm ") + tetMeshName + string(".trainingset");
  system(valuename.c_str());
  system(vectorname.c_str());
  system(surfacename.c_str());
  system(trainingset.c_str());

  Real meshScale           = configFile.getFloat("mesh scale", 1.0);
  
  cout << "Got mesh scale " << meshScale << endl;

  cout << "==================================================" << endl;
  cout << " Reading mesh file " << tetMeshName.c_str() << endl;
  cout << "==================================================" << endl;
  tetMesh = new SUBSPACE_TET_MESH(tetMeshName.c_str(), materials, totalMaterials,
                                  false, NULL, NULL, meshScale);
                                  // Jeff: these aren't in the constructor anymore?
                                  //false, NULL, NULL, meshScale,
                                  //useVolumeScale, useDensity);
  /*
  CUBATURE_GENERATOR generator(tetMesh);
  generator.maxKeyTets() = maxKeyTets;
  generator.candidatesPerTry() = candidatesPerTry;
  generator.errorTolerance() = errorTolerance;
  generator.randSeed(randSeed);
  */

  // stomp any previous basis
  if (tetMesh->basisLoaded())
  {
    tetMesh->U().resizeAndWipe(0,0);
    tetMesh->eigenvalues().resizeAndWipe(0);
  }

#if 0
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " MODIFYING MASS FOR THE TREE" << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  
  SPARSE_MATRIX& masses = tetMesh->massMatrix();
  int unconstrained = tetMesh->unconstrainedNodes();
  vector<VEC3F>& vertices = tetMesh->restPose();

  // modify the masses in the original meshes
  cout << " Computing new masses " << endl;
  for (unsigned int x = 0; x < unconstrained; x++)
  {
    int index = 3 * x;
    VEC3F& vertex = vertices[x];

    // is it inside the trunk?
    if (vertex[1] < 0.347) continue;

    Real xComponent = vertex[0] - 0.529623;
    Real zComponent = vertex[2] - 0.525765;

    Real radius = sqrt(xComponent * xComponent + zComponent * zComponent);
    if (radius < 0.03) continue;

    float massFactor = 100.0;

    masses(index, index) *= massFactor;
    masses(index + 1,index + 1) *= massFactor; 
    masses(index + 2, index + 2) *= massFactor;
  }
#endif

  // if we're doing this for a skinned mesh, it's for basis enrichment, so do a bunch of them
  if (skinned)
    rank *= 4;

  cout << "==================================================" << endl;
  cout << " Generating basis of rank = " << rank << endl;
  cout << "==================================================" << endl;
  cout << " Generating stiffness matrix ... ";
  flush(cout);
  SPARSE_MATRIX& stiffness = tetMesh->generateSparseStiffnessMatrix(true);
  cout << "done " << endl;

  // create the output directory
  string mkdir("mkdir ");
  mkdir = mkdir + outputPath;
  system(mkdir.c_str());

  //generateLMABasisArpackTaucsLU(rank);
  //generateLMABasisArpackTaucsCholesky(rank);
  //generateLMABasisArpackTaucsLDLT(rank);
  generateLMABasisArpackUmfpack(rank);

  /*  
  // solve for twice as many as we need to
  int totalValid = generateLMABasisArpackUmfpack(2 * rank);

  // only keep the valid ones
  int duds = 2 * rank - totalValid;

  cout << " Found " << duds << " dud eigenvalues. " << endl;

  // delete the dud eigenvalues
  MATRIX& U = tetMesh->U();
  VECTOR& eigs = tetMesh->eigenvalues();

  MATRIX newU(U.rows(), rank);
  VECTOR newEigs(rank);

  // build the dud-less matrices
  for (int x = 0; x < rank; x++)
  {
    newEigs[x] = eigs[x + duds];

    VECTOR column = U.getColumn(x + duds);
    newU.setColumn(column, x);
  }
  //cout << " New (culled) eigenvalues: " << newEigs << endl;

  valuename = tetMeshName + string(".eigenvalues.vector");
  vectorname = tetMeshName + string(".eigenvectors.matrix");
  cout << " Writing eigenvalues to file " << valuename.c_str() << endl;
  cout << " Writing eigenvectors to file " << vectorname.c_str() << endl;
  newEigs.write(valuename.c_str());
  newU.write(vectorname.c_str());
  */

  cout << " Skipping cubature training. Call CubatureGenerator." << endl;
  cout << " Make sure to delete any existing training sets. " << endl;

  delete tetMesh;
  delete materials;
}
