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
// PARTITIONED_CUBATURE_GENERATOR.h: interface for the PARTITIONED_CUBATURE_GENERATOR class.
//
//////////////////////////////////////////////////////////////////////

#include "PARTITIONED_CUBATURE_GENERATOR.h"
#include <float.h>
#ifdef USING_OPENMP
#include <omp.h>
#endif

//////////////////////////////////////////////////////////////////////
// Constructor for tetrahedra
//////////////////////////////////////////////////////////////////////
PARTITIONED_CUBATURE_GENERATOR::PARTITIONED_CUBATURE_GENERATOR(PARTITIONED_SUBSPACE_TET_MESH* tetMesh,
    const char* filename) :
  _partitionedMesh(tetMesh),
  // init Mersenne Twister with deterministic seed
  // so that runs are reproducible
  _trainingTwister(123456),
  _tolerance(1e-7),
  _maxKeyTets(1000),
  _candidatesPerTry(100),
  _errorTolerance(0.01f),
  _maxIterationsNNLS(3)
{
  // allocate training sample arrays for each partition
  int partitions = _partitionedMesh->partitions();
  _partitions = partitions;
  _trainingForcesPartitioned = new vector<VECTOR>[partitions];
  _trainingQsPartitioned = new vector<VECTOR>[partitions];
  _trainingColumnsPartitioned = new vector<VECTOR>[partitions];
  _trainingConstraintsPartitioned = new vector<VECTOR>[partitions];

  for (int x = 0; x < _partitions; x++)
    _forceMagnitudesPartitioned.push_back(VECTOR());

  // try reading in a precomputed training set
  if (filename != NULL)
  {
    _filename = string(filename);
    read(filename);
    cout << " Precomputed training set found: " << _filename << endl;
  }
  else 
  {
    // try probing for the default filename
    _filename = tetMesh->filename() + string(".trainingset");
    FILE* file = fopen(_filename.c_str(), "rb");
    bool exists = (file != NULL);
    if (file != NULL) fclose(file);

    if (exists)
    {
      cout << " Precomputed training set found: " << _filename << endl;
      read(_filename.c_str());
    }
    else
    {
      cout << " No precomputed training set (even tried " << _filename.c_str() << ") " << endl;
    }
  }

  // allocate key tet and weight vectors for each partition
  _keyTets = new vector<TET*>[partitions];
  _keyWeights = new vector<Real>[partitions];

  // allocate a bunch of independant RNGs
  _twisters = new MERSENNETWISTER[partitions];
}

PARTITIONED_CUBATURE_GENERATOR::~PARTITIONED_CUBATURE_GENERATOR()
{
  delete[] _trainingForcesPartitioned;
  delete[] _trainingQsPartitioned;
  delete[] _trainingColumnsPartitioned;
  delete[] _trainingConstraintsPartitioned;
  delete[] _keyTets;
  delete[] _keyWeights;
  delete[] _twisters;
}

//////////////////////////////////////////////////////////////////////
// Reinitialize for a different tet mesh
//////////////////////////////////////////////////////////////////////
void PARTITIONED_CUBATURE_GENERATOR::reinitialize(PARTITIONED_SUBSPACE_TET_MESH* tetMesh)
{
  delete[] _trainingForcesPartitioned;
  delete[] _trainingQsPartitioned;
  if (_trainingColumnsPartitioned)
  {
    delete[] _trainingColumnsPartitioned;
    _trainingConstraintsPartitioned = NULL;
  }
  if (_trainingConstraintsPartitioned)
  {
    delete[] _trainingConstraintsPartitioned;
    _trainingConstraintsPartitioned = NULL;
  }
  delete[] _keyTets;
  delete[] _keyWeights;
  delete[] _twisters;

  _trainingTwister = 123456;
  _partitionedMesh = tetMesh;
  _maxKeyTets = 1000;
  _candidatesPerTry = 100;
  _errorTolerance = 0.01f;
  _maxIterationsNNLS = 3;
  _tolerance = 1e-7;

  // allocate training sample arrays for each partition
  int partitions = _partitionedMesh->partitions();
  _partitions = partitions;
  _trainingForcesPartitioned = new vector<VECTOR>[partitions];
  _trainingQsPartitioned = new vector<VECTOR>[partitions];
  _trainingColumnsPartitioned = new vector<VECTOR>[partitions];
  _trainingConstraintsPartitioned = new vector<VECTOR>[partitions];

  for (int x = 0; x < _partitions; x++)
    _forceMagnitudesPartitioned.push_back(VECTOR());

  // try probing for the default filename
  _filename = tetMesh->filename() + string(".trainingset");
  FILE* file = fopen(_filename.c_str(), "rb");
  bool exists = (file != NULL);
  if (file != NULL) fclose(file);

  if (exists)
  {
    cout << " Precomputed training set found: " << _filename << endl;
    read(_filename.c_str());
  }
  else
  {
    cout << " No precomputed training set (even tried " << _filename.c_str() << ") " << endl;
  }

  // allocate key tet and weight vectors for each partition
  _keyTets = new vector<TET*>[partitions];
  _keyWeights = new vector<Real>[partitions];

  // allocate a bunch of independant RNGs
  _twisters = new MERSENNETWISTER[partitions];
}

//////////////////////////////////////////////////////////////////////
// Generate a basis for the mesh using linear modal analysis (LMA)
//////////////////////////////////////////////////////////////////////
void PARTITIONED_CUBATURE_GENERATOR::generateLMABasis(int totalModes)
{
  SUBSPACE_TET_MESH* originalMesh = (SUBSPACE_TET_MESH*)(_partitionedMesh->originalMesh());
  
  // retrieve the raw data
  MATRIX& UBasis = originalMesh->U();
  VECTOR& eigenvalues = originalMesh->eigenvalues();
  string meshFilename = originalMesh->filename();
  
  // make sure the stiffness matrix has been allocated
  cout << " Building stiffness matrix ... ";
  flush(cout);
  originalMesh->generateSparseStiffnessMatrix();
  cout << "done." << endl;
  
  PetscInt N = originalMesh->stiffnessMatrix().rows();

  Mat         	 A, B;
  EPS         	 eps;
  PetscReal   	 error, tol;
  int         	 maxit, i, its;

  // resize the result arrays
  eigenvalues.resizeAndWipe(totalModes);
  UBasis.resizeAndWipe(N, totalModes);

  SlepcInitialize(NULL,NULL,NULL,NULL);

  // get preallocation hints for PetSc
  vector<int> nonZerosPerRow = originalMesh->stiffnessMatrix().nonZerosPerRow();
  PetscInt* nonZeros = new PetscInt[nonZerosPerRow.size()];
  int maxNonZeros = 0;
  for (unsigned int x = 0; x < nonZerosPerRow.size(); x++)
  {
    if (nonZerosPerRow[x] > maxNonZeros)
      maxNonZeros = nonZerosPerRow[x];
    nonZeros[x] = nonZerosPerRow[x];
  }

  // create A
  MatCreateSeqAIJ(PETSC_COMM_WORLD,N,N,maxNonZeros,nonZeros, &A);
  MatSetFromOptions(A);
 
  // create B
  MatCreateSeqAIJ(PETSC_COMM_WORLD,N,N,1, PETSC_NULL, &B);
  MatSetFromOptions(B);

  // delete preallocation hints for PetSc
  delete[] nonZeros;

  // copy entries into stiffness matrix
  // Do it a generic but memory-hogging way
  cout << " Setting up stiffness eigensystem ... ";
  vector<int> rows;
  vector<int> cols;
  vector<Real> values;
  originalMesh->stiffnessMatrix().entries(rows, cols, values);
  for (unsigned int x = 0; x < rows.size(); x++)
  {
    // format is:
    //
    // MatSetValues(matrix, 
    //  number of rows, row index, 
    //  number of cols, col index, 
    //  actual value, INSERT_VALUES);
    PetscScalar 	 v = values[x];
    MatSetValues(A,1,&(rows[x]),1,&(cols[x]),&v,INSERT_VALUES);
  }
  cout << " done. " << endl;
  flush(cout);

  // copy entries into mass matrix
  // Do it a generic but memory-hogging way
  cout << " Setting up mass eigensystem ... ";
  originalMesh->massMatrix().entries(rows, cols, values);
  for (unsigned int x = 0; x < rows.size(); x++)
  {
    PetscScalar 	 v = values[x];
    MatSetValues(B,1,&(rows[x]),1,&(cols[x]),&v,INSERT_VALUES);
  }
  cout << " done. " << endl;

  // finalize matrices
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);

  // Create eigensolver context
  EPSCreate(PETSC_COMM_WORLD,&eps);

  // Set operators. In this case, it is a generalized eigenvalue problem
  EPSSetOperators(eps,A,B);
  EPSSetProblemType(eps,EPS_GHEP);
  EPSSetWhichEigenpairs(eps, EPS_SMALLEST_MAGNITUDE);

  // total number of eigenvalues to search for --
  // pad by 6 for the rigid modes if the mesh is
  // unconstrained
  int totalEigenvalues = totalModes;
  if (!originalMesh->constrained())
  {
    cout << " Solving for 6 extra modes because the mesh is unconstrained." << endl;
    totalEigenvalues += 6;
  }

  // Slepc recommends using at least twice as many eigenvectors
  // for the search space
  //int totalEigenvectors = totalEigenvalues * 2 + 1;
  //EPSSetDimensions(eps, totalEigenvalues, totalEigenvectors);
  // SLEPc 3.0 interface change
  EPSSetDimensions(eps, totalEigenvalues, PETSC_DECIDE, PETSC_DECIDE);

  PetscTruth hermitian = PETSC_TRUE;
  EPSIsHermitian(eps, &hermitian);

  // set the tolerances
  int maxIterations = 10000;
  EPSSetTolerances(eps, _tolerance, maxIterations);

  // Send Slepc the solver parameters
  EPSSetFromOptions(eps);

  // Solve the eigensystem
  TIMER eigenTimer;
  cout << " Solving eigensystem ... ";
  EPSSolve(eps);
  EPSGetIterationNumber(eps, &its);
  EPSGetTolerances(eps,&tol,&maxit);
  cout << "done in " << its << " of " << maxit << " max iterations " << endl;
  cout << "time: " << eigenTimer.timing() << endl;

  // Display solution and clean up
  int totalFound = 0;
  EPSGetConverged(eps, &totalFound);

  if (totalFound < totalModes)
    cout << __FILE__ << " " << __LINE__ << " : Not enough eigenvalues were found! " << endl;

  Vec realVector;
  VecCreateSeq(PETSC_COMM_SELF, N, &realVector);

  // make sure no imaginary components were found
  bool imaginaryFound = false;

  if (totalFound > 0) {
    int currentMode = 0;
    for (i = 0; i < totalEigenvalues; i++) {
      PetscScalar realValue, imagValue;

      // Get converged eigenpairs: i-th eigenvalue is stored in 
      // kr (real part) and ki (imaginary part)
      EPSGetEigenpair(eps, i, &realValue, &imagValue, realVector, NULL);

      // Compute the relative error associated to each eigenpair
      EPSComputeRelativeError(eps,i,&error);

      Real re = realValue;
      Real im = imagValue;

      // make sure the eigenvalue is non-zero within the
      // solver tolerance before storing, and also make
      // sure we haven't already stored enough eigenvectors
      if (fabs(re) > _tolerance && currentMode < totalModes)
      {
        // store the eigenvalue
        eigenvalues(currentMode) = re;

        if (im!=0.0) 
          imaginaryFound = true;

        PetscScalar* a;
        VecGetArray(realVector, &a);

        // store the eigenvector
        for (int x = 0; x < N; x++)
          UBasis(x, currentMode) = a[x];
        VecRestoreArray(realVector, &a);

        currentMode++;
      }
    }

    // check if we found enough non-near-zero modes
    if (currentMode < totalModes)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " 
           << " Not enough usable eigenvectors were found!" << endl
           << " Try using a finer mesh or reducing the rank." << endl;
    }
  }
  VecDestroy(realVector);
  
  // Free work space
  EPSDestroy(eps);
  MatDestroy(A);

  if (imaginaryFound)
    cout << __FILE__ << " " << __LINE__ << " : Eigenvalue solver found an imaginary component! There must be a problem with the stiffness matrix. " << endl;

  cout << " Eigenvalues : " << endl << eigenvalues;

  string valuename = meshFilename + string(".eigenvalues.vector");
  string vectorname = meshFilename + string(".eigenvectors.matrix");
  cout << " Writing eigenvalues to file " << valuename.c_str() << endl;
  cout << " Writing eigenvectors to file " << vectorname.c_str() << endl;
  eigenvalues.write(valuename.c_str());
  UBasis.write(vectorname.c_str());
}

//////////////////////////////////////////////////////////////////////
// Generate cubature training samples
//////////////////////////////////////////////////////////////////////
void PARTITIONED_CUBATURE_GENERATOR::generateTrainingSamples(int partition, int totalSamples, Real magnitude)
{
  SUBSPACE_TET_MESH* mesh = (SUBSPACE_TET_MESH*)(_partitionedMesh->mesh(partition));
  VECTOR& eigenvalues = mesh->eigenvalues();
  VECTOR q(eigenvalues.size());

  // the user wants the first mode to have the given magnitude
  _magnitude = magnitude * sqrt(eigenvalues(0));

  int attempts = 0;
  int successes = 0;

  // clear old samples just in case
  _trainingForcesPartitioned[partition].clear();
  _trainingQsPartitioned[partition].clear();

  // track if a Nan has been seen before so that the error message
  // doesn't get output too many times
  bool firstNan = true;

  cout << " Generating training set for partition " << partition << " ...";
  for(int i = 0; i < totalSamples; i++)
  {
    int uninvertedTries = 0;
    bool tetInverted = true;

    // try really hard to generate an uninverted sample
    while (tetInverted && uninvertedTries < 100)
    {
      // create random sample
      for (int x = 0; x < q.size(); x++)
      {
        // scale by inverse freq squared (ie. divide by eig val)
        // 4.0 because.. gauss is usually between -4 and 4
        q(x) = _magnitude / sqrt(eigenvalues(x)) / 4.0 * _trainingTwister.randNorm();

        // if the eigenvalue is small enough to produce a nan
        // (this has happened before) try to do something graceful
        if (isnan(q(x))) 
        {
          q(x) = 0.0;
          if (firstNan)
          {
            firstNan = false;
            cout << __FILE__ << " " << __LINE__ << " : Very small eigenvalue detected!" << endl
                 << " You should reduce your rank if it is very large "<< endl
                 << " or increase the resolution of your mesh if it is very coarse." << endl;
          }
        }
      }

      mesh->q() = q;
      mesh->updateFullMesh();

      // check for inversion
      tetInverted = mesh->inverted();

      // DEBUG
      tetInverted = false;

      uninvertedTries++;
      attempts++;
    }

    // throw a flag if the tet is still inverted, but forge ahead
    if (tetInverted)
    {
      cout << __FILE__ << " " << __LINE__ << " : Inverted Tet Mesh! "
           << " You should try again with a smaller magnitude" << endl;
    }
    else
      successes++;

    // store q and the force it produces
    VECTOR forceSample = mesh->projectedInternalForce();
    _trainingForcesPartitioned[partition].push_back(forceSample);
    _trainingQsPartitioned[partition].push_back(q);

    // DEBUG
    if (isnan(forceSample[0]))
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " force sample: " << forceSample << endl;
      cout << " q: " << q << endl;

      bool found = false;
      MATRIX& U = mesh->U();
      for (int x = 0; x < U.rows(); x++)
        for (int y = 0; y < U.cols(); y++)
          if (isnan(U(x,y)))
          {
            cout << " Nan in U at " << x << ", " << y << ": " << U(x,y) << endl;
            found = true;
          }
      if (found)
        cout << " Nan found in U" << endl;
      else
        cout << " No Nan found in U " << endl;

      found = false;
      vector<VEC3F>& forces = mesh->forceVectors();
      for (unsigned int x = 0; x < forces.size(); x++)
        if (isnan(forces[x][0]) ||
            isnan(forces[x][1]) ||
            isnan(forces[x][2]))
          found = true;
      if (found)
        cout << " Nan found in R" << endl;
      else
        cout << " No Nan found in R" << endl;
      
      exit(0);
    }

    // output status every 10%
    if (i % (int)(totalSamples / 10) == 0)
    {
      cout << 100 * ((Real)i / totalSamples) << "% ";
      flush(cout);
    }
  }
  cout << " done. " << endl;
  cout << " " << successes << " of " << attempts << " attempts successful (" 
        << 100.0f * (Real)successes / attempts << "%)" << endl;
}

//////////////////////////////////////////////////////////////////////
// Generate cubature training samples
//////////////////////////////////////////////////////////////////////
void PARTITIONED_CUBATURE_GENERATOR::generateTrainingSamples(int totalSamples, Real magnitude)
{
  for (int x = 0; x < _partitions; x++)
    generateTrainingSamples(x, totalSamples, magnitude);

  char buffer[256];
  sprintf(buffer, "%i", _partitions);

  // cache results to the default filename
  string trainingFilename = _partitionedMesh->filename();
  trainingFilename = trainingFilename + string(".partitions.") + string(buffer) + string(".trainingset");
  write(trainingFilename.c_str());
}

//////////////////////////////////////////////////////////////////////
// write out a training set
//////////////////////////////////////////////////////////////////////
void PARTITIONED_CUBATURE_GENERATOR::write(const char* filename)
{
  FILE* file;
  file = fopen(filename, "wb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __LINE__ << " Couldn't write file: " << filename << endl;
    return;
  }
  cout << " Writing training set: " << filename << endl;

  // write dimensions
  for (int y = 0; y < _partitions; y++)
  {
    int size = _trainingForcesPartitioned[y].size();
    fwrite((void*)&size, sizeof(int), 1, file);

    // write out the forces
    for (int x = 0; x < size; x++)
      _trainingForcesPartitioned[y][x].write(file);

    // write out the qs
    for (int x = 0; x < size; x++)
      _trainingQsPartitioned[y][x].write(file);
  }

  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// read in a training set
//////////////////////////////////////////////////////////////////////
void PARTITIONED_CUBATURE_GENERATOR::read(const char* filename)
{
  FILE* file;
  file = fopen(filename, "rb");
  cout << " Reading training set: " << filename << endl;

  if (file == NULL)
  {
    cout << " Could not open file " << filename << " !!! " << endl;
    exit(0);
  }

  // read dimensions
  for (int y = 0; y < _partitions; y++)
  {
    int size;
    fread((void*)&size, sizeof(int), 1, file);
    cout << "   Reading " << size << " samples in partition " << y << endl; 
    _trainingForcesPartitioned[y].clear();
    _trainingQsPartitioned[y].clear();

    // read in the forces
    for (int x = 0; x < size; x++)
    {
      VECTOR sample(file);
      _trainingForcesPartitioned[y].push_back(sample);
    }

    // read in the qs
    for (int x = 0; x < size; x++)
    {
      VECTOR sample(file);
      _trainingQsPartitioned[y].push_back(sample);
    }
  }

  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// pick candidates for key tet selection
//////////////////////////////////////////////////////////////////////
vector<TET*> PARTITIONED_CUBATURE_GENERATOR::pickCandidates(int partition, int candidatesPerTry)
{
  TET_MESH* tetMesh = _partitionedMesh->mesh(partition);
  vector<TET*>& keyTets = _keyTets[partition];

  vector<TET>& tets = tetMesh->tets();
  int totalTets = tetMesh->totalTets();

  // build a lookup table for the already selected key tets
  map<TET*, bool> alreadyUsed;
  for (unsigned int x = 0; x < keyTets.size(); x++)
    alreadyUsed[keyTets[x]] = true;

  // construct a list of all non-key tets
  list<TET*> nonKeyTets;
  for (int x = 0; x < totalTets; x++)
  {
    // check that it is not already a key tet
    if (alreadyUsed.find(&(tets[x])) == alreadyUsed.end())
      nonKeyTets.push_back(&(tets[x]));
  }

  // cache list size since it is a linear time op
  // for an STL list
  int listSize = nonKeyTets.size();

  // cache which twister we want to use
  MERSENNETWISTER& twister = _twisters[partition];
  
  // construct the final vector of candidate tets
  vector<TET*> candidates;
  for (int x = 0; x < candidatesPerTry; x++)
  {
    // pick a random tet
    int index = twister.randInt(listSize - 1);

    // get an iterator to that tet
    list<TET*>::iterator i = nonKeyTets.begin();
    advance(i, index);

    // record it as a candidate tet
    candidates.push_back(*i);

    // remove it from the list
    nonKeyTets.erase(i);

    // record the new size of the list
    listSize--;
  }

  return candidates;
}

//////////////////////////////////////////////////////////////////////
// Generate the training column corresponding to a tet
//////////////////////////////////////////////////////////////////////
VECTOR PARTITIONED_CUBATURE_GENERATOR::getTrainingColumn(int partition, TET* tet)
{
  SUBSPACE_TET_MESH* tetMesh = (SUBSPACE_TET_MESH*)(_partitionedMesh->mesh(partition));
  vector<VECTOR>& trainingForces = _trainingForcesPartitioned[partition];
  VECTOR& forceMagnitudes = _forceMagnitudesPartitioned[partition];

  int sampleSize = trainingForces[0].size();
  int totalRows = trainingForces.size() * sampleSize;
  VECTOR gColumn(totalRows);

  // submatrix of _Ubasis corresponding to this tet
  MATRIX tetU = tetMesh->tetSubBasis(tet);

  // generate column corresponding to stacked force samples
  for (unsigned int x = 0; x < trainingForces.size(); x++)
  {
    // generate the displacement from the forcing, but only for
    // the vertices associated with this tet
    VECTOR displacement = tetU * _trainingQsPartitioned[partition][x];
 
    // copy the rest pose of the tet
    VEC3F vertices[4];
    for (int y = 0; y < 4; y++)
      vertices[y] = *(tet->vertices[y]);

    // create a tet from the vertices -- this must be done before adding
    // the displacement so that the material inverse is calculated
    // correctly
    TET newTet(&vertices[0], &vertices[1], &vertices[2], &vertices[3]);

    // add the displacement
    for (int y = 0; y < 4; y++)
    {
      vertices[y][0] += displacement(3 * y);
      vertices[y][1] += displacement(3 * y + 1);
      vertices[y][2] += displacement(3 * y + 2);
    }

    // evaluate force density with respect to u
    VECTOR forceDensity(9);
    MATERIAL* material = tetMesh->materials()[tet->materialIndex()];
    MATRIX3 F3 = newTet.F();
    VECTOR F = TET::flattenF(F3);
    material->forceDensity(F.data(), forceDensity.data());
   
    // convert force density to du
    MATRIX pFpu(9,12);
    material->computePFPu(newTet, pFpu);
    VECTOR densityVector = forceDensity * pFpu;

    // get reduced version of force density
    VECTOR forceReduced = tetU ^ densityVector;

    // copy this reduced force vector into the correct
    // place of the overall training column, making sure to 
    // normalize by the same factor as in generateCubature()
    for (int y = 0; y < sampleSize; y++)
      gColumn(x * sampleSize + y) = forceReduced(y) * forceMagnitudes(x);
  }

  return gColumn;
}

//////////////////////////////////////////////////////////////////////
// Generate the training column corresponding to a tet
//////////////////////////////////////////////////////////////////////
VECTOR PARTITIONED_CUBATURE_GENERATOR::getTrainingColumnPCA(int partition, TET* tet)
{
  vector<VECTOR>& trainingColumns = _trainingColumnsPartitioned[partition];
  vector<VECTOR>& trainingForces = _trainingForcesPartitioned[partition];
  vector<VECTOR>& trainingQs = _trainingQsPartitioned[partition];
  vector<VECTOR>& trainingConstraints = _trainingConstraintsPartitioned[partition];
  SUBSPACE_TET_MESH* tetMesh = (SUBSPACE_TET_MESH*)(_partitionedMesh->mesh(partition));
  VECTOR& forceMagnitudes = _forceMagnitudesPartitioned[partition];

  // DEBUG: turning this off to see why we can't generate it on demand
  // if the training forces were explicitly stored, just return those
  if (trainingColumns.size() > 0)
    return retrieveTrainingColumn(partition, tet);

  int sampleSize = trainingForces[0].size();
  int totalRows = trainingForces.size() * sampleSize;
  VECTOR gColumn(totalRows);

  // submatrix of _Ubasis corresponding to this tet
  MATRIX tetU = tetMesh->tetSubBasis(tet);

  // generate column corresponding to stacked force samples
  for (unsigned int x = 0; x < trainingForces.size(); x++)
  {
    // generate the displacement from the forcing, but only for
    // the vertices associated with this tet
    VECTOR displacement = tetU * trainingQs[x];

    // copy the rest pose of the tet
    VEC3F vertices[4];
    for (int y = 0; y < 4; y++)
    {
      int vertexID = tetMesh->vertexID(tet->vertices[y]);
      vertices[y] = *(tetMesh->restVertices(vertexID));
    }

    // create a tet from the vertices -- this must be done before adding
    // the displacement so that the material inverse is calculated
    // correctly
    TET newTet(&vertices[0], &vertices[1], &vertices[2], &vertices[3]);

    // add the displacement
    for (int y = 0; y < 4; y++)
    {
      vertices[y][0] += displacement(3 * y);
      vertices[y][1] += displacement(3 * y + 1);
      vertices[y][2] += displacement(3 * y + 2);
    }
   
    // check if any of the vertices are constrained, and add those
    // displacements as well
    for (int y = 0; y < 4; y++)
      if (tetMesh->isConstrained(tet->vertices[y]))
      {
        int vertexID = tetMesh->vertexID(tet->vertices[y]);
        vertexID -= tetMesh->unconstrainedNodes();

        vertices[y][0] += trainingConstraints[x][3 * vertexID];
        vertices[y][1] += trainingConstraints[x][3 * vertexID + 1];
        vertices[y][2] += trainingConstraints[x][3 * vertexID + 2];
      }

    // evaluate force density with respect to u
    VECTOR forceDensity(9);
    MATERIAL* material = tetMesh->materials()[tet->materialIndex()];
    MATRIX3 F3 = newTet.F();
    VECTOR F = TET::flattenF(F3);
    material->forceDensity(F.data(), forceDensity.data());
  
    // convert force density to du
    MATRIX pFpu(9,12);
    material->computePFPu(newTet, pFpu);
    VECTOR densityVector = forceDensity * pFpu;

    // get reduced version of force density
    VECTOR forceReduced = tetU ^ densityVector;

    // copy this reduced force vector into the correct
    // place of the overall training column, making sure to 
    // normalize by the same factor as in generateCubature()
    for (int y = 0; y < sampleSize; y++)
      gColumn(x * sampleSize + y) = forceReduced(y) * forceMagnitudes(x);
  }

  return gColumn;
}

//////////////////////////////////////////////////////////////////////
// pick candidates for key tet selection
//////////////////////////////////////////////////////////////////////
vector<TET*> PARTITIONED_CUBATURE_GENERATOR::pickCandidatesPCA(int partition, int candidatesPerTry)
{
  TET_MESH* tetMesh = _partitionedMesh->mesh(partition);
  vector<TET>& tets = tetMesh->tets();
  vector<TET*>& keyTets = _keyTets[partition];

  // cache which twister we want to use
  MERSENNETWISTER& twister = _twisters[partition];

  int totalTets = tetMesh->totalTets();

  // build a lookup table for the already selected key tets
  map<TET*, bool> alreadyUsed;
  for (unsigned int x = 0; x < keyTets.size(); x++)
    alreadyUsed[keyTets[x]] = true;

  // construct the final vector of candidate tets
  vector<TET*> candidates;
  while (candidates.size() != (unsigned int)candidatesPerTry)
  {
    // pick a random tet
    int index = twister.randInt(totalTets - 1);

    // if it hasn't already been used, add it
    if (alreadyUsed.find(&(tets[index])) == alreadyUsed.end())
    {
      candidates.push_back(&tets[index]);
      alreadyUsed[&tets[index]] = true;
    }
  }

  return candidates;
}

//////////////////////////////////////////////////////////////////////
// Find the best fitting tet from the list of candidates
//////////////////////////////////////////////////////////////////////
int PARTITIONED_CUBATURE_GENERATOR::findBestTet(int partition, vector<TET*>& candidates, VECTOR& residual)
{
  Real rNorm = residual.norm2();

  // candidate dot products;
  vector<Real> dots;
  dots.resize(candidates.size());

  // calculate all the dot products
  for (unsigned int x = 0; x < candidates.size(); x++)
  {
    // get the 'g' column corresponding to this tet
    VECTOR gColumn = getTrainingColumn(partition, candidates[x]);
    Real gNorm = gColumn.norm2();

    // I wonder if this check is necessary, not very robust...
    if (fabs(gNorm) < 1e-7)
      dots[x] = -1e9;
    else
      dots[x] = gColumn * residual / (gNorm * rNorm);
  }

  int maxDot = 0;
  for (unsigned int x = 0; x < candidates.size(); x++)
    maxDot = (dots[x] > dots[maxDot]) ? x : maxDot;

  return maxDot;
}

//////////////////////////////////////////////////////////////////////
// Find the best fitting tet from the list of candidates
//////////////////////////////////////////////////////////////////////
int PARTITIONED_CUBATURE_GENERATOR::findBestTetPCA(int partition, vector<TET*>& candidates, VECTOR& residual)
{
  Real rNorm = residual.norm2();

  // candidate dot products;
  vector<Real> dots;
  dots.resize(candidates.size());

  // calculate all the dot products
  for (unsigned int x = 0; x < candidates.size(); x++)
  {
    // get the 'g' column corresponding to this tet
    VECTOR gColumn = getTrainingColumnPCA(partition, candidates[x]);
    Real gNorm = gColumn.norm2();

    // I wonder if this check is necessary, not very robust...
    if (fabs(gNorm) < 1e-7)
      dots[x] = -1e9;
    else
      dots[x] = gColumn * residual / (gNorm * rNorm);
  }

  int maxDot = 0;
  Real maxFound = 0;
  for (unsigned int x = 0; x < candidates.size(); x++)
  {
    if (dots[x] > dots[maxDot])
      maxFound = dots[x];
    
    maxDot = (dots[x] > dots[maxDot]) ? x : maxDot;
  }
  return maxDot;
}

//////////////////////////////////////////////////////////////////////
// generate cubature for all the partitions
//////////////////////////////////////////////////////////////////////
void PARTITIONED_CUBATURE_GENERATOR::generateCubature()
{
  // generate the cubature (finally)
  for (int x = 0; x < _partitions; x++)
  {
    cout << " Generating cubature for partition " << x << endl;
    generateCubature(x);
  }

  // write out the results
  cout << " Writing out generated cubatures " << endl;
  writePartitionedCubatures();
}

//////////////////////////////////////////////////////////////////////
// generate cubature for all the partitions using PCA data
//////////////////////////////////////////////////////////////////////
void PARTITIONED_CUBATURE_GENERATOR::generatePCACubature()
{
  // generate the cubature (finally)
  for (int x = 0; x < _partitions; x++)
  {
    cout << "=============================================" << endl;
    cout << " Generating PCA cubature for partition " << x << endl;
    cout << "=============================================" << endl;
    generatePCACubature(x);
  }

  // write out the results
  cout << " Writing out generated cubatures " << endl;
  writePartitionedCubatures();
}

//////////////////////////////////////////////////////////////////////
// Generate the cubature points
//////////////////////////////////////////////////////////////////////
void PARTITIONED_CUBATURE_GENERATOR::generateCubature(int partition)
{
  // retrieve the specific mesh and forces for this partition
  SUBSPACE_TET_MESH* tetMesh = (SUBSPACE_TET_MESH*)(_partitionedMesh->mesh(partition));
  vector<VECTOR>& trainingForces = _trainingForcesPartitioned[partition];
  VECTOR& forceMagnitudes = _forceMagnitudesPartitioned[partition];

  // see if a reasonable number of key tets were requested
  if (_maxKeyTets > tetMesh->totalTets())
    cout << __FILE__ << " " << __LINE__ << " : " 
         << "You've asked for more key tets than there are in the mesh! " 
         << "This could potentially loop forever." << endl;

  // reset the tet mesh to the rest pose since we will be
  // using its rest pose in getTrainingColumn()
  tetMesh->resetToRestPose();

  // size of a single force vector
  int forceSize = trainingForces[0].size();
  
  // total number of force samples
  int totalSamples = trainingForces.size();

  // number of rows in a column where all the samples are flattened
  // into a big vector
  int totalRows = forceSize * totalSamples;

  // calculate the squared magnitudes of the samples
  forceMagnitudes.resizeAndWipe(totalSamples);
  for (int x = 0; x < totalSamples; x++)
    for (int y = 0; y < forceSize; y++)
    {
      // bracket paren notation []() is a little confusing, I'll try to 
      // keep it to a minimum
      forceMagnitudes(x) += trainingForces[x](y) * trainingForces[x](y);
    }

  // calculate the inverse sqrt of the squared magnitudes
  for (int x = 0; x < totalSamples; x++)
    // CAREFUL WITH SINGLE PRECISION HERE
    if (fabs(forceMagnitudes(x)) > 1e-7)
       forceMagnitudes(x) = 1.0f / sqrtf(forceMagnitudes(x));
    else
       forceMagnitudes(x) = 1.0f;

  // stack normalized force samples into one big 'b' vector
  VECTOR b(totalRows);
  int index = 0;
  for (int x = 0; x < totalSamples; x++)
    for (int y = 0; y < forceSize; y++, index++)
      b(index) = trainingForces[x](y) * forceMagnitudes(x);

  // init the residual to b, as is traditional
  VECTOR residual(b);

  // "A" matrix to send to the the NNLS solver
  double* A = new double[totalRows * _maxKeyTets];

  // the Non-Negative Least Squares (NNLS) solver
  NNLS_SOLVER nnls(totalRows, _maxKeyTets);
  nnls.maxIter() = _maxIterationsNNLS;

  // scratch space for the NNLS solver
  double* bNNLS = new double[totalRows];
  double* weightsNNLS = new double[_maxKeyTets];
  for (int x = 0; x < _maxKeyTets; x++)
    weightsNNLS[x] = 0.0;

  // error metrics
  Real bNorm = b.norm2();
  double rNorm = bNorm;
  double relativeError = 1.0f;

  // retrieve the specific key tet and weight vectors for this partitions
  vector<TET*>& keyTets = _keyTets[partition];
  vector<Real>& keyWeights = _keyWeights[partition];
  
  // tracking key tets
  keyTets.clear();
  keyWeights.clear();
  int totalKeyTets = 0;

  // start training
  int totalIterations = 0;
  while (relativeError > _errorTolerance && 
         totalKeyTets < _maxKeyTets &&
         totalKeyTets < tetMesh->totalTets())
  {
    int candidatesPerTry = _candidatesPerTry;
    
    // get new set of randomly selected candidates
    if (tetMesh->totalTets() - totalKeyTets < _candidatesPerTry)
      candidatesPerTry = tetMesh->totalTets() - totalKeyTets;
    vector<TET*> candidates = pickCandidates(partition, candidatesPerTry);

    // pick the best tet out of the list of candidates
    int bestTet = findBestTet(partition, candidates, residual);

    // store the new key tet
    keyTets.push_back(candidates[bestTet]);
    totalKeyTets++;

    // recompute its column (better than storing them all)
    VECTOR bestTetColumn = getTrainingColumn(partition, candidates[bestTet]);

    // copy the column to "A" matrix --
    // note that A is column major, so it appears to be copying to a row
    for (int x = 0; x < totalRows; x++)
      A[(totalKeyTets - 1) * totalRows + x] = bestTetColumn(x);

    // copy b to double precision, also since the NNLS_SOLVER is going
    // to clobber it
    for (int x = 0; x < totalRows; x++)
      bNNLS[x] = b(x);
   
    // do the NNLS fit - note that both weightsNNLS and rNorm
    // both return values
    bool converged = nnls.solve(A, totalKeyTets, bNNLS, weightsNNLS, rNorm);

    // update residual of the fit converged
    if (converged)
    {
      relativeError = rNorm / bNorm;

      // calculate the residual by doing an explicit A*x, 
      // but since A is column major, it is really (A^T) * x
      residual = b;
      for (int i = 0; i < totalKeyTets; i++)
        for (int j = 0; j < totalRows; j++)
        {
          int index = i * totalRows + j;
          residual(j) -= weightsNNLS[i] * A[index];
        }
    }
    cout << "Partition " << partition << " iteration " << totalIterations << ": \tWeight: " << weightsNNLS[totalKeyTets-1] 
         << "\tError: " << relativeError << endl;
    totalIterations++;

    if (isnan(relativeError))
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " NaN error found!!!!" << endl;
      cout << " rNorm: " << rNorm << endl;
      cout << " bNorm: " << bNorm << endl;
      cout << " mesh rank: " << tetMesh->rank() << endl;
      cout << " total key tets: " << totalKeyTets << endl;
      cout << " force magnitudes: " << forceMagnitudes << endl;

      cout << " training samples: " << endl;
      for (int i = 0; i < totalSamples; i++)
        cout << trainingForces[i] << endl;
      exit(0);
    }
  }

  // copy the key tet weights to single precision
  keyWeights.clear();
  for (int x = 0; x < totalKeyTets; x++)
    keyWeights.push_back(weightsNNLS[x]);

  cout << " Total key tets before filtering: " << keyTets.size() << endl;
  filterKeyTets(partition);
  cout << " Total key tets after filtering: " << keyTets.size() << endl;
  
  // clean up scratch space
  delete[] bNNLS;
  delete[] weightsNNLS;
  delete[] A;
}

//////////////////////////////////////////////////////////////////////
// Filter out the dud tets
//////////////////////////////////////////////////////////////////////
void PARTITIONED_CUBATURE_GENERATOR::filterKeyTets(int partition)
{
  vector<TET*>& keyTets = _keyTets[partition];
  vector<Real>& keyWeights = _keyWeights[partition];
  
  vector<TET*> filteredTets;
  vector<Real> filteredWeights;

  // filter the tets
  for (unsigned int x = 0; x < keyTets.size(); x++)
    if (keyWeights[x] > 1e-12)
    {
      filteredTets.push_back(keyTets[x]);
      filteredWeights.push_back(keyWeights[x]);
    }

  // copy them all back
  keyTets.resize(filteredTets.size());
  keyWeights.resize(filteredTets.size());
  for (unsigned int x = 0; x < filteredTets.size(); x++)
  {
    keyTets[x] = filteredTets[x];
    keyWeights[x] = filteredWeights[x];
  }
}

//////////////////////////////////////////////////////////////////////
// write out the cubature for a single partition
//////////////////////////////////////////////////////////////////////
void PARTITIONED_CUBATURE_GENERATOR::writePartitionedCubature(int i, const char* filename)
{
  string finalfile("");
  SUBSPACE_TET_MESH* tetMesh = (SUBSPACE_TET_MESH*)(_partitionedMesh->mesh(i));
  vector<TET*>& keyTets = _keyTets[i];
  vector<Real>& keyWeights = _keyWeights[i];
  
  // if no filename is provided, generate one based on the mesh name
  if (filename == NULL)
  {
    finalfile.append(tetMesh->filename());
    finalfile.append(".cubature");
  }
  else
    finalfile.append(filename);

  // check that a cubature even exists to write out
  if (keyTets.size() == 0)
  {
    cout << __FILE__ << " " << __LINE__ << " : " 
         << " Cubature hasn't been generated yet!" 
         << " Did you remember to call generateCubature()?" 
         << endl;
    return;
  }
 
  // write out the cubature file
  FILE* file;
  file = fopen(finalfile.c_str(), "wb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Failed to open cubature file " << finalfile.c_str() << " for writing! Does the directory exist?" << endl;
  }

  // write dimensions
  int size = keyTets.size();
  fwrite((void*)&size, sizeof(int), 1, file);

  // write out the key tet indices
  for (int x = 0; x < size; x++)
  {
    int tetID = tetMesh->tetID(keyTets[x]);
    fwrite((void*)&tetID, sizeof(int), 1, file);
  }

  // write out the key tet weights
  for (int x = 0; x < size; x++)
  {
    double weight = keyWeights[x];
    fwrite((void*)&weight, sizeof(double), 1, file);
  }
  
  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// write out the cubature for all the different partitions
//////////////////////////////////////////////////////////////////////
void PARTITIONED_CUBATURE_GENERATOR::writePartitionedCubatures(const char* filename)
{
  for (int i = 0; i < _partitions; i++)
    writePartitionedCubature(i, filename);
}

//////////////////////////////////////////////////////////////////////////////
// Perform PCA on the data matrix, returning only the 1st "rank" columns
//////////////////////////////////////////////////////////////////////////////
void PARTITIONED_CUBATURE_GENERATOR::pcaLapack(MATRIX& data, int rank, MATRIX& components, VECTOR& values, bool shift)
{
  // subtract out the means
  int rows = data.rows();
  int cols = data.cols();
  VECTOR meanColumn(data.rows());
  if (shift)
  {
    for (int x = 0; x < rows; x++)
    {
      Real mean = 0.0;
      for (int y = 0; y < cols; y++)
        mean += data(x,y);
      mean /= cols;

      for (int y = 0; y < cols; y++)
        data(x,y) -= mean;

      meanColumn[x] = mean;
    }

    // divide through by sqrt(N - 1)
    data *= 1.0 / sqrt(cols - 1);
  }

  VECTOR::printVertical = false;

  MATRIX U;
  MATRIX VT;
  VECTOR S;

  // do a full SVD on the whole shebang
  data.SVD(U,S,VT);

  // only store the first "rank" vectors
  components.resizeAndWipe(rows, rank);
  values.resizeAndWipe(rank);

  for (int x = 0; x < rank; x++)
  {
    for (int y = 0; y < rows; y++)
      components(y,x) = U(y,x);
    values[x] = S[x] * S[x];
  }
}

//////////////////////////////////////////////////////////////////////////////
// Perform PCA on the data matrix, returning only the 1st "rank" columns
//////////////////////////////////////////////////////////////////////////////
void PARTITIONED_CUBATURE_GENERATOR::pcaSlepc(MATRIX& data, int rank, MATRIX& components, VECTOR& values, bool shift)
{
  // subtract out the means
  int rows = data.rows();
  int cols = data.cols();
  if (shift)
  {
    for (int x = 0; x < rows; x++)
    {
      Real mean = 0.0;
      for (int y = 0; y < cols; y++)
        mean += data(x,y);
      mean /= cols;

      for (int y = 0; y < cols; y++)
        data(x,y) -= mean;
    }

    // divide through by sqrt(N - 1)
    data *= 1.0 / sqrt(cols - 1);
  }


  // copy everything to a SPARSE_PETSC_MATRIX
  _A.equals(data);

  // finalize the current matrix
  _A.finalize();

  static bool first = true;
  
  // Initialize the solver
  if (first)
  {
    SlepcInitialize(NULL,NULL,NULL,NULL);
    first = false;
  }

  // create the SVD solver
  SVD svd;
  SVDCreate(PETSC_COMM_WORLD, &svd);

  // set it to solve the current matrix
  SVDSetOperator(svd, _A.petscData());

  // set tolerances
  int maxIterations = 1000000;
  SVDSetTolerances(svd, 1e-8, maxIterations);

  // set number of singular values to solve for
  SVDSetDimensions(svd, rank, PETSC_DECIDE, PETSC_DECIDE);
  SVDSetFromOptions(svd);
  SVDSolve(svd);

  int its, maxit;
  PetscReal tol;
  SVDGetIterationNumber(svd, &its);
  SVDGetTolerances(svd,&tol,&maxit);

  // see how many triplets converged
  PetscInt numConverged;
  SVDGetConverged(svd, &numConverged);

  // see how many we want to keep
  int totalKept = (numConverged < rank) ? numConverged : rank;

  // resize the results containers
  components.resizeAndWipe(rows, totalKept);
  values.resizeAndWipe(totalKept);

  // copy the results into the results containers
  Vec u,v;
  VecCreateSeq(PETSC_COMM_SELF, rows, &u);
  VecCreateSeq(PETSC_COMM_SELF, cols, &v);
  PetscReal sigma;
  for (int x = 0; x < totalKept; x++)
  {
    PetscInt i = x;
    SVDGetSingularTriplet(svd, i, &sigma, u, v);

    // copy the V into the container
    PetscScalar* uArray;
    VecGetArray(u, &uArray);
    for (int y = 0; y < rows; y++)
      components(y,x) = uArray[y];
    VecRestoreArray(u, &uArray);
    
    // principal components are the singular values squared
    values[x] = sigma * sigma;
  }
  // clean up
  SVDDestroy(svd);
  VecDestroy(u);
  VecDestroy(v);
}

//////////////////////////////////////////////////////////////////////
// Distribute the _UBasis of _originalMesh to the partitions
//////////////////////////////////////////////////////////////////////
void PARTITIONED_CUBATURE_GENERATOR::distributeSubbases(MATRIX& modalDerivatives, VECTOR& modalEigenvalues, vector<int> maxRanks)
{
  // get the U basis
  SUBSPACE_TET_MESH* originalMesh = (SUBSPACE_TET_MESH*)_partitionedMesh->originalMesh();
  MATRIX& UBasis = originalMesh->U();
  VECTOR& eigenvalues = originalMesh->eigenvalues();

  // create a vector of all the combined eigenvalues
  VECTOR combinedEigenvalues(eigenvalues.size() + modalEigenvalues.size());
  for (int x = 0; x < eigenvalues.size(); x++)
    combinedEigenvalues[x] = eigenvalues[x];
  for (int x = 0; x < modalEigenvalues.size(); x++)
    combinedEigenvalues[eigenvalues.size() + x] = modalEigenvalues[x];

  // for each partition, extract the rows corresponding to each
  // vertex and concatenate them into a basis
  for (int x = 0; x < _partitions; x++)
  {
    // get the basis of the partition
    SUBSPACE_TET_MESH* mesh = (SUBSPACE_TET_MESH*)_partitionedMesh->mesh(x);
    MATRIX& subbasis = mesh->U();
    VECTOR& subvalues = mesh->eigenvalues();

    // resize it appropriately
    int partitionNodes = mesh->unconstrainedNodes();
    subbasis.resizeAndWipe(partitionNodes * 3, UBasis.cols() + modalDerivatives.cols());

    // for each vertex in the partition
    for (int y = 0; y < partitionNodes; y++)
    {
      // get the vertex index in the original mesh
      int originalID = _partitionedMesh->originalID(x, y);

      // get the submatrix corresponding to that vertex
      SUBMATRIX vertexBasis(UBasis, originalID * 3, 3);

      // copy the vertex basis into the subbasis
      vertexBasis.copiesInto(subbasis, 3 * y);

      // get the submatrix from the modal derivatives
      SUBMATRIX derivativeBasis(modalDerivatives, originalID * 3, 3);
      
      // copy the modal derivatives into the subbasis
      derivativeBasis.copiesInto(subbasis, 3 * y, UBasis.cols());
    }

    // copy all the eigenvalues
    subvalues = combinedEigenvalues;

    // stomp the q size so the resize doesn't die
    mesh->q().resizeAndWipe(0);
    mesh->qOld().resizeAndWipe(0);
  }

  pcaSubbases(maxRanks);
}

//////////////////////////////////////////////////////////////////////
// Distribute the _UBasis of _originalMesh to the partitions
//////////////////////////////////////////////////////////////////////
void PARTITIONED_CUBATURE_GENERATOR::distributeSubbases(vector<int> maxRanks)
{
  // get the U basis
  SUBSPACE_TET_MESH* originalMesh = (SUBSPACE_TET_MESH*)_partitionedMesh->originalMesh();
  MATRIX& UBasis = originalMesh->U();
  VECTOR& eigenvalues = originalMesh->eigenvalues();

  // for each partition, extract the rows corresponding to each
  // vertex and concatenate them into a basis
  for (int x = 0; x < _partitions; x++)
  {
    // get the basis of the partition
    SUBSPACE_TET_MESH* mesh = (SUBSPACE_TET_MESH*)_partitionedMesh->mesh(x);
    MATRIX& subbasis = mesh->U();
    VECTOR& subvalues = mesh->eigenvalues();

    // resize it appropriately
    int partitionNodes = mesh->unconstrainedNodes();
    subbasis.resizeAndWipe(partitionNodes * 3, UBasis.cols());

    // for each vertex in the partition
    for (int y = 0; y < partitionNodes; y++)
    {
      // get the vertex index in the original mesh
      int originalID = _partitionedMesh->originalID(x, y);

      // get the submatrix corresponding to that vertex
      SUBMATRIX vertexBasis(UBasis, originalID * 3, 3);

      // copy the vertex basis into the subbasis
      vertexBasis.copiesInto(subbasis, 3 * y);
    }

    // copy all the eigenvalues
    subvalues = eigenvalues;
  }

  pcaSubbases(maxRanks);
}

//////////////////////////////////////////////////////////////////////
// Assuming the raw basis has been set for each partition,
// run PCA on them
//////////////////////////////////////////////////////////////////////
void PARTITIONED_CUBATURE_GENERATOR::pcaSubbasis(vector<int> maxRanks, int partition)
{
  int x = partition;

  // get the basis of the partition
  SUBSPACE_TET_MESH* mesh = (SUBSPACE_TET_MESH*)_partitionedMesh->mesh(x);
  MATRIX subbasis = mesh->U();
  VECTOR weights = mesh->eigenvalues();

  MATRIX testMatrix(subbasis);
  VECTOR rowMeans = subbasis.rowMeans();

  VECTOR pcaValues;
  MATRIX pcaComponents;
  pcaLapack(subbasis, subbasis.cols(), pcaComponents, pcaValues);
  VECTOR values(pcaValues);
  MATRIX components(pcaComponents);

  // shouldn't need to do this -- the true data mean is the origin, so the
  // snapshots are all pure variance
  VECTOR::printVertical = false;
  cout << " PCA results: " << values << endl;

  // Relative PCA cutoff
  // ORIGINAL VALUE IS 1e-5
  //Real relativeCutoff = values[0] * 1e-5;
  //Real relativeCutoff = values[0] * 1e-7;
  //Real relativeCutoff = 1.0;
  //Real relativeCutoff = 1e-5;
  //Real relativeCutoff = 1e-7;
  //Real relativeCutoff = 1e-12;
  Real relativeCutoff = 1e-16;

  int significant = 0;
  vector<int> significantIndices;
  for (int y = 0; y < values.size(); y++)
    if (((values(y) > relativeCutoff) || (values(y) > 10.0)) && !isnan(values(y)))
    {
      significant++;
      significantIndices.push_back(y);
    }
  cout << " Partition " << x << " has " << significant << " of " << values.size() << " significant components" << endl;
  if (significant > maxRanks[x])
    significant = maxRanks[x];
  cout << " Using " << significant << " of these components" << endl;

  // build the new basis
  MATRIX newBasis(subbasis.rows(), significant);
  for (int j = 0; j < subbasis.rows(); j++)
    for (int i = 0; i < significant; i++)
      newBasis(j,i) = components(j,significantIndices[i]);

  newBasis.orthogonalize();

  // build the new principal values
  VECTOR newValues(significant);
  for (int i = 0; i < significant; i++)
    newValues[i] = values[significantIndices[i]];

  mesh->updateBasis(newBasis);
  mesh->eigenvalues() = newValues;

  cout << " new eigenvalues: " << newValues << endl;

  // check the projection error
  Real meanDiffNorm = 0;
  Real maxDiffNorm = 0;
  Real meanRelativeDiffNorm = 0;
  Real maxRelativeDiffNorm = 0;
  vector<Real> diffNorms;
  vector<Real> relativeDiffNorms;

  Real maxRelativeActualNorm;
  for (int i = 0; i < testMatrix.cols(); i++)
  {
    VECTOR testColumn = testMatrix.getColumn(i);

    VECTOR testProjection = newBasis * (newBasis ^ testColumn);
    VECTOR diff = testColumn - testProjection;
    
    diffNorms.push_back(diff.norm2());

    meanDiffNorm += diff.norm2();

    Real relativeError = diff.norm2() / testColumn.norm2();
    meanRelativeDiffNorm += relativeError;

    relativeDiffNorms.push_back(relativeError);

    if (diff.norm2() > maxDiffNorm)
      maxDiffNorm = diff.norm2();

    if (relativeError > maxRelativeDiffNorm)
    {
      maxRelativeDiffNorm = relativeError;
      maxRelativeActualNorm = maxDiffNorm;
    }
  }
  cout << "======================================================" << endl;
  cout << " PCA test results " << endl;
  cout << "======================================================" << endl;
  cout << " Mean absolute projection error: " << meanDiffNorm / testMatrix.cols() << endl;
  cout << " Max absolute projection error: " << maxDiffNorm << endl ;
  cout << " Mean relative projection error: " << meanRelativeDiffNorm / testMatrix.cols() << endl;
  cout << " Max relative projection error: " << maxRelativeDiffNorm << endl;
  cout << " Actual norm of max relative error sample: " << maxRelativeActualNorm << endl;
  cout << "======================================================" << endl;
}

//////////////////////////////////////////////////////////////////////
// Assuming the raw basis has been set for each partition,
// run PCA on them
//////////////////////////////////////////////////////////////////////
void PARTITIONED_CUBATURE_GENERATOR::pcaSubbases(vector<int> maxRanks)
{
  // PCA each of the subbases, see if they are degenerate
  for (int x = 0; x < _partitions; x++)
    pcaSubbasis(maxRanks, x);
}

//////////////////////////////////////////////////////////////////////
// Generate the cubature points using PCA data
//////////////////////////////////////////////////////////////////////
void PARTITIONED_CUBATURE_GENERATOR::generatePCACubatureDebug(int partition)
{
  // retrieve the specific mesh and forces for this partition
  SUBSPACE_TET_MESH* tetMesh = (SUBSPACE_TET_MESH*)(_partitionedMesh->mesh(partition));
  vector<VECTOR>& trainingForces = _trainingForcesPartitioned[partition];
  VECTOR& forceMagnitudes = _forceMagnitudesPartitioned[partition];

  if (trainingForces.size() == 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << "No training forces found! Did you call ONLINE_SUBSPACE_INTEGRATOR::populateTrainingSet()?" << endl;
  }
  
  // size of a single force vector
  int forceSize = trainingForces[0].size();
  
  // total number of force samples
  int totalSamples = trainingForces.size();

  // number of rows in a column where all the samples are flattened
  // into a big vector
  int totalRows = forceSize * totalSamples;

  // calculate the squared magnitudes of the samples
  forceMagnitudes.resizeAndWipe(totalSamples);
  for (int x = 0; x < totalSamples; x++)
    forceMagnitudes(x) = trainingForces[x] * trainingForces[x];

  // calculate the inverse sqrt of the squared magnitudes
  for (int x = 0; x < totalSamples; x++)
  {
    // CAREFUL WITH SINGLE PRECISION HERE
    if (fabs(forceMagnitudes(x)) > 1e-7)
       forceMagnitudes(x) = 1.0 / sqrt(forceMagnitudes(x));
    else
       forceMagnitudes(x) = 1.0;
  }

  // stack normalized force samples into one big 'b' vector
  VECTOR b(totalRows);
  int index = 0;
  for (int x = 0; x < totalSamples; x++)
  {
    for (int y = 0; y < forceSize; y++, index++)
      b(index) = trainingForces[x](y) * forceMagnitudes(x);
  }

  // init the residual to b, as is traditional
  VECTOR residual(b);

  // "A" matrix to send to the the NNLS solver
  double* A = new double[totalRows * _maxKeyTets];

  // the Non-Negative Least Squares (NNLS) solver
  NNLS_SOLVER nnls(totalRows, _maxKeyTets);
  nnls.maxIter() = 10;

  // scratch space for the NNLS solver
  double* bNNLS = new double[totalRows];
  double* weightsNNLS = new double[_maxKeyTets];
  for (int x = 0; x < _maxKeyTets; x++)
    weightsNNLS[x] = 0.0;

  // error metrics
  Real bNorm = b.norm2();
  double rNorm = bNorm;
  double relativeError = 1.0;

  // retrieve the specific key tet and weight vectors for this partitions
  vector<TET*>& keyTets = _keyTets[partition];
  vector<Real>& keyWeights = _keyWeights[partition];

  // incorporate previous key tets
  int totalKeyTets = keyTets.size();
  cout << " Total starting key tets: " << totalKeyTets << endl;

  // start training
  int totalIterations = totalKeyTets;
  while (relativeError > _errorTolerance && 
         totalKeyTets < _maxKeyTets &&
         totalKeyTets < tetMesh->totalTets())
  {
    // get new set of randomly selected candidates
    if (tetMesh->totalTets() - totalKeyTets < _candidatesPerTry)
      _candidatesPerTry = tetMesh->totalTets() - totalKeyTets;

    vector<TET*> candidates = pickCandidatesPCA(partition, _candidatesPerTry);

    // pick the best tet out of the list of candidates
    int bestTet = findBestTetPCA(partition, candidates, residual);

    // store the new key tet
    keyTets.push_back(candidates[bestTet]);
    totalKeyTets++;

    VECTOR bestTetColumn = getTrainingColumnPCA(partition, candidates[bestTet]);

    // copy the column to "A" matrix --
    // note that A is column major, so it appears to be copying to a row
    for (int x = 0; x < totalRows; x++)
      A[(totalKeyTets - 1) * totalRows + x] = bestTetColumn(x);

    // copy b to double precision, also since the NNLS_SOLVER is going
    // to clobber it
    for (int x = 0; x < totalRows; x++)
      bNNLS[x] = b(x);

    // do the NNLS fit - note that both weightsNNLS and rNorm
    // both return values
    nnls.solve(A, totalKeyTets, bNNLS, weightsNNLS, rNorm);

    // update residual
    relativeError = rNorm / bNorm;

    // calculate the residual by doing an explicit A*x, 
    // but since A is column major, it is really (A^T) * x
    residual = b;
    for (int i = 0; i < totalKeyTets; i++)
      for (int j = 0; j < totalRows; j++)
      {
        int index = i * totalRows + j;
        residual(j) -= weightsNNLS[i] * A[index];
      }

    VECTOR currentForce(forceSize);
    for (int i = 0; i < totalKeyTets; i++)
      for (int j = 0; j < forceSize; j++)
      {
        int index = i * totalRows + j;
        currentForce(j) += weightsNNLS[i] * A[index];
      }
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " estimate: " << (1.0 / forceMagnitudes[0]) * currentForce << endl;
    cout << " training force: " << trainingForces[0] << endl;

    totalIterations++;

    cout << "Key Tet " << totalIterations << ": \tWeight: " << weightsNNLS[totalKeyTets-1] 
         << "\tError: " << relativeError << endl;
  }

  // copy the key tet weights to single precision
  keyWeights.clear();
  for (int x = 0; x < totalKeyTets; x++)
    keyWeights.push_back(weightsNNLS[x]);

  // delete the key tets with weight 0
  cout << " Total key tets before filtering: " << keyTets.size() << endl;
  filterKeyTets(partition);
  cout << " Total key tets after filtering: " << keyTets.size() << endl;

  // clean up scratch space
  delete[] bNNLS;
  delete[] weightsNNLS;
  delete[] A;
  cout << "=============================================" << endl;
}

//////////////////////////////////////////////////////////////////////
// Generate the cubature points using PCA data
//////////////////////////////////////////////////////////////////////
void PARTITIONED_CUBATURE_GENERATOR::generatePCACubature(int partition)
{
  // retrieve the specific mesh and forces for this partition
  SUBSPACE_TET_MESH* tetMesh = (SUBSPACE_TET_MESH*)(_partitionedMesh->mesh(partition));
  vector<VECTOR>& trainingForces = _trainingForcesPartitioned[partition];
  VECTOR& forceMagnitudes = _forceMagnitudesPartitioned[partition];

  if (trainingForces.size() == 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << "No training forces found! Did you call ONLINE_SUBSPACE_INTEGRATOR::populateTrainingSet()?" << endl;
  }
  
  // size of a single force vector
  int forceSize = trainingForces[0].size();
  
  // total number of force samples
  int totalSamples = trainingForces.size();

  // number of rows in a column where all the samples are flattened
  // into a big vector
  int totalRows = forceSize * totalSamples;

  // calculate the squared magnitudes of the samples
  forceMagnitudes.resizeAndWipe(totalSamples);
  for (int x = 0; x < totalSamples; x++)
    forceMagnitudes(x) = trainingForces[x] * trainingForces[x];

  // calculate the inverse sqrt of the squared magnitudes
  for (int x = 0; x < totalSamples; x++)
  {
    // CAREFUL WITH SINGLE PRECISION HERE
    if (fabs(forceMagnitudes(x)) > 1e-7)
       forceMagnitudes(x) = 1.0 / sqrt(forceMagnitudes(x));
    else
       forceMagnitudes(x) = 1.0;
  }

  // stack normalized force samples into one big 'b' vector
  VECTOR b(totalRows);
  int index = 0;
  for (int x = 0; x < totalSamples; x++)
    for (int y = 0; y < forceSize; y++, index++)
      b(index) = trainingForces[x](y) * forceMagnitudes(x);

  // init the residual to b, as is traditional
  VECTOR residual(b);

  // "A" matrix to send to the the NNLS solver
  double* A = new double[totalRows * _maxKeyTets];

  // the Non-Negative Least Squares (NNLS) solver
  NNLS_SOLVER nnls(totalRows, _maxKeyTets);
  nnls.maxIter() = 10;

  // scratch space for the NNLS solver
  double* bNNLS = new double[totalRows];
  double* weightsNNLS = new double[_maxKeyTets];
  for (int x = 0; x < _maxKeyTets; x++)
    weightsNNLS[x] = 0.0;

  // error metrics
  Real bNorm = b.norm2();
  double rNorm = bNorm;
  double relativeError = 1.0;

  // retrieve the specific key tet and weight vectors for this partitions
  vector<TET*>& keyTets = _keyTets[partition];
  vector<Real>& keyWeights = _keyWeights[partition];

  // incorporate previous key tets
  int totalKeyTets = keyTets.size();
  cout << " Total starting key tets: " << totalKeyTets << endl;

  // start training
  int totalIterations = totalKeyTets;
  while (relativeError > _errorTolerance && 
         totalKeyTets < _maxKeyTets &&
         totalKeyTets < tetMesh->totalTets())
  {
    // get new set of randomly selected candidates
    if (tetMesh->totalTets() - totalKeyTets < _candidatesPerTry)
      _candidatesPerTry = tetMesh->totalTets() - totalKeyTets;

    vector<TET*> candidates = pickCandidatesPCA(partition, _candidatesPerTry);

    // pick the best tet out of the list of candidates
    int bestTet = findBestTetPCA(partition, candidates, residual);

    // store the new key tet
    keyTets.push_back(candidates[bestTet]);
    totalKeyTets++;

    VECTOR bestTetColumn = getTrainingColumnPCA(partition, candidates[bestTet]);

    // copy the column to "A" matrix --
    // note that A is column major, so it appears to be copying to a row
    for (int x = 0; x < totalRows; x++)
      A[(totalKeyTets - 1) * totalRows + x] = bestTetColumn(x);

    // copy b to double precision, also since the NNLS_SOLVER is going
    // to clobber it
    for (int x = 0; x < totalRows; x++)
      bNNLS[x] = b(x);

    // do the NNLS fit - note that both weightsNNLS and rNorm
    // both return values
    nnls.solve(A, totalKeyTets, bNNLS, weightsNNLS, rNorm);

    // update residual
    relativeError = rNorm / bNorm;

    // calculate the residual by doing an explicit A*x, 
    // but since A is column major, it is really (A^T) * x
    residual = b;
    for (int i = 0; i < totalKeyTets; i++)
      for (int j = 0; j < totalRows; j++)
      {
        int index = i * totalRows + j;
        residual(j) -= weightsNNLS[i] * A[index];
      }
    totalIterations++;

    cout << "Key Tet " << totalIterations << ": \tWeight: " << weightsNNLS[totalKeyTets-1] 
         << "\tError: " << relativeError << endl;
  }
  if (totalKeyTets == tetMesh->totalTets())
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Ran out of tets to add!!! " << endl;
    cout << " If error is not zero at this point, there's a bug. " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  }

  // copy the key tet weights to single precision
  keyWeights.clear();
  for (int x = 0; x < totalKeyTets; x++)
    keyWeights.push_back(weightsNNLS[x]);

  // delete the key tets with weight 0
  cout << " Total key tets before filtering: " << keyTets.size() << endl;
  filterKeyTets(partition);
  cout << " Total key tets after filtering: " << keyTets.size() << endl;

  // clean up scratch space
  delete[] bNNLS;
  delete[] weightsNNLS;
  delete[] A;
  cout << "=============================================" << endl;
}

//////////////////////////////////////////////////////////////////////
// Retrieve the explicitly stored training column
//////////////////////////////////////////////////////////////////////
VECTOR PARTITIONED_CUBATURE_GENERATOR::retrieveTrainingColumn(int partition, TET* tet)
{
  SUBSPACE_TET_MESH* tetMesh = (SUBSPACE_TET_MESH*)(_partitionedMesh->mesh(partition));
  vector<VECTOR>& trainingForces = _trainingForcesPartitioned[partition];
  vector<VECTOR>& trainingColumns = _trainingColumnsPartitioned[partition];
  VECTOR& forceMagnitudes = _forceMagnitudesPartitioned[partition];

  int sampleSize = trainingForces[0].size();
  int totalRows = trainingForces.size() * sampleSize;
  VECTOR gColumn(totalRows);

  // submatrix of _Ubasis corresponding to this tet
  MATRIX tetU = tetMesh->tetSubBasis(tet);

  // generate column corresponding to stacked force samples
  for (unsigned int x = 0; x < trainingForces.size(); x++)
  {
    // get the density vector
    VECTOR& forceColumn = trainingColumns[x];

    // get the tet index in the vector
    int tetID = tetMesh->tetID(tet);

    // copy the force into a 12 vector
    int index = 12 * tetID;
    VECTOR force(12);
    for (int y = 0; y < 12; y++)
      force[y] = forceColumn[index + y];
 
    // get reduced version of force density
    VECTOR forceReduced = tetU ^ force;

    // copy this reduced force vector into the correct
    // place of the overall training column, making sure to 
    // normalize by the same factor as in generateCubature()
    for (int y = 0; y < sampleSize; y++)
      gColumn(x * sampleSize + y) = forceReduced(y) * forceMagnitudes(x);
  }

  return gColumn;
}

//////////////////////////////////////////////////////////////////////////////
// discretely sample the rigid rotations and factor them out
//////////////////////////////////////////////////////////////////////////////
void PARTITIONED_CUBATURE_GENERATOR::projectOutRotations(int partition)
{
  if (_partitionedMesh->constrained(partition))
    return;

  cout << " Projecting rotations out of partition " << partition << endl;

  // wipe the mesh to the rest state
  UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_partitionedMesh->mesh(partition);
  mesh->q() *= 0;
  mesh->rigidTranslation() *= 0;
  mesh->rotationQuaternion() = QUATERNION();
  mesh->updateFullMesh();

  // which axes to sample
  vector<VEC3F> axes;
  axes.push_back(VEC3F(1,0,0));
  axes.push_back(VEC3F(0,1,0));
  axes.push_back(VEC3F(0,0,1));

  // how many samples per axis
  int totalSlices = 100;

  // rest pose positions
  vector<VEC3F>& restPose = mesh->restPose();

  // project out the vectors
  VECTOR toProject(restPose.size() * 3);
  for (int axis = 0; axis < axes.size(); axis++)
  {
    for (int thetaSlices = 0; thetaSlices < totalSlices; thetaSlices++)
    {
      // build the vector to project out
      Real theta = ((Real)thetaSlices) / totalSlices * 2.0 * M_PI;
      MATRIX3 rotation = MATRIX3::rotation(axes[axis], theta);
      for (int x = 0; x < restPose.size(); x++)
      {
        VEC3F rotated = rotation * restPose[x];

        toProject[3 * x] = rotated[0];
        toProject[3 * x + 1] = rotated[1];
        toProject[3 * x + 2] = rotated[2];
      }

      // get the projection
      VECTOR projection = mesh->U() ^ toProject;

      // subtract out from the basis
      for (int x = 0; x < mesh->U().cols(); x++)
      {
        VECTOR rotationComponent = projection[x] * toProject;
        rotationComponent *= -1.0;
        mesh->U().addColumn(rotationComponent, x);
      }

      // remember to re-orthogonalize
      mesh->U().orthogonalize();

      if (thetaSlices % (totalSlices / 10) == 0)
      {
        cout << (100.0 * thetaSlices / totalSlices) << "% "; flush(cout);
      }
    }
    cout << " Done projecting axis " << axes[axis] << endl;
  }
}
