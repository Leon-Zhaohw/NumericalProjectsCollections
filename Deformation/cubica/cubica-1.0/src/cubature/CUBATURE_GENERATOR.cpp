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
// CUBATURE_GENERATOR.h: interface for the CUBATURE_GENERATOR class.
//
//////////////////////////////////////////////////////////////////////

#include "CUBATURE_GENERATOR.h"
#include <float.h>

//////////////////////////////////////////////////////////////////////
// Constructor for tetrahedra
//////////////////////////////////////////////////////////////////////
CUBATURE_GENERATOR::CUBATURE_GENERATOR(SUBSPACE_TET_MESH* tetMesh,
    const char* filename) :
  // init Mersenne Twister with deterministic seed
  // so that runs are reproducible
  _tetMesh(tetMesh),
  _twister(123456),
  _tolerance(1e-15),
  _maxKeyTets(1000),
  _candidatesPerTry(100),
  _errorTolerance(0.01f),
  _maxIterationsNNLS(4)
{
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
  }
}

//////////////////////////////////////////////////////////////////////
// Generate a basis for the mesh using linear modal analysis (LMA)
//////////////////////////////////////////////////////////////////////
void CUBATURE_GENERATOR::generateLMABasis(int totalModes, bool finalize)
{
  // retrieve the raw data
  MATRIX& UBasis = _tetMesh->U();
  VECTOR& eigenvalues = _tetMesh->eigenvalues();
  string meshFilename = _tetMesh->filename();
  PetscInt N = _tetMesh->stiffnessMatrix().rows();

  Mat         	 A, B;
  EPS         	 eps;
  PetscReal   	 error, tol;
  int         	 maxit, i, its;

  // resize the result arrays
  eigenvalues.resizeAndWipe(totalModes);
  UBasis.resizeAndWipe(N, totalModes);

  SlepcInitialize(NULL,NULL,NULL,NULL);

  // get preallocation hints for PetSc
  vector<int> nonZerosPerRow = _tetMesh->stiffnessMatrix().nonZerosPerRow();
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
  _tetMesh->stiffnessMatrix().entries(rows, cols, values);
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
  cout << " Stiffness size: " << _tetMesh->stiffnessMatrix().rows() << " x " << _tetMesh->stiffnessMatrix().cols() << endl;
  flush(cout);

  // copy entries into mass matrix
  // Do it a generic but memory-hogging way
  cout << " Setting up mass eigensystem ... ";
  _tetMesh->massMatrix().entries(rows, cols, values);
  for (unsigned int x = 0; x < rows.size(); x++)
  {
    PetscScalar 	 v = values[x];

    if (isnan(values[x]))
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " NAN encountered! " << endl;
    }
    MatSetValues(B,1,&(rows[x]),1,&(cols[x]),&v,INSERT_VALUES);
  }
  cout << " done. " << endl;
  cout << " Mass size: " << _tetMesh->massMatrix().rows() << " x " << _tetMesh->massMatrix().cols() << endl;

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

  // total number of eigenvalues to search for --
  // pad by 6 for the rigid modes if the mesh is
  // unconstrained
  int totalEigenvalues = totalModes;
  if (!_tetMesh->constrained())
  {
    cout << " Solving for 6 extra modes because the mesh is unconstrained." << endl;
    totalEigenvalues += 6;
  }

  // set the tolerances
  _tolerance = 1e-2;
  int maxIterations = 1000000;
  EPSSetTolerances(eps, _tolerance, maxIterations);

  // Slepc recommends using at least twice as many eigenvectors
  // for the search space
  EPSSetDimensions(eps, totalEigenvalues, PETSC_DECIDE, PETSC_DECIDE);
  EPSSetWhichEigenpairs(eps, EPS_SMALLEST_MAGNITUDE);

  PetscTruth hermitian = PETSC_TRUE;
  EPSIsHermitian(eps, &hermitian);

  // Send Slepc the solver parameters
  EPSSetFromOptions(eps);

  // Solve the eigensystem
  TIMER eigenTimer;
  cout << " Solving eigensystem ... ";
  EPSSolve(eps);
  EPSGetIterationNumber(eps, &its);
  EPSGetTolerances(eps,&tol,&maxit);
  cout << "done in " << its << " of " << maxit << " max iterations " << endl;
  cout << "took " << eigenTimer.timing() / 60.0 << " minutes" << endl;

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
  if (finalize)
    SlepcFinalize();

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
void CUBATURE_GENERATOR::generateTrainingSamples(int totalSamples, Real magnitude)
{
  VECTOR& eigenvalues = _tetMesh->eigenvalues();
  VECTOR q(eigenvalues.size());

  // filter the eigenvalues (in case we're using OGDEN)
  int firstNonZero = 0;
  for (int x = 0; x < eigenvalues.size(); x++)
    if (eigenvalues[x] > 0.0)
    {
      firstNonZero = x;
      break;
    }
  for (int x = 0; x < firstNonZero; x++)
    eigenvalues[x] = eigenvalues[firstNonZero];  

  // the user wants the first mode to have the given magnitude
  _magnitude = magnitude * sqrt(eigenvalues(0));

  int attempts = 0;
  int successes = 0;

  // clear old samples just in case
  _trainingForces.clear();
  _trainingQs.clear();

  // track if a Nan has been seen before so that the error message
  // doesn't get output too many times
  bool firstNan = true;

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
        q(x) = _magnitude / sqrt(eigenvalues(x)) / 4.0 * _twister.randNorm();

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

      _tetMesh->q() = q;
      _tetMesh->updateFullMesh();

      // check for inversion
      tetInverted = _tetMesh->inverted();

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
    VECTOR forceSample = _tetMesh->projectedInternalForce();
    _trainingForces.push_back(forceSample);
    _trainingQs.push_back(q);

    for (int x = 0; x < forceSample.size(); x++)
      if (isnan(forceSample[x]))
      {
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
        cout << " NAN found in force sample! " << endl;
        cout << forceSample << endl;
        exit(1);
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

  // cache results to the default filename
  write(_filename.c_str());
}

//////////////////////////////////////////////////////////////////////
// write out a training set
//////////////////////////////////////////////////////////////////////
void CUBATURE_GENERATOR::write(const char* filename)
{
  FILE* file;
  file = fopen(filename, "wb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __LINE__ << " Couldn't write file: " << filename << endl;
    return;
  }

  // write dimensions
  int size = _trainingForces.size();
  fwrite((void*)&size, sizeof(int), 1, file);

  // write out the forces
  for (int x = 0; x < size; x++)
    _trainingForces[x].write(file);

  // write out the qs
  for (int x = 0; x < size; x++)
    _trainingQs[x].write(file);

  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// read in a training set
//////////////////////////////////////////////////////////////////////
void CUBATURE_GENERATOR::read(const char* filename)
{
  FILE* file;
  file = fopen(filename, "rb");

  // read dimensions
  int size;
  fread((void*)&size, sizeof(int), 1, file);

  // read in the forces
  for (int x = 0; x < size; x++)
  {
    VECTOR sample(file);
    _trainingForces.push_back(sample);
  }

  // read in the qs
  for (int x = 0; x < size; x++)
  {
    VECTOR sample(file);
    _trainingQs.push_back(sample);
  }

  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// write out the generated cubature
//////////////////////////////////////////////////////////////////////
void CUBATURE_GENERATOR::writeCubature(const char* filename)
{
  string finalfile;
  
  // if no filename is provided, generate one based on the mesh name
  if (filename == NULL)
    finalfile = _tetMesh->filename() + string(".cubature");
  else
    finalfile = string(filename);

  // check that a cubature even exists to write out
  if (_keyTets.size() == 0)
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

  // write dimensions
  int size = _keyTets.size();
  fwrite((void*)&size, sizeof(int), 1, file);

  // write out the key tet indices
  for (int x = 0; x < size; x++)
  {
    int tetID = _tetMesh->tetID(_keyTets[x]);
    fwrite((void*)&tetID, sizeof(int), 1, file);
  }

  // write out the key tet weights
  for (int x = 0; x < size; x++)
  {
    double weight = _keyWeights[x];
    fwrite((void*)&weight, sizeof(double), 1, file);
  }
  
  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// pick candidates for key tet selection
//////////////////////////////////////////////////////////////////////
vector<TET*> CUBATURE_GENERATOR::pickCandidates(int _candidatesPerTry)
{
  vector<TET>& tets = _tetMesh->tets();
  int totalTets = _tetMesh->totalTets();

  // build a lookup table for the already selected key tets
  map<TET*, bool> alreadyUsed;
  for (unsigned int x = 0; x < _keyTets.size(); x++)
    alreadyUsed[_keyTets[x]] = true;

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

  // construct the final vector of candidate tets
  vector<TET*> candidates;
  for (int x = 0; x < _candidatesPerTry; x++)
  {
    // pick a random tet
    int index = _twister.randInt(listSize - 1);

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
VECTOR CUBATURE_GENERATOR::getTrainingColumn(TET* tet)
{
  int sampleSize = _trainingForces[0].size();
  int totalRows = _trainingForces.size() * sampleSize;
  VECTOR gColumn(totalRows);

  // submatrix of _Ubasis corresponding to this tet
  MATRIX tetU = _tetMesh->tetSubBasis(tet);

  // generate column corresponding to stacked force samples
  for (unsigned int x = 0; x < _trainingForces.size(); x++)
  {
    // generate the displacement from the forcing, but only for
    // the vertices associated with this tet
    VECTOR displacement = tetU * _trainingQs[x];
 
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
    MATERIAL* material = _tetMesh->materials()[tet->materialIndex()];
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
      gColumn(x * sampleSize + y) = forceReduced(y) * _forceMagnitudes(x);
  }

  return gColumn;
}

//////////////////////////////////////////////////////////////////////
// Find the best fitting tet from the list of candidates
//////////////////////////////////////////////////////////////////////
int CUBATURE_GENERATOR::findBestTet(vector<TET*>& candidates, VECTOR& residual)
{
  Real rNorm = residual.norm2();

  // candidate dot products;
  vector<Real> dots;
  dots.resize(candidates.size());

  // calculate all the dot products
  for (unsigned int x = 0; x < candidates.size(); x++)
  {
    // get the 'g' column corresponding to this tet
    VECTOR gColumn = getTrainingColumn(candidates[x]);
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
// Generate the cubature points
//////////////////////////////////////////////////////////////////////
void CUBATURE_GENERATOR::generateCubature()
{
  // UNCOMMENT FOR LEGACY REGRESSION TESTS
  //_twister.seed(123456);

  // see if a reasonable number of key tets were requested
  if (_maxKeyTets > _tetMesh->totalTets())
    cout << __FILE__ << " " << __LINE__ << " : " 
         << "You've asked for more key tets than there are in the mesh! " 
         << "This will probably loop forever." << endl;

  // reset the tet mesh to the rest pose since we will be
  // using its rest pose in getTrainingColumn()
  _tetMesh->resetToRestPose();

  // size of a single force vector
  int forceSize = _trainingForces[0].size();
  
  // total number of force samples
  int totalSamples = _trainingForces.size();

  // number of rows in a column where all the samples are flattened
  // into a big vector
  int totalRows = forceSize * totalSamples;

  // calculate the squared magnitudes of the samples
  _forceMagnitudes.resizeAndWipe(totalSamples);
  for (int x = 0; x < totalSamples; x++)
  {
    for (int y = 0; y < forceSize; y++)
    {
      // bracket paren notation []() is a little confusing, I'll try to 
      // keep it to a minimum
      _forceMagnitudes(x) += _trainingForces[x](y) * _trainingForces[x](y);
    }
  }

  // calculate the inverse sqrt of the squared magnitudes
  for (int x = 0; x < totalSamples; x++)
    // CAREFUL WITH SINGLE PRECISION HERE
    if (fabs(_forceMagnitudes(x)) > 1e-7)
       _forceMagnitudes(x) = 1.0f / sqrtf(_forceMagnitudes(x));
    else
       _forceMagnitudes(x) = 1.0f;

  // stack normalized force samples into one big 'b' vector
  VECTOR b(totalRows);
  int index = 0;
  for (int x = 0; x < totalSamples; x++)
    for (int y = 0; y < forceSize; y++, index++)
      b(index) = _trainingForces[x](y) * _forceMagnitudes(x);

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

  // tracking key tets
  _keyTets.clear();
  _keyWeights.clear();
  int totalKeyTets = 0;

  // start training
  int totalIterations = 0;
  while (relativeError > _errorTolerance && 
         totalKeyTets < _maxKeyTets &&
         totalKeyTets < _tetMesh->totalTets())
  {
    // get new set of randomly selected candidates
    if (_tetMesh->totalTets() - totalKeyTets < _candidatesPerTry)
      _candidatesPerTry = _tetMesh->totalTets() - totalKeyTets;
    vector<TET*> candidates = pickCandidates(_candidatesPerTry);

    // pick the best tet out of the list of candidates
    int bestTet = findBestTet(candidates, residual);

    // store the new key tet
    _keyTets.push_back(candidates[bestTet]);
    totalKeyTets++;

    // recompute its column (better than storing them all)
    VECTOR bestTetColumn = getTrainingColumn(candidates[bestTet]);

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

      cout << "Iteration " << totalIterations << ": \tWeight: " << weightsNNLS[totalKeyTets-1] 
           << "\tError: " << relativeError << endl;
    }
    else
    {
      cout << " NNLS did not converge! " << endl;

      totalKeyTets--;
      _keyTets.resize( totalKeyTets );
    }

    if (isnan(relativeError))
      relativeError = 1.0;
    totalIterations++;
  }

  // copy the key tet weights to single precision
  _keyWeights.clear();
  for (int x = 0; x < totalKeyTets; x++)
    _keyWeights.push_back(weightsNNLS[x]);

  cout << " Total key tets before filtering: " << _keyTets.size() << endl;
  filterKeyTets();
  cout << " Total key tets after filtering: " << _keyTets.size() << endl;

  // clean up scratch space
  delete[] bNNLS;
  delete[] weightsNNLS;
}

//////////////////////////////////////////////////////////////////////
// Filter out the dud tets
//////////////////////////////////////////////////////////////////////
void CUBATURE_GENERATOR::filterKeyTets()
{
  vector<TET*>& keyTets = _keyTets;
  vector<Real>& keyWeights = _keyWeights;
  
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
// Assuming the raw basis has been set in _tetMesh, run PCA on it
//////////////////////////////////////////////////////////////////////
void CUBATURE_GENERATOR::pcaBasis(int maxRank)
{
  // get the basis of the partition
  MATRIX subbasis = _tetMesh->U();
  VECTOR weights = _tetMesh->eigenvalues();

  MATRIX testMatrix(subbasis);
  VECTOR rowMeans = subbasis.rowMeans();

  VECTOR pcaValues;
  MATRIX pcaComponents;
  pcaLapack(subbasis, subbasis.cols(), pcaComponents, pcaValues);
  VECTOR values(pcaValues);
  MATRIX components(pcaComponents);

  VECTOR::printVertical = false;
  cout << " PCA results: " << values << endl;

  // Relative PCA cutoff
  Real relativeCutoff = 1e-5;

  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " DIALED PCA DISCARD BACK TO 1e-6" << endl;
  int significant = 0;
  vector<int> significantIndices;
  for (int y = 0; y < values.size(); y++)
    if (((values(y) > relativeCutoff) || (values(y) > 10.0)) && !isnan(values(y)))
    {
      significant++;
      significantIndices.push_back(y);
    }
  cout << " Mesh has " << significant << " of " << values.size() << " significant components" << endl;
  if (significant > maxRank)
    significant = maxRank;
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

  _tetMesh->updateBasis(newBasis);
  _tetMesh->eigenvalues() = newValues;

  cout << " new eigenvalues: " << newValues << endl;

  // check the projection error
  Real meanDiffNorm = 0;
  Real maxDiffNorm = 0;
  Real meanRelativeDiffNorm = 0;
  Real maxRelativeDiffNorm = 0;
  vector<Real> diffNorms;
  for (int i = 0; i < testMatrix.cols(); i++)
  {
    VECTOR testColumn = testMatrix.getColumn(i);

    VECTOR testProjection = newBasis * (newBasis ^ testColumn);
    VECTOR diff = testColumn - testProjection;
    
    diffNorms.push_back(diff.norm2());

    meanDiffNorm += diff.norm2();

    Real relativeError = diff.norm2() / testColumn.norm2();
    meanRelativeDiffNorm += relativeError;

    if (diff.norm2() > maxDiffNorm)
      maxDiffNorm = diff.norm2();

    if (relativeError > maxRelativeDiffNorm)
      maxRelativeDiffNorm = relativeError;
  }
  cout << "======================================================" << endl;
  cout << " PCA test results " << endl;
  cout << "======================================================" << endl;
  cout << " Mean absolute projection error: " << meanDiffNorm / testMatrix.cols() << endl;
  cout << " Max absolute projection error: " << maxDiffNorm << endl ;
  cout << " Mean relative projection error: " << meanRelativeDiffNorm / testMatrix.cols() << endl;
  cout << " Max relative projection error: " << maxRelativeDiffNorm << endl;
  cout << " diffs: " << VECTOR(diffNorms) << endl;
  cout << "======================================================" << endl;

  string meshFilename = _tetMesh->filename();
  string valuename = meshFilename + string(".eigenvalues.vector");
  string vectorname = meshFilename + string(".eigenvectors.matrix");
  cout << " Writing eigenvalues to file " << valuename.c_str() << endl;
  cout << " Writing eigenvectors to file " << vectorname.c_str() << endl;
  _tetMesh->eigenvalues().write(valuename.c_str());
  _tetMesh->U().write(vectorname.c_str());
}

//////////////////////////////////////////////////////////////////////////////
// Perform PCA on the data matrix, returning only the 1st "rank" columns
//////////////////////////////////////////////////////////////////////////////
void CUBATURE_GENERATOR::pcaLapack(MATRIX& data, int rank, MATRIX& components, VECTOR& values, bool shift)
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
