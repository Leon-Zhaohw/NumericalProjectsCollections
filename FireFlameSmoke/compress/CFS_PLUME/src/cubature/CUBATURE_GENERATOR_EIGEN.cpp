/*
 This file is part of SSFR (Zephyr).
 
 Zephyr is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 Zephyr is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Zephyr.  If not, see <http://www.gnu.org/licenses/>.
 
 Copyright 2013 Theodore Kim
 */
// CUBATURE_GENERATOR_EIGEN.h: interface for the CUBATURE_GENERATOR_EIGEN class.
//
//////////////////////////////////////////////////////////////////////

#include "CUBATURE_GENERATOR_EIGEN.h"
#include <float.h>

//////////////////////////////////////////////////////////////////////
// Constructor for tetrahedra
//////////////////////////////////////////////////////////////////////
CUBATURE_GENERATOR_EIGEN::CUBATURE_GENERATOR_EIGEN(SUBSPACE_FLUID_3D_EIGEN* fluid,
    const char* filename) :
  // init Mersenne Twister with deterministic seed
  // so that runs are reproducible
  _fluid(fluid),
  _twister(123456),
  _maxKeyCells(10000),
  _candidatesPerTry(10000),
  _errorTolerance(0.01f),
  _tolerance(1e-15),
  _maxIterationsNNLS(4),
  _grabAllCandidates(false),
  _importanceSamples(17000)
{
  cout << " Cubature generator received a fluid with dims: " << fluid->xPeeled() << " " << fluid->yPeeled() << " " << fluid->zPeeled() << endl;
  cout << " Cubature generator received a fluid with dims: " << fluid->xRes() << " " << fluid->yRes() << " " << fluid->zRes() << endl;
}

//////////////////////////////////////////////////////////////////////
// write out the generated cubature
//////////////////////////////////////////////////////////////////////
void CUBATURE_GENERATOR_EIGEN::writeCubature(const string& filename)
{
  cout << " Writing cubature to file " << filename.c_str() << " ... "; flush(cout);
  string finalfile;
  finalfile = string(filename);

  // check that a cubature even exists to write out
  if (_keyCells.size() == 0)
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
  int size = _keyCells.size();
  fwrite((void*)&size, sizeof(int), 1, file);

  // write out the key tet indices
  for (int x = 0; x < size; x++)
  {
    int cellIndex = _keyCells[x];
    fwrite((void*)&cellIndex, sizeof(int), 1, file);
  }

  // write out the key tet weights
  for (int x = 0; x < size; x++)
  {
    double weight = _keyWeights[x];
    fwrite((void*)&weight, sizeof(double), 1, file);
  }
  
  fclose(file);

  cout << " done." << endl;
}

//////////////////////////////////////////////////////////////////////
// pick candidates for key tet selection
//////////////////////////////////////////////////////////////////////
void CUBATURE_GENERATOR_EIGEN::pickImportanceSampledCandidatesOMP(const int pass, map<int, bool> alreadyUsed, const VectorXd& residual, vector<int>& candidates, int totalCandidates)
{
  TIMER functionTimer(__FUNCTION__);
  cout << " Picking importance sampled candidates ... " << endl;

  // construct a list of all non-key tets, but only iterate through
  // peeled center cells
  vector<int> nonKeyCells;
  int xRes = _fluid->xRes();
  int yRes = _fluid->yRes();
  int zRes = _fluid->zRes();

  int xResPeeled = xRes - 2;
  int yResPeeled = yRes - 2;
  int slabPeeled = xResPeeled * yResPeeled;

  for (int z = 0; z < zRes - 2; z++)
    for (int y = 0; y < yRes - 2; y++)
      for (int x = 0; x < xRes - 2; x++)
      {
        int index = x + y * xResPeeled + z * slabPeeled;

        // check that it is not already a key cell
        if (alreadyUsed.find(index) == alreadyUsed.end())
          nonKeyCells.push_back(index);
      }

  // record what key cells have already been picked
  vector<bool> usedAlready(nonKeyCells.size());
  for (unsigned int x = 0; x < usedAlready.size(); x++)
    usedAlready[x] = false;

  Real residualDot = residual.dot(residual);
  int rejected = 0;
  const int totalPossible = usedAlready.size();

  // get the total number of candidates each thread needs to generate
  int totalThreads = 1; 
  int candidatesPerThread = totalCandidates / totalThreads;
  vector<vector<int> > perThreadCandidates(totalThreads);

#pragma omp parallel
#pragma omp for  schedule(static)
  for (int x = 0; x < totalThreads; x++)
  {
    int threadID = 0; 
    MERSENNETWISTER twister(threadID + pass * 314159);

    int candidatesFound = 0;
    int tries = 0;
    while (candidatesFound < candidatesPerThread)
    {
      // pick a random tet
      int index = twister.randInt(nonKeyCells.size() - 1);
      tries++;

      // in case we pick the same one again by accident
      if (usedAlready[index])
        continue;

      // pick it based on a weighted probability
      VectorXd gColumn = getStagedTrainingColumn(nonKeyCells[index]);
    
      // div by just rNorm, not gNorm as well, since we just care 
      // how large it is relative to the total residual 
      Real dot = (gColumn.dot(residual)) / (residualDot / totalPossible);
      dot = fabs(dot);

      if (tries % 1000 == 0)
      {
#pragma omp critical
        {
        cout << " Thread " << x << " picked " << candidatesFound << " of " << candidatesPerThread << " candidates with " << rejected << " rejections " << endl;
        }
      }

      Real pick = _twister.rand();
      if (pick >= dot)
      {
        rejected++;
        continue;
      }

      // record it as a candidate tet
      perThreadCandidates[x].push_back(nonKeyCells[index]);
      candidatesFound++;

      // make sure it's not picked again, but as this is shared across
      // all threads, make it a critical region
#pragma omp critical
      {
        usedAlready[index] = true;
      }
    }
  }

  // merge all the candidates into one
  for (int x = 0; x < totalThreads; x++)
  {
    for (unsigned int y = 0; y < perThreadCandidates[x].size(); y++)
    {
      int toAdd = perThreadCandidates[x][y];

      // just to be safe, make sure some other thread didn't add this
      // at some point as well
      if (alreadyUsed.find(toAdd) == alreadyUsed.end())
      {
        candidates.push_back(toAdd);
        alreadyUsed[toAdd] = true;
      }
    }
  }
  cout << " Found " << candidates.size() << " candidates of the requested " << totalCandidates << endl;
}

//////////////////////////////////////////////////////////////////////
// pick candidates for key tet selection
//////////////////////////////////////////////////////////////////////
void CUBATURE_GENERATOR_EIGEN::pickImportanceSampledCandidatesOMPMemory(const int pass, map<int, bool> alreadyUsed, const VectorXd& residual, vector<int>& candidates, int totalCandidates)
{
  TIMER functionTimer(__FUNCTION__);
  cout << " Picking importance sampled candidates ... " << endl;

  // construct a list of all non-key tets, but only iterate through
  // peeled center cells
  vector<int> nonKeyCells;
  int xRes = _fluid->xRes();
  int yRes = _fluid->yRes();
  int zRes = _fluid->zRes();

  int xResPeeled = xRes - 2;
  int yResPeeled = yRes - 2;
  int slabPeeled = xResPeeled * yResPeeled;

  for (int z = 0; z < zRes - 2; z++)
    for (int y = 0; y < yRes - 2; y++)
      for (int x = 0; x < xRes - 2; x++)
      {
        int index = x + y * xResPeeled + z * slabPeeled;

        // check that it is not already a key cell
        if (alreadyUsed.find(index) == alreadyUsed.end())
          nonKeyCells.push_back(index);
      }

  // record what key cells have already been picked
  vector<bool> usedAlready(nonKeyCells.size());
  for (unsigned int x = 0; x < usedAlready.size(); x++)
    usedAlready[x] = false;

  Real residualDot = residual.dot(residual);
  int rejected = 0;
  const int totalPossible = usedAlready.size();

  // get the total number of candidates each thread needs to generate
  int totalThreads = 1; 
  int candidatesPerThread = totalCandidates / totalThreads;
  vector<vector<int> > perThreadCandidates(totalThreads);

#pragma omp parallel
#pragma omp for  schedule(static)
  for (int x = 0; x < totalThreads; x++)
  {
    int threadID = 0; 
    MERSENNETWISTER twister(threadID + pass * 314159);

    int candidatesFound = 0;
    int tries = 0;
    while (candidatesFound < candidatesPerThread)
    {
      // pick a random tet
      int index = twister.randInt(nonKeyCells.size() - 1);
      tries++;

      // in case we pick the same one again by accident
      if (usedAlready[index])
        continue;

      // pick it based on a weighted probability
      VectorXd gColumn = getStagedTrainingColumnMemory(nonKeyCells[index]);
    
      // div by just rNorm, not gNorm as well, since we just care 
      // how large it is relative to the total residual 
      Real dot = (gColumn.dot(residual)) / (residualDot / totalPossible);
      dot = fabs(dot);

      if (tries % 1000 == 0)
      {
#pragma omp critical
        {
        cout << " Thread " << x << " picked " << candidatesFound << " of " << candidatesPerThread << " candidates with " << rejected << " rejections " << endl;
        }
      }

      Real pick = _twister.rand();
      if (pick >= dot)
      {
        rejected++;
        continue;
      }

      // record it as a candidate tet
      perThreadCandidates[x].push_back(nonKeyCells[index]);
      candidatesFound++;

      // make sure it's not picked again, but as this is shared across
      // all threads, make it a critical region
#pragma omp critical
      {
        usedAlready[index] = true;
      }
    }
  }

  // merge all the candidates into one
  for (int x = 0; x < totalThreads; x++)
  {
    for (unsigned int y = 0; y < perThreadCandidates[x].size(); y++)
    {
      int toAdd = perThreadCandidates[x][y];

      // just to be safe, make sure some other thread didn't add this
      // at some point as well
      if (alreadyUsed.find(toAdd) == alreadyUsed.end())
      {
        candidates.push_back(toAdd);
        alreadyUsed[toAdd] = true;
      }
    }
  }
  cout << " Found " << candidates.size() << " candidates of the requested " << totalCandidates << endl;
}

//////////////////////////////////////////////////////////////////////
// Generate the training column corresponding to a tet
//////////////////////////////////////////////////////////////////////
VectorXd CUBATURE_GENERATOR_EIGEN::getStagedTrainingColumn(int index)
{
  TIMER functionTimer(__FUNCTION__);
  // DEBUG
  // puts("Inside vanilla getStagedTrainingColumn!");

  int xPeeled = _fluid->xRes() - 2;
  int yPeeled = _fluid->yRes() - 2;
  int slabPeeled = xPeeled * yPeeled;

  // decompose into x,y,z
  int decompose = index;
  int z = decompose / slabPeeled;
  decompose -= z * slabPeeled;
  int y = decompose / xPeeled;
  decompose -= y * xPeeled;

  Real dt0 = _fluid->dt() / _fluid->dx();

  int sampleSize = _trainingPostadvection[0].size();
  int totalRows = _trainingPostadvection.size() * sampleSize;
  VectorXd gColumn(totalRows);

  const MatrixXd& preU = _fluid->preadvectU();
  const MatrixXd& postU = _fluid->prediffuseU();

  // submatrix of _Ubasis corresponding to this tet
  MatrixXd cellPreU = _fluid->cellBasisPeeled(preU, index);
  MatrixXd cellPostU = _fluid->cellBasisPeeled(postU, index);

  // generate column corresponding to stacked force samples
  for (unsigned int x = 0; x < _trainingPostadvection.size(); x++)
  {
    // do an advection on just this cell, for the velocity field pair
    // trainingPreadvection[x], trainingPostadvection[x]
    VectorXd postVelocity = _fluid->advectCellStamPeeled(preU, dt0, _trainingPreadvection[x], index);

    // project the compute nonlinear velocity
    VectorXd projectedPostVelocity = cellPostU.transpose() * postVelocity;

    // copy this reduced force vector into the correct
    // place of the overall training column, making sure to 
    // normalize by the same factor as in generateCubature()
    for (int y = 0; y < sampleSize; y++)
      gColumn(x * sampleSize + y) = projectedPostVelocity(y) * _postadvectionMagnitudes(x);
  }

  return gColumn;
}


//////////////////////////////////////////////////////////////////////
// Generate the training column corresponding to a tet
//////////////////////////////////////////////////////////////////////
VectorXd CUBATURE_GENERATOR_EIGEN::getStagedTrainingColumnMemory(int index)
{
  TIMER functionTimer(__FUNCTION__);
  int xPeeled = _fluid->xRes() - 2;
  int yPeeled = _fluid->yRes() - 2;
  int slabPeeled = xPeeled * yPeeled;

  // decompose into x,y,z
  int decompose = index;
  int z = decompose / slabPeeled;
  decompose -= z * slabPeeled;
  int y = decompose / xPeeled;
  decompose -= y * xPeeled;

  Real dt0 = _fluid->dt() / _fluid->dx();

  int sampleSize = _trainingPostadvection[0].size();
  int totalRows = _trainingPostadvection.size() * sampleSize;
  VectorXd gColumn(totalRows);

  const MatrixXd& preU = _fluid->preadvectU();

  // submatrix of _Ubasis corresponding to this tet
  MatrixXd cellPreU = _fluid->cellBasisPeeled(preU, index);

  // this fills the member variable submatrix _postU_T

  _fluid->cellBasisPeeledMemory(_prediffuseU_rows, _prediffuseU_cols, _prediffuseFile,
      this->postU_T_data(), index);


  // generate column corresponding to stacked force samples
  for (unsigned int x = 0; x < _trainingPostadvection.size(); x++)
  {
    // do an advection on just this cell, for the velocity field pair
    // trainingPreadvection[x], trainingPostadvection[x]
    VectorXd postVelocity = _fluid->advectCellStamPeeled(preU, dt0, _trainingPreadvection[x], index);

    // project the compute nonlinear velocity
    // ADJ: we *don't* use transpose since we read it in already transposed!
    VectorXd projectedPostVelocity = _postU_T * postVelocity;

    // copy this reduced force vector into the correct
    // place of the overall training column, making sure to 
    // normalize by the same factor as in generateCubature()
    for (int y = 0; y < sampleSize; y++)
      gColumn(x * sampleSize + y) = projectedPostVelocity(y) * _postadvectionMagnitudes(x);
  }

  return gColumn;
}

//////////////////////////////////////////////////////////////////////////////
// call Matlab to solve the non-negative least squares problem using the 
// Bro and Jong fast modification
//////////////////////////////////////////////////////////////////////////////
VectorXd CUBATURE_GENERATOR_EIGEN::matlabFNNLS(const MatrixXd& A, const VectorXd& b, bool verbose)
{
  TIMER functionTimer(__FUNCTION__);
  cout << " Calling Matlab fast non-negative least squares ... " << flush;
  string matlabBin("/Applications/MATLAB_R2011a.app/bin/matlab");
  int cols = A.cols();
  int rows = A.rows();
  MATRIX inputMatrix(rows, cols);

  for (int y = 0; y < cols; y++)
    for (int x = 0; x < rows; x++)
      inputMatrix(x,y) = A(x,y);

  VECTOR inputVector = EIGEN::convert(b);

  // write it all out
  string inputFile = string("LS.matrix");
  inputMatrix.write(inputFile.c_str());

  inputFile = string("LS.vector");
  inputVector.write(inputFile.c_str());

  // copy the script to the top directory
  string copy("cp ./src/matlab/fnnls*.m .");
  system(copy.c_str());

  // run matlab script
  string matlabCall = matlabBin + string(" -nodisplay -nosplash -r fnnlsFromFile");
  if (!verbose)
    matlabCall = matlabCall + string(" >& matlabSplash.txt");
  system(matlabCall.c_str());

  // read result in
  string outputFile = string("LS_result.vector");
  VECTOR result(outputFile.c_str());

  // clean up the files
  string rm("rm LS.matrix LS.vector LS_result.vector fnnls*.m matlabSplash.txt");
  system(rm.c_str());

  cout << " done. " << endl;

  return EIGEN::convert(result);
}

//////////////////////////////////////////////////////////////////////////////
// call Matlab to solve the non-negative least squares problem using the 
// Bro and Jong fast modification
//////////////////////////////////////////////////////////////////////////////
VectorXd CUBATURE_GENERATOR_EIGEN::matlabFNNLS(const vector<VectorXd>& columns, const VectorXd& b, bool verbose)
{
  TIMER functionTimer(__FUNCTION__);
  cout << " Calling Matlab fast non-negative least squares ... " << flush;
  string matlabBin("/Applications/MATLAB_R2011a.app/bin/matlab");
  int cols = columns.size();
  int rows = columns[0].size();
  MATRIX inputMatrix(rows, cols);

  for (int y = 0; y < cols; y++)
    for (int x = 0; x < rows; x++)
      inputMatrix(x,y) = columns[y][x];

  VECTOR inputVector = EIGEN::convert(b);

  // write it all out
  string inputFile = string("LS.matrix");
  inputMatrix.write(inputFile.c_str());

  inputFile = string("LS.vector");
  inputVector.write(inputFile.c_str());

  // copy the script to the top directory
  string copy("cp ./src/matlab/fnnls*.m .");
  system(copy.c_str());

  // run matlab script
  string matlabCall = matlabBin + string(" -nodisplay -nosplash -r fnnlsFromFile");
  if (!verbose)
    matlabCall = matlabCall + string(" >& matlabSplash.txt");
  system(matlabCall.c_str());

  // read result in
  string outputFile = string("LS_result.vector");
  VECTOR result(outputFile.c_str());

  // clean up the files
  string rm("rm LS.matrix LS.vector LS_result.vector fnnls*.m matlabSplash.txt");
  system(rm.c_str());

  //cout << "NNLS result: " << result << endl;
  cout << " done. " << endl;

  return EIGEN::convert(result);
}

//////////////////////////////////////////////////////////////////////
// refine the key cell set
//////////////////////////////////////////////////////////////////////
void CUBATURE_GENERATOR_EIGEN::cullCubatureZeros(VectorXd& weights, vector<int>& candidates)
{
  int totalInitialCandidates = candidates.size();

  vector<int> refinedCandidates;
  vector<Real> refinedWeights;
  for (unsigned int x = 0; x < candidates.size(); x++)
  {
    if (fabs(weights[x]) < 1e-8) 
      continue;

    // keep this cell for the next fit
    refinedCandidates.push_back(candidates[x]);
    refinedWeights.push_back(weights[x]);
  }
  candidates = refinedCandidates;

  weights.resize(candidates.size());
  for (unsigned int x = 0; x < candidates.size(); x++)
    weights[x] = refinedWeights[x];

  cout << " By culling zeros, " << totalInitialCandidates << " points reduced to " << candidates.size() << " points" << endl;
}

//////////////////////////////////////////////////////////////////////
// Compute a cubature scheme using divide and conquer
//////////////////////////////////////////////////////////////////////
void CUBATURE_GENERATOR_EIGEN::generateImportanceSampledCubature()
{
  TIMER functionTimer(__FUNCTION__);

  // see if a reasonable number of key tets were requested
  if (_maxKeyCells > _fluid->size())
    cout << __FILE__ << " " << __LINE__ << " : " 
         << "You've asked for more key cells than there are in the volume! " 
         << "This will probably loop forever." << endl;

  if (_candidatesPerTry > _fluid->size())
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Too many candidate per try, " << _candidatesPerTry << "!" << endl;
    cout << " Fluid only contains " << _fluid->size() << " cells! " << endl;
    exit(0);
  }

  // size of a single advection vector
  int advectionSize = _trainingPostadvection[0].size();
  cout << "advectionSize is: " << advectionSize << endl;
  
  // total number of force samples
  int totalSamples = _trainingPostadvection.size();

  // number of rows in a column where all the samples are flattened
  // into a big vector
  const int totalRows = advectionSize * totalSamples;

  // calculate the squared magnitudes of the samples
  _postadvectionMagnitudes.resize(totalSamples);
  _postadvectionMagnitudes.setZero();

  for (int x = 0; x < totalSamples; x++)
    _postadvectionMagnitudes(x) = _trainingPostadvection[x].dot(_trainingPostadvection[x]);

  // calculate the inverse sqrt of the squared magnitudes
  for (int x = 0; x < totalSamples; x++)
    // CAREFUL WITH SINGLE PRECISION HERE
    if (fabs(_postadvectionMagnitudes(x)) > 1e-7)
       _postadvectionMagnitudes(x) = 1.0f / sqrtf(_postadvectionMagnitudes(x));
    else
       _postadvectionMagnitudes(x) = 1.0f;

  // stack normalized force samples into one big 'b' vector
  // ADJ: this could create potential memory problems if totalRows is too large
  // We also may need to declare the variable totalRows and indices into it as type *long*
  VectorXd b(totalRows);
  int index = 0;
  for (int x = 0; x < totalSamples; x++)
    for (int y = 0; y < advectionSize; y++, index++)
      b(index) = _trainingPostadvection[x][y] * _postadvectionMagnitudes(x);

  VectorXd residual = b;

  // error metrics
  //Real bNorm = b.norm();
  double relativeError = 1.0f;

  // generate a single candidate partition
  int totalCandidates = 5000;
  vector<int> candidates;
  VectorXd weights;
  int totalPasses = 6;
  map<int, bool> alreadyUsed;

  for (int x = 0; x < totalPasses; x++)
  {
    pickImportanceSampledCandidatesOMP(x, alreadyUsed, residual, candidates, totalCandidates);

    //totalCandidates = candidates.size();
    cout << " Using " << totalCandidates << " candidates " << endl;
    TIMER solverTimer("Solving partitions");

    // build the system for the fit
    cout << " Building " << candidates.size() << " columns ... " << flush;
    MatrixXd A(totalRows, candidates.size());

    // TODO: as each column is built independently, this could actually be parallelized
    for (unsigned int y = 0; y < candidates.size(); y++)
    {
      VectorXd column = getStagedTrainingColumn(candidates[y]);
      for (int i = 0; i < totalRows; i++)
        A(i,y) = column[i];
    }
    cout << " build done." << endl;

#if 1
    cout << " Solving NNLS (this would be a lot faster if you used Matlab) ... " << flush; 
    NNLS_SOLVER nnls;
    weights = nnls.solve(A, b);
    cout << " solve done." << endl;
#else
    cout << " Solving FNNLS ... " << flush; 
    weights = matlabFNNLS(A, b, false);
    cout << " solve done." << endl;
#endif
    
    // look at the error for this partition
    residual = b - A * weights;
    relativeError = residual.norm() / b.norm();
    cout << " Error after solve : " << relativeError << endl;
    solverTimer.stop();

    // now hack back the cubature set -- only cull non-zero entries
    // if we're looking at the final results
    cullCubatureZeros(weights, candidates);

    if (relativeError < 0.01)
      break;

    TIMER::printTimings();
  }

  _keyCells = candidates;
  _keyWeights.resize(candidates.size());

  for (unsigned int x = 0; x < candidates.size(); x++)
    _keyWeights[x] = weights[x];

  TIMER::printTimings();
}


//////////////////////////////////////////////////////////////////////
// Compute a cubature scheme using divide and conquer
//////////////////////////////////////////////////////////////////////
void CUBATURE_GENERATOR_EIGEN::generateImportanceSampledCubatureMemory()
{
  TIMER functionTimer(__FUNCTION__);

  // see if a reasonable number of key tets were requested
  if (_maxKeyCells > _fluid->size())
    cout << __FILE__ << " " << __LINE__ << " : " 
         << "You've asked for more key cells than there are in the volume! " 
         << "This will probably loop forever." << endl;

  if (_candidatesPerTry > _fluid->size())
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Too many candidate per try, " << _candidatesPerTry << "!" << endl;
    cout << " Fluid only contains " << _fluid->size() << " cells! " << endl;
    exit(0);
  }

  // size of a single advection vector
  int advectionSize = _trainingPostadvection[0].size();
  cout << "advectionSize is: " << advectionSize << endl;
  
  // total number of force samples
  int totalSamples = _trainingPostadvection.size();

  // number of rows in a column where all the samples are flattened
  // into a big vector
  const int totalRows = advectionSize * totalSamples;

  // calculate the squared magnitudes of the samples
  _postadvectionMagnitudes.resize(totalSamples);
  _postadvectionMagnitudes.setZero();

  for (int x = 0; x < totalSamples; x++)
    _postadvectionMagnitudes(x) = _trainingPostadvection[x].dot(_trainingPostadvection[x]);

  // calculate the inverse sqrt of the squared magnitudes
  for (int x = 0; x < totalSamples; x++)
    // CAREFUL WITH SINGLE PRECISION HERE
    if (fabs(_postadvectionMagnitudes(x)) > 1e-7)
       _postadvectionMagnitudes(x) = 1.0f / sqrtf(_postadvectionMagnitudes(x));
    else
       _postadvectionMagnitudes(x) = 1.0f;

  // stack normalized force samples into one big 'b' vector
  // ADJ: this could create potential memory problems if totalRows is too large
  // We also may need to declare the variable totalRows and indices into it as type *long*
  VectorXd b(totalRows);
  int index = 0;
  for (int x = 0; x < totalSamples; x++)
    for (int y = 0; y < advectionSize; y++, index++)
      b(index) = _trainingPostadvection[x][y] * _postadvectionMagnitudes(x);

  VectorXd residual = b;

  // error metrics
  //Real bNorm = b.norm();
  double relativeError = 1.0f;

  // generate a single candidate partition
  int totalCandidates = 5000;
  vector<int> candidates;
  VectorXd weights;
  int totalPasses = 6;
  map<int, bool> alreadyUsed;

  for (int x = 0; x < totalPasses; x++)
  {
    pickImportanceSampledCandidatesOMPMemory(x, alreadyUsed, residual, candidates, totalCandidates);

    //totalCandidates = candidates.size();
    cout << " Using " << totalCandidates << " candidates " << endl;
    TIMER solverTimer("Solving partitions");

    // build the system for the fit
    cout << " Building " << candidates.size() << " columns ... " << flush;
    MatrixXd A(totalRows, candidates.size());

    // TODO: as each column is built independently, this could actually be parallelized
    for (unsigned int y = 0; y < candidates.size(); y++)
    {
      VectorXd column = getStagedTrainingColumnMemory(candidates[y]);
      for (int i = 0; i < totalRows; i++)
        A(i,y) = column[i];
    }
    cout << " build done." << endl;

#if 1
    cout << " Solving NNLS (this would be a lot faster if you used Matlab) ... " << flush; 
    NNLS_SOLVER nnls;
    weights = nnls.solve(A, b);
    cout << " solve done." << endl;
#else
    cout << " Solving FNNLS ... " << flush; 
    weights = matlabFNNLS(A, b, false);
    cout << " solve done." << endl;
#endif
    
    // look at the error for this partition
    residual = b - A * weights;
    relativeError = residual.norm() / b.norm();
    cout << " Error after solve : " << relativeError << endl;
    solverTimer.stop();

    // now hack back the cubature set -- only cull non-zero entries
    // if we're looking at the final results
    cullCubatureZeros(weights, candidates);

    if (relativeError < 0.01)
      break;

    TIMER::printTimings();
  }

  _keyCells = candidates;
  _keyWeights.resize(candidates.size());

  for (unsigned int x = 0; x < candidates.size(); x++)
    _keyWeights[x] = weights[x];

  TIMER::printTimings();
}

void CUBATURE_GENERATOR_EIGEN::preprocess_prediffuse(const string& prediffusePath)
{
  _prediffuseFile = fopen(prediffusePath.c_str(), "rb");
  if (_prediffuseFile==NULL) {perror("Error opening file!"); exit(EXIT_FAILURE);}

  int rows, cols;
  fread(&rows, sizeof(int), 1, _prediffuseFile);
  fread(&cols, sizeof(int), 1, _prediffuseFile);
  printf("Parsing prediffuse matrix (rows, cols) as (%i, %i).\n", rows, cols);
  _prediffuseU_rows = rows;
  _prediffuseU_cols = cols;
}

