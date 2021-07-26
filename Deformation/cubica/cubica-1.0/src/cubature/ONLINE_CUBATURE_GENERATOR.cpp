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
// ONLINE_CUBATURE_GENERATOR.h: interface for the ONLINE_CUBATURE_GENERATOR class.
//
//////////////////////////////////////////////////////////////////////

#include "ONLINE_CUBATURE_GENERATOR.h"
#ifdef USING_OPENMP
#include <omp.h>
#endif

//////////////////////////////////////////////////////////////////////
// Constructor for tetrahedra
//////////////////////////////////////////////////////////////////////
ONLINE_CUBATURE_GENERATOR::ONLINE_CUBATURE_GENERATOR(SUBSPACE_TET_MESH* tetMesh) :
  // init Mersenne Twister with deterministic seed
  // so that runs are reproducible
  _tetMesh(tetMesh),
  _twister(123456),
  _tolerance(1e-7),
  _maxKeyTets(1000),
  _candidatesPerTry(100),
  _errorTolerance(0.01f),
  _relativeError(0.0),
  _skeleton(NULL),
  _nnlsMaxIter(1),
  _bailed(false)
{
}

//////////////////////////////////////////////////////////////////////
// pick candidates for key tet selection
//////////////////////////////////////////////////////////////////////
vector<TET*> ONLINE_CUBATURE_GENERATOR::pickCandidates(int _candidatesPerTry)
{
  vector<TET>& tets = _tetMesh->tets();
  int totalTets = _tetMesh->totalTets();

  // build a lookup table for the already selected key tets
  map<TET*, bool> alreadyUsed;
  for (unsigned int x = 0; x < _keyTets.size(); x++)
    alreadyUsed[_keyTets[x]] = true;

  // construct the final vector of candidate tets
  vector<TET*> candidates;
  while (candidates.size() != (unsigned int)_candidatesPerTry)
  {
    // pick a random tet
    int index = _twister.randInt(totalTets - 1);

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
// pick all candidates for key tet selection
//////////////////////////////////////////////////////////////////////
vector<TET*> ONLINE_CUBATURE_GENERATOR::pickCandidatesAll()
{
  vector<TET>& tets = _tetMesh->tets();
  int totalTets = _tetMesh->totalTets();

  // build a lookup table for the already selected key tets
  map<TET*, bool> alreadyUsed;
  for (unsigned int x = 0; x < _keyTets.size(); x++)
    alreadyUsed[_keyTets[x]] = true;

  // construct a list of all non-key tets
  vector<TET*> candidates;
  for (int x = 0; x < totalTets; x++)
  {
    // check that it is not already a key tet
    if (alreadyUsed.find(&(tets[x])) == alreadyUsed.end())
      candidates.push_back(&(tets[x]));
  }

  return candidates;
}

//////////////////////////////////////////////////////////////////////
// Retrieve the explicitly stored training column
//////////////////////////////////////////////////////////////////////
VECTOR ONLINE_CUBATURE_GENERATOR::retrieveTrainingColumn(TET* tet)
{
  int sampleSize = _trainingForces[0].size();
  int totalRows = _trainingForces.size() * sampleSize;
  VECTOR gColumn(totalRows);

  // submatrix of _Ubasis corresponding to this tet
  MATRIX tetU = _tetMesh->tetSubBasis(tet);

  // generate column corresponding to stacked force samples
  for (unsigned int x = 0; x < _trainingForces.size(); x++)
  {
    // get the density vector
    VECTOR& forceColumn = _trainingColumns[x];

    // get the tet index in the vector
    int tetID = _tetMesh->tetID(tet);

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
      gColumn(x * sampleSize + y) = forceReduced(y) * _forceMagnitudes(x);
  }

  return gColumn;
}

//////////////////////////////////////////////////////////////////////
// Generate the training column corresponding to a tet
//////////////////////////////////////////////////////////////////////
VECTOR ONLINE_CUBATURE_GENERATOR::getTrainingColumn(TET* tet)
{
  // if the training forces were explicitly stored, just return those
  if (_trainingColumns.size() > 0)
    return retrieveTrainingColumn(tet);

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
    {
      int vertexID = _tetMesh->vertexID(tet->vertices[y]);
      vertices[y] = *(_tetMesh->restVertices(vertexID));
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
      if (_tetMesh->isConstrained(tet->vertices[y]))
      {
        int vertexID = _tetMesh->vertexID(tet->vertices[y]);
        vertexID -= _tetMesh->unconstrainedNodes();

        vertices[y][0] += _trainingConstraints[x][3 * vertexID];
        vertices[y][1] += _trainingConstraints[x][3 * vertexID + 1];
        vertices[y][2] += _trainingConstraints[x][3 * vertexID + 2];
      }

    // if there is a skinning, we must fix up the skinning matrices for the tet
    MATRIX skinningMatrices[4];
    if (_skeleton)
    {
      for (int y = 0; y < 4; y++)
      {
        int vertexID = _tetMesh->vertexID(tet->vertices[y]);
        skinningMatrices[y] = _skeleton->recreateHeadSkinning(vertexID, _snapshotTimes[x]);
        newTet.skinningMatrix(y) = &(skinningMatrices[y]);
      }
    }

    // if there is a skinning, we must fix up the deformed to rest hash
    map<VEC3F*, VEC3F*> deformedToRest;
    if (_skeleton)
      for (int y = 0; y < 4; y++)
      {
        deformedToRest[newTet.vertices[y]] = (*(tet->deformedToRest))[tet->vertices[y]];
      }
    newTet.deformedToRest = &deformedToRest;

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
int ONLINE_CUBATURE_GENERATOR::findBestTet(vector<TET*>& candidates, VECTOR& residual)
{
  Real rNorm = residual.norm2();

  // candidate dot products;
  vector<Real> dots;
  dots.resize(candidates.size());

#if USING_OPENMP  
  int size = candidates.size();
#pragma omp parallel 
  {
	  const int num = omp_get_num_threads();
    const int id = omp_get_thread_num();
    int begin = (size / num) * id;
    int end   = (id == num-1) ? size : (size / num) * (id+1);
    for (int x = begin; x < end; x++)
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
  } 
#else
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
#endif

  int maxDot = 0;
  Real maxFound = dots[maxDot];
  for (unsigned int x = 0; x < candidates.size(); x++)
  {
    if (dots[x] > dots[maxDot])
      maxFound = dots[x];
    
    maxDot = (dots[x] > dots[maxDot]) ? x : maxDot;
  }
  return maxDot;
}

//////////////////////////////////////////////////////////////////////
// Generate the cubature points
//////////////////////////////////////////////////////////////////////
void ONLINE_CUBATURE_GENERATOR::generateCubature(map<string, double>& timings)
{
  TIMER preamble;
  if (_trainingForces.size() == 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << "No training forces found! Did you call ONLINE_SUBSPACE_INTEGRATOR::populateTrainingSet()?" << endl;
  }

  cout << "=============================================" << endl;
  cout << " Creating cubature" << endl;
  cout << "=============================================" << endl;
  
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
    _forceMagnitudes(x) = _trainingForces[x] * _trainingForces[x];

  // calculate the inverse sqrt of the squared magnitudes
  for (int x = 0; x < totalSamples; x++)
  {
    // CAREFUL WITH SINGLE PRECISION HERE
    if (fabs(_forceMagnitudes(x)) > 1e-7)
       _forceMagnitudes(x) = 1.0 / sqrt(_forceMagnitudes(x));
    else
       _forceMagnitudes(x) = 1.0;
  }

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
  timings["Cub. Preamble"] += preamble.timing();

  TIMER previousTimer;
  // incorporate previous key tets
  int totalKeyTets = _keyTets.size();
  for (int x = 0; x < totalKeyTets; x++)
  {
    // get training column for this tet
    VECTOR bestTetColumn = getTrainingColumn(_keyTets[x]);

    // copy this column into A
    for (int y = 0; y < totalRows; y++)
      A[x * totalRows + y] = bestTetColumn(y);
  }
  timings["Cub. Prev. Columns"] += previousTimer.timing();

  // see how many times the last mod 100 coefficient was zero
  int zeroCoefficient = 0;

  // we didn't bail until we find out otherwise
  _bailed = false;

  // start training
  int totalIterations = totalKeyTets;
  while (relativeError > _errorTolerance && 
         totalKeyTets < _maxKeyTets &&
         totalKeyTets < _tetMesh->totalTets())
  {
    // get new set of randomly selected candidates
    if (_tetMesh->totalTets() - totalKeyTets < _candidatesPerTry)
      _candidatesPerTry = _tetMesh->totalTets() - totalKeyTets;

    TIMER pickTimer;
    vector<TET*> candidates = pickCandidates(_candidatesPerTry);
    timings["Cub. Pick Candidates"] += pickTimer.timing();

    TIMER bestTetTimer;
    // pick the best tet out of the list of candidates
    int bestTet = findBestTet(candidates, residual);
    timings["Cub. Best Tet"] += bestTetTimer.timing();

    // store the new key tet
    _keyTets.push_back(candidates[bestTet]);
    totalKeyTets++;

    TIMER columnTimer;
    // recompute its column (better than storing them all)
    VECTOR bestTetColumn = getTrainingColumn(candidates[bestTet]);
    timings["Cub. Compute Column"] += columnTimer.timing();

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
    TIMER nnlsTimer;
    nnls.solve(A, totalKeyTets, bNNLS, weightsNNLS, rNorm);
    timings["Cub. NNLS"] += nnlsTimer.timing();

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

    if (weightsNNLS[totalKeyTets-1] <= 0.0)
      zeroCoefficient++;
    else
      zeroCoefficient = 0;

    // bail if it happened five times in a row
    if (zeroCoefficient == 100)
    {
      cout << " Too many zero coefficients found! Bailing ... " << endl;
      break;
    }

    cout << "Key Tet " << totalIterations << ": \tWeight: " << weightsNNLS[totalKeyTets-1] 
         << "\tError: " << relativeError << endl;
  }
  _priorWeights.push_back(0.0);

  // copy the key tet weights to single precision
  _keyWeights.clear();
  for (int x = 0; x < totalKeyTets; x++)
    _keyWeights.push_back(weightsNNLS[x]);

  // delete the key tets with weight 0
  cout << " Total key tets before filtering: " << _keyTets.size() << endl;
  vector<int> validTets = filterKeyTets();
  cout << " Total key tets after filtering: " << _keyTets.size() << endl;
  _totalTetsAfterUpdates.push_back(_keyTets.size());
  cout << " Relative error: " << relativeError << endl;
  _relativeError = relativeError;

  if (relativeError > _errorTolerance)
    _bailed = true;

  // clean up scratch space
  delete[] bNNLS;
  delete[] weightsNNLS;
  delete[] A;
  cout << "=============================================" << endl;
}

//////////////////////////////////////////////////////////////////////
// train the Gauss-Kronrod pairs
//////////////////////////////////////////////////////////////////////
void ONLINE_CUBATURE_GENERATOR::trainGaussKronrod(double* A, VECTOR& b, double* weightsNNLS, NNLS_SOLVER& nnls, vector<int>& validTets)
{
  cout << __FILE__ << " " << __LINE__ << " : " << endl;
  cout << " Gauss-Kronrod training disabled! " << endl;
}

//////////////////////////////////////////////////////////////////////
// Update the cubature points
// all priors are assigned a single weight version
//////////////////////////////////////////////////////////////////////
void ONLINE_CUBATURE_GENERATOR::updateCubature(map<string, double>& timings)
{
  if (_trainingForces.size() == 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << "No training forces found! Did you call ONLINE_SUBSPACE_INTEGRATOR::populateTrainingSet()?" << endl;
  }

  // see if a reasonable number of key tets were requested
  if (_maxKeyTets > _tetMesh->totalTets())
    cout << __FILE__ << " " << __LINE__ << " : " 
         << "You've asked for more key tets than there are in the mesh! " 
         << "This will probably loop forever." << endl;

  cout << "=============================================" << endl;
  cout << " Updating cubature" << endl;
  cout << "=============================================" << endl;

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
    _forceMagnitudes(x) = _trainingForces[x] * _trainingForces[x];

  // calculate the inverse sqrt of the squared magnitudes
  for (int x = 0; x < totalSamples; x++)
  {
    // CAREFUL WITH SINGLE PRECISION HERE
    if (fabs(_forceMagnitudes(x)) > 1e-7)
       _forceMagnitudes(x) = 1.0 / sqrt(_forceMagnitudes(x));
    else
       _forceMagnitudes(x) = 1.0;
  }

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

  // scratch space for the NNLS solver
  double* bNNLS = new double[totalRows];
  double* weightsNNLS = new double[_maxKeyTets];
  for (int x = 0; x < _maxKeyTets; x++)
    weightsNNLS[x] = 0.0;

  // error metrics
  Real bNorm = b.norm2();
  double rNorm = bNorm;
  double relativeError = 1.0;

  // incorporate previous key tets
  int totalKeyTets = _keyTets.size();
  int newKeyTets = 1;

  // lump all the priors into one column
  VECTOR priorColumn(totalRows);
  for (int x = 0; x < totalKeyTets; x++)
  {
    VECTOR column = getTrainingColumn(_keyTets[x]);
    priorColumn += _keyWeights[x] * column;
  }

  // copy this column into the first column of A
  for (int x = 0; x < totalRows; x++)
    A[x] = priorColumn(x);

  // start training
  int totalIterations = totalKeyTets;
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
    newKeyTets++;

    // recompute its column (better than storing them all)
    VECTOR bestTetColumn = getTrainingColumn(candidates[bestTet]);

    // copy the column to "A" matrix --
    // note that A is column major, so it appears to be copying to a row
    //
    // Skip the first column, since that is where all the old tets have
    // been piled
    for (int x = 0; x < totalRows; x++)
      A[newKeyTets * totalRows + x] = bestTetColumn(x);

    // copy b to double precision, also since the NNLS_SOLVER is going
    // to clobber it
    for (int x = 0; x < totalRows; x++)
      bNNLS[x] = b(x);
   
    // do the NNLS fit - note that both weightsNNLS and rNorm
    // both return values
    bool converged = nnls.solve(A, newKeyTets, bNNLS, weightsNNLS, rNorm);

    // update residual if the fit converged
    if (converged)
    {
      relativeError = rNorm / bNorm;

      // calculate the residual by doing an explicit A*x, 
      // but since A is column major, it is really (A^T) * x
      residual = b;
      for (int i = 0; i < newKeyTets; i++)
        for (int j = 0; j < totalRows; j++)
        {
          int index = i * totalRows + j;
          residual(j) -= weightsNNLS[i] * A[index];
        }
    }
    if (totalIterations != 0 && totalIterations % 100 == 0)
      cout << "Key tet " << _keyTets.size() << ": \t Priors weight: " << weightsNNLS[0] << "\tWeight: " << weightsNNLS[newKeyTets-1] 
           << "\tError: " << relativeError << endl;

    totalIterations++;
  }
  cout << " Priors weighting: " << weightsNNLS[0] << endl;
  _priorWeights.push_back(weightsNNLS[0]);

  // multiply the old weights by the new constant
  for (unsigned int x = 0; x < _keyWeights.size(); x++)
    _keyWeights[x] *= weightsNNLS[0];

  // record all the new tets
  for (int x = 1; x < newKeyTets; x++)
    _keyWeights.push_back(weightsNNLS[x]);

  // delete the key tets with weight 0
  vector<int> validTets = filterKeyTets();
  cout << " Total key tets: " << _keyTets.size() << endl;
  _totalTetsAfterUpdates.push_back(_keyTets.size());

  // train the Gauss-Kronrod rule
  trainGaussKronrod(A, b, weightsNNLS, nnls, validTets);

  // clean up scratch space
  delete[] bNNLS;
  delete[] weightsNNLS;
  delete[] A;
  
  cout << "=============================================" << endl;
}

//////////////////////////////////////////////////////////////////////
// Filter out the dud tets
//////////////////////////////////////////////////////////////////////
vector<int> ONLINE_CUBATURE_GENERATOR::filterKeyTets()
{
  int beforeFiltering = _keyTets.size();
  
  vector<TET*> filteredTets;
  vector<Real> filteredWeights;
  vector<int> valid;

  // filter the tets
  for (unsigned int x = 0; x < _keyTets.size(); x++)
    if (_keyWeights[x] > 1e-12)
    {
      filteredTets.push_back(_keyTets[x]);
      filteredWeights.push_back(_keyWeights[x]);
      valid.push_back(x);
    }

  // copy them all back
  _keyTets.resize(filteredTets.size());
  _keyWeights.resize(filteredTets.size());
  for (unsigned int x = 0; x < filteredTets.size(); x++)
  {
    _keyTets[x] = filteredTets[x];
    _keyWeights[x] = filteredWeights[x];
  }

  // record how many duds we found
  _totalDuds.push_back(beforeFiltering - _keyTets.size());

  return valid;
}

//////////////////////////////////////////////////////////////////////
// Reset everything
//////////////////////////////////////////////////////////////////////
void ONLINE_CUBATURE_GENERATOR::reset()
{
  _trainingForces.resize(0);
  _trainingQs.resize(0);
  _trainingConstraints.resize(0);
  _forceMagnitudes.resizeAndWipe(0);
  
  _keyTets.resize(0);
  _keyWeights.resize(0);
  _kronrodWeights.resize(0);
  _totalTetsAfterUpdates.resize(0);
  _totalDuds.resize(0);
  _priorWeights.resize(0);
}

//////////////////////////////////////////////////////////////////////
// Write out the state
//////////////////////////////////////////////////////////////////////
void ONLINE_CUBATURE_GENERATOR::writeState(string filename)
{
  cout << " Dumping state to file " << filename.c_str() << endl;
  FILE* file = fopen(filename.c_str(), "wb");
  int size;

  size = _trainingForces.size();
  fwrite((void*)&size, sizeof(int), 1, file);
  for (unsigned int x = 0; x < _trainingForces.size(); x++)
    _trainingForces[x].write(file);

  size = _trainingQs.size();
  fwrite((void*)&size, sizeof(int), 1, file);
  for (unsigned int x = 0; x < _trainingQs.size(); x++)
    _trainingQs[x].write(file);

  size = _trainingConstraints.size();
  fwrite((void*)&size, sizeof(int), 1, file);
  for (unsigned int x = 0; x < _trainingConstraints.size(); x++)
    _trainingConstraints[x].write(file);

  size = _snapshotTimes.size();
  fwrite((void*)&size, sizeof(int), 1, file);
  for (unsigned int x = 0; x < _snapshotTimes.size(); x++)
  {
    Real time = _snapshotTimes[x];
    fwrite((void*)&time, sizeof(Real), 1, file);
  }

  // previous key tets?
  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// Read in the state
//////////////////////////////////////////////////////////////////////
void ONLINE_CUBATURE_GENERATOR::readState(string filename)
{
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " NOT IMPLEMENTED " << endl;
}
