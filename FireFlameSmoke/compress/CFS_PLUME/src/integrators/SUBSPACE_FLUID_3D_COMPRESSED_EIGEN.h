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

//////////////////////////////////////////////////////////////////////
// SUBSPACE_FLUID_3D_COMPRESSED_EIGEN.h: interface for the SUBSPACE_FLUID_3D_COMPRESSED_EIGEN class.
//////////////////////////////////////////////////////////////////////

#ifndef SUBSPACE_FLUID_3D_COMPRESSED_EIGEN_H
#define SUBSPACE_FLUID_3D_COMPRESSED_EIGEN_H

// Eigen really doesn't like to be included after anything else
#include "EIGEN.h"

#include "FLUID_3D_MIC.h"
#include "SPARSE_MATRIX_ARRAY.h"
#include "COMPRESSION_DATA.h"
#include "MATRIX_COMPRESSION_DATA.h"
#include "COMPRESSION.h"

using namespace std;
using namespace Eigen;

class SUBSPACE_FLUID_3D_COMPRESSED_EIGEN : public FLUID_3D_MIC
{
public:
  // constructor for out-of-core SVD results
  SUBSPACE_FLUID_3D_COMPRESSED_EIGEN(int xRes, int yRes, int zRes, const string& reducedPath, unsigned int* boundaries = NULL, bool usingIOP = false, bool loadNothing = false);
  virtual ~SUBSPACE_FLUID_3D_COMPRESSED_EIGEN();

  const int xPeeled() const { return _xPeeled; };
  const int yPeeled() const { return _yPeeled; };
  const int zPeeled() const { return _zPeeled; };
  Real& discardThreshold() { return _discardThreshold; };
  unsigned int& domainBcFront()  { return _domainBcFront; };
  unsigned int& domainBcBack()   { return _domainBcBack; };
  unsigned int& domainBcLeft()   { return _domainBcLeft; };
  unsigned int& domainBcRight()  { return _domainBcRight; };
  unsigned int& domainBcTop()    { return _domainBcTop; };
  unsigned int& domainBcBottom() { return _domainBcBottom; };

  const VectorXd& l2Error() { return _l2Error; };

  void setCompressionPath(const string& path) { _compressionPath = path; };

  void stepReorderedCubatureStam();
  void stepPlume();
  void stepWithObstacle();
  void stepMovingObstacle(BOX* box);
  void stepMovingObstacleDebug(BOX* box);

  // const MatrixXd& U() const { return _U; };
  const MatrixXd& preadvectU() const { return _preadvectU; };

  MATRIX_COMPRESSION_DATA& U_final_data() { return _U_final_data; }
  MATRIX_COMPRESSION_DATA& U_preadvect_data() { return _U_preadvect_data; }

  const MatrixXd& prediffuseU() const { return _prediffuseU; };
  const MatrixXd& prevorticityU() const { return _prevorticityU; };
  string& fullRankPath() { return _fullRankPath; };
  vector<Real>& velocityErrorAbs() { return _velocityErrorAbs; };
  vector<Real>& densityErrorAbs() { return _densityErrorAbs; };
  vector<Real>& velocityErrorRelative() { return _velocityErrorRelative; };
  vector<Real>& densityErrorRelative() { return _densityErrorRelative; };

  // get the sub-basis of a cell from a specific basis  
  MatrixXd cellBasisPeeled(const MatrixXd& U, const int index);
  // get the sub-basis of a cell from a specific basis, compressed version
  MatrixXd cellBasisCompressedPeeled(MATRIX_COMPRESSION_DATA& U_data, const int index);
  // advect a single cell using Semi-Lagrangian,
  // assuming that "index" is a peeled index, not a full grid index


  VectorXd advectCellStamPeeled(MATRIX_COMPRESSION_DATA& U_data, const MatrixXd& cellU, 
      Real dt, const VectorXd& qDot, int index, MatrixXd* submatrix);

  void reducedSetMovingBox(BOX* box);

  // write the dims of the subspace error matrix to its file
  void writeCompressedErrorMatrixDims(int simulationSnapshots);

  // write out the most recent compressed subspace vector to the compressed subspace matrix file
  void appendCompressedSubspaceVectors();

private:
  struct CUBATURE_DATA {
    int index;
    int ixxx;
    Real wxxx;
  };

public:
  void accumAdvectRequests(const MatrixXd& cellU, 
      const Real dt, const VectorXd& qDot, const int cubatureIndex, const int index, 
      const VEC3I& dims, multimap<int, CUBATURE_DATA>& requestedBlocks);

  VectorXd advectCellStamNoProject(const Real& dt, const int index);

  // stomp the other matrices and load the ones needed for cubature training
  void loadCubatureTrainingBases();
  
  // stomp all loaded bases
  void stompAllBases();

  // initialize the matrix compression data
  void initCompressionData();

  // load the bases needed for cubature runtime
  void loadReducedRuntimeBases(string path = string(""));
  void loadReducedIOP(string path = string(""));
  void loadReducedIOPAll(string path = string(""));

  // load the qDotMatrix for error comparison
  void loadSubspaceComparisonMatrix();

  // compare the current compressed qDot against the uncompressed subspace reference
  void compareSubspace(int step);

  // read in a cubature scheme
  void readAdvectionCubature();

  // debugging compressed-space vs full space
  void diffTruth(const VECTOR3_FIELD_3D& testVelocity, const FIELD_3D& testDensity);

  // build matrices assuming that a limited number of matrices fit in memory
  // void buildOutOfCoreMatrices();

protected: 
  MatrixXd _U;

  BOX _box;
  int _totalReducedSteps;

  MATRIX_COMPRESSION_DATA _U_final_data;
  MATRIX_COMPRESSION_DATA _U_preadvect_data;
  MATRIX_COMPRESSION_DATA _U_preproject_data;

  // boundary bases -- ordering is:
  // 0 - left
  // 1 - right
  // 2 - bottom
  // 3 - top
  // 4 - near
  // 5 - far
  MatrixXd _Ub[6];

  // middle velocity vector
  VectorXd _qDot;
  
  // boundary velocity vectors
  // ordering is the same as Ub
  VectorXd _qbDot[6];

  // file paths
  string _reducedPath;
  string _compressionPath;
  
  // cache the damping matrix
  MatrixXd _dampingMatrixReduced;

  // build a peeled damping matrix -- boundaries have been peeled off
  SPARSE_MATRIX _peeledDampingMatrix;

  // build a peeled damping matrix for each boundary slab
  vector<SPARSE_MATRIX> _peeledBoundaryDampingMatrix;

  // build a matrix that translates the boundary damping values to the 
  // correct positions in the peeled core vector
  vector<SPARSE_MATRIX> _peeledBoundaryTranslationMatrix;
 
  // precompute once the product of _peeledBoundaryDampingMatrix and
  // _peeledBoundaryTranslationMatrix 
  vector<SPARSE_MATRIX> _peeledBoundaryMatrix;

  // reduced boundary damping matrices
  vector<MatrixXd> _boundaryDampingMatrixReduced;

  // matrix to derive divergence from velocity
  SPARSE_MATRIX _velocityToDivergence;

  // reduced version of _velocityToDivergence
  MatrixXd _reducedVelocityToDivergence;

  // matrix to project pressure out of velocity
  SPARSE_MATRIX _pressureToVelocity;

  // reduced version of _pressureToVelocity
  MatrixXd _reducedPressureToVelocity;

  // matrix to project the full IOP matrix into the subspace
  MatrixXd _projectionIOP;

  // precomputed _prejectionIOP' * _preprojectU for moving obstacle
  MatrixXd _projectionIOP_T_preprojectU;

  // the complement matrix to the neumann IOP matrix, flipping 0's and 1's on the
  // diagonal
  SPARSE_MATRIX _neumannIOPcomplement;

  // the homogeneous vector n that is appended as the last column of the IOP matrix
  VectorXd _neumannVector;

  // the reduced space version of _neumannVector
  VectorXd _reducedNeumannVector;

  // reduced matrix version of stomping the interior of an obstacle for IOP
  MatrixXd _reducedIOP;

  // a pressure basis -- _velocityToDivergence * _U
  MatrixXd _pressureU;

  // the Poisson matrix
  SPARSE_MATRIX _Asparse;

  // the projected Poisson matrix
  MatrixXd _reducedA;

  // factored projected Poisson matrix
  LDLT<MatrixXd> _factoredReducedA;

  // currently selected key cells
  vector<int> _keyAdvectionCells;
  vector<int> _keyVorticityCells;

  // key tet weights
  vector<Real> _keyAdvectionWeights;
  vector<Real> _keyVorticityWeights;

  // full rank snapshots path
  string _fullRankPath;

  // keep the peeled dimensions around
  int _xPeeled;
  int _yPeeled;
  int _zPeeled;
  int _slabPeeled;

  // threshold at which to start discarding PCA columns
  Real _discardThreshold;

  // for output purposes, keep track of the difference against ground truth
  vector<Real> _velocityErrorAbs;
  vector<Real> _velocityErrorRelative;
  vector<Real> _densityErrorAbs;
  vector<Real> _densityErrorRelative;

  // all the different bases -- eventually these will never get loaded, just created
  // once to project the matrix
  MatrixXd _preprojectU;
  // TODO: use compressed version here?
  MatrixXd _preadvectU;
  MatrixXd _prediffuseU;
  MatrixXd _prevorticityU;

  // projection step matrices
  MatrixXd _preprojectToPreadvect;
  MatrixXd _preprojectToFinal;
  MatrixXd _inverseProduct;

  // subspace error matrix to compare against
  MatrixXd _qDotMatrix;
  // list of the relative l2 errors against _qDotMatrix
  VectorXd _l2Error;

  // domain boundary conditions
  unsigned int _domainBcFront;
  unsigned int _domainBcBack;
  unsigned int _domainBcBottom;
  unsigned int _domainBcTop;
  unsigned int _domainBcLeft;
  unsigned int _domainBcRight;

  // cached advection cubature matrices
  vector<MatrixXd> _advectionCubatureBefore;
  vector<MatrixXd> _advectionCubatureAfter;
  
  // cached vorticity cubature matrices
  vector<MatrixXd> _vorticityCubatureBefore;
  vector<MatrixXd> _vorticityCubatureAfter;

  // a basis for after the IOP projection
  MatrixXd _iopU;

  // are we using IOP for boundaries?
  bool _usingIOP;

  // projected Neumann IOP matrices
  vector<MatrixXd> _projectedNeumann;

  // perform reduced order diffusion with separate boundary slabs
  void reducedPeeledDiffusion();

  
  // initialize the staged version, assuming the Us were computed out of core and will not
  // all fit in memory
  void initOutOfCore();
 
  // initialize the staged version for IOP 
  void initOutOfCoreIOP();
  
  // build a peeled version of the damping matrix
  void buildPeeledDampingMatrix();
  void buildPeeledDampingMatrixFlat(SPARSE_MATRIX_ARRAY& final);
  
  // get the projection error of a vector with respect to a basis
  Real projectionError(const MatrixXd& basis, const VectorXd& v);

  // compute the velocity-to-divergence matrix
  void computeVelocityToDivergence();
 
  // compute pressure-to-velocity matrix
  void computePressureToVelocity();

  // build a sparse version of the Poisson matrix
  void buildFlatA(SPARSE_MATRIX_ARRAY& sparseA, unsigned char* skip);

  // do reduced stomping of the obstacle interior
  void reducedSetZeroSphere();

  // do a staged reduced order pressure projection
  void reducedStagedProject();
  
  // do a staged reduced order pressure projection for IOP
  void reducedStagedProjectIOP();

  // do a full-rank advection of heat and density using semi-Lagrangian
  void advectHeatAndDensityStam();
  
  // diff the current sim results against ground truth
  void diffGroundTruth();
  
  // do Stam-style adveiction using cubautre
  void reducedAdvectStagedStamFast();

  // do Stam-style adveiction using cubautre
  void reducedAdvectCompressionFriendly();

  // check of a file exists
  bool fileExists(const string& filename);

  // purge inactive memory on OSX
  void purge() { cout << " Purging inactive OSX memory" << endl; system("purge"); };

  // reduced order IOP -- both orthogonal and pressure projection
  void reducedIOP();
  
  // add a new orthogonalized column to the basis
  void addNewColumn(const VectorXd& newColumn, MatrixXd& U);


};

#endif

