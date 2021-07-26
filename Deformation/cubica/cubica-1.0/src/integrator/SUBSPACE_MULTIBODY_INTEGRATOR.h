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
// SUBSPACE_MULTIBODY_INTEGRATOR.h: interface for the SUBSPACE_MULTIBODY_INTEGRATOR class.
//
//////////////////////////////////////////////////////////////////////

#ifndef SUBSPACE_MULTIBODY_INTEGRATOR_H
#define SUBSPACE_MULTIBODY_INTEGRATOR_H

#include <PARTITIONED_SUBSPACE_TET_MESH.h>
#include <SUBSPACE_INTEGRATOR.h>
#include <string>
#include <MERSENNETWISTER.h>
#include <SANDWICH_TRANSFORM.h>
#include <UNCONSTRAINED_SUBSPACE_INTEGRATOR.h>
#include <FULLSPACE_INTEGRATOR.h>
#include <TIMER.h>
#include <BLOCK_VECTOR.h>
#include <GNUPLOT.h>
#include <SPARSE_PETSC_MATRIX.h>

#define RECENTER_ROTATION_FORCES 1
#define USING_IMPLICIT_QUADRATICS 0
#define USING_EXPLICIT_QUADRATICS 0
#define USING_IMPROVED_DAMPING 1
#define USING_IMPROVED_DAMPING_RESIDUAL 1

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
class SUBSPACE_MULTIBODY_INTEGRATOR {

public:
  SUBSPACE_MULTIBODY_INTEGRATOR(PARTITIONED_SUBSPACE_TET_MESH* tetMesh, Real dt = 1.0 / 60.0, Real springConst = 10.0, Real dampingConst = 0.01, Real alpha = 0.01, Real beta = 0.01);
  virtual ~SUBSPACE_MULTIBODY_INTEGRATOR();

  void stepReducedSkinnedDynamic(bool verbose = true);
  void stepReducedSkinnedDynamicOMP(bool verbose = true);
  void stepReducedSkinnedDynamicWithCollisionsOMP(vector<vector<pair<SURFACE*, int> > >& collisions, bool verbose);

  // mouse support
  void click(VEC3F& point, int forcePartition = -1);
  void drag(VEC3F& point);
  void unclick();
  void drawForceVector();

  // accessors
  SUBSPACE_INTEGRATOR* integrator(int x) { return _integrators[x]; };
  int partitions() { return _partitions; };
  const int totalSteps() { return _totalSteps; };
  const Real dt() { return _dt; };
  Real& solverEps() { return _solverEps; };
  const int totalNewtonStepsSeen() { return _totalNewtonStepsSeen; };
  const int maxNewtonStepsSeen() { return _maxNewtonStepsSeen; };
  int& maxNewtonSteps() { return _maxNewtonSteps; };
  int& newtonStalls() { return _newtonStalls; };
  Real& forceMultiplier() { return _forceMultiplier; };
  Real& totalTime() { return _totalTime; };
  
  void addGravity();
  void changeGravity(VEC3F gravityDown, Real gravityMagnitude);
  void changeGravity(int partition, VEC3F gravityDown, Real gravityMagnitude);
  void addBodyForce(VEC3F bodyForce);

  // dump out timing info
  void printTimingBreakdown();

  // reset simulation
  void reset();

  // compute ODE versions of the sandwiches -- the rest poses change
  void precomputeOdeSandwiches(string dataPath);

protected:
  SUBSPACE_INTEGRATOR** _integrators;
  PARTITIONED_SUBSPACE_TET_MESH* _partitionedMesh;

  // total number of blocks, taking into account rigid components
  int _totalBlocks;
  int _partitions;
  Real _dt;
  int _totalSteps;
  Real _solverEps;
  Real _rayleighAlpha;
  Real _rayleighBeta;
  map<string, double> _timingBreakdown;
  Real _totalTime;

  // Newmark consts
  Real _beta;
  Real _gamma;

  // penalty spring consts
  Real _springConst;
  Real _dampingConst;

  // force multiplier for all partitions
  Real _forceMultiplier;

  // per-interface spring consts
  map<pair<int, int>, Real> _interfaceSpringConsts;

  // maximum Newton steps allowed in a single solve
  int _maxNewtonSteps;

  // total Newton steps seen over entire lifetime
  int _totalNewtonStepsSeen;

  // maximum Newton steps seen over entire lifetime
  int _maxNewtonStepsSeen;

  // total number of Newton stalls seen
  int _newtonStalls;

  // i->x translation
  vector<int> _unconstrainedIDs;

  // x->i translation
  map<int, int> _inverseUnconstrainedIDs;

  // force vectors, including rigid terms
  vector<BLOCK_VECTOR> _forceVectors;

  // mouse support
  int _clickedPartition;

  // block number that partitions begin at
  vector<int> _startingBlock;

  // block system to solve
  BLOCK_MATRIX _blockA;
  BLOCK_VECTOR _blockB;

  // thread specific copies of _blockA
  vector<BLOCK_MATRIX> _blockAs;

  //////////////////////////////////////////////////////////////////////
  // Fast transform sandwiches 
  //
  // The naming convention is <left matrix>_<T or R>_<right matrix>
  // where <left matrix> and <right matrix> try to describe the 
  // constant matrix portions of the sandwich somehow, and T and R
  // denote whether the transform expects a regular rotation (R)
  // or a tensor (T).
  //////////////////////////////////////////////////////////////////////
  // used to compute translation partial angular
  map<pair<int, int>, SANDWICH_TRANSFORM> _I3x3n_T_restU;
  // used to compute translation partial defo and translation partial angular
  map<pair<int, int>, SANDWICH_TRANSFORM> _I3x3n_R_Uij;
  map<pair<int, int>, SANDWICH_TRANSFORM> _I3x3n_R_Ukj;

  // used to compute angular partial translation, angular, and defo
  map<pair<int, int>, SANDWICH_TRANSFORM> _tildeUo_R_I3nx3;
  map<pair<int, int>, SANDWICH_TRANSFORM> _tildeUo_R_uiBarO;
  map<pair<int, int>, SANDWICH_TRANSFORM> _tildeUo_R_Uij;
  map<pair<int, int>, SANDWICH_TRANSFORM> _tildeUo_R_ukBarO;
  map<pair<int, int>, SANDWICH_TRANSFORM> _tildeUo_R_Ukj;
  map<pair<int, int>, SANDWICH_TRANSFORM> _I3x3n_R_tildeUo;
  map<pair<int, int>, SANDWICH_TRANSFORM> _I3x3n_R_uiBarO;
  map<pair<int, int>, SANDWICH_TRANSFORM> _I3x3n_R_ukBarO;
  map<pair<int, int>, TENSOR3> _UijTildeU;
  map<pair<int, int>, MATRIX> _tildeUoTUij;
  map<pair<int, int>, VECTOR> _tildeUo_uiBarO;
  map<pair<int, int>, MATRIX> _tildeUTuiBarO;
  map<pair<int, int>, TENSOR3> _tildeUTUij;

  map<pair<int, int>, SANDWICH_TRANSFORM> _Ui_R_I3nx3;
  map<pair<int, int>, SANDWICH_TRANSFORM> _Ui_R_uiBarO;
  map<pair<int, int>, SANDWICH_TRANSFORM> _Ui_R_Ui;
  map<pair<int, int>, SANDWICH_TRANSFORM> _Ui_R_ukBarO;
  map<pair<int, int>, SANDWICH_TRANSFORM> _Ui_R_Uk;
  map<int, MATRIX> _Ui_Ui;
  map<int, MATRIX> _spring_Ui_Ui;
  map<pair<int, int>, MATRIX> _Ui_I3nx3;
  map<pair<int, int>, MATRIX> _I3nx3_Uk;
  map<pair<int, int>, TENSOR3> _tildeU_Ui;

  // constrained only -- sandwich reduced to a matrix
  map<pair<int, int>, VECTOR> _Ui_uiBarO;
  map<pair<int, int>, VECTOR> _Ui_ukBarO;
  map<pair<int, int>, MATRIX> _Ui_Uk;

  // workspace tensors
  map<pair<int, int>, TENSOR3> _rightProduct;
  map<pair<int, int>, TENSOR3> _workspace3xRxL;
  map<pair<int, int>, TENSOR3> _workspace3xRx3;
  map<pair<int, int>, TENSOR3> _workspaceLxRx3;
  map<pair<int, int>, TENSOR3> _workspaceRxLx3;
  map<pair<int, int>, MATRIX> _workspace3xR;
  map<pair<int, int>, MATRIX> _workspaceRx3;
  map<pair<int, int>, MATRIX> _workspaceLxR;
  map<pair<int, int>, TENSOR3> _workspaceRx1x3;

  map<int, TENSOR3> _workspace3xLx3;
  map<int, TENSOR3> _workspace3x3xL;
  map<int, TENSOR3> _workspace3x1xL;
  map<int, TENSOR3> _workspace3x1x3;
  map<int, TENSOR3> _workspaceLx3x3;
  map<int, TENSOR3> _workspaceLx1x3;
  map<int, TENSOR3> _workspaceLxLx3;
  map<int, TENSOR3> _workspace3xLxL;
  map<int, MATRIX> _workspace3x3;
  map<int, MATRIX> _workspace3xL;
  map<int, MATRIX> _workspaceLx3;
 
  // sandwiches for external force computation
  map<int, SANDWICH_TRANSFORM> _fullTildeUo_Ui;

  // presummed sandwiches
  map<int, SANDWICH_TRANSFORM> _summed_I3x3n_R_tildeUo;
  map<int, SANDWICH_TRANSFORM> _summed_tildeUo_R_I3nx3;
  map<int, SANDWICH_TRANSFORM> _summed_tildeUo_R_uiBarO;
  map<int, SANDWICH_TRANSFORM> _summed_tildeUo_R_Uij;
  map<int, MATRIX> _summed_tildeUoTUij;
  map<int, TENSOR3> _summed_UijTildeU;
  map<int, MATRIX> _summed_tildeUTuiBarO;
  map<int, TENSOR3> _summed_tildeUTUij;
  map<int, SANDWICH_TRANSFORM> _summed_I3x3n_T_restU;
  map<int, SANDWICH_TRANSFORM> _summed_I3x3n_R_Uij;
  map<int, SANDWICH_TRANSFORM> _summed_Ui_R_I3nx3;
  map<int, SANDWICH_TRANSFORM> _summed_Ui_R_Ui;
  map<int, SANDWICH_TRANSFORM> _summed_Ui_R_uiBarO;
  map<int, SANDWICH_TRANSFORM> _summed_I3x3n_R_uiBarO;
  map<int, VECTOR> _summed_Ui_uiBarO;
  map<int, VECTOR> _summed_tildeUo_uiBarO;
  map<int, TENSOR3> _summed_tildeU_Ui;

  // spring const weighted presummed sandwiches
  map<int, SANDWICH_TRANSFORM> _spring_summed_I3x3n_R_tildeUo;
  map<int, SANDWICH_TRANSFORM> _spring_summed_tildeUo_R_I3nx3;
  map<int, SANDWICH_TRANSFORM> _spring_summed_tildeUo_R_uiBarO;
  map<int, SANDWICH_TRANSFORM> _spring_summed_tildeUo_R_Uij;
  map<int, MATRIX> _spring_summed_tildeUoTUij;
  map<int, TENSOR3> _spring_summed_UijTildeU;
  map<int, MATRIX> _spring_summed_tildeUTuiBarO;
  map<int, TENSOR3> _spring_summed_tildeUTUij;
  map<int, SANDWICH_TRANSFORM> _spring_summed_I3x3n_T_restU;
  map<int, SANDWICH_TRANSFORM> _spring_summed_I3x3n_R_Uij;
  map<int, SANDWICH_TRANSFORM> _spring_summed_Ui_R_I3nx3;
  map<int, SANDWICH_TRANSFORM> _spring_summed_Ui_R_Ui;
  map<int, SANDWICH_TRANSFORM> _spring_summed_Ui_R_uiBarO;
  map<int, SANDWICH_TRANSFORM> _spring_summed_I3x3n_R_uiBarO;
  map<int, VECTOR> _spring_summed_Ui_uiBarO;
  map<int, VECTOR> _spring_summed_tildeUo_uiBarO;
  map<int, TENSOR3> _spring_summed_tildeU_Ui;

  vector<Real> _clonesTimesSpringConsts;

  // compute the non-linear mass matrix
  void computeMassMatrix(int partition, BLOCK_MATRIX& massMatrix);

  // compute the external force vectors
  void computeExternalForces(int partition, BLOCK_VECTOR& forceVector, VECTOR& externalForces);
  
  // add the inter-partition spring forces for a rigid body
  void computeDeformableSpringForces(int partition, BLOCK_VECTOR& forceVector, bool debug = false);
  void computeDeformableSpringForcesDefoOnly(int partition, VECTOR& forceVector, bool debug = false);
  void computeDeformableSpringForcesReduced(int partition, BLOCK_VECTOR& forceVector);
  void computeConstrainedSpringForces(int partition, VECTOR& forceVector, bool debug = false);
  void computeConstrainedSpringForcesReduced(int partition, VECTOR& forceVector);

  // compute the quadratic velocity vectors, explicit in time
  void computeQuadraticsExplicit(int partition, BLOCK_VECTOR& quadratics);

  // compute the quadratic velocity vectors, implicit in time
  // PRELIMINARY -- DO NOT USE
  void computeQuadraticsImplicit(int partition, BLOCK_VECTOR& quadratics);

  // compute the jacobians of the implicit quadratic velocities
  // PRELIMINARY -- DO NOT USE
  void computeQuadraticsJacobians(int partition, BLOCK_MATRIX& jacobians);

  // compute the kinetic energy of a partition
  Real kineticEnergy(BLOCK_MATRIX& mass, UNCONSTRAINED_SUBSPACE_INTEGRATOR* integrator);

  // helper for verifyGlobalSpringJacobian -- compute the current spring force
  BLOCK_VECTOR computeCurrentSpringForce();

  // precompute the fast sandwich transforms
  void precomputeSandwiches();

  // presum any sandwiches possible
  void presumSandwiches();

  // read/write out the precomputed sandwiches to disk
  void writeSandwichCache();
  bool readSandwichCache();

  // allocate the block system matrices
  void resizeBlockSystem();

  // compute per-interface spring consts
  void computeSpringConsts();

  // recompute spring-const dependent matrices
  void computeSpringMatrices();

  //////////////////////////////////////////////////////////////////////
  // Skinned versions of the above functions
  //////////////////////////////////////////////////////////////////////
  void computeSkinnedDiagonal(int partition, MATRIX& system, VECTOR& residual, bool dynamic = false);
  void computeSkinnedSpringForces(int partition, VECTOR& forceVector);
  void computeSkinnedSpringJacobians(int partition, MATRIX& springMatrix, bool debug = false, bool dynamic = false);
  void computeSkinnedCouplingJacobian(int unconstrained, int constrained, MATRIX& springMatrix);

  void computeSkinnedDiagonalReduced(int partition, MATRIX& system, VECTOR& residual, bool dynamic = true);
  void computeSkinnedSpringForcesReduced(int partition, VECTOR& forceVector);
  void computeSkinnedSpringJacobiansReduced(int partition, MATRIX& springMatrix, bool dynamic = true);
  void computeSkinnedCouplingJacobianReduced(int unconstrained, int constrained, MATRIX& springMatrix);

  //////////////////////////////////////////////////////////////////////
  // Constraint projection variables and functions start here
  //////////////////////////////////////////////////////////////////////
  BLOCK_MATRIX _constraintJacobian;
  vector<int> _startingConstraintColumn;
  vector<pair<int, int> > _constraintPairs;
  map<pair<int, int>, int> _constraintIndex;
  int _constraintStartingBlock;
  
  // block lagrange multipliers
  BLOCK_VECTOR _lambdas;

  // workspace for constraint Jacobians 
  map<int, MATRIX> _constraintWorkspace3xL;
  map<int, MATRIX> _constraintWorkspaceLx3;
  map<int, MATRIX> _constraintWorkspace3xR;
  map<int, MATRIX> _constraintWorkspaceRx3;
  map<int, MATRIX> _constraintWorkspace3x3;
  map<int, TENSOR3> _constraintWorkspace3xLx3;
  map<int, TENSOR3> _constraintWorkspace3xRx3;

  //////////////////////////////////////////////////////////////////////
  // Try solving the block system using a variety of solvers
  //////////////////////////////////////////////////////////////////////
  int* _Ap;
  int* _Ai;
  double* _Ax;
  vector<Real> _blockValues;
  SPARSE_PETSC_MATRIX _sparseA;

  // compute the spring energy of a partition
  Real springEnergy(int partition);

  // DEBUG: store the full state of the springs in the system
  VECTOR _debugFullSprings;

  // cache the spring forces so that they can be examined later
  vector<VECTOR> _rotationSpringForces;
  vector<VECTOR> _defoSpringForces;

  Real _defoSpringFactor;
};

#endif
