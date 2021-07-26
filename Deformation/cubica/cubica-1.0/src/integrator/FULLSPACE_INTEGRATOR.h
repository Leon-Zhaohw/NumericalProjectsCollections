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
// FULLSPACE_INTEGRATOR.h: interface for the FULLSPACE_INTEGRATOR class.
//
//////////////////////////////////////////////////////////////////////

#ifndef FULLSPACE_INTEGRATOR_H
#define FULLSPACE_INTEGRATOR_H

#include <SETTINGS.h>
#include <TET_MESH.h>
#include <SPARSE_PCG_MATRIX.h>
#if !defined(_WIN32) && !defined(IGNORE_PETSC)
#include <SPARSE_PETSC_MATRIX.h>
#endif
#include <PCG_MATRIX.h>
#include <INCOMPLETE_CHOLESKY.h>
#include <DIAGONAL.h>
#include <COLLISION_RESPONSE.h>
#include <PLANE_COLLISION.h>

//////////////////////////////////////////////////////////////////////
// Timestepping class for the tet mesh, follows the notation
// of [Barbic and James], Appendix C. The variables names
// adhere as closely as possible to those therein.
//
// There are lots of instances where a function call to _tetMesh is
// prepended with TET_MESH:: -- this is so that a pointer of type
// SUBSPACE_TET_MESH* can be passed to the constructor, but the
// fullspace integration can proceed. The TET_MESH:: forces the pointer
// to call the superclass version, not the reduced version.
// 
//////////////////////////////////////////////////////////////////////
class FULLSPACE_INTEGRATOR {

public:
  FULLSPACE_INTEGRATOR(TET_MESH* tetMesh, Real dt = 1.0 / 60.0,
                       Real alpha = 0.01, Real beta = 0.01, Real gravity = 9.8);
  virtual ~FULLSPACE_INTEGRATOR() { delete _preconditioner; };

  // various solver options are available
  virtual void stepFullImplicit();
  virtual void stepFullExplicit();
  virtual void stepSparseImplicit(bool computeExternals = true);
  virtual void stepSparseExplicit();
  virtual void stepSparseQuasistatic();
  virtual bool stepSparseQuasistaticInvertible();
  virtual bool stepSparseImplicitInvertible(bool computeExternals = true,
                                            bool substepping = false);
  virtual bool stepSkinnedImplicit();
  virtual bool stepSkinnedQuasistatic();
  virtual bool stepSkinnedQuasistaticWithCollisions(vector<pair<SURFACE*, int> > collisionPairs);
  virtual bool stepSparseImplicitWithCollisions(vector<pair<SURFACE*, int> > collisionPairs);

  // same as stepSparseImplicitAcceleration, but with acceleration as the primary variable
  virtual void stepImplicitAcceleration();
  void computeAlphas();
  
  void debugStiffness();

  double FPS() { return _totalSteps / _totalTime; };
  double SPF() { return _totalTime / _totalSteps; };
  int totalSteps() { return _totalSteps; };
  Real& forceMultiplier() { return _forceMultiplier; };

  // print out a detailed timing breakdown
  virtual void printTimingBreakdown() { printTimingBreakdown(string("FULLSPACE_INTEGRATOR"), _timingBreakdown, _totalTime); };
  void printSolverTimingBreakdown() { printTimingBreakdown(string("PCG SOLVER"), _ASparse.timingBreakdown(), _ASparse.totalTime()); };
  map<string, double>& timingBreakdown() { return _timingBreakdown; };
  map<string, double>& solverTimingBreakdown() { return _ASparse.timingBreakdown(); };

  // mouse support
  void click(VEC3F& point);
  void drag(VEC3F& point) { _draggedPosition = point; };
  void unclick() { _clickedNode = NULL; };
 
  // drawing options
  void drawClickedNode();
  void drawForceVector();

  // poke the mesh in a deterministic way for debugging
  void poke();

  // add a gravity force
  void addGravity(VEC3F down, float gravity = 1.0f);
  //void addGravity() { addGravity(_gravityDown, _gravityMagnitude); };
  virtual void addGravity() { _externalForces += _gravity; };
  void changeGravity(VEC3F gravityDown, Real gravityMagnitude);
  //void addGravity() { addGravity(VEC3F(0,0,1), -9.8); };
  //void addGravity() { addGravity(VEC3F(0,1,0), -0.1); };

  // collision detection and response for a floor (hacky -- debugging only)
  virtual bool collideWithFloor(VEC3F down, float floorPosition = 0.0);

  // compute rigid translation, rotation, and velocites
  virtual void printRigidComponents();

  // clear all forces
  virtual void clearForces() { _forceVectors.clear(); _forceNodes.clear(); };

  // add an arbitrary force
  void addForce(VEC3F* point, VEC3F& forceVector) { _forceNodes.push_back(point); _forceVectors.push_back(forceVector); };
  
  // accessors
  VECTOR& position()              { return _position; };
  VECTOR& positionOld()           { return _positionOld; };
  VECTOR& velocity()              { return _velocity; };
  VECTOR& velocityOld()           { return _velocityOld; };
  VECTOR& acceleration()          { return _acceleration; };
  VECTOR& accelerationOld()       { return _accelerationOld; };
  Real&   dt()                    { return _dt; };
  Real&   time()                  { return _time; };
  VECTOR& externalForces()        { return _externalForces; };
  VECTOR& externalForcesOld()     { return _externalForcesOld; };
  VECTOR& residual()              { return _residual; };
  VECTOR& temp()                  { return _temp; };
  Real&   solverEps()             { return _solverEps; };
  PRECONDITIONER*& preconditioner() { return _preconditioner; };
  int& maxNewtonSteps()           { return _maxNewtonSteps; };
  int& PCGDigits()                { return _PCGDigits; };
#if defined(_WIN32) || defined(IGNORE_PETSC)
  SPARSE_PCG_MATRIX& solver()     { return _ASparse; };
#else
  SPARSE_PETSC_MATRIX& solver()   { return _ASparse; };
#endif
  vector<VEC3F>& forceVectors()   { return _forceVectors; };
  vector<VEC3F*>& forceNodes()    { return _forceNodes; };
  TET_MESH* tetMesh()             { return _tetMesh; };

  Real& rayleighAlpha()           { return _rayleighAlpha; };
  Real& rayleighBeta()            { return _rayleighBeta; };
  Real* alpha()                   { return &_alpha[0]; };
  Real* accelerationAlpha()       { return &_accelerationAlpha[0]; };
  vector<VEC3F*>& frozenVertices() { return _frozenVertices; };
  bool& substeppingEnabled()      { return _substeppingEnabled; };
  VECTOR& gravity()               { return _gravity; };
  VECTOR& totalMeshForces()       { return _totalMeshForces; };
  VEC3F* clickedNode()            { return _clickedNode; }; 
  VEC3F& clickedPosition()        { return _clickedPosition; };
  VEC3F& draggedPosition()        { return _draggedPosition; };

  VECTOR& positionSample()        { return _positionSample; };
  //SPARSE_PARDISO_MATRIX& dampingMatrix() { return _CDampingSparse; };
  SPARSE_MATRIX& dampingMatrix()  { return _CDampingSparse; };
  SPARSE_MATRIX& A()              { return _tetMesh->stiffnessMatrix(); };
  VECTOR& b()                     { return _residual; };
  SPARSE_MATRIX& systemMatrix()   { return _systemMatrix; };
  int rank()                      { return _rank; };
  
  // cache the static damping matrix
  void cacheStaticDampingFull();
  void cacheStaticDampingSparse();

  // dump state to a filename
  void writeState(string filename);
  void readState(string filename);

  // stream state to already open file
  void writeState(FILE* file);
  void readState(FILE* file);

  // compute external force vector
  virtual void computeExternalForces();

  // compute the current force residual;
  //void forceResidual(VECTOR& position, VECTOR& velocity, VECTOR& acceleration, VECTOR& residual);
  VECTOR& computeForceResidual(bool invertible);

  // add a collision response handler
  void addCollisionResponse(COLLISION_RESPONSE* response) { _collisionResponses.push_back(response); };

  // freeze this vertex
  void freezeVertex(VEC3F* vertex) { _frozenVertices.push_back(vertex); };
  
  // unfreeze all
  void unfreezeAll() { _frozenVertices.clear(); };

  // set everything to zero
  void reset();

  // trying this out for the cubature generator
  VECTOR& computeInternalForces(bool invertible);

  // needed when calling the interface stiffness auto-tuner
  void generateImplicitMatrices(SPARSE_MATRIX& ASparse, bool invertible);
  void generateImplicitAccelerationMatrices(map<string,double>& timingBreakdown);
  virtual void generateQuasistaticMatrices(map<string,double>& timingBreakdown);

  virtual void initializeQuasistaticStep();
  virtual void finalizeQuasistaticStep();

  virtual void initializeImplicitStep();
  virtual void finalizeImplicitAccelerationStep();

  // update the current state using acceleration -- does not update olds however
  virtual void updateState();
  virtual void updateStateUsingAcceleration();

  // build list of unique collisions
  void buildCollisionList(vector<pair<SURFACE*, int> >& collisions,
                          vector<pair<VEC3F*, SURFACE*> >& collisionVertices,
                          vector<SURFACE*>& collisionSurfaces);
  
  // compute collision residual based on collision pairs
  VECTOR computeCollisionResidual(vector<pair<VEC3F*, SURFACE*> >& collisionVertices);

  // compute collision jacobian based on collision pairs
  void addCollisionForceJacobians(vector<pair<VEC3F*, SURFACE*> >& collisionVertices);

protected:
  TET_MESH* _tetMesh;

  // max number of Newton steps
  int _maxNewtonSteps;

  // system size being solved
  int _rank;
  
  // integration constants
  Real _alpha[6];
  Real _accelerationAlpha[5];
  Real _beta;
  Real _gamma;
  Real _dt;
  Real _rayleighAlpha;
  Real _rayleighBeta;
  
  // current simulation time
  Real _time;

  // variables to solve for
  VECTOR _position;
  VECTOR _velocity;
  VECTOR _acceleration;
  VECTOR _positionOld;
  VECTOR _velocityOld;
  VECTOR _accelerationOld;

  // auxiliary solve variables
  VECTOR _residual;
  VECTOR _temp;

  // damping matrix -- gcc doesn't like just '_C'
  MATRIX _CDamping;

  // system matrix
  PCG_MATRIX _A;

  // sparse versions of system and damping matrix
  // if this is not used then virtually no space is allocated,
  // so keeping this around in addition to _A is not a big deal
#if defined(_WIN32) || defined(IGNORE_PETSC)
  SPARSE_PCG_MATRIX _ASparse;
#else
  SPARSE_PETSC_MATRIX _ASparse;
#endif
  //SPARSE_PARDISO_MATRIX _CDampingSparse;
  SPARSE_MATRIX _CDampingSparse;
  PRECONDITIONER* _preconditioner;

  // when we want a cached version of the system matrix that is non PetSc
  SPARSE_MATRIX _systemMatrix;

  // solution vector for sparse matrix -- Pardiso does not let
  // you use the same 'b' and 'x' vector like LAPACK
  VECTOR _solution;

  // timings
  map<string, double> _timingBreakdown;
  double _totalTime;

  // current simulation step
  int _totalSteps;

  // digits that the PCG solver needs to obtain
  int _PCGDigits;

  // mouse support
  VEC3F* _clickedNode;
  VEC3F _clickedPosition;
  VEC3F _draggedPosition;

  // external force variables
  VECTOR _externalForces;
  VECTOR _externalForcesOld;
  vector<VEC3F> _forceVectors;
  vector<VEC3F*> _forceNodes;
  Real _forceMultiplier;
  Real _gravityMagnitude;
  VEC3F _gravityDown;

  VECTOR _gravity;

public:  
  // cached implicit invertible element matrices
  vector<MATRIX3> _Us;  // U from F diagonalization
  vector<MATRIX3> _Vs;  // V from F diagonalization
  vector<MATRIX3> _Fhats;  // Fhat from F diagonalization
  vector<MATRIX> _stiffnesses;  // diagonalized and clamped stiffness

protected:  
  // precision of the Newton-Raphson and CG solvers
  Real _solverEps;
 
  // collision response handling
  vector<COLLISION_RESPONSE*> _collisionResponses;

  // collision response force vector
  VECTOR _collisionForces;
  
  // friction force vector
  VECTOR _frictionForces;

  // collision response force jacobian
  SPARSE_MATRIX _collisionForceJacobians;
  
  // keep the position around for online cubature
  VECTOR _positionSample;

  // is substepping enabled?
  bool _substeppingEnabled;

  // total newton steps seen so far
  int _totalNewtonStepsSeen;

  // solve using the invertible finite element stiffness matrix
  // from "Robust Quasisatic Finite Elements and Flesh Simulation"
  // Teran, Sifakis, Irving, Fedkiw, SCA 2005
  void solveInvertibleStiffness(VECTOR& x, VECTOR& b, bool quasistatic);
  
  // Do the matrix-vector multiply for solveInvertibleQuasistatic()
  void invertibleStiffnessMultiply(VECTOR& direction, VECTOR& answer);

  // keep this one around just in case (does not use diagonalization)
  void stiffnessMultiply(VECTOR& direction, VECTOR& answer);

  VECTOR _totalMeshForces;

public:  
  // cache the deformation gradient diagonalizations
  void cacheDiagonalizations();

protected:
  // print out some arbitrary timing breakdown
  void printTimingBreakdown(string name, map<string, double>& timingBreakdown, double totalTime);

  // compute the friction forces and store them in _frictionForces
  VECTOR& computeFrictionForces();

  // compute the collision forces and store them in _collisionForces
  VECTOR& computeCollisionForces();

  // compute the collision jacobians and store them in _collisionForceJacobians
  SPARSE_MATRIX& computeCollisionForceJacobians();
  
  // update the velocity from the current position --
  // this is necessary for the proper collision forces to be computed
  void updateVelocities();

  // freeze current vertex positions
  vector<VEC3F*> _frozenVertices;

  // apply constraints, ie "frozen vertices"
  void applyResidualConstraints();
  void applyMatrixConstraints();
};

#endif
