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
// SUBSPACE_INTEGRATOR.h: interface for the SUBSPACE_INTEGRATOR class.
//
//////////////////////////////////////////////////////////////////////

#ifndef SUBSPACE_INTEGRATOR_H
#define SUBSPACE_INTEGRATOR_H

#include <SETTINGS.h>
#include <SUBSPACE_TET_MESH.h>
#include <TIMER.h>
#include <PCG_MATRIX.h>
#include <FULLSPACE_INTEGRATOR.h>
#include <SURFACE.h>

#include <ThreadSpecificData.h>

//////////////////////////////////////////////////////////////////////
// Timestepping class for the tet mesh, follows the notation
// of [Barbic and James], Appendix C. The variables names
// adhere as closely as possible to those therein.
//
// This does not explicitly subclass FULLSPACE_INTEGRATOR, it just
// mirrors the functionality. The two classes are sufficiently
// different that subclassing would just introduce compatibility 
// headaches. You wouldn't swap one with the other at runtime anyway,
// so there's no reason to keep their pointers interchangeable.
//////////////////////////////////////////////////////////////////////
class SUBSPACE_INTEGRATOR {

public:
  // The last two parameters control whether or not we
  // use inexact stiffness matrices in the newton solve
  // and the tolerance for adding new key tet contributions
  // to the stiffness matrix.
  SUBSPACE_INTEGRATOR(SUBSPACE_TET_MESH* tetMesh, Real dt = 1.0 / 60.0,
                      Real alpha = 0.01, Real beta = 0.01,
                      Real gravity = 9.8,
                      bool inexactNewton = false,
                      Real deformationTolerance = -1.0);
  virtual ~SUBSPACE_INTEGRATOR() {};

  void stepQuasistatic();
  virtual void stepInvertibleQuasistatic(bool directSolve = true);
  void stepImplicit(bool finalize = true);
  virtual void stepInvertibleImplicit(bool finalize = true, bool debugPrint = false);
  void stepInvertibleSemiImplicit(bool directSolve = true);

  // same as stepInvertibleImplicit, but using acceleration as the primary variable
  // returns the number of Newton steps taken
  virtual int stepImplicitAcceleration();
  virtual int stepImplicitWithExplicitCollisions(vector<pair<SURFACE*, int> >& collisions);
  virtual int stepImplicitWithImplicitCollisions(vector<pair<SURFACE*, int> >& collisions);
  void stepExplicit();
  Real& forceMultiplier() { return _forceMultiplier; };
  void finalizeStep() { _tetMesh->updateSurfaceMesh(); _externalForces.clear(); };

  // same as stepInvertibleImplicit, but using velocity as the primary variable
  virtual int stepImplicitVelocity();

  // mouse support
  virtual bool click(VEC3F& point, Real maxRadius = -1.0);
  void drag(VEC3F& point) { _draggedPosition = point; };
  void unclick() { _clickedNode = NULL; _clickedID = -1; };
 
  // drawing options
  virtual void drawClickedNode();
  virtual void drawForceVector();

  // poke the mesh in a deterministic way for debugging
  void poke();

  // gradually pull the mesh in a deterministic way for debugging
  void pull();

  // accessors
  VECTOR& externalForces()  { return _externalForces; };
  VECTOR& position()        { return _position; };
  VECTOR& positionOld()     { return _positionOld; };
  VECTOR& velocity()        { return _velocity; };
  VECTOR& velocityOld()     { return _velocityOld; };
  VECTOR& acceleration()    { return _acceleration; };
  VECTOR& accelerationOld() { return _accelerationOld; };
  MATRIX& dampingMatrix()   { return _CDamping; };
  Real&   dt()              { return _dt; };
  Real&   solverEps()       { return _solverEps; };
  vector<VEC3F>& forceVectors() { return _forceVectors; };
  vector<VEC3F*>& forceNodes() { return _forceNodes; };
  vector<int>& forceNodeIDs() { return _forceNodeIDs; };
  Real* alpha() { return _alpha; };
  Real* accelerationAlpha() { return _accelerationAlpha; };
  Real* velocityAlpha() { return _velocityAlpha; };
  int& maxNewtonSteps()     { return _maxNewtonSteps; };
  bool& useKryslStiffness() { return _useKryslStiffness; };
  int& totalSteps()         { return _totalSteps; };
  int rank()                { return _rank; };
  const Real beta()         { return _beta; };
  const Real gamma()        { return _gamma; };
  VEC3F* clickedNode()      { return _clickedNode; };
  vector<MATRIX3>& Us()     { return _Us; };
  vector<MATRIX3>& Vs()     { return _Vs; };
  vector<MATRIX3>& Fhats()  { return _Fhats; };
  vector<MATRIX>& stiffnesses() { return _stiffnesses; };
  MATRIX& reducedMass()     { return _tetMesh->reducedMass(); };

  // print out a detailed timing breakdown
  virtual void printTimingBreakdown();
  map<string, double>& timingBreakdown() { return _timingBreakdown; };
 
  // partitioned mesh support function -- 
  // form A matrix (stiffness, damping, etc.) and b (residual
  // so that it can be passed out and then solved by the partitioned 
  // integrator
  virtual void initializeImplicitStep();

  // Generate the matrices for an implicit step -- should be called 
  // right after initializeImplicitStep
  //
  // This has been forked out into its own function in case we want
  // to do more than one Newton step
  virtual void generateImplicitMatricesInvertible(map<string,double>& timingBreakdown);
  virtual void generateMultibodyMatrices();
  virtual void generateImplicitMatrices(map<string,double>& timingBreakdown, bool useKrysl = false);
  virtual void generateImplicitAccelerationMatrices(map<string,double>& timingBreakdown, bool useKrysl = false);
  virtual void generateImplicitVelocityMatrices(map<string,double>& timingBreakdown);
  void generateImplicitMatricesDebug();

  // compute velocity, acceleration and position using final values
  virtual void finalizeImplicitStep();

  // compute position and velocity using final acceleration values
  virtual void finalizeImplicitAccelerationStep();
  
  // compute position and velocity using final velocity values
  virtual void finalizeImplicitVelocityStep();
  
  // same as above, but for quasistatic solve
  virtual void initializeQuasistaticStep();
  virtual void generateQuasistaticMatrices(map<string,double>& timingBreakdown);
  virtual void finalizeQuasistaticStep();
  MATRIX* A() { return &_A; };
  VECTOR& b() { return _residual; };
  VECTOR& workspace() { return _temp; };
  VECTOR& reducedGravity() { return _gravity; };
  VECTOR& totalMeshForces() { return _totalMeshForces; }; 

  // add precomputed gravity;
  virtual void addGravity() { _externalForces += _gravity; };
  virtual void addBodyForce(VEC3F bodyForce) { _externalForces += _bodyForceU * bodyForce.toVector(); };
  virtual void changeGravity(VEC3F gravityDown, Real gravityMagnitude);
  VECTOR& gravity() { return _gravity; };
  
  // clear all forces
  virtual void clearForces()
  {
    _forceVectors.clear();
    _forceNodes.clear();
    _forceNodeIDs.clear();
    //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    //cout << " EXTERNAL CLEAR DISABLED " << endl;
    _externalForces.clear();
  }

  // update quantities after the tet mesh basis has been altered
  void updateBasis(MATRIX& oldBasis, MATRIX& newBasis);

  // update quantities after the cubature points have been altered
  void updateCubature(bool invertible);

  // stomp everything related to the reduced model
  void resetBasis();

  // read and write state
  void readState(string filename);
  void writeState(string filename);
 
  VEC3F& clickedPosition() { return _clickedPosition; };
  VEC3F& draggedPosition() { return _draggedPosition; };
  SUBSPACE_TET_MESH* tetMesh() { return _tetMesh; };
  int clickedID() { return _clickedID; };

  // reset everything to zero
  virtual void reset();

  // compute and return full coordinate velocity
  virtual VECTOR fullCoordinateVelocity();

  // compute velocity using the current _tetMesh q
  virtual VECTOR computeNewestReducedVelocity();

  // Checks the cached deformation state of the given key tet
  // against the current deformation.  If the state has changed
  // by a given threshold, then mark the tet as dirty, and
  // update its deformation state.
  bool updateTetDeformation( int i, Real threshold,
                             MATRIX3 &F,
                             MATRIX3 &U, MATRIX3 &Fhat, MATRIX3 &V );

  // update the current state using acceleration -- does not update olds however
  virtual void updateStateUsingAcceleration();
  virtual void updateStateUsingVelocity();
  virtual void updateState();

  // compute the kinetic energy of the entire mesh
  virtual Real kineticEnergy();

  // get the velocity of a specific vertex
  VEC3F vertexVelocity(int vertexID);
  VEC3F vertexVelocityOld(int vertexID);

  // set the timestep size to something else (i.e. recompute the Newmark consts)
  void setTimestep(Real dt);

public:
  // This is a bad, lazy way of doing things.  I will
  // fix it later.
  // FIXME
  // The maximum external force magnitude that we can
  // apply to a node on the surface.
  static Real MAX_EXTERNAL_FORCE;

protected:
  SUBSPACE_TET_MESH* _tetMesh;

  // max number of Newton steps
  int _maxNewtonSteps;

  // system size being solved
  int _rank;
  
  // integration constants
  Real _alpha[6];
  Real _accelerationAlpha[5]; // alphas when performing an acceleration level update
  Real _velocityAlpha[5]; // alphas when performing a velocity level update
  Real _beta;
  Real _gamma;
  Real _dt;
  Real _rayleighAlpha;
  Real _rayleighBeta;
  
  // current simulation time
  Real _time;

  // precomputed gravity vector
  VECTOR _gravity;
  Real _gravityMagnitude;
  VEC3F _gravityDown;
  
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
  VECTOR _totalMeshForces;

  // damping matrix -- gcc doesn't like just '_C'
  MATRIX _CDamping;

  // system matrix
  PCG_MATRIX _A;

  // mouse support
  VEC3F* _clickedNode;
  int _clickedID;
  VEC3F _clickedPosition;
  VEC3F _draggedPosition;

  // external force variables
  VECTOR _externalForces;
  vector<VEC3F> _forceVectors;
  vector<VEC3F*> _forceNodes;
  vector<int> _forceNodeIDs;
  Real _forceMultiplier;

  // A workspace for storing vertex subbases extracted
  // from the main basis.
  ThreadSpecificData<MATRIX> _subBasisWorkspace;
 
  // cached implicit invertible element matrices
  vector<MATRIX3> _Us;    // U from F diagonalization
  vector<MATRIX3> _Vs;    // V from F diagonalization
  vector<MATRIX3> _Fhats; // Fhats from F diagonalization
  vector<MATRIX> _stiffnesses;  // diagonalized and clamped stiffness

  // Whether or not to use inexact stiffness matrices in
  // Newton solves
  bool _inexactNewton;

  // Validity flags for key tets.  This is used for inexact
  // Newton solves when we wish to only perform partial
  // updates to the reduced stiffness matrix
  vector<bool> _updateKeyTetDeformation;

  // Threshold for deciding whether or not a tet has deformed
  // enough from its previous state.
  Real _deformationTolerance;

  // Diagnostic information for inexact newton stuff
  Real _totalSolveSteps;
  Real _totalSkipped;

  // timings
  map<string, double> _timingBreakdown;
  double _totalTime;

  // precision of Newton-Raphson and CG solvers
  Real _solverEps;
  
  // current simulation step
  int _totalSteps;

  // do a Krysl-style reduced stiffness matrix instead of using
  // cubature
  bool _useKryslStiffness;

  // node being pulled by pull()
  VEC3F* _pulledNode;

  // fast matrix to add a body force in any direction
  MATRIX _bodyForceU;

  // external force functions
  virtual void computeExternalForces();
  
  // solve using the invertible finite element stiffness matrix
  // from "Robust Quasisatic Finite Elements and Flesh Simulation"
  // Teran, Sifakis, Irving, Fedkiw, SCA 2005
  void solveInvertibleStiffness(VECTOR& x, VECTOR& b, bool quasistatic);
  
  // cache the deformation gradient diagonalizations.
  // If checkTetState == true, then we will use the updateTetDeformation
  // function to check whether or not a key tet needs to have its
  // stiffness matrix contribution updated.  If not, we will not
  // bother caching its stiffness density.
  virtual void cacheDiagonalizations( bool checkTetState = false );
  
  // if doing Krysl-style computation, cache the diagonalized
  // deformation gradients for all the tets
  virtual void cacheKryslDiagonalizations();
  
  // Do the matrix-vector multiply for solveInvertibleQuasistatic()
  void invertibleStiffnessMultiply(VECTOR& direction, VECTOR& answer, int verbose = 0);

  // debug the stiffness matrix multiply
  void debugStiffness();

  // add explicit collision forces
  virtual void addExplicitCollisions(vector<pair<SURFACE*, int> >& collisions);

  // add implicit collision forces directly to the given residual
  virtual void addImplicitCollisionResiduals(vector<pair<SURFACE*, int> >& collisions);

  // add implicit collision forces directly to the Jacobian
  void addImplicitCollisionJacobians(vector<pair<SURFACE*, int> >& collisions);
};

#endif
