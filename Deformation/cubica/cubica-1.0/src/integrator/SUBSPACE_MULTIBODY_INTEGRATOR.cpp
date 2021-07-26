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

#include "SUBSPACE_MULTIBODY_INTEGRATOR.h"
#include "PARTITIONED_SKINNED_SUBSPACE_TET_MESH.h"

// uses OpenMP in SUBSPACE_MULTIBODY_INTEGRATOR
#define USING_MULTIBODY_OPENMP 1

#ifdef USING_MULTIBODY_OPENMP
#include <omp.h>
#endif

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
SUBSPACE_MULTIBODY_INTEGRATOR::SUBSPACE_MULTIBODY_INTEGRATOR(PARTITIONED_SUBSPACE_TET_MESH* tetMesh, Real dt, Real springConst, Real dampingConst, Real alpha, Real beta) :
  _partitionedMesh(tetMesh),
  _dt(dt),
  _totalSteps(0),
  _solverEps(1e-3),
  _rayleighAlpha(alpha),
  _rayleighBeta(beta),
  _totalTime(0.0),
  _springConst(springConst),
  _dampingConst(dampingConst),
  _forceMultiplier(1.0),
  _maxNewtonSteps(20),
  _totalNewtonStepsSeen(0),
  _maxNewtonStepsSeen(0),
  _newtonStalls(0),
  _clickedPartition(-1),
  _constraintStartingBlock(-1),
  _Ap(NULL),
  _Ai(NULL),
  _Ax(NULL)
{
  _partitions = tetMesh->partitions();
  assert(_partitions > 0);

  if (_partitions <= 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Zero partitions found! Skeleton path is probably incorrect! " << endl;
    exit(0);
  }

  // if there's no sandwich cache, there is a bad inertia cache
  string sandwichFile = _partitionedMesh->partitionPath() + string(".sandwiches");
  cout << " Reading from sandwich cache " << sandwichFile.c_str() << " ... "; flush(cout);

  // allocate the integrators
  _integrators = new SUBSPACE_INTEGRATOR*[_partitions];
  for (int x = 0; x < _partitions; x++)
  {
    //if (tetMesh->mesh(x)->unconstrained())
    if (tetMesh->unconstrained(x))
    {
      cout << " Created unconstrained integrator for partition " << x << endl;
      _integrators[x] = new UNCONSTRAINED_SUBSPACE_INTEGRATOR((SUBSPACE_TET_MESH*)(tetMesh->mesh(x)), _dt, alpha, beta);
    }
    else
    {
      cout << " Created constrained integrator for partition " << x << endl;
      _integrators[x] = new SUBSPACE_INTEGRATOR((SUBSPACE_TET_MESH*)(tetMesh->mesh(x)), _dt, alpha, beta);
    }
  }

  // use whatever Newmark consts the partition integrators are using
  _beta = _integrators[0]->beta();
  _gamma = _integrators[0]->gamma();

  // get unconstrained partition IDs
  _unconstrainedIDs.clear();
  _partitionedMesh->getUnconstrainedIDs(_unconstrainedIDs);

  for (unsigned int x = 0; x < _unconstrainedIDs.size(); x++)
    _inverseUnconstrainedIDs[_unconstrainedIDs[x]] = x;

  // compute how many blocks will be needed, including rigid blocks
  _totalBlocks = 0;
  for (int x = 0; x < _partitions; x++)
    if (_partitionedMesh->unconstrained(x))
      _totalBlocks += 3;
    else
      _totalBlocks++;

  // compute the block number that each partition starts at
  for (int x = 0, i = 0; x < _partitions; x++)
  {
    _startingBlock.push_back(i);
    i += (_partitionedMesh->unconstrained(x)) ? 3 : 1;
  }

  // allocate the necessary blocks
  _forceVectors.resize(_partitions);
  for (int x = 0; x < _partitions; x++)
  {
    int rank = _partitionedMesh->rank(x);
    if (_partitionedMesh->unconstrained(x))
    {
      _forceVectors[x].resizeAndWipe(3);
      _forceVectors[x].resizeAndWipeBlock(0, 3);
      _forceVectors[x].resizeAndWipeBlock(1, 3);
      _forceVectors[x].resizeAndWipeBlock(2, rank);
    }
    else
    {
      _forceVectors[x].resizeAndWipe(1);
      _forceVectors[x].resizeAndWipeBlock(0, rank);
    }
  }

  // preallocate the block matrices
  resizeBlockSystem();

  // look for a sandwich cache, and compute them if nothing
  // is found
  if (readSandwichCache() == false)
    precomputeSandwiches();

  // compute per-interface spring consts
  computeSpringConsts();

  // sum together whatever sandwiches you can
  presumSandwiches();

  // recompute spring-const dependent matrices
  computeSpringMatrices();

  // cache the spring forces so that they can be examined later
  _rotationSpringForces.resize(_partitions);
  _defoSpringForces.resize(_partitions);

  _defoSpringFactor = 1;
}

SUBSPACE_MULTIBODY_INTEGRATOR::~SUBSPACE_MULTIBODY_INTEGRATOR()
{
  for (int x = 0; x < _partitions; x++)
    delete _integrators[x];
  delete[] _integrators;
}

//////////////////////////////////////////////////////////////////////
// Add gravity
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::addGravity()
{
  for (int x = 0; x < _partitions; x++)
    _integrators[x]->addGravity();
}

//////////////////////////////////////////////////////////////////////
// Add a specified body force
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::addBodyForce(VEC3F bodyForce)
{
  for (int x = 0; x < _partitions; x++)
    _integrators[x]->addBodyForce(bodyForce);
}

//////////////////////////////////////////////////////////////////////
// Change gravity direction
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::changeGravity(int partition, VEC3F gravityDown, Real gravityMagnitude)
{
  // first call the original so that the vectors are set correctly
  _integrators[partition]->changeGravity(gravityDown, gravityMagnitude);
}

//////////////////////////////////////////////////////////////////////
// Change gravity direction
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::changeGravity(VEC3F gravityDown, Real gravityMagnitude)
{
  Real totalMass = 0;
  for (int x = 0; x < _partitions; x++)
  {
    // get the mass of this partition
    //Real mass = _partitionedMesh->mesh(x)->totalMass();
    Real mass = _partitionedMesh->mesh(x)->totalMass();
    totalMass += mass;

    // first call the original so that the vectors are set correctly
    _integrators[x]->changeGravity(gravityDown, mass * gravityMagnitude);
  }
  cout << " Total mass: " << totalMass << endl;
}

//////////////////////////////////////////////////////////////////////
// compute the non-linear mass matrix
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::computeMassMatrix(int partition, BLOCK_MATRIX& massMatrix)
{
  massMatrix.clear();

  UNCONSTRAINED_SUBSPACE_INTEGRATOR* integrator = (UNCONSTRAINED_SUBSPACE_INTEGRATOR*)_integrators[partition];
  UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_partitionedMesh->mesh(partition);

  // get the interia tensor in anticipation of computing rotation A
  const MATRIX& Ithetatheta = mesh->inertiaTensor();

  /////////////////////////////////////////////////////////////////////////
  // copy in rotation A
  /////////////////////////////////////////////////////////////////////////
  
  // shouldn't assume that the quadratics are implicit -- make a separate function
  // to add these in if necessary
  MATRIX rotationA = Ithetatheta;
  massMatrix.add(rotationA, 1,1);

  /////////////////////////////////////////////////////////////////////////
  // Assemble translation A
  /////////////////////////////////////////////////////////////////////////

  // count the clones
  int totalClones = 0;
  for (int x = 0; x < _partitions; x++)
  {
    vector<pair<int,int> > clones = _partitionedMesh->clonedVertices(partition, x);
    if (clones.size() == 0) continue;
    totalClones += clones.size();
  }

  // get the translation mass matrix
  MATRIX translationA = mesh->totalMass() * MATRIX3::I();

  // finalize the mass matrix
  massMatrix.add(translationA, 0,0);

  /////////////////////////////////////////////////////////////////////////
  // Assemble rigid-rotation coupling terms
  /////////////////////////////////////////////////////////////////////////
  mesh->refreshSitBar();
  MATRIX3 SitBarTilde = MATRIX3::cross(mesh->SitBar());

  // get the current rotation
  MATRIX3 A = integrator->rotation().toExplicitMatrix3x3();

  MATRIX AST = A * SitBarTilde.transpose();
  massMatrix.add(AST, 0, 1);

  MATRIX SAT = AST.transpose();
  massMatrix.add(SAT, 1, 0);

  // if it's rigid only, we're done
  if (massMatrix.blockRows() == 2) return;
 
  /////////////////////////////////////////////////////////////////////////
  // Assemble deformable-translation couplings
  /////////////////////////////////////////////////////////////////////////
  MATRIX& SiBar = mesh->SiBar();
  MATRIX ASiBar = A * SiBar;
  massMatrix.add(ASiBar, 0, 2);
  MATRIX SiBarTAT = ASiBar.transpose();
  massMatrix.add(SiBarTAT, 2, 0);

  /////////////////////////////////////////////////////////////////////////
  // Assemble deformable-rotation coupling
  /////////////////////////////////////////////////////////////////////////
  MATRIX Ithetaf = mesh->rotationDefoTensor();

  massMatrix.add(Ithetaf, 1, 2);
  MATRIX IthetafT = Ithetaf.transpose();
  massMatrix.add(IthetafT, 2, 1);

  /////////////////////////////////////////////////////////////////////////
  // add deformable components
  /////////////////////////////////////////////////////////////////////////
  massMatrix.add(mesh->reducedMass(), 2,2);
}

//////////////////////////////////////////////////////////////////////
// compute the generalized forces
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::computeExternalForces(int partition, BLOCK_VECTOR& forceVector, VECTOR& externalForces)
{
  UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_partitionedMesh->mesh(partition);
  UNCONSTRAINED_SUBSPACE_INTEGRATOR* integrator = (UNCONSTRAINED_SUBSPACE_INTEGRATOR*)_integrators[partition];

  VECTOR fullExternals = mesh->U() * externalForces;

  // compute the rigid translation component
  VECTOR QR(3);
  for (int x = 0; x < fullExternals.size() / 3; x++)
  {
    QR[0] += fullExternals[3 * x];
    QR[1] += fullExternals[3 * x + 1];
    QR[2] += fullExternals[3 * x + 2];
  }
  QR[0] += integrator->translationForce()[0];
  QR[1] += integrator->translationForce()[1];
  QR[2] += integrator->translationForce()[2];

  // TODO: O(N)
  // compute rotational force
  mesh->SUBSPACE_TET_MESH::updateFullMesh();
  vector<VEC3F>& vertices = mesh->vertices();
  QUATERNION& quaternion = mesh->rotationQuaternion();
  MATRIX A = quaternion.toExplicitMatrix3x3();
  VECTOR QTheta(3);
  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    VEC3F& vertex = vertices[x];
    MATRIX uBar = MATRIX::cross(vertex);

    MATRIX AuBar = A * uBar;
    MATRIX AuBarT = AuBar.transpose();

    VECTOR F(3);
    F(0) = fullExternals[3 * x];
    F(1) = fullExternals[3 * x + 1];
    F(2) = fullExternals[3 * x + 2];

    QTheta += -1.0 * AuBarT * F;
  }

  // populate the entry
  forceVector[0] += QR;
  forceVector[1] += QTheta;
}

//////////////////////////////////////////////////////////////////////
// add the inter-partition spring forces
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::computeConstrainedSpringForces(int partition, VECTOR& forceVector, bool debug)
{
  SUBSPACE_INTEGRATOR* integrator = _integrators[partition];
  SUBSPACE_TET_MESH* leftMesh = (SUBSPACE_TET_MESH*)_partitionedMesh->mesh(partition);
  leftMesh->updateFullMesh();
  vector<VEC3F>& leftVertices = leftMesh->vertices();

  int totalClones = 0;
  VECTOR fullSprings(leftMesh->U().rows());
  VECTOR fullDampers(leftMesh->U().rows());
  VECTOR fullVelocity = leftMesh->U() * integrator->velocity();
  for (int x = 0; x < _partitions; x++)
  {
    vector<pair<int,int> > clones = _partitionedMesh->clonedVertices(partition, x);
    if (clones.size() == 0) continue;

    SUBSPACE_TET_MESH* rightMesh = (SUBSPACE_TET_MESH*)_partitionedMesh->mesh(x);
    rightMesh->updateFullMesh();
    vector<VEC3F>& rightVertices = rightMesh->vertices();

    MATRIX3 rightRotation(MATRIX3::I());
    VEC3F rightTranslation;
    rightRotation = _partitionedMesh->quaternionRotation(x).toExplicitMatrix3x3();
    rightTranslation = _partitionedMesh->rigidTranslation(x);

    VECTOR rightVelocity = rightMesh->U() * _integrators[x]->velocity();

    for (unsigned int y = 0; y < clones.size(); y++)
    {
      int leftIndex = clones[y].first;
      int rightIndex = clones[y].second;

      VEC3F leftVertex = leftVertices[leftIndex];
      VEC3F rightVertex = rightVertices[rightIndex];

      // no transform necessary for left vertex
      rightVertex = rightRotation * rightVertex + rightTranslation;

      pair<int, int> xy(partition, x);
      VEC3F force = -_defoSpringFactor * _interfaceSpringConsts[xy] * (leftVertex - rightVertex);

      fullSprings[3 * leftIndex] += force[0];
      fullSprings[3 * leftIndex + 1] += force[1];
      fullSprings[3 * leftIndex + 2] += force[2];

      // rotate the velocity
      VEC3F rightVelocityRotated;
      rightVelocityRotated[0] = rightVelocity[3 * rightIndex];
      rightVelocityRotated[1] = rightVelocity[3 * rightIndex + 1];
      rightVelocityRotated[2] = rightVelocity[3 * rightIndex + 2];
      rightVelocityRotated = rightRotation * rightVelocityRotated;

      fullDampers[3 * leftIndex] += _dampingConst * (fullVelocity[3 * leftIndex] - rightVelocityRotated[0]);
      fullDampers[3 * leftIndex + 1] += _dampingConst * (fullVelocity[3 * leftIndex + 1] - rightVelocityRotated[1]);
      fullDampers[3 * leftIndex + 2] += _dampingConst * (fullVelocity[3 * leftIndex + 2] - rightVelocityRotated[2]);

      totalClones++;
    }
  }
  // DEBUG: store full springs to peek at later
  _debugFullSprings = fullSprings;

  // compute the deformable component
  forceVector = leftMesh->U() ^ fullSprings;

  VECTOR deformableDamping = leftMesh->U() ^ fullDampers;
  forceVector -= deformableDamping;
}

//////////////////////////////////////////////////////////////////////
// add the inter-partition spring forces
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::computeConstrainedSpringForcesReduced(int partition, VECTOR& forceVector)
{
  SUBSPACE_INTEGRATOR* integrator = _integrators[partition];
  SUBSPACE_TET_MESH* leftMesh = (SUBSPACE_TET_MESH*)_partitionedMesh->mesh(partition);

  forceVector.resizeAndWipe(leftMesh->rank()); 
  for (int x = 0; x < _partitions; x++)
  {
    if (!_partitionedMesh->neighbors(partition, x)) continue;
    vector<pair<int,int> > clones = _partitionedMesh->clonedVertices(partition, x);
   
    pair<int, int> xy(partition, x);
    Real springConst = _interfaceSpringConsts[xy];
    forceVector.axpy(-springConst, _Ui_uiBarO[xy]);

    SUBSPACE_TET_MESH* rightMesh = (SUBSPACE_TET_MESH*)_partitionedMesh->mesh(x);

    if (_partitionedMesh->constrained(x))
    {
      forceVector.axpy(springConst, _Ui_ukBarO[xy]);
      _Ui_Uk[xy].gemvInplace(springConst, rightMesh->q(), forceVector);
    }
    else
    {
      MATRIX3 rightRotation = _partitionedMesh->quaternionRotation(x).toExplicitMatrix3x3();
      VECTOR rightTranslation = _partitionedMesh->rigidTranslation(x).toVector();

      _Ui_I3nx3[xy].gemvInplace(springConst, rightTranslation, forceVector);
      _Ui_R_ukBarO[xy].vectorTransformAxpy(springConst, rightRotation, forceVector);

      MATRIX workspaceLxR(_partitionedMesh->rank(partition), _partitionedMesh->rank(x));
      _Ui_R_Uk[xy].transformInplace(rightRotation, workspaceLxR);
      workspaceLxR.gemvInplace(springConst, rightMesh->q(), forceVector);
    }
  }

  // spring const is baked into _spring_Ui_Ui
  _spring_Ui_Ui[partition].gemvInplace(-1.0, leftMesh->q(), forceVector);

  _Ui_Ui[partition].gemvInplace(-_dampingConst, integrator->velocity(), forceVector);
}

//////////////////////////////////////////////////////////////////////
// add the inter-partition spring forces
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::computeSkinnedSpringForcesReduced(int partition, VECTOR& forceVector)
{
  UNCONSTRAINED_SUBSPACE_INTEGRATOR* integrator = (UNCONSTRAINED_SUBSPACE_INTEGRATOR*)_integrators[partition];
  UNCONSTRAINED_SUBSPACE_TET_MESH* leftMesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_partitionedMesh->mesh(partition);

  MATRIX3 leftRotation = integrator->rotation().toExplicitMatrix3x3();
  MATRIX3 leftRotationT = leftRotation.transpose();
  VECTOR leftTranslation = integrator->translation().toVector();

  // compute the total clones
  int totalClones = 0;
  for (int x = 0; x < _partitions; x++)
  {
    if (!_partitionedMesh->neighbors(partition, x)) continue;
    vector<pair<int,int> > clones = _partitionedMesh->clonedVertices(partition, x);
    totalClones += clones.size();
  }

  // deformable component
  VECTOR deformableForce(leftMesh->rank());
  VECTOR tempForce(leftMesh->rank());
  VECTOR deformableDamping(leftMesh->rank());
  for (int x = 0; x < _partitions; x++)
  {
    if (!_partitionedMesh->neighbors(partition, x)) continue;
    vector<pair<int,int> > clones = _partitionedMesh->clonedVertices(partition, x);
    pair<int, int> xy(partition, x);
    Real springConst = _interfaceSpringConsts[xy];

    SUBSPACE_TET_MESH* rightMesh = (SUBSPACE_TET_MESH*)_partitionedMesh->mesh(x);
    UNCONSTRAINED_SUBSPACE_INTEGRATOR* rightIntegrator = (UNCONSTRAINED_SUBSPACE_INTEGRATOR*)_integrators[x];
    MATRIX3 rightRotation = _partitionedMesh->quaternionRotation(x).toExplicitMatrix3x3();
    VECTOR rightTranslation = _partitionedMesh->rigidTranslation(x).toVector();
    MATRIX3 RTR = leftRotationT * rightRotation;

    MATRIX& workspaceLx3 = _workspaceLx3[partition];
    MATRIX& workspaceLxR = _workspaceLxR[xy];

    // add the rotated right translation
    _Ui_R_I3nx3[xy].transformInplace(leftRotationT, workspaceLx3);
    workspaceLx3.gemvInplace(springConst, rightTranslation, deformableForce);

    // add the rotated right rest pose
    _Ui_R_ukBarO[xy].vectorTransformAxpy(springConst, RTR, deformableForce);

    // add the rotated right deformation 
    _Ui_R_Uk[xy].transformInplace(RTR, workspaceLxR);
    workspaceLxR.gemvInplace(springConst, rightMesh->q(), deformableForce);
    
    //deformableDamping += _dampingConst * _partitionedMesh->interfaceU(partition, x) * integrator->velocity();
    workspaceLxR.gemvInplace(-_dampingConst, rightIntegrator->velocity(), deformableDamping);
  }

  // add the left translation
  MATRIX& workspaceLx3 = _workspaceLx3[partition];
  _spring_summed_Ui_R_I3nx3[partition].transformInplace(leftRotationT, workspaceLx3);
  workspaceLx3.gemvInplace(-1.0, leftTranslation, deformableForce);

  // add the left rest pose
  deformableForce.axpy(-1.0, _spring_summed_Ui_uiBarO[partition]);

  // add the left deformation
  _spring_Ui_Ui[partition].gemvInplace(-1.0, leftMesh->q(), deformableForce);

  // new damping model, add this at the end
  deformableDamping += _dampingConst * _Ui_Ui[partition] * integrator->velocity();

  forceVector = deformableForce - deformableDamping;
}

//////////////////////////////////////////////////////////////////////
// add the inter-partition spring forces
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::computeSkinnedSpringForces(int partition, VECTOR& forceVector)
{
  UNCONSTRAINED_SUBSPACE_TET_MESH* leftMesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_partitionedMesh->mesh(partition);
  UNCONSTRAINED_SUBSPACE_INTEGRATOR* integrator = (UNCONSTRAINED_SUBSPACE_INTEGRATOR*)_integrators[partition];
  vector<VEC3F>& leftVertices = leftMesh->vertices();
  leftMesh->updateFullMesh();

  QUATERNION rotationQuaternion = leftMesh->rotationQuaternion();
  MATRIX3 rotation = rotationQuaternion.toExplicitMatrix3x3();
  MATRIX3 rotationTranspose = rotation.transpose();

  // compute new translation quantities
  VEC3F translation = leftMesh->rigidTranslation();

  // get the current velocity
  VECTOR qDot = integrator->velocity();

  int totalClones = 0;
  VECTOR fullSprings(leftMesh->U().rows());
  VECTOR fullDampers(leftMesh->U().rows());
  VECTOR fullVelocity = leftMesh->U() * qDot;

  // DEBUG
  VECTOR fullTemp = fullVelocity;
  fullTemp *= 0;

  for (int x = 0; x < _partitions; x++)
  {
    vector<pair<int,int> > clones = _partitionedMesh->clonedVertices(partition, x);
    if (clones.size() == 0) continue;

    SUBSPACE_TET_MESH* rightMesh = (SUBSPACE_TET_MESH*)_partitionedMesh->mesh(x);
    vector<VEC3F>& rightVertices = rightMesh->vertices();
    rightMesh->updateFullMesh();

    MATRIX3 rightRotation(MATRIX3::I());
    VEC3F rightTranslation;
    rightRotation = _partitionedMesh->quaternionRotation(x).toExplicitMatrix3x3();
    rightTranslation = _partitionedMesh->rigidTranslation(x);

    VECTOR rightVelocity = rightMesh->U() * _integrators[x]->velocity();

    pair<int, int> xy(partition, x);
    
    for (unsigned int y = 0; y < clones.size(); y++)
    {
      int leftIndex = clones[y].first;
      int rightIndex = clones[y].second;

      VEC3F leftVertex = leftVertices[leftIndex];
      VEC3F rightVertex = rightVertices[rightIndex];

      leftVertex = rotation * leftVertex + translation;
      rightVertex = rightRotation * rightVertex + rightTranslation;

      pair<int, int> xy(partition, x);
      VEC3F force = -_interfaceSpringConsts[xy] * (leftVertex - rightVertex);

      VEC3F rightRest = rightMesh->restPose()[rightIndex];
      VEC3F leftRest = leftMesh->restPose()[leftIndex];

      VEC3F temp = -_interfaceSpringConsts[xy] * (leftVertices[leftIndex] - 
                                                  rotation.transpose() * rightRotation * rightVertices[rightIndex] + 
                                                  rotation.transpose() * (translation - rightTranslation));
      fullTemp[3 * leftIndex] += temp[0];
      fullTemp[3 * leftIndex + 1] += temp[1];
      fullTemp[3 * leftIndex + 2] += temp[2];

      fullSprings[3 * leftIndex] += force[0];
      fullSprings[3 * leftIndex + 1] += force[1];
      fullSprings[3 * leftIndex + 2] += force[2];
      
      VEC3F rightVelocityRotated;
      rightVelocityRotated[0] = rightVelocity[3 * rightIndex];
      rightVelocityRotated[1] = rightVelocity[3 * rightIndex + 1];
      rightVelocityRotated[2] = rightVelocity[3 * rightIndex + 2];
      rightVelocityRotated = rotation.transpose() * rightRotation * rightVelocityRotated;

      fullDampers[3 * leftIndex] += _dampingConst * (fullVelocity[3 * leftIndex] - rightVelocityRotated[0]);
      fullDampers[3 * leftIndex + 1] += _dampingConst * (fullVelocity[3 * leftIndex + 1] - rightVelocityRotated[1]);
      fullDampers[3 * leftIndex + 2] += _dampingConst * (fullVelocity[3 * leftIndex + 2] - rightVelocityRotated[2]);

      totalClones++;
    }
  }

  // compute the deformable component
  VECTOR rotatedSprings(fullSprings.size());
  for (int x = 0; x < rotatedSprings.size() / 3; x++)
  {
    VEC3F force;
    force[0] = fullSprings[3 * x];
    force[1] = fullSprings[3 * x + 1];
    force[2] = fullSprings[3 * x + 2];
    force = rotation.transpose() * force;
    rotatedSprings[3 * x] = force[0];
    rotatedSprings[3 * x + 1] = force[1];
    rotatedSprings[3 * x + 2] = force[2];
  }
  VECTOR deformableForce = leftMesh->U() ^ rotatedSprings;

  // new damping model
  VECTOR deformableDamping = leftMesh->U() ^ fullDampers;

  forceVector = deformableForce - deformableDamping;
}

//////////////////////////////////////////////////////////////////////
// add the inter-partition spring forces
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::computeDeformableSpringForcesDefoOnly(int partition, VECTOR& forceVector, bool debug)
{
  UNCONSTRAINED_SUBSPACE_INTEGRATOR* integrator = (UNCONSTRAINED_SUBSPACE_INTEGRATOR*)_integrators[partition];
  UNCONSTRAINED_SUBSPACE_TET_MESH* leftMesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_partitionedMesh->mesh(partition);
  vector<VEC3F>& leftVertices = leftMesh->vertices();
  leftMesh->updateFullMesh();

  // get the updated vars
  VEC3F& translationDotDot = integrator->translationAcceleration();
  VEC3F& angularAcceleration = integrator->angularAcceleration();

  // get the old vars
  VEC3F& translationOld = integrator->translationOld();
  VEC3F& translationDotOld = integrator->translationVelocityOld();
  VEC3F& translationDotDotOld = integrator->translationAccelerationOld();

  QUATERNION& rotationOld = integrator->rotationOld();
  VEC3F& angularVelocityOld = integrator->angularVelocityOld();
  VEC3F& angularAccelerationOld = integrator->angularAccelerationOld();

  VECTOR& qOld = integrator->positionOld();
  VECTOR& qDotOld = integrator->velocityOld();
  VECTOR& qDotDotOld = integrator->accelerationOld();
  VECTOR& qDotDot = integrator->acceleration();

  // compute new angular quantities
  Real* alpha = integrator->accelerationAlpha();
  VEC3F angularDelta = alpha[2] * angularVelocityOld + 
                       alpha[3] * angularAccelerationOld + 
                       alpha[4] * angularAcceleration; 
  VEC3F angularVelocity = angularVelocityOld + 
                          alpha[0] * angularAccelerationOld + 
                          alpha[1] * angularAcceleration;
  QUATERNION update = QUATERNION::fromAxisAngle(angularDelta);
  QUATERNION rotationQuaternion = update * rotationOld;
  MATRIX3 rotation = rotationQuaternion.toExplicitMatrix3x3();
  MATRIX3 rotationTranspose = rotation.transpose();

  // compute new translation quantities
  VEC3F translation = translationOld + alpha[2] * translationDotOld + alpha[3] * translationDotDotOld + alpha[4] * translationDotDot;
  VEC3F translationVelocity = translationDotOld + alpha[0] * translationDotDotOld + alpha[1] * translationDotDot;

  // compute new deformation quantities
  VECTOR qDot = qDotOld + alpha[0] * qDotDotOld + alpha[1] * qDotDot;

  int totalClones = 0;
  VECTOR fullSprings(leftMesh->U().rows());
  VECTOR fullVelocity = leftMesh->U() * qDot;
  for (int x = 0; x < _partitions; x++)
  {
    vector<pair<int,int> > clones = _partitionedMesh->clonedVertices(partition, x);
    if (clones.size() == 0) continue;

    SUBSPACE_TET_MESH* rightMesh = (SUBSPACE_TET_MESH*)_partitionedMesh->mesh(x);
    vector<VEC3F>& rightVertices = rightMesh->vertices();
    rightMesh->updateFullMesh();

    MATRIX3 rightRotation(MATRIX3::I());
    VEC3F rightTranslation;
    VEC3F rightTranslationVelocity;
    VEC3F rightAngularVelocity;
    if (_partitionedMesh->unconstrained(x))
    {
      rightRotation = _partitionedMesh->quaternionRotation(x).toExplicitMatrix3x3();
      rightTranslation = _partitionedMesh->rigidTranslation(x);

      UNCONSTRAINED_SUBSPACE_INTEGRATOR* rightIntegrator = (UNCONSTRAINED_SUBSPACE_INTEGRATOR*)_integrators[x];
      rightIntegrator->updateRigids();
      rightTranslationVelocity = rightIntegrator->translationVelocity();
      rightAngularVelocity = rightIntegrator->angularVelocity();
    }
    VECTOR rightVelocity = rightMesh->U() * _integrators[x]->velocity();

    for (unsigned int y = 0; y < clones.size(); y++)
    {
      int leftIndex = clones[y].first;
      int rightIndex = clones[y].second;

      VEC3F leftLocal = leftVertices[leftIndex];
      VEC3F rightLocal= rightVertices[rightIndex];

      VEC3F leftGlobal = rotation * leftLocal + translation;
      VEC3F rightGlobal = rightRotation * rightLocal + rightTranslation;

      pair<int, int> xy(partition, x);
      VEC3F force = -_defoSpringFactor * _interfaceSpringConsts[xy] * (leftGlobal - rightGlobal);

      fullSprings[3 * leftIndex] += force[0];
      fullSprings[3 * leftIndex + 1] += force[1];
      fullSprings[3 * leftIndex + 2] += force[2];

      // get the right deformation velocity
      VEC3F rightVelocityRotated;
      rightVelocityRotated[0] = rightVelocity[3 * rightIndex];
      rightVelocityRotated[1] = rightVelocity[3 * rightIndex + 1];
      rightVelocityRotated[2] = rightVelocity[3 * rightIndex + 2];
      
      // add the right rotational velocity
      rightVelocityRotated += cross(rightAngularVelocity, rightLocal);

      // rotate it
      rightVelocityRotated = rightRotation * rightVelocityRotated;

      // add the right translation velocity
      rightVelocityRotated += rightTranslationVelocity;

      // get the left deformation velocity
      VEC3F leftVelocityRotated;
      leftVelocityRotated[0] = fullVelocity[3 * leftIndex];
      leftVelocityRotated[1] = fullVelocity[3 * leftIndex + 1];
      leftVelocityRotated[2] = fullVelocity[3 * leftIndex + 2];
      
      // add the left rotational velocity
      leftVelocityRotated += cross(angularVelocity, leftLocal);

      // rotate it
      leftVelocityRotated = rotation * leftVelocityRotated;

      // add the right translation velocity
      leftVelocityRotated += translationVelocity;

      // get the final damping force
      VEC3F dampingForce = -_dampingConst * (leftVelocityRotated - rightVelocityRotated);
      fullSprings[3 * leftIndex] += dampingForce[0];
      fullSprings[3 * leftIndex + 1] += dampingForce[1];
      fullSprings[3 * leftIndex + 2] += dampingForce[2];

      totalClones++;
    }
  }

  // compute the deformable component
  VECTOR rotatedSprings(fullSprings.size());
  for (int x = 0; x < rotatedSprings.size() / 3; x++)
  {
    VEC3F force;
    force[0] = fullSprings[3 * x];
    force[1] = fullSprings[3 * x + 1];
    force[2] = fullSprings[3 * x + 2];
    force = rotation.transpose() * force;
    rotatedSprings[3 * x] = force[0];
    rotatedSprings[3 * x + 1] = force[1];
    rotatedSprings[3 * x + 2] = force[2];
  }
  VECTOR deformableForce = leftMesh->U() ^ rotatedSprings;

  forceVector += deformableForce;
}

//////////////////////////////////////////////////////////////////////
// add the inter-partition spring forces
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::computeDeformableSpringForces(int partition, BLOCK_VECTOR& forceVector, bool debug)
{
  UNCONSTRAINED_SUBSPACE_INTEGRATOR* integrator = (UNCONSTRAINED_SUBSPACE_INTEGRATOR*)_integrators[partition];
  UNCONSTRAINED_SUBSPACE_TET_MESH* leftMesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_partitionedMesh->mesh(partition);
  vector<VEC3F>& leftVertices = leftMesh->vertices();
  leftMesh->updateFullMesh();

  // get the updated vars
  VEC3F& translationDotDot = integrator->translationAcceleration();
  VEC3F& angularAcceleration = integrator->angularAcceleration();

  // get the old vars
  VEC3F& translationOld = integrator->translationOld();
  VEC3F& translationDotOld = integrator->translationVelocityOld();
  VEC3F& translationDotDotOld = integrator->translationAccelerationOld();

  QUATERNION& rotationOld = integrator->rotationOld();
  VEC3F& angularVelocityOld = integrator->angularVelocityOld();
  VEC3F& angularAccelerationOld = integrator->angularAccelerationOld();

  //VECTOR& qOld = integrator->positionOld();
  VECTOR& qDotOld = integrator->velocityOld();
  VECTOR& qDotDotOld = integrator->accelerationOld();
  VECTOR& qDotDot = integrator->acceleration();

  // compute new angular quantities
  Real* alpha = integrator->accelerationAlpha();
  VEC3F angularDelta = alpha[2] * angularVelocityOld + 
                       alpha[3] * angularAccelerationOld + 
                       alpha[4] * angularAcceleration; 
  VEC3F angularVelocity = angularVelocityOld + 
                          alpha[0] * angularAccelerationOld + 
                          alpha[1] * angularAcceleration;
  QUATERNION update = QUATERNION::fromAxisAngle(angularDelta);
  QUATERNION rotationQuaternion = update * rotationOld;
  MATRIX3 rotation = rotationQuaternion.toExplicitMatrix3x3();
  MATRIX3 rotationTranspose = rotation.transpose();

  // compute new translation quantities
  VEC3F translation = translationOld + alpha[2] * translationDotOld + alpha[3] * translationDotDotOld + alpha[4] * translationDotDot;
  VEC3F translationVelocity = translationDotOld + alpha[0] * translationDotDotOld + alpha[1] * translationDotDot;

  // compute new deformation quantities
  VECTOR qDot = qDotOld + alpha[0] * qDotDotOld + alpha[1] * qDotDot;

  int totalClones = 0;
  VECTOR fullSprings(leftMesh->U().rows());
  VECTOR fullDampers(leftMesh->U().rows());
  VECTOR fullVelocity = leftMesh->U() * qDot;
  for (int x = 0; x < _partitions; x++)
  {
    vector<pair<int,int> > clones = _partitionedMesh->clonedVertices(partition, x);
    if (clones.size() == 0) continue;

    SUBSPACE_TET_MESH* rightMesh = (SUBSPACE_TET_MESH*)_partitionedMesh->mesh(x);
    vector<VEC3F>& rightVertices = rightMesh->vertices();
    rightMesh->updateFullMesh();

    MATRIX3 rightRotation(MATRIX3::I());
    VEC3F rightTranslation;
    VEC3F rightTranslationVelocity;
    VEC3F rightAngularVelocity;
    if (_partitionedMesh->unconstrained(x))
    {
      rightRotation = _partitionedMesh->quaternionRotation(x).toExplicitMatrix3x3();
      rightTranslation = _partitionedMesh->rigidTranslation(x);

      UNCONSTRAINED_SUBSPACE_INTEGRATOR* rightIntegrator = (UNCONSTRAINED_SUBSPACE_INTEGRATOR*)_integrators[x];
      rightIntegrator->updateRigids();
      rightTranslationVelocity = rightIntegrator->translationVelocity();
      rightAngularVelocity = rightIntegrator->angularVelocity();
    }
    VECTOR rightVelocity = rightMesh->U() * _integrators[x]->velocity();

    for (unsigned int y = 0; y < clones.size(); y++)
    {
      int leftIndex = clones[y].first;
      int rightIndex = clones[y].second;

      VEC3F leftLocal = leftVertices[leftIndex];
      VEC3F rightLocal= rightVertices[rightIndex];

      VEC3F leftGlobal = rotation * leftLocal + translation;
      VEC3F rightGlobal = rightRotation * rightLocal + rightTranslation;

      pair<int, int> xy(partition, x);
      VEC3F force = -_interfaceSpringConsts[xy] * (leftGlobal - rightGlobal);

      fullSprings[3 * leftIndex] += force[0];
      fullSprings[3 * leftIndex + 1] += force[1];
      fullSprings[3 * leftIndex + 2] += force[2];

      // get the right deformation velocity
      VEC3F rightVelocityRotated;
      rightVelocityRotated[0] = rightVelocity[3 * rightIndex];
      rightVelocityRotated[1] = rightVelocity[3 * rightIndex + 1];
      rightVelocityRotated[2] = rightVelocity[3 * rightIndex + 2];
      
      // add the right rotational velocity
      rightVelocityRotated += cross(rightAngularVelocity, rightLocal);

      // rotate it
      rightVelocityRotated = rightRotation * rightVelocityRotated;

      // add the right translation velocity
      rightVelocityRotated += rightTranslationVelocity;

      // get the left deformation velocity
      VEC3F leftVelocityRotated;
      leftVelocityRotated[0] = fullVelocity[3 * leftIndex];
      leftVelocityRotated[1] = fullVelocity[3 * leftIndex + 1];
      leftVelocityRotated[2] = fullVelocity[3 * leftIndex + 2];
      
      // add the left rotational velocity
      leftVelocityRotated += cross(angularVelocity, leftLocal);

      // rotate it
      leftVelocityRotated = rotation * leftVelocityRotated;

      // add the right translation velocity
      leftVelocityRotated += translationVelocity;

      // get the final damping force
      VEC3F dampingForce = -_dampingConst * (leftVelocityRotated - rightVelocityRotated);
      fullSprings[3 * leftIndex] += dampingForce[0];
      fullSprings[3 * leftIndex + 1] += dampingForce[1];
      fullSprings[3 * leftIndex + 2] += dampingForce[2];

      totalClones++;
    }
  }

  // compute the deformable component
  VECTOR rotatedSprings(fullSprings.size());
  for (int x = 0; x < rotatedSprings.size() / 3; x++)
  {
    VEC3F force;
    force[0] = fullSprings[3 * x];
    force[1] = fullSprings[3 * x + 1];
    force[2] = fullSprings[3 * x + 2];
    force = rotation.transpose() * force;
    rotatedSprings[3 * x] = force[0];
    rotatedSprings[3 * x + 1] = force[1];
    rotatedSprings[3 * x + 2] = force[2];
  }
  VECTOR deformableForce = leftMesh->U() ^ rotatedSprings;

  // compute the translation component
  VECTOR meanForce(3);
  for (int x = 0; x < fullSprings.size() / 3; x++)
  {
    meanForce[0] += fullSprings[3 * x];
    meanForce[1] += fullSprings[3 * x + 1];
    meanForce[2] += fullSprings[3 * x + 2];
  }

  // compute the rotation component
  VEC3F SitBar = leftMesh->SitBar();

  VECTOR rotationForce(3);
  for (unsigned int x = 0; x < leftVertices.size(); x++)
  {
    VEC3F& vertex = leftVertices[x];
    VEC3F centered = vertex - SitBar;
    MATRIX uBarTilde = MATRIX::cross(centered);

    MATRIX AuBar = rotation * uBarTilde;
    MATRIX AuBarT = AuBar.transpose();

    VECTOR F(3);
    F(0) = fullSprings[3 * x];
    F(1) = fullSprings[3 * x + 1];
    F(2) = fullSprings[3 * x + 2];

    rotationForce += -1.0 * AuBarT * F;
  }
  
  forceVector[0] += meanForce;
  forceVector[1] += rotationForce;
  forceVector[2] += deformableForce;

  _rotationSpringForces[partition] = rotationForce;
}

//////////////////////////////////////////////////////////////////////
// compute the coupling jacobian between a constrained and
// unconstrained partition
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::computeSkinnedCouplingJacobianReduced(int left, int right, MATRIX& springMatrix)
{
  UNCONSTRAINED_SUBSPACE_INTEGRATOR* leftIntegrator = (UNCONSTRAINED_SUBSPACE_INTEGRATOR*)_integrators[left];
  UNCONSTRAINED_SUBSPACE_TET_MESH* leftMesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_partitionedMesh->mesh(left);
  UNCONSTRAINED_SUBSPACE_INTEGRATOR* rightIntegrator = (UNCONSTRAINED_SUBSPACE_INTEGRATOR*)_integrators[right];
  UNCONSTRAINED_SUBSPACE_TET_MESH* rightMesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_partitionedMesh->mesh(right);

  pair<int, int> xy(left, right);
  int leftRank = leftMesh->rank();
  int rightRank = rightMesh->rank();
  Real springConst = _interfaceSpringConsts[xy];

  MATRIX3 leftRotation = leftIntegrator->rotation().toExplicitMatrix3x3();
  MATRIX3 leftRotationT = leftRotation.transpose();
  MATRIX3 rightRotation = rightIntegrator->rotation().toExplicitMatrix3x3();
  
  MATRIX RTR = leftRotationT * rightRotation;
  MATRIX reducedDefoPartialDefo(leftRank, rightRank);
  _Ui_R_Uk[xy].transformInplace(RTR, reducedDefoPartialDefo);
  Real* alpha = leftIntegrator->alpha();
  reducedDefoPartialDefo *= (springConst + _dampingConst * alpha[3]);
  springMatrix = reducedDefoPartialDefo;
}

//////////////////////////////////////////////////////////////////////
// compute the coupling jacobian between a constrained and
// unconstrained partition
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::computeSkinnedCouplingJacobian(int left, int right, MATRIX& springMatrix)
{
  UNCONSTRAINED_SUBSPACE_INTEGRATOR* leftIntegrator = (UNCONSTRAINED_SUBSPACE_INTEGRATOR*)_integrators[left];
  UNCONSTRAINED_SUBSPACE_TET_MESH* leftMesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_partitionedMesh->mesh(left);
  UNCONSTRAINED_SUBSPACE_INTEGRATOR* rightIntegrator = (UNCONSTRAINED_SUBSPACE_INTEGRATOR*)_integrators[right];
  UNCONSTRAINED_SUBSPACE_TET_MESH* rightMesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_partitionedMesh->mesh(right);

  MATRIX3 leftRotation = leftIntegrator->rotation().toExplicitMatrix3x3();
  MATRIX3 rightRotation = rightIntegrator->rotation().toExplicitMatrix3x3();

  int leftRank = leftMesh->rank();
  int rightRank = rightMesh->rank();

  vector<VEC3F>& leftVertices = leftMesh->vertices();
  vector<VEC3F>& rightVertices = rightMesh->vertices();

  MATRIX defoPartialDefo(leftRank, rightRank);

  // get both partials at the same time
  vector<pair<int,int> > clones = _partitionedMesh->clonedVertices(left, right);

  // get all the spring forces
  for (unsigned int y = 0; y < clones.size(); y++)
  {
    int leftIndex = clones[y].first;
    int rightIndex = clones[y].second;

    VEC3F leftVertex = leftVertices[leftIndex];
    VEC3F rightVertex = rightVertices[rightIndex];

    pair<int, int> xy(left, right);

    SUBMATRIX Ukj(rightMesh->U(), 3 * rightIndex, 3);

    // defo terms
    SUBMATRIX Uij(leftMesh->U(), 3 * leftIndex, 3);
    MATRIX rotatedLeftBasis = leftRotation * Uij;
    MATRIX rotatedLeftBasisT = rotatedLeftBasis.transpose();

    // new damping model
    Real* alpha = leftIntegrator->alpha();
    defoPartialDefo += (_interfaceSpringConsts[xy] + _dampingConst * alpha[3]) * rotatedLeftBasisT * (rightRotation * Ukj);
  }

  springMatrix = defoPartialDefo;
}

//////////////////////////////////////////////////////////////////////
// Get the jacobian of the spring forces for all of the blocks
// in the rigid body mass matrix
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::computeSkinnedSpringJacobiansReduced(int partition, MATRIX& systemMatrix, bool dynamic)
{
  UNCONSTRAINED_SUBSPACE_INTEGRATOR* integrator = (UNCONSTRAINED_SUBSPACE_INTEGRATOR*)_integrators[partition];
  UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_partitionedMesh->mesh(partition);

  Real* alpha = integrator->alpha();
  int rank = mesh->rank();

  // deformation terms
  MATRIX reducedDefoPartialDefo(rank, rank);

  reducedDefoPartialDefo = -1.0 * _spring_Ui_Ui[partition];
  if (dynamic)
    reducedDefoPartialDefo.axpy((-_dampingConst * alpha[3]), _Ui_Ui[partition]);

  systemMatrix = reducedDefoPartialDefo;
  systemMatrix *= -1;
}

//////////////////////////////////////////////////////////////////////
// Get the jacobian of the spring forces for all of the blocks
// in the rigid body mass matrix
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::computeSkinnedSpringJacobians(int partition, MATRIX& systemMatrix, bool debug, bool dynamic)
{
  UNCONSTRAINED_SUBSPACE_INTEGRATOR* integrator = (UNCONSTRAINED_SUBSPACE_INTEGRATOR*)_integrators[partition];
  UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_partitionedMesh->mesh(partition);

  VECTOR qDot = integrator->velocity();
  VECTOR fullVelocity = mesh->U() * qDot;
  MATRIX3 leftRotation = integrator->rotation().toExplicitMatrix3x3();

  int rank = mesh->rank();
  MATRIX defoPartialDefo(rank,rank);

  // get both partials at the same time
  int totalClones = 0;
  for (int x = 0; x < _partitions; x++)
  {
    vector<pair<int,int> > clones = _partitionedMesh->clonedVertices(partition, x);
    if (clones.size() == 0) continue;

    MATRIX3 rightRotation(MATRIX3::I());
    VEC3F rightTranslation;
    if (_partitionedMesh->unconstrained(x))
    {
      rightRotation = _partitionedMesh->quaternionRotation(x).toExplicitMatrix3x3();
      rightTranslation = _partitionedMesh->rigidTranslation(x);
    }

    // get all the spring forces
    for (unsigned int y = 0; y < clones.size(); y++)
    {
      int leftIndex = clones[y].first;

      pair<int, int> xy(partition, x);

      // compute translation partial defo, and the transpose
      SUBMATRIX Uij(mesh->U(), 3 * leftIndex, 3);
      if (dynamic)
      {
        Real* alpha = integrator->alpha();
        defoPartialDefo += (-_interfaceSpringConsts[xy] - _dampingConst * alpha[3]) * Uij ^ Uij;
      }
      else
        defoPartialDefo += (-_interfaceSpringConsts[xy]) * Uij ^ Uij;

      totalClones++;
    }
  }

  systemMatrix = defoPartialDefo;
  systemMatrix *= -1;
}

//////////////////////////////////////////////////////////////////////
// process a mouse click event
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::click(VEC3F& point, int forcePartition)
{
  // find the closest partition
  if (forcePartition == -1)
    _clickedPartition = _partitionedMesh->closestPartition(point);
  else
    _clickedPartition = forcePartition;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " Clicked partition: " << _clickedPartition << endl;
  cout << " Point: " << point << endl;

  // send the click to that partition
  _integrators[_clickedPartition]->click(point);
}

//////////////////////////////////////////////////////////////////////
// Unclick mouse support
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::unclick()
{ 
  if (_clickedPartition == -1) return;
  _integrators[_clickedPartition]->unclick();
  _clickedPartition = -1; 
}

//////////////////////////////////////////////////////////////////////
// Drag mouse support
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::drag(VEC3F& point)
{
  _integrators[_clickedPartition]->drag(point);

  // if the partition is unconstrained, add the rigid transform
  if (_partitionedMesh->unconstrained(_clickedPartition))
  {
    VEC3F& draggedPosition = _integrators[_clickedPartition]->draggedPosition();
    MATRIX3 rotation = ((UNCONSTRAINED_SUBSPACE_TET_MESH*)_partitionedMesh->mesh(_clickedPartition))->rotationQuaternion().toExplicitMatrix3x3();
    VEC3F translation = ((UNCONSTRAINED_SUBSPACE_TET_MESH*)_partitionedMesh->mesh(_clickedPartition))->rigidTranslation();

    draggedPosition -= translation;
    draggedPosition = rotation.transpose() * draggedPosition;
  }
}

//////////////////////////////////////////////////////////////////////
// draw the force vector -- if the click and unclick were done
// properly, a simple passthrough should work
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::drawForceVector()
{
  if (_clickedPartition != -1)
    if (!_partitionedMesh->unconstrained(_clickedPartition))
      _integrators[_clickedPartition]->drawForceVector();
    else
    {
      // apply the rigid component
      VEC3F& clickedPosition = _integrators[_clickedPartition]->clickedPosition();
      VEC3F& dragPosition = _integrators[_clickedPartition]->draggedPosition();
      VEC3F clickBackup = clickedPosition;
      VEC3F dragBackup = dragPosition;
      clickedPosition = _partitionedMesh->rigidRotation(_clickedPartition) * clickedPosition + 
                        _partitionedMesh->rigidTranslation(_clickedPartition);
      dragPosition = _partitionedMesh->rigidRotation(_clickedPartition) * dragPosition + 
                     _partitionedMesh->rigidTranslation(_clickedPartition);
      
      // draw the vector
      _integrators[_clickedPartition]->drawForceVector();

      // undo the rigid component;
      clickedPosition = clickBackup;
      dragPosition = dragBackup;
    }
}

//////////////////////////////////////////////////////////////////////
// compute the quadratic velocity vectors in implicit form
//
// make sure updateFullMesh(qOld) was not called previously
//
// PRELIMINARY -- DO NOT USE
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::computeQuadraticsImplicit(int partition, BLOCK_VECTOR& quadratics)
{
  UNCONSTRAINED_SUBSPACE_INTEGRATOR* integrator = (UNCONSTRAINED_SUBSPACE_INTEGRATOR*)_integrators[partition];
  UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_partitionedMesh->mesh(partition);
  VEC3F& angularVelocityOld = integrator->angularVelocityOld();
  VEC3F& angularAccelerationOld = integrator->angularAccelerationOld();
  VEC3F& angularAcceleration = integrator->angularAcceleration();

  VECTOR& qDotDot = integrator->acceleration();
  VECTOR& qDotOld = integrator->velocityOld();
  VECTOR& qDotDotOld = integrator->accelerationOld();

  Real* alpha = integrator->accelerationAlpha();
  VEC3F angularVelocity = angularVelocityOld + 
                          alpha[0] * angularAcceleration + 
                          alpha[1] * angularAccelerationOld;
  VECTOR qDot = qDotOld + 
                alpha[0] * qDotDot + 
                alpha[1] * qDotDotOld;

  /////////////////////////////////////////////////////////////////////////
  // Add angular quadratics
  /////////////////////////////////////////////////////////////////////////

  // First term: -\tilde{\bar{\omega}}^i (\bar \II^i_{\theta \theta} \bar\omega^i)
  // get the quadratic term (seems to be near zero for now) -- implicit version
  MATRIX3 omegaTilde = MATRIX3::cross(angularVelocity);
  VECTOR quadratic = omegaTilde * mesh->inertiaTensor() * angularVelocity.toVector();
  VECTOR QTheta = -1.0 * quadratic;

  // second term: {\dot{\bar{\II}}^i_{\theta \theta}} \bar\omega^i
  quadratic = mesh->inertiaTensorDt() * angularVelocity.toVector();
  QTheta -= quadratic;

  // third term: \tilde{\bar \omega}({\bar\II}^i_{\theta f} \dot \qq^i_f)
  MATRIX& U = mesh->U();
  vector<VEC3F>& vertices = mesh->vertices();
  MATRIX Ithetaf(3, U.cols());
  for (int x = 0; x < U.rows() / 3; x++)
  {
    Real mass = mesh->mass(x);
    MATRIX uBarTilde = MATRIX::cross(vertices[x]);
    SUBMATRIX submatrix(U, 3 * x, 3);

    Ithetaf += mass * uBarTilde * submatrix;
  }
  quadratic = omegaTilde * Ithetaf * qDot;
  QTheta -= quadratic;

  quadratics.clear();
  quadratics.add(QTheta, 1);

  /////////////////////////////////////////////////////////////////////////
  // Add translation quadratics
  /////////////////////////////////////////////////////////////////////////
  VEC3F SitBar = mesh->SitBar();
  MATRIX& SiBar = mesh->SiBar();
  MATRIX3 omegaTilde2 = omegaTilde * omegaTilde;

  VEC3F translationFinal = omegaTilde2 * SitBar;
  translationFinal += 2.0 * omegaTilde * SiBar * qDot;
  MATRIX3 A = integrator->rotation().toExplicitMatrix3x3();
  translationFinal = -A * translationFinal;
  quadratics.add(translationFinal.toVector(), 0);

  /////////////////////////////////////////////////////////////////////////
  // Add deformation quadratics
  /////////////////////////////////////////////////////////////////////////
  VECTOR fullVelocity = U * qDot;
  VECTOR deformationFinal(U.cols());
  for (int x = 0; x < U.rows() / 3; x++)
  {
    Real mass = mesh->mass(x);
    VEC3F uBar = vertices[x];
    VEC3F uDotBar;
    uDotBar[0] = fullVelocity[3 * x];
    uDotBar[1] = fullVelocity[3 * x + 1];
    uDotBar[2] = fullVelocity[3 * x + 2];

    VEC3F sum = omegaTilde2 * uBar + 2.0 * omegaTilde2  * uDotBar;
    SUBMATRIX submatrix(U, 3 * x, 3);
    deformationFinal -= mass * submatrix ^ sum.toVector();
  }
  quadratics.add(deformationFinal, 2);
}

//////////////////////////////////////////////////////////////////////
// compute the quadratic velocity vectors
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::computeQuadraticsExplicit(int partition, BLOCK_VECTOR& quadratics)
{
  UNCONSTRAINED_SUBSPACE_INTEGRATOR* integrator = (UNCONSTRAINED_SUBSPACE_INTEGRATOR*)_integrators[partition];
  UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_partitionedMesh->mesh(partition);
  VEC3F& angularVelocityOld = integrator->angularVelocityOld();

  VECTOR& qDotOld = integrator->velocityOld();

  /////////////////////////////////////////////////////////////////////////
  // Add angular quadratics
  /////////////////////////////////////////////////////////////////////////

  // First term: -\tilde{\bar{\omega}}^i (\bar \II^i_{\theta \theta} \bar\omega^i)
  MATRIX3 omegaTilde = MATRIX3::cross(angularVelocityOld);
  VECTOR quadratic = omegaTilde * mesh->inertiaTensorOld() * angularVelocityOld.toVector();
  VECTOR QTheta = -1.0 * quadratic;

  // second term: {\dot{\bar{\II}}^i_{\theta \theta}} \bar\omega^i 
  quadratic = mesh->inertiaTensorDtOld() * angularVelocityOld.toVector();
  QTheta -= quadratic;

  // third term: \tilde{\bar \omega}({\bar\II}^i_{\theta f} \dot \qq^i_f)
  MATRIX& U = mesh->U();
  vector<VEC3F>& vertices = mesh->vertices();
  MATRIX Ithetaf(3, U.cols());
  for (int x = 0; x < U.rows() / 3; x++)
  {
    Real mass = mesh->mass(x);
    MATRIX uBarTilde = MATRIX::cross(vertices[x]);
    SUBMATRIX submatrix(U, 3 * x, 3);

    Ithetaf += mass * uBarTilde * submatrix;
  }
  quadratic = omegaTilde * Ithetaf * qDotOld;
  QTheta -= quadratic;

  quadratics.clear();
  quadratics.add(QTheta, 1);

  /////////////////////////////////////////////////////////////////////////
  // Add translation quadratics
  /////////////////////////////////////////////////////////////////////////
  VEC3F SitBar = mesh->SitBar();
  MATRIX& SiBar = mesh->SiBar();
  MATRIX3 omegaTilde2 = omegaTilde * omegaTilde;

  VEC3F translationFinal = omegaTilde2 * SitBar;
  translationFinal += 2.0 * omegaTilde * SiBar * qDotOld;
  MATRIX3 A = integrator->rotationOld().toExplicitMatrix3x3();
  translationFinal = -A * translationFinal;
  quadratics.add(translationFinal.toVector(), 0);

  /////////////////////////////////////////////////////////////////////////
  // Add deformation quadratics
  /////////////////////////////////////////////////////////////////////////
  VECTOR fullVelocity = U * qDotOld;
  VECTOR deformationFinal(U.cols());
  for (int x = 0; x < U.rows() / 3; x++)
  {
    Real mass = mesh->mass(x);
    VEC3F uBar = vertices[x];
    VEC3F uDotBar;
    uDotBar[0] = fullVelocity[3 * x];
    uDotBar[1] = fullVelocity[3 * x + 1];
    uDotBar[2] = fullVelocity[3 * x + 2];

    VEC3F sum = omegaTilde2 * uBar + 2.0 * omegaTilde2  * uDotBar;
    SUBMATRIX submatrix(U, 3 * x, 3);
    deformationFinal -= mass * submatrix ^ sum.toVector();
  }
  quadratics.add(deformationFinal, 2);
}

//////////////////////////////////////////////////////////////////////
// compute the quadratic velocity vectors in implicit form
// PRELIMINARY -- DO NOT USE
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::computeQuadraticsJacobians(int partition, BLOCK_MATRIX& jacobian)
{
  UNCONSTRAINED_SUBSPACE_INTEGRATOR* integrator = (UNCONSTRAINED_SUBSPACE_INTEGRATOR*)_integrators[partition];
  UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_partitionedMesh->mesh(partition);
  VEC3F& angularVelocityOld = integrator->angularVelocityOld();
  VEC3F& angularAccelerationOld = integrator->angularAccelerationOld();
  VEC3F& angularAcceleration = integrator->angularAcceleration();

  VECTOR& qDotDot = integrator->acceleration();
  VECTOR& qDotOld = integrator->velocityOld();
  VECTOR& qDotDotOld = integrator->accelerationOld();

  Real* alpha = integrator->accelerationAlpha();
  VEC3F angularVelocity = angularVelocityOld + 
                          alpha[0] * angularAcceleration + 
                          alpha[1] * angularAccelerationOld;
  VECTOR qDot = qDotOld + 
                alpha[0] * qDotDot + 
                alpha[1] * qDotDotOld;

  // get the interia tensor in anticipation of computing rotation A
  const MATRIX& Ithetatheta = mesh->inertiaTensor();

  /////////////////////////////////////////////////////////////////////////
  // Compute rotation jacobians (rigid only!)
  /////////////////////////////////////////////////////////////////////////
  VEC3F omegaP = angularVelocityOld + _dt * _gamma * angularAccelerationOld;
  MATRIX omegaPTilde = MATRIX::cross(omegaP);
  VEC3F IomegaP = Ithetatheta * omegaP.toVector();
  MATRIX IomegaPTilde = MATRIX::cross(IomegaP);

  VEC3F Ialpha = Ithetatheta * angularAcceleration.toVector();
  MATRIX alphaTilde = MATRIX::cross(angularAcceleration);
  MATRIX IalphaTilde = MATRIX::cross(Ialpha);

  MATRIX term1 = (omegaPTilde * Ithetatheta - IomegaPTilde);
  term1 *= 0.5 * _dt;

  // typo in Krysl and Endres? Should work out to 0.25, not 0.5.
  MATRIX term2 = (alphaTilde * Ithetatheta - IalphaTilde);
  term2 *= 0.25 * _dt * _dt;

  jacobian.add(term1, 1,1);
  jacobian.add(term2, 1,1);

  /////////////////////////////////////////////////////////////////////////
  // Compute translation jacobians (rigid only!)
  /////////////////////////////////////////////////////////////////////////
  MATRIX3 omegaTilde = MATRIX3::cross(angularVelocity);
  // build the matrix derivatives
  vector<MATRIX> partials;
  for (int x = 0; x < 3; x++)
  {
    // depending on the component, the derivative of the
    // cross product (tilde) operator
    MATRIX partial(3,3);
    Real halfDt2 = alpha[4];
    
    if (x == 0)
    {
      partial(1,2) = -halfDt2;
      partial(2,1) = halfDt2;
    }
    if (x == 1)
    {
      partial(0,2) = halfDt2;
      partial(2,0) = -halfDt2;
    }
    if (x == 2)
    {
      partial(0,1) = -halfDt2;
      partial(1,0) = halfDt2;
    }
    partials.push_back(partial);
  }

  // build the derivative of the matrix exponential
  VEC3F angularDelta = alpha[2] * angularVelocityOld + 
                       alpha[3] * angularAccelerationOld + 
                       alpha[4] * angularAcceleration;
  QUATERNION update = QUATERNION::fromAxisAngle(angularDelta);
  MATRIX deltaR = update.toExplicitMatrix3x3();
  TENSOR3 partialDelta(deltaR, partials);

  QUATERNION oldRotation = integrator->rotationOld();
  QUATERNION newRotation = update * integrator->rotationOld();
  VEC3F SitBar = mesh->SitBar();

  VEC3F eta(alpha[1], alpha[1], alpha[1]);
  MATRIX etaTilde = MATRIX::cross(eta);

  // first partial of exponential term
  VEC3F temp = oldRotation.toExplicitMatrix3x3() * omegaTilde * omegaTilde * SitBar;
  VECTOR RHS = temp.toVector();
  MATRIX final = partialDelta.modeOneProduct(RHS);

  // compute partial of cross product terms
  TENSOR3 omegaPartialAlpha(3,3,3);
  omegaPartialAlpha(1,2,0) = -alpha[1];
  omegaPartialAlpha(2,1,0) = alpha[1];
  omegaPartialAlpha(0,2,1) = alpha[1];
  omegaPartialAlpha(2,0,1) = -alpha[1];
  omegaPartialAlpha(0,1,2) = -alpha[1];
  omegaPartialAlpha(1,0,2) = alpha[1];

  // compute with partial first
  temp = omegaTilde * SitBar;
  RHS = temp.toVector();
  MATRIX matrixRHS = omegaPartialAlpha.modeOneProduct(RHS);
  final += newRotation.toExplicitMatrix3x3() * matrixRHS;

  // compute with partial second
  VECTOR vSitBar = SitBar.toVector();
  matrixRHS = omegaPartialAlpha.modeOneProduct(vSitBar);
  final += newRotation.toExplicitMatrix3x3() * omegaTilde * matrixRHS;

  final *= -1.0;
  jacobian.add(final, 0,1);
}

//////////////////////////////////////////////////////////////////////
// compute the current kinetic energy of a body
//////////////////////////////////////////////////////////////////////
Real SUBSPACE_MULTIBODY_INTEGRATOR::kineticEnergy(BLOCK_MATRIX& mass, UNCONSTRAINED_SUBSPACE_INTEGRATOR* integrator)
{
  VEC3F translationDot = integrator->translationVelocity();
  VEC3F angularVelocity = integrator->angularVelocity();
  VECTOR deformationVelocity = integrator->velocity();

  BLOCK_VECTOR velocity(3);
  velocity.add(translationDot.toVector(), 0);
  velocity.add(angularVelocity.toVector(), 1);
  velocity.add(deformationVelocity, 2);

  BLOCK_VECTOR product = mass * velocity;
  Real final = 0.5 * (velocity ^ product);

  return final;
}

//////////////////////////////////////////////////////////////////////
// compute the system matrix and residual for a skinned unconstrained 
// partition
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::computeSkinnedDiagonalReduced(int partition, MATRIX& system, VECTOR& residual, bool dynamic)
{
  UNCONSTRAINED_SUBSPACE_INTEGRATOR* currentIntegrator = (UNCONSTRAINED_SUBSPACE_INTEGRATOR*)_integrators[partition];

  /////////////////////////////////////////////////////////////////////////
  // Assemble rotation b
  /////////////////////////////////////////////////////////////////////////

  // get the spring forces
  VECTOR forceVector;

  // make a skinned version
  computeSkinnedSpringForcesReduced(partition, forceVector);

  /////////////////////////////////////////////////////////////////////////
  // Assemble deformation b
  /////////////////////////////////////////////////////////////////////////
  residual = currentIntegrator->b();
  residual -= forceVector;

  if (!dynamic) return;

  // add in the quadratic forces
  BLOCK_VECTOR& quadratics = currentIntegrator->quadraticForces();
  if (quadratics.totalBlocks() == 3 && quadratics[2].size() == residual.size())
    residual -= quadratics[2];
}

//////////////////////////////////////////////////////////////////////
// compute the system matrix and residual for a skinned unconstrained 
// partition
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::computeSkinnedDiagonal(int partition, MATRIX& system, VECTOR& residual, bool dynamic)
{
  UNCONSTRAINED_SUBSPACE_INTEGRATOR* currentIntegrator = (UNCONSTRAINED_SUBSPACE_INTEGRATOR*)_integrators[partition];

  /////////////////////////////////////////////////////////////////////////
  // Assemble rotation b
  /////////////////////////////////////////////////////////////////////////

  // get the spring forces
  VECTOR forceVector;

  // DONE: TODO: make a skinned version
  computeSkinnedSpringForces(partition, forceVector);

  /////////////////////////////////////////////////////////////////////////
  // Assemble deformation b
  /////////////////////////////////////////////////////////////////////////
  residual = currentIntegrator->b();
  residual -= forceVector;
}

//////////////////////////////////////////////////////////////////////
// helper for verifyGlobalSpringJacobian -- compute the current 
// spring force
//////////////////////////////////////////////////////////////////////
BLOCK_VECTOR SUBSPACE_MULTIBODY_INTEGRATOR::computeCurrentSpringForce()
{
  int totalBlocks = 0;
  for (int x = 0; x < _partitions; x++)
    if (_partitionedMesh->constrained(x))
      totalBlocks += 1;
    else
      totalBlocks += 3;

  // build a list of which blocks start where
  BLOCK_VECTOR springForce(totalBlocks);
  for (int x = 0, i = 0; x < _partitions; x++)
  {
    if (_partitionedMesh->constrained(x))
    {
      VECTOR springForces;
      computeConstrainedSpringForces(x, springForces);
      springForce.add(springForces, i);
      i++;
    }
    else
    {
      vector<int> blockSizes;
      blockSizes.push_back(3);
      blockSizes.push_back(3);
      blockSizes.push_back(_partitionedMesh->rank(x));
      BLOCK_VECTOR forceVectors(blockSizes);
      computeDeformableSpringForces(x, forceVectors);
      springForce.add(forceVectors[0], i++);
      springForce.add(forceVectors[1], i++);
      springForce.add(forceVectors[2], i++);
    }
  }
  return springForce;
}

//////////////////////////////////////////////////////////////////////
// Print out a detailed timing breakdown
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::printTimingBreakdown()
{
  // create an inverse map so that it will sort by time
  map<double, string> inverseMap;
  map<string, double>::iterator forwardIter;
  for (forwardIter = _timingBreakdown.begin(); forwardIter != _timingBreakdown.end(); forwardIter++)
    inverseMap[forwardIter->second] = forwardIter->first;

  // print the map out backwards since it sorts from least to greatest
  cout << "SUBSPACE_MULTIBODY_INTEGRATOR TIMING BREAKDOWN: " << endl;
  cout << "==============================================================================================" << endl;
  map<double,string>::reverse_iterator backwardIter;
  double totalSeen = 0.0;
  char buffer[256];
  for (backwardIter = inverseMap.rbegin(); backwardIter != inverseMap.rend(); backwardIter++)
  {
    string name = (*backwardIter).second;
    name = name.substr(0,45);
    while (name.size() < 45)
      name = name + string(" ");

    sprintf(buffer, "%f%%", (*backwardIter).first / _totalTime * 100.0);
    string percent(buffer);
    while (percent.size() < 12)
      percent = percent + string(" ");

    //cout << "[" << (*backwardIter).first / _totalTime * 100.0 << "%\t]: "
    cout << "[" << percent.c_str() << "]: "
         << name.c_str() << "\t" << (*backwardIter).first / _totalSteps << "s / frame" << endl;
    totalSeen += (*backwardIter).first;
  }
  Real misc = (_totalTime - totalSeen) / _totalTime * 100.0;
  sprintf(buffer, "%f%%", misc);
  string percent(buffer);
  while (percent.size() < 12)
    percent = percent + string(" ");
  cout << "[" << percent.c_str() << "]: " << "Misc. " << endl;
  cout << "==============================================================================================" << endl;
  cout << " Current FPS: " << _totalSteps / _totalTime << endl;
  cout << " Newton steps per second: " << _totalNewtonStepsSeen / _totalTime << endl;
  cout << " Mean Newton steps: " << (Real)_totalNewtonStepsSeen / (Real)_totalSteps << endl;
  cout << " Total Newton stalls: " << _newtonStalls << endl;
  cout << "==============================================================================================" << endl;
  if (misc < 0.0)
  {
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cout << " BREAKDOWN ADDS UP TO MORE THAN 100! TIMERS ARE OVERLAPPING! " << endl;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  }
}

//////////////////////////////////////////////////////////////////////
// precompute the fast sandwich transforms
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::precomputeSandwiches()
{
  cout << " Precomputing sandwiches ... "; flush(cout);

  // compute sandwiches for diagonal blocks of unconstrained partitions
  for (int x = 0; x < _partitions; x++)
  {
    // compute the common unconstrained and constrained stuff
    int rank = _partitionedMesh->rank(x);
    MATRIX UiUi(rank, rank);
    for (int y = 0; y < _partitions; y++)
    {
      vector<pair<int,int> > clones = _partitionedMesh->clonedVertices(x,y);
      if (clones.size() == 0) continue;

      MATRIX& Ui = _partitionedMesh->interfaceU(x,y);
      MATRIX& Uk = _partitionedMesh->interfaceU(y,x);
      MATRIX I3nx3 = MATRIX::columnOfIdentities(clones.size());
      UiUi += Ui ^ Ui;

      pair<int,int> xy(x,y);

      _Ui_I3nx3[xy] = Ui ^ I3nx3;

      SANDWICH_TRANSFORM I3x3n_R_Ukj(I3nx3, Uk);
      SANDWICH_TRANSFORM I3x3n_R_Uij(I3nx3, Ui);
      VECTOR uiBarO = _partitionedMesh->interfaceRests(x,y);
      VECTOR ukBarO = _partitionedMesh->interfaceRests(y,x);
      SANDWICH_TRANSFORM I3x3n_R_uiBarO(I3nx3, uiBarO);
      SANDWICH_TRANSFORM I3x3n_R_ukBarO(I3nx3, ukBarO);

      _I3x3n_R_uiBarO[xy] = I3x3n_R_uiBarO;
      _I3x3n_R_ukBarO[xy] = I3x3n_R_ukBarO;
      _I3x3n_R_Ukj[xy] = I3x3n_R_Ukj;
      _I3x3n_R_Uij[xy] = I3x3n_R_Uij;
      if (_partitionedMesh->constrained(x))
      {
        _Ui_uiBarO[xy] = Ui ^ uiBarO;

        if (_partitionedMesh->constrained(y))
        {
          _Ui_Uk[xy] = Ui ^ Uk;
          _Ui_ukBarO[xy] = Ui ^ ukBarO;
        }
        else
        {
          SANDWICH_TRANSFORM Ui_R_ukBarO(Ui, ukBarO);
          SANDWICH_TRANSFORM Ui_R_Uk(Ui, Uk);
          _Ui_R_ukBarO[xy] = Ui_R_ukBarO;
          _Ui_R_Uk[xy] = Ui_R_Uk;
        }
      }
    }
    _Ui_Ui[x] = UiUi;

    // do unconstrained only stuff
    if (_partitionedMesh->constrained(x)) continue;

    for (int y = 0; y < _partitions; y++)
    {
      vector<pair<int,int> > clones = _partitionedMesh->clonedVertices(x,y);
      if (clones.size() == 0) continue;

      // everybody uses these
      MATRIX I3x3n = MATRIX::columnOfIdentities(clones.size());
      MATRIX& Uij = _partitionedMesh->interfaceU(x,y);

      // build the rest pose term of translation partial angular
      VECTOR uRest = _partitionedMesh->interfaceRestVertices(x,y);
      SANDWICH_TRANSFORM I3x3n_T_restU(I3x3n, uRest);

      pair<int,int> xy(x,y);
      _I3x3n_T_restU[xy] = I3x3n_T_restU;

      // build translation partial defo sandwich
      MATRIX tildeUo = _partitionedMesh->interfaceRestTildeU(x,y);
      MATRIX I3nx3 = MATRIX::columnOfIdentities(clones.size());
      VECTOR uiBarO = _partitionedMesh->interfaceRests(x,y);
      VECTOR ukBarO = _partitionedMesh->interfaceRests(y,x);
      MATRIX& Ukj = _partitionedMesh->interfaceU(y,x);
      TENSOR3 tildeU = _partitionedMesh->interfaceTildeU(x,y);

      // angular terms
      SANDWICH_TRANSFORM tildeUo_R_I3nx3(tildeUo, I3nx3);
      SANDWICH_TRANSFORM tildeUo_R_uiBarO(tildeUo, uiBarO);
      SANDWICH_TRANSFORM tildeUo_R_Uij(tildeUo, Uij);
      SANDWICH_TRANSFORM tildeUo_R_ukBarO(tildeUo, ukBarO);
      SANDWICH_TRANSFORM tildeUo_R_Ukj(tildeUo, Ukj);
      SANDWICH_TRANSFORM I3x3n_R_tildeUo(I3x3n, tildeUo);

      _tildeUo_R_I3nx3[xy] = tildeUo_R_I3nx3;
      _tildeUo_R_uiBarO[xy] = tildeUo_R_uiBarO;
      _tildeUo_R_Uij[xy] = tildeUo_R_Uij;
      _tildeUo_R_ukBarO[xy] = tildeUo_R_ukBarO;
      _tildeUo_R_Ukj[xy] = tildeUo_R_Ukj;

      _I3x3n_R_tildeUo[xy] = I3x3n_R_tildeUo;
    
      _UijTildeU[xy] = tildeU.transpose() * Uij;
      _tildeUoTUij[xy] = tildeUo.transpose() * Uij;
      _tildeUTuiBarO[xy] = tildeU.transpose().modeOneProduct(uiBarO);
      _tildeUTUij[xy] = tildeU.transpose() * Uij;

      // deformation terms
      MATRIX& Ui = Uij;
      MATRIX& Uk = Ukj;
      SANDWICH_TRANSFORM Ui_R_I3nx3(Ui, I3nx3);
      SANDWICH_TRANSFORM Ui_R_uiBarO(Ui, uiBarO);
      SANDWICH_TRANSFORM Ui_R_Ui(Ui, Ui);
      SANDWICH_TRANSFORM Ui_R_ukBarO(Ui, ukBarO);
      SANDWICH_TRANSFORM Ui_R_Uk(Ui, Uk);

      _Ui_R_I3nx3[xy] = Ui_R_I3nx3;
      _Ui_R_uiBarO[xy]= Ui_R_uiBarO;
      _Ui_R_Ui[xy] = Ui_R_Ui;
      _Ui_R_ukBarO[xy] = Ui_R_ukBarO;
      _Ui_R_Uk[xy] = Ui_R_Uk;
      _Ui_uiBarO[xy] = Ui ^ uiBarO;

      // force terms
      VECTOR tildeUo_uiBarO = tildeUo ^ uiBarO;
      _tildeUo_uiBarO[xy] = tildeUo_uiBarO;
      _tildeU_Ui[xy] = tildeU.transpose() * Uij;

      if (_partitionedMesh->constrained(y))
      {
        _Ui_Uk[xy] = Ui ^ Uk;
        _Ui_ukBarO[xy] = Ui ^ ukBarO;
      }
    }
  }

  // do full sandwiches for external computation
  for (int x = 0; x < _partitions; x++)
  {
    cout << " Computing full sandwich " << x << " of " << _partitions << endl;
    MATRIX& Ui = _partitionedMesh->U(x);

    MATRIX tildeUo = _partitionedMesh->restTildeU(x);
    TENSOR3 tildeUi = _partitionedMesh->tildeU(x);

    SANDWICH_TRANSFORM fullTildeUo_Ui(tildeUo, Ui);

    _fullTildeUo_Ui[x] = fullTildeUo_Ui;
  }

  cout << " done." << endl;

  cout << " Writing sandwich cache to disk ... "; flush(cout);
  writeSandwichCache();
  cout << " done. " << endl;
}

//////////////////////////////////////////////////////////////////////
// write out the precomputed sandwiches to disk
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::writeSandwichCache()
{
  string sandwichFile = _partitionedMesh->partitionPath() + string(".sandwiches");
  cout << " Writing out sandwich cache " << sandwichFile.c_str() << " ... ";
  FILE* file = fopen(sandwichFile.c_str(), "wb");
  if (file == NULL)
  {
    cout << " Could not write out sandwich cache! " << endl;
    return;
  }

  // compute sandwiches for diagonal blocks of unconstrained partitions
  for (int x = 0; x < _partitions; x++)
  {
    // compute the common unconstrained and constrained stuff
    for (int y = 0; y < _partitions; y++)
    {
      vector<pair<int,int> > clones = _partitionedMesh->clonedVertices(x,y);
      if (clones.size() == 0) continue;

      pair<int,int> xy(x,y);

      _Ui_I3nx3[xy].write(file);

      _I3x3n_R_uiBarO[xy].write(file);
      _I3x3n_R_ukBarO[xy].write(file);
      _I3x3n_R_Ukj[xy].write(file);
      _I3x3n_R_Uij[xy].write(file);
      if (_partitionedMesh->constrained(x))
      {
        _Ui_uiBarO[xy].write(file);

        if (_partitionedMesh->constrained(y))
        {
          _Ui_Uk[xy].write(file);
          _Ui_ukBarO[xy].write(file);
        }
        else
        {
          _Ui_R_ukBarO[xy].write(file);
          _Ui_R_Uk[xy].write(file);
        }
      }
    }
    _Ui_Ui[x].write(file);

    // do unconstrained only stuff
    if (_partitionedMesh->constrained(x)) continue;

    for (int y = 0; y < _partitions; y++)
    {
      vector<pair<int,int> > clones = _partitionedMesh->clonedVertices(x,y);
      if (clones.size() == 0) continue;

      pair<int,int> xy(x,y);
      _I3x3n_T_restU[xy].write(file);

      _tildeUo_R_I3nx3[xy].write(file);
      _tildeUo_R_uiBarO[xy].write(file);
      _tildeUo_R_Uij[xy].write(file);
      _tildeUo_R_ukBarO[xy].write(file);
      _tildeUo_R_Ukj[xy].write(file);

      _I3x3n_R_tildeUo[xy].write(file);
    
      _UijTildeU[xy].write(file);
      _tildeUoTUij[xy].write(file);
      _tildeUTuiBarO[xy].write(file);
      _tildeUTUij[xy].write(file);

      _Ui_R_I3nx3[xy].write(file);
      _Ui_R_uiBarO[xy].write(file);
      _Ui_R_Ui[xy].write(file);
      _Ui_R_ukBarO[xy].write(file);
      _Ui_R_Uk[xy].write(file);
      _Ui_uiBarO[xy].write(file);

      _tildeUo_uiBarO[xy].write(file);
      _tildeU_Ui[xy].write(file);

      if (_partitionedMesh->constrained(y))
      {
        _Ui_Uk[xy].write(file);
        _Ui_ukBarO[xy].write(file);
      }
    }
  }

  // do full sandwiches for external for computation
  for (int x = 0; x < _partitions; x++)
  {
    _fullTildeUo_Ui[x].write(file);
  }

  fclose(file);
  cout << " done." << endl;
}

//////////////////////////////////////////////////////////////////////
// read in the precomputed sandwiches from disk
//////////////////////////////////////////////////////////////////////
bool SUBSPACE_MULTIBODY_INTEGRATOR::readSandwichCache()
{
  string sandwichFile = _partitionedMesh->partitionPath() + string(".sandwiches");
  cout << " Reading from sandwich cache " << sandwichFile.c_str() << " ... "; flush(cout);
  FILE* file = fopen(sandwichFile.c_str(), "rb");
  if (file == NULL)
  {
    cout << " did not find a cache. " << endl;
    return false;
  }

  // compute sandwiches for diagonal blocks of unconstrained partitions
  for (int x = 0; x < _partitions; x++)
  {
    // compute the common unconstrained and constrained stuff
    for (int y = 0; y < _partitions; y++)
    {
      vector<pair<int,int> > clones = _partitionedMesh->clonedVertices(x,y);
      if (clones.size() == 0) continue;

      pair<int,int> xy(x,y);

      _Ui_I3nx3[xy].read(file);
      _I3x3n_R_uiBarO[xy].read(file);
      _I3x3n_R_ukBarO[xy].read(file);
      _I3x3n_R_Ukj[xy].read(file);
      _I3x3n_R_Uij[xy].read(file);
      if (_partitionedMesh->constrained(x))
      {
        _Ui_uiBarO[xy].read(file);

        if (_partitionedMesh->constrained(y))
        {
          _Ui_Uk[xy].read(file);
          _Ui_ukBarO[xy].read(file);
        }
        else
        {
          _Ui_R_ukBarO[xy].read(file);
          _Ui_R_Uk[xy].read(file);
        }
      }
    }
    _Ui_Ui[x].read(file);

    // do unconstrained only stuff
    if (_partitionedMesh->constrained(x)) continue;

    for (int y = 0; y < _partitions; y++)
    {
      vector<pair<int,int> > clones = _partitionedMesh->clonedVertices(x,y);
      if (clones.size() == 0) continue;

      pair<int,int> xy(x,y);
      _I3x3n_T_restU[xy].read(file);

      _tildeUo_R_I3nx3[xy].read(file);
      _tildeUo_R_uiBarO[xy].read(file);
      _tildeUo_R_Uij[xy].read(file);
      _tildeUo_R_ukBarO[xy].read(file);
      _tildeUo_R_Ukj[xy].read(file);

      _I3x3n_R_tildeUo[xy].read(file);
    
      _UijTildeU[xy].read(file);
      _tildeUoTUij[xy].read(file);
      _tildeUTuiBarO[xy].read(file);
      _tildeUTUij[xy].read(file);

      _Ui_R_I3nx3[xy].read(file);
      _Ui_R_uiBarO[xy].read(file);
      _Ui_R_Ui[xy].read(file);
      _Ui_R_ukBarO[xy].read(file);
      _Ui_R_Uk[xy].read(file);
      _Ui_uiBarO[xy].read(file);

      _tildeUo_uiBarO[xy].read(file);
      _tildeU_Ui[xy].read(file);

      if (_partitionedMesh->constrained(y))
      {
        _Ui_Uk[xy].read(file);
        _Ui_ukBarO[xy].read(file);
      }
    }
  }

  // do full sandwiches for external for computation
  for (int x = 0; x < _partitions; x++)
  {
    _fullTildeUo_Ui[x].read(file);
  }

  fclose(file);
  cout << " done." << endl;
  return true;
}

//////////////////////////////////////////////////////////////////////
// Presum any possible sandwiches
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::presumSandwiches()
{
  for (int x = 0; x < _partitions; x++)
  {
    if (_partitionedMesh->unconstrained(x))
    {
      bool first = true;
      SUBSPACE_TET_MESH* leftMesh = (SUBSPACE_TET_MESH*)(_partitionedMesh->mesh(x)); 
      
      for (int y = 0; y < _partitions; y++)
      {
        vector<pair<int,int> > clones = _partitionedMesh->clonedVertices(x,y);
        if (clones.size() == 0) continue;
        pair<int,int> xy(x,y);

        if (first)
        {
          _summed_I3x3n_R_tildeUo[x] = _I3x3n_R_tildeUo[xy];
          _summed_tildeUo_R_I3nx3[x] = _tildeUo_R_I3nx3[xy];
          _summed_tildeUo_R_uiBarO[x] = _tildeUo_R_uiBarO[xy];
          _summed_tildeUo_R_Uij[x] = _tildeUo_R_Uij[xy];
          _summed_tildeUoTUij[x] = _tildeUoTUij[xy];
          _summed_UijTildeU[x] = _UijTildeU[xy];
          _summed_tildeUTuiBarO[x] = _tildeUTuiBarO[xy];
          _summed_tildeUTUij[x] = _tildeUTUij[xy];
          _summed_I3x3n_T_restU[x] = _I3x3n_T_restU[xy];
          _summed_I3x3n_R_Uij[x] = _I3x3n_R_Uij[xy];
          _summed_Ui_R_I3nx3[x] = _Ui_R_I3nx3[xy];
          _summed_Ui_R_Ui[x] = _Ui_R_Ui[xy];
          _summed_Ui_R_uiBarO[x] = _Ui_R_uiBarO[xy];
          _summed_I3x3n_R_uiBarO[x] = _I3x3n_R_uiBarO[xy];
          _summed_Ui_uiBarO[x] = _Ui_uiBarO[xy];
          _summed_tildeUo_uiBarO[x] = _tildeUo_uiBarO[xy];
          _summed_tildeU_Ui[x] = _tildeU_Ui[xy];

          first = false;
        }
        else
        {
          _summed_I3x3n_R_tildeUo[x] += _I3x3n_R_tildeUo[xy];
          _summed_tildeUo_R_I3nx3[x] += _tildeUo_R_I3nx3[xy];
          _summed_tildeUo_R_uiBarO[x] += _tildeUo_R_uiBarO[xy];
          _summed_tildeUo_R_Uij[x] += _tildeUo_R_Uij[xy];
          _summed_tildeUoTUij[x] += _tildeUoTUij[xy];
          _summed_UijTildeU[x] += _UijTildeU[xy];
          _summed_tildeUTuiBarO[x] += _tildeUTuiBarO[xy];
          _summed_tildeUTUij[x] += _tildeUTUij[xy];
          _summed_I3x3n_T_restU[x] += _I3x3n_T_restU[xy];
          _summed_I3x3n_R_Uij[x] += _I3x3n_R_Uij[xy];
          _summed_Ui_R_I3nx3[x] += _Ui_R_I3nx3[xy];
          _summed_Ui_R_Ui[x] += _Ui_R_Ui[xy];
          _summed_Ui_R_uiBarO[x] += _Ui_R_uiBarO[xy];
          _summed_I3x3n_R_uiBarO[x] += _I3x3n_R_uiBarO[xy];
          _summed_Ui_uiBarO[x] += _Ui_uiBarO[xy];
          _summed_tildeUo_uiBarO[x] += _tildeUo_uiBarO[xy];
          _summed_tildeU_Ui[x] += _tildeU_Ui[xy];
        }

        // allocate workspace as well
        SUBSPACE_TET_MESH* rightMesh = (SUBSPACE_TET_MESH*)(_partitionedMesh->mesh(y)); 
        _workspace3xRxL[xy] = TENSOR3(3, rightMesh->rank(), leftMesh->rank());
        _workspace3xRx3[xy] = TENSOR3(3, rightMesh->rank(), 3);
        _workspaceLxRx3[xy] = TENSOR3(leftMesh->rank(), rightMesh->rank(), 3);
        _workspaceRxLx3[xy] = TENSOR3(rightMesh->rank(), leftMesh->rank(), 3);
        _workspace3xR[xy] = MATRIX(3, rightMesh->rank());
        _workspaceRx3[xy] = MATRIX(rightMesh->rank(), 3);
        _workspaceLxR[xy] = MATRIX(leftMesh->rank(), rightMesh->rank());
        _workspaceRx1x3[xy] = TENSOR3(rightMesh->rank(),1,3);
      }
      int rank = leftMesh->rank();
      _workspace3xLx3[x] = TENSOR3(3,rank,3);
      _workspace3x3xL[x] = TENSOR3(3,3,rank);
      _workspace3x1xL[x] = TENSOR3(3,1,rank);
      _workspace3x1x3[x] = TENSOR3(3,1,3);
      _workspaceLx3x3[x] = TENSOR3(rank,3,3);
      _workspaceLx1x3[x] = TENSOR3(rank,1,3);
      _workspaceLxLx3[x] = TENSOR3(rank,rank,3);
      _workspace3xLxL[x] = TENSOR3(3,rank,rank);
      _workspace3x3[x] = MATRIX(3,3);
      _workspace3xL[x] = MATRIX(3,rank);
      _workspaceLx3[x] = MATRIX(rank,3);
    }
  }
}

//////////////////////////////////////////////////////////////////////
// allocate the block system matrices
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::resizeBlockSystem()
{
  // allocate the block matrix and vector
  _blockA.resizeAndWipe(_totalBlocks, _totalBlocks);
  _blockB.resizeAndWipe(_totalBlocks);

  // create thread-specific block As
  for (int x = 0; x < _partitions; x++)
    _blockAs.push_back(_blockA);

  // set the diagonal size first
  for (int x = 0; x < _partitions; x++)
  {
    int currentRow = _startingBlock[x];
    int leftRank = _partitionedMesh->rank(x);
    if (_partitionedMesh->unconstrained(x))
    {
      _blockA.resizeAndWipeBlock(currentRow, currentRow, 3,3);
      _blockA.resizeAndWipeBlock(currentRow + 1, currentRow + 1, 3,3);
      _blockA.resizeAndWipeBlock(currentRow + 2, currentRow + 2, leftRank, leftRank);
      
      _blockAs[x].resizeAndWipeBlock(currentRow, currentRow, 3,3);
      _blockAs[x].resizeAndWipeBlock(currentRow + 1, currentRow + 1, 3,3);
      _blockAs[x].resizeAndWipeBlock(currentRow + 2, currentRow + 2, leftRank, leftRank);

      _blockB.resizeAndWipeBlock(currentRow, 3);
      _blockB.resizeAndWipeBlock(currentRow + 1, 3);
      _blockB.resizeAndWipeBlock(currentRow + 2, leftRank);
    }
    else
    {
      _blockA.resizeAndWipeBlock(currentRow, currentRow, leftRank, leftRank);
      _blockAs[x].resizeAndWipeBlock(currentRow, currentRow, leftRank, leftRank);
      _blockB.resizeAndWipeBlock(currentRow, leftRank);
    }
  }

  // set the off-diagonal sizes
  for (int x = 0; x < _partitions; x++)
  {
    // cache the index of the block we're currently on
    int currentBlock = _startingBlock[x];
    int leftRank = _partitionedMesh->rank(x);

    // constrained case
    if (_partitionedMesh->constrained(x))
    {
      for (int y = 0; y < _partitions; y++)
      {
        if (!_partitionedMesh->neighbors(x,y)) continue;

        int yBlock = _startingBlock[y];
        int rightRank = _partitionedMesh->rank(y);
        if (_partitionedMesh->constrained(y))
        {
          _blockA.resizeAndWipeBlock(currentBlock, yBlock, leftRank, rightRank);
          _blockAs[x].resizeAndWipeBlock(currentBlock, yBlock, leftRank, rightRank);
        }
        else
        {
          // populate both the current row and the unconstrained row
          _blockA.resizeAndWipeBlock(yBlock, currentBlock, 3, leftRank);
          _blockA.resizeAndWipeBlock(yBlock + 1, currentBlock, 3, leftRank);
          _blockA.resizeAndWipeBlock(yBlock + 2, currentBlock, rightRank, leftRank);
         
          _blockA.resizeAndWipeBlock(currentBlock, yBlock, leftRank, 3);
          _blockA.resizeAndWipeBlock(currentBlock, yBlock + 1, leftRank, 3);
          _blockA.resizeAndWipeBlock(currentBlock, yBlock + 2, leftRank, rightRank);

          _blockAs[x].resizeAndWipeBlock(yBlock, currentBlock, 3, leftRank);
          _blockAs[x].resizeAndWipeBlock(yBlock + 1, currentBlock, 3, leftRank);
          _blockAs[x].resizeAndWipeBlock(yBlock + 2, currentBlock, rightRank, leftRank);
         
          _blockAs[x].resizeAndWipeBlock(currentBlock, yBlock, leftRank, 3);
          _blockAs[x].resizeAndWipeBlock(currentBlock, yBlock + 1, leftRank, 3);
          _blockAs[x].resizeAndWipeBlock(currentBlock, yBlock + 2, leftRank, rightRank);
        }
      }
    }
    else
    {
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
          int rowSize = (i < 2) ? 3 : leftRank;
          int colSize = (j < 2) ? 3 : leftRank;
          _blockA.resizeAndWipeBlock(i + currentBlock, j + currentBlock, rowSize, colSize);
          _blockAs[x].resizeAndWipeBlock(i + currentBlock, j + currentBlock, rowSize, colSize);
        }

      // resize couplings to other unconstrained partitions
      for (int y = 0; y < _partitions; y++)
      {
        if (!_partitionedMesh->neighbors(x,y)) continue;
        if (_partitionedMesh->constrained(y)) continue;
        int currentRow = _startingBlock[x];
        int currentCol = _startingBlock[y];
        int rightRank = _partitionedMesh->rank(y);

        _blockA.resizeAndWipeBlock(currentRow, currentCol, 3, 3);
        _blockA.resizeAndWipeBlock(currentRow + 1, currentCol, 3, 3);
        _blockA.resizeAndWipeBlock(currentRow + 2, currentCol, leftRank, 3);

        _blockA.resizeAndWipeBlock(currentRow, currentCol + 1, 3, 3);
        _blockA.resizeAndWipeBlock(currentRow + 1, currentCol + 1, 3, 3);
        _blockA.resizeAndWipeBlock(currentRow + 2, currentCol + 1, leftRank, 3);

        _blockA.resizeAndWipeBlock(currentRow, currentCol + 2, 3, rightRank);
        _blockA.resizeAndWipeBlock(currentRow + 1, currentCol + 2, 3, rightRank);
        _blockA.resizeAndWipeBlock(currentRow + 2, currentCol + 2, leftRank, rightRank);

        _blockAs[x].resizeAndWipeBlock(currentRow, currentCol, 3, 3);
        _blockAs[x].resizeAndWipeBlock(currentRow + 1, currentCol, 3, 3);
        _blockAs[x].resizeAndWipeBlock(currentRow + 2, currentCol, leftRank, 3);

        _blockAs[x].resizeAndWipeBlock(currentRow, currentCol + 1, 3, 3);
        _blockAs[x].resizeAndWipeBlock(currentRow + 1, currentCol + 1, 3, 3);
        _blockAs[x].resizeAndWipeBlock(currentRow + 2, currentCol + 1, leftRank, 3);

        _blockAs[x].resizeAndWipeBlock(currentRow, currentCol + 2, 3, rightRank);
        _blockAs[x].resizeAndWipeBlock(currentRow + 1, currentCol + 2, 3, rightRank);
        _blockAs[x].resizeAndWipeBlock(currentRow + 2, currentCol + 2, leftRank, rightRank);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
// compute per-interface spring consts
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::computeSpringConsts()
{
  Real meanSpringConst = 0;
  Real maxSpringConst = 0;
  Real minSpringConst = 0;
  bool first = true;
  int totalInterfaces = 0;

  for (int x = 0; x < _partitions; x++)
    for (int y = 0; y < _partitions; y++)
    {
      pair<int, int> xy(x,y);
      int totalClones = _partitionedMesh->totalClones(x,y);

      if (totalClones == 0)
      {
        _interfaceSpringConsts[xy] = 0;
        continue;
      }

      // detect if Pinocchio weightings need to be used
      if (!_partitionedMesh->hasSkeleton())
      {
        _interfaceSpringConsts[xy] = _springConst * _partitionedMesh->interfaceArea(x,y);
        _interfaceSpringConsts[xy] *= 1.0 / totalClones;
        if (x == 0)
        {
          cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
          cout << " interface area " << x << ", " << y << ": " << _partitionedMesh->interfaceArea(x,y) << endl;
          cout << " total clones: " << totalClones << endl;
        }
      }
      else
        _interfaceSpringConsts[xy] = _springConst;

      if (first)
      {
        first = false;
        maxSpringConst = _interfaceSpringConsts[xy];
        minSpringConst = _interfaceSpringConsts[xy];
      }
      else
      {
        if (_interfaceSpringConsts[xy] < minSpringConst)
          minSpringConst = _interfaceSpringConsts[xy];
        if (_interfaceSpringConsts[xy] > maxSpringConst)
          maxSpringConst = _interfaceSpringConsts[xy];
      }
      meanSpringConst += _interfaceSpringConsts[xy];
      totalInterfaces++;
    }
  meanSpringConst *= 1.0 / totalInterfaces;
  cout << " Min spring const: " << minSpringConst << endl;
  cout << " Max spring const: " << maxSpringConst << endl;
  cout << " Mean spring const: " << meanSpringConst << endl;
}

//////////////////////////////////////////////////////////////////////
// recompute spring-const dependent matrices
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::computeSpringMatrices()
{
  cout << " Computing spring matrices... "; flush(cout);

  // Do the Ui-Ui matrix
  for (int x = 0; x < _partitions; x++)
  {
    // compute the common unconstrained and constrained stuff
    int rank = _partitionedMesh->rank(x);
    MATRIX UiUi(rank, rank);
    for (int y = 0; y < _partitions; y++)
    {
      vector<pair<int,int> > clones = _partitionedMesh->clonedVertices(x,y);
      if (clones.size() == 0) continue;
      pair<int,int> xy(x,y);

      MATRIX& Ui = _partitionedMesh->interfaceU(x,y);
      UiUi += _interfaceSpringConsts[xy] * (Ui ^ Ui);
    }
    _spring_Ui_Ui[x] = UiUi;
  }
  cout << " done." << endl;

  // do the summed sandwiches
  for (int x = 0; x < _partitions; x++)
  {
    if (_partitionedMesh->unconstrained(x))
    {
      bool first = true;
      
      for (int y = 0; y < _partitions; y++)
      {
        if (!_partitionedMesh->neighbors(x, y)) continue;
        pair<int,int> xy(x,y);

        Real springConst = _interfaceSpringConsts[xy];

        if (first)
        {
          _spring_summed_I3x3n_R_tildeUo[x] = springConst * _I3x3n_R_tildeUo[xy];
          _spring_summed_tildeUo_R_I3nx3[x] = springConst * _tildeUo_R_I3nx3[xy];
          _spring_summed_tildeUo_R_uiBarO[x] = springConst * _tildeUo_R_uiBarO[xy];
          _spring_summed_tildeUo_R_Uij[x] = springConst * _tildeUo_R_Uij[xy];
          _spring_summed_tildeUoTUij[x] = springConst * _tildeUoTUij[xy];
          _spring_summed_UijTildeU[x] = springConst * _UijTildeU[xy];
          _spring_summed_tildeUTuiBarO[x] = springConst * _tildeUTuiBarO[xy];
          _spring_summed_tildeUTUij[x] = springConst * _tildeUTUij[xy];
          _spring_summed_I3x3n_T_restU[x] = springConst * _I3x3n_T_restU[xy];
          _spring_summed_I3x3n_R_Uij[x] = springConst * _I3x3n_R_Uij[xy];
          _spring_summed_Ui_R_I3nx3[x] = springConst * _Ui_R_I3nx3[xy];
          _spring_summed_Ui_R_Ui[x] = springConst * _Ui_R_Ui[xy];
          _spring_summed_Ui_R_uiBarO[x] = springConst * _Ui_R_uiBarO[xy];
          _spring_summed_I3x3n_R_uiBarO[x] = springConst * _I3x3n_R_uiBarO[xy];
          _spring_summed_Ui_uiBarO[x] = springConst * _Ui_uiBarO[xy];
          _spring_summed_tildeUo_uiBarO[x] = springConst * _tildeUo_uiBarO[xy];
          _spring_summed_tildeU_Ui[x] = springConst * _tildeU_Ui[xy];

          first = false;
        }
        else
        {
          _spring_summed_I3x3n_R_tildeUo[x] += springConst * _I3x3n_R_tildeUo[xy];
          _spring_summed_tildeUo_R_I3nx3[x] += springConst * _tildeUo_R_I3nx3[xy];
          _spring_summed_tildeUo_R_uiBarO[x] += springConst * _tildeUo_R_uiBarO[xy];
          _spring_summed_tildeUo_R_Uij[x] += springConst * _tildeUo_R_Uij[xy];
          _spring_summed_tildeUoTUij[x] += springConst * _tildeUoTUij[xy];
          _spring_summed_UijTildeU[x] += springConst * _UijTildeU[xy];
          _spring_summed_tildeUTuiBarO[x] += springConst * _tildeUTuiBarO[xy];
          _spring_summed_tildeUTUij[x] += springConst * _tildeUTUij[xy];
          _spring_summed_I3x3n_T_restU[x] += springConst * _I3x3n_T_restU[xy];
          _spring_summed_I3x3n_R_Uij[x] += springConst * _I3x3n_R_Uij[xy];
          _spring_summed_Ui_R_I3nx3[x] += springConst * _Ui_R_I3nx3[xy];
          _spring_summed_Ui_R_Ui[x] += springConst * _Ui_R_Ui[xy];
          _spring_summed_Ui_R_uiBarO[x] += springConst * _Ui_R_uiBarO[xy];
          _spring_summed_I3x3n_R_uiBarO[x] += springConst * _I3x3n_R_uiBarO[xy];
          _spring_summed_Ui_uiBarO[x] += springConst * _Ui_uiBarO[xy];
          _spring_summed_tildeUo_uiBarO[x] += springConst * _tildeUo_uiBarO[xy];
          _spring_summed_tildeU_Ui[x] += springConst * _tildeU_Ui[xy];
        }
      }
    }
  }

  // do the sum of the spring consts over all interfaces
  _clonesTimesSpringConsts.resize(_partitions);
  for (int x = 0; x < _partitions; x++)
    for (int y = 0; y < _partitions; y++)
    {
      vector<pair<int,int> > clones = _partitionedMesh->clonedVertices(x,y);
      if (clones.size() == 0) continue;
      pair<int,int> xy(x,y);
      Real springConst = _interfaceSpringConsts[xy];

      _clonesTimesSpringConsts[x] += springConst * clones.size();
    }
}

//////////////////////////////////////////////////////////////////////
// Do a single Newton solve where spring forces are computed in
// full coordinates, and the mesh is assumed to be skinned
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::stepReducedSkinnedDynamic(bool verbose)
{
  TIMER totalTimer;
  if (verbose)
  {
    cout << "==================================================" << endl;
    cout << " Beginning full skinned dynamic spring solve " << _totalSteps << endl;
    cout << "==================================================" << endl;
  }

  // if there are any new bone transforms, they are in the tet meshes, so
  // copy their values to the integrators
  for (int x = 0; x < _partitions; x++)
  {
    UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_partitionedMesh->mesh(x);
    UNCONSTRAINED_SUBSPACE_INTEGRATOR* integrator = (UNCONSTRAINED_SUBSPACE_INTEGRATOR*)_integrators[x];
    integrator->translation() = mesh->rigidTranslation();
    integrator->rotation() = mesh->rotationQuaternion();
  }

  // initialize all integrator (copy current state to olds)
  TIMER initImplicitTimer;
  for (int x = 0; x < _partitions; x++)
    _integrators[x]->initializeImplicitStep();
  _timingBreakdown["Init integrators"] += initImplicitTimer.timing();

  // if we're using explicit quadratics, compute them here before anything else
  TIMER quadraticsTimer;
  for (int x = 0; x < _partitions; x++)
    if (_partitionedMesh->unconstrained(x))
      ((UNCONSTRAINED_SUBSPACE_INTEGRATOR*)_integrators[x])->computeSkinnedQuadraticForcesReduced();
  _timingBreakdown["Compute quadratics"] += quadraticsTimer.timing();

  /////////////////////////////////////////////////////////////////////////
  // Newton-Raphson loop
  /////////////////////////////////////////////////////////////////////////
  // get rotation forces
  Real newtonEps = _solverEps;
  int newtonIterations = 0;
  //int maxNewtonSteps = 10;
  int maxNewtonSteps = _maxNewtonSteps;
  while (newtonIterations < maxNewtonSteps)
  {
    TIMER updateTimer;
    vector<MATRIX> masses;
    vector<VECTOR> residuals;

    for (int x = 0; x < _partitions; x++)
      _integrators[x]->updateState();

    // generate all system matrices first, as this also computes
    // the new q
    for (int x = 0; x < _partitions; x++)
      _integrators[x]->SUBSPACE_INTEGRATOR::generateImplicitMatrices(_timingBreakdown);

    // update all unconstrained meshes
    for (int x = 0; x < _partitions; x++)
    {
      TIMER unconstrainedDiagonalTimer;
      // get the diagonal entry and the residual
      MATRIX mass;
      VECTOR residual;
      computeSkinnedDiagonalReduced(x, mass, residual, true);

      // store the result
      masses.push_back(mass);
      residuals.push_back(residual);
      _timingBreakdown["Reduced unconstrained Diagonals"] += unconstrainedDiagonalTimer.timing();
    }

    /////////////////////////////////////////////////////////////////////////
    // Peek at total residual
    /////////////////////////////////////////////////////////////////////////
    TIMER residualTimer;
    Real residualNorm = 0;
    for (int x = 0; x < _partitions; x++)
      residualNorm += residuals[x].sum2();
    residualNorm = sqrt(residualNorm);
    if (verbose)
      cout << " residual norm: " << residualNorm << endl;
    _timingBreakdown["Compute residual"] += residualTimer.timing();
    if (residualNorm < newtonEps) break;

    /////////////////////////////////////////////////////////////////////////
    // Assemble and solve block version
    /////////////////////////////////////////////////////////////////////////

    // final system matrix
    BLOCK_MATRIX blockA(_partitions, _partitions);

    for (int x = 0; x < _partitions; x++)
    { 
      // add the spring-related diagonal and off-diagonal terms
      TIMER unconstrainedJacobianTimer;
      MATRIX springs;
      computeSkinnedSpringJacobiansReduced(x, springs, true);
      _timingBreakdown["Diagonal unconstrained Jacobian"] += unconstrainedJacobianTimer.timing();

      MATRIX mass = springs;
     
      // add deformation system matrix in
      MATRIX deformationA = *_integrators[x]->A();
      mass += deformationA;

      // build the block A
      blockA.add(mass, x,x);

      // compute coupling Jacobian to other unconstrained partitions
      for (int y = 0; y < _partitions; y++)
      {
        if (!_partitionedMesh->neighbors(x,y)) continue;

        // constrained blocks are populated by constrained partitions
        if (_partitionedMesh->constrained(y)) continue;
        if (x == y) continue;

        TIMER unconstrainedUnconstrainedTimer;
        MATRIX couplingJacobian;

        computeSkinnedCouplingJacobianReduced(x, y, couplingJacobian);
        blockA.subtract(couplingJacobian, x, y);
        _timingBreakdown["Unconstrained-unconstrained Jacobian"] += unconstrainedUnconstrainedTimer.timing();
      }
    }

    // build block b
    TIMER finalSystemTimer;
    BLOCK_VECTOR blockB(_partitions);
    for (int x = 0; x < _partitions; x++)
      blockB.set(residuals[x], x);

    MATRIX fullA = blockA.full();
    VECTOR fullB = blockB.full();
    _timingBreakdown["Build final matrix"] += finalSystemTimer.timing();

    // do the solve
    TIMER solverTimer;
    fullA.factorLU();
    fullA.solveLU(fullB);
    _timingBreakdown["LU solve time"] += solverTimer.timing();

    // TODO: update velocity and accel as well
    // commit the solve results
    TIMER accelUpdateTimer;
    for (int x = 0, i = 0; x < _partitions; x++)
    {
      int rank = _partitionedMesh->rank(x);
      for (int y = 0; y < rank; y++)
        _integrators[x]->position()[y] -= fullB[i + y];
      i += rank;
    }
    _timingBreakdown["Acceleration update"] += accelUpdateTimer.timing();

    newtonIterations++;
  }
  _totalNewtonStepsSeen += newtonIterations;
  if (newtonIterations > _maxNewtonStepsSeen)
    _maxNewtonStepsSeen = newtonIterations;

  TIMER finalizeTimer;
  // finalize everything
  for (int x = 0; x < _partitions; x++)
  {
    *(_partitionedMesh->q(x)) = _integrators[x]->position();
    _integrators[x]->SUBSPACE_INTEGRATOR::finalizeImplicitStep();
  }

  if (verbose)
  {
    cout << " Newton steps: " << newtonIterations << endl;
    if (newtonIterations == maxNewtonSteps)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " NEWTON STEPS MAXXED OUT!!! " << endl;
      cout << " timestep: " << _totalSteps << endl;
      //system("read -p \"Press any key to continue\""); 
    }
  }

  // wipe for next step
  for (int x = 0; x < _partitions; x++)
    _integrators[x]->clearForces();
  
  for (int x = 0; x < _partitions; x++)
    _forceVectors[x].clear();

  _totalSteps++;
  _timingBreakdown["Finalize"] += finalizeTimer.timing();
  _totalTime += totalTimer.timing();
}

//////////////////////////////////////////////////////////////////////
// Do a single Newton solve where spring forces are computed in
// full coordinates, and the mesh is assumed to be skinned
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::stepReducedSkinnedDynamicWithCollisionsOMP(vector<vector<pair<SURFACE*, int> > >& collisions, bool verbose)
{
  TIMER totalTimer;
  if (verbose)
  {
    cout << "==================================================" << endl;
    cout << " Beginning OMP full skinned dynamic spring solve with collisions " << _totalSteps << endl;
    cout << "==================================================" << endl;
  }

  // if there are any new bone transforms, they are in the tet meshes, so
  // copy their values to the integrators
  for (int x = 0; x < _partitions; x++)
  {
    UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_partitionedMesh->mesh(x);
    UNCONSTRAINED_SUBSPACE_INTEGRATOR* integrator = (UNCONSTRAINED_SUBSPACE_INTEGRATOR*)_integrators[x];
    integrator->translation() = mesh->rigidTranslation();
    integrator->rotation() = mesh->rotationQuaternion();
  }

  // NEW: build collision lists for all unconstrained partitions
  TIMER collisionListTimer;
  for (int x = 0; x < _partitions; x++)
    if (_partitionedMesh->unconstrained(x))
      ((UNCONSTRAINED_SUBSPACE_INTEGRATOR*)_integrators[x])->buildCollisionList(collisions[x]);
  _timingBreakdown["Collision Lists"] += collisionListTimer.timing();

  // initialize all integrator (copy current state to olds)
  TIMER initImplicitTimer;
  for (int x = 0; x < _partitions; x++)
    _integrators[x]->initializeImplicitStep();
  _timingBreakdown["Init integrators"] += initImplicitTimer.timing();

  // if we're using explicit quadratics, compute them here before anything else
  TIMER quadraticsTimer;
  for (int x = 0; x < _partitions; x++)
    if (_partitionedMesh->unconstrained(x))
      ((UNCONSTRAINED_SUBSPACE_INTEGRATOR*)_integrators[x])->computeSkinnedQuadraticForcesReduced();
  _timingBreakdown["Compute quadratics"] += quadraticsTimer.timing();

  vector<MATRIX> masses;
  vector<VECTOR> residuals;
  masses.resize(_partitions);
  residuals.resize(_partitions);

  /////////////////////////////////////////////////////////////////////////
  // Newton-Raphson loop
  /////////////////////////////////////////////////////////////////////////
  // get rotation forces
  Real newtonEps = _solverEps;
  int newtonIterations = 0;
  //int maxNewtonSteps = 10;
  int maxNewtonSteps = _maxNewtonSteps;
  Real residualNorm = 0;
  while (newtonIterations < maxNewtonSteps)
  {
    newtonIterations++;
    //for (int x = 0; x < _partitions; x++)
    //  _integrators[x]->updateState();

    TIMER stiffnessTimer;
#if USING_MULTIBODY_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
    for (int x = 0; x < _partitions; x++)
    {
      _integrators[x]->updateState();

      map<string, double> dummyTimers;
      _integrators[x]->SUBSPACE_INTEGRATOR::generateImplicitMatrices(dummyTimers);
    }
    _timingBreakdown["OMP Stiffness construction"] += stiffnessTimer.timing();

    for (int x = 0; x < _partitions; x++)
    {
      TIMER unconstrainedDiagonalTimer;
      computeSkinnedDiagonalReduced(x, masses[x], residuals[x], true);
      _timingBreakdown["Reduced unconstrained Diagonals"] += unconstrainedDiagonalTimer.timing();
      
      // NEW: adding collision residual
      UNCONSTRAINED_SUBSPACE_INTEGRATOR* integrator = (UNCONSTRAINED_SUBSPACE_INTEGRATOR*)_integrators[x];
      BLOCK_VECTOR collisionResidual = integrator->blockImplicitCollisionResiduals();

      //cout << " collision residual " << x << ": " << collisionResidual(2)->norm2() << endl;
      residuals[x] += *(collisionResidual(2));
    }

    /////////////////////////////////////////////////////////////////////////
    // Peek at total residual
    /////////////////////////////////////////////////////////////////////////
    TIMER residualTimer;
    residualNorm = 0;
    for (int x = 0; x < _partitions; x++)
      residualNorm += residuals[x].sum2();
    residualNorm = sqrt(residualNorm);
    if (verbose)
      cout << " residual norm: " << residualNorm << endl;
    _timingBreakdown["Compute residual"] += residualTimer.timing();
    if (residualNorm < newtonEps) break;

    /////////////////////////////////////////////////////////////////////////
    // Assemble and solve block version
    /////////////////////////////////////////////////////////////////////////

    // final system matrix
    BLOCK_MATRIX blockA(_partitions, _partitions);

    for (int x = 0; x < _partitions; x++)
    { 
      // add the spring-related diagonal and off-diagonal terms
      TIMER unconstrainedJacobianTimer;
      MATRIX springs;
      computeSkinnedSpringJacobiansReduced(x, springs, true);
      _timingBreakdown["Diagonal unconstrained Jacobian"] += unconstrainedJacobianTimer.timing();

      MATRIX mass = springs;
      MATRIX deformationA = *_integrators[x]->A();
      mass += deformationA;
      blockA.add(mass, x,x);

      // NEW: adding collision jacobian
      UNCONSTRAINED_SUBSPACE_INTEGRATOR* integrator = (UNCONSTRAINED_SUBSPACE_INTEGRATOR*)_integrators[x];
      BLOCK_MATRIX collisionJacobian;
      map<string, double> dummy;
      integrator->blockImplicitCollisionJacobians(collisionJacobian, dummy);
      blockA.add(*collisionJacobian(2,2), x,x);

      // compute coupling Jacobian to other unconstrained partitions
      for (int y = 0; y < _partitions; y++)
      {
        if (!_partitionedMesh->neighbors(x,y)) continue;

        // constrained blocks are populated by constrained partitions
        if (_partitionedMesh->constrained(y)) continue;
        if (x == y) continue;

        TIMER unconstrainedUnconstrainedTimer;
        MATRIX couplingJacobian;
        computeSkinnedCouplingJacobianReduced(x, y, couplingJacobian);
        blockA.subtract(couplingJacobian, x, y);
        _timingBreakdown["Unconstrained-unconstrained Jacobian"] += unconstrainedUnconstrainedTimer.timing();
      }
    }

    // build block b
    TIMER finalSystemTimer;
    BLOCK_VECTOR blockB(_partitions);
    for (int x = 0; x < _partitions; x++)
      blockB.set(residuals[x], x);

    MATRIX fullA = blockA.full();
    VECTOR fullB = blockB.full();
    _timingBreakdown["Build final matrix"] += finalSystemTimer.timing();

    // do the solve
    TIMER solverTimer;
    fullA.factorLU();
    fullA.solveLU(fullB);
    _timingBreakdown["LU solve time"] += solverTimer.timing();

    // TODO: update velocity and accel as well
    // commit the solve results
    TIMER accelUpdateTimer;
    for (int x = 0, i = 0; x < _partitions; x++)
    {
      int rank = _partitionedMesh->rank(x);
      for (int y = 0; y < rank; y++)
        _integrators[x]->position()[y] -= fullB[i + y];
      i += rank;
    }
    _timingBreakdown["Acceleration update"] += accelUpdateTimer.timing();
  }
  _totalNewtonStepsSeen += newtonIterations;
  if (newtonIterations > _maxNewtonStepsSeen)
    _maxNewtonStepsSeen = newtonIterations;

  TIMER finalizeTimer;
  // finalize everything
  for (int x = 0; x < _partitions; x++)
  {
    *(_partitionedMesh->q(x)) = _integrators[x]->position();
    _integrators[x]->SUBSPACE_INTEGRATOR::finalizeImplicitStep();
  }

  if (verbose)
  {
    cout << " Newton steps: " << newtonIterations << endl;
    if (newtonIterations == maxNewtonSteps)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " NEWTON STEPS MAXXED OUT!!! " << endl;
      cout << " timestep: " << _totalSteps << endl;
      //system("read -p \"Press any key to continue\"");
    }
  }

  // wipe for next step
  for (int x = 0; x < _partitions; x++)
    _integrators[x]->clearForces();
  
  for (int x = 0; x < _partitions; x++)
    _forceVectors[x].clear();

  _totalSteps++;
  _timingBreakdown["Finalize"] += finalizeTimer.timing();
  _totalTime += totalTimer.timing();
}
//////////////////////////////////////////////////////////////////////
// Do a single Newton solve where spring forces are computed in
// full coordinates, and the mesh is assumed to be skinned
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::stepReducedSkinnedDynamicOMP(bool verbose)
{
  TIMER totalTimer;
  if (verbose)
  {
    cout << "==================================================" << endl;
    cout << " Beginning OMP full skinned dynamic spring solve " << _totalSteps << endl;
    cout << "==================================================" << endl;
  }

  // if there are any new bone transforms, they are in the tet meshes, so
  // copy their values to the integrators
  for (int x = 0; x < _partitions; x++)
  {
    UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_partitionedMesh->mesh(x);
    UNCONSTRAINED_SUBSPACE_INTEGRATOR* integrator = (UNCONSTRAINED_SUBSPACE_INTEGRATOR*)_integrators[x];
    integrator->translation() = mesh->rigidTranslation();
    integrator->rotation() = mesh->rotationQuaternion();
  }

  // initialize all integrator (copy current state to olds)
  TIMER initImplicitTimer;
  for (int x = 0; x < _partitions; x++)
    _integrators[x]->initializeImplicitStep();
  _timingBreakdown["Init integrators"] += initImplicitTimer.timing();

  // if we're using explicit quadratics, compute them here before anything else
  TIMER quadraticsTimer;
  for (int x = 0; x < _partitions; x++)
    if (_partitionedMesh->unconstrained(x))
      ((UNCONSTRAINED_SUBSPACE_INTEGRATOR*)_integrators[x])->computeSkinnedQuadraticForcesReduced();
  _timingBreakdown["Compute quadratics"] += quadraticsTimer.timing();

  vector<MATRIX> masses;
  vector<VECTOR> residuals;
  masses.resize(_partitions);
  residuals.resize(_partitions);

  /////////////////////////////////////////////////////////////////////////
  // Newton-Raphson loop
  /////////////////////////////////////////////////////////////////////////
  // get rotation forces
  Real newtonEps = _solverEps;
  int newtonIterations = 0;
  //int maxNewtonSteps = 10;
  int maxNewtonSteps = _maxNewtonSteps;
  Real residualNorm = 0;
  while (newtonIterations < maxNewtonSteps)
  {
    newtonIterations++;
    //for (int x = 0; x < _partitions; x++)
    //  _integrators[x]->updateState();

    TIMER stiffnessTimer;
#if USING_MULTIBODY_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
    for (int x = 0; x < _partitions; x++)
    {
      _integrators[x]->updateState();

      map<string, double> dummyTimers;
      _integrators[x]->SUBSPACE_INTEGRATOR::generateImplicitMatrices(dummyTimers);
    }
    _timingBreakdown["OMP Stiffness construction"] += stiffnessTimer.timing();

    for (int x = 0; x < _partitions; x++)
    {
      TIMER unconstrainedDiagonalTimer;
      computeSkinnedDiagonalReduced(x, masses[x], residuals[x], true);
      //computeSkinnedDiagonal(x, masses[x], residuals[x], true);
      _timingBreakdown["Reduced unconstrained Diagonals"] += unconstrainedDiagonalTimer.timing();
    }

    /////////////////////////////////////////////////////////////////////////
    // Peek at total residual
    /////////////////////////////////////////////////////////////////////////
    TIMER residualTimer;
    residualNorm = 0;
    for (int x = 0; x < _partitions; x++)
      residualNorm += residuals[x].sum2();
    residualNorm = sqrt(residualNorm);
    if (verbose)
      cout << " residual norm: " << residualNorm << endl;
    _timingBreakdown["Compute residual"] += residualTimer.timing();
    if (residualNorm < newtonEps) break;

    /////////////////////////////////////////////////////////////////////////
    // Assemble and solve block version
    /////////////////////////////////////////////////////////////////////////

    // final system matrix
    BLOCK_MATRIX blockA(_partitions, _partitions);

    for (int x = 0; x < _partitions; x++)
    { 
      // add the spring-related diagonal and off-diagonal terms
      TIMER unconstrainedJacobianTimer;
      MATRIX springs;
      computeSkinnedSpringJacobiansReduced(x, springs, true);
      //computeSkinnedSpringJacobians(x, springs, true);
      _timingBreakdown["Diagonal unconstrained Jacobian"] += unconstrainedJacobianTimer.timing();

      MATRIX mass = springs;
      MATRIX deformationA = *_integrators[x]->A();
      mass += deformationA;
      blockA.add(mass, x,x);

      // compute coupling Jacobian to other unconstrained partitions
      for (int y = 0; y < _partitions; y++)
      {
        if (!_partitionedMesh->neighbors(x,y)) continue;

        // constrained blocks are populated by constrained partitions
        if (_partitionedMesh->constrained(y)) continue;
        if (x == y) continue;

        TIMER unconstrainedUnconstrainedTimer;
        MATRIX couplingJacobian;
        computeSkinnedCouplingJacobianReduced(x, y, couplingJacobian);
        //computeSkinnedCouplingJacobian(x, y, couplingJacobian);
        blockA.subtract(couplingJacobian, x, y);
        _timingBreakdown["Unconstrained-unconstrained Jacobian"] += unconstrainedUnconstrainedTimer.timing();
      }
    }

    // build block b
    TIMER finalSystemTimer;
    BLOCK_VECTOR blockB(_partitions);
    for (int x = 0; x < _partitions; x++)
      blockB.set(residuals[x], x);

    MATRIX fullA = blockA.full();
    VECTOR fullB = blockB.full();
    _timingBreakdown["Build final matrix"] += finalSystemTimer.timing();

    // do the solve
    TIMER solverTimer;
    fullA.factorLU();
    fullA.solveLU(fullB);
    _timingBreakdown["LU solve time"] += solverTimer.timing();

    // TODO: update velocity and accel as well
    // commit the solve results
    TIMER accelUpdateTimer;
    for (int x = 0, i = 0; x < _partitions; x++)
    {
      int rank = _partitionedMesh->rank(x);
      for (int y = 0; y < rank; y++)
        _integrators[x]->position()[y] -= fullB[i + y];
      i += rank;
    }
    _timingBreakdown["Acceleration update"] += accelUpdateTimer.timing();
  }
  _totalNewtonStepsSeen += newtonIterations;
  if (newtonIterations > _maxNewtonStepsSeen)
    _maxNewtonStepsSeen = newtonIterations;

  TIMER finalizeTimer;
  // finalize everything
  for (int x = 0; x < _partitions; x++)
  {
    *(_partitionedMesh->q(x)) = _integrators[x]->position();
    _integrators[x]->SUBSPACE_INTEGRATOR::finalizeImplicitStep();
  }

  if (verbose)
  {
    cout << " Newton steps: " << newtonIterations << endl;
    if (newtonIterations == maxNewtonSteps)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " NEWTON STEPS MAXXED OUT!!! " << endl;
      cout << " timestep: " << _totalSteps << endl;
      //system("read -p \"Press any key to continue\"");
    }
  }

  // wipe for next step
  for (int x = 0; x < _partitions; x++)
    _integrators[x]->clearForces();
  
  for (int x = 0; x < _partitions; x++)
    _forceVectors[x].clear();

  _totalSteps++;
  _timingBreakdown["Finalize"] += finalizeTimer.timing();
  _totalTime += totalTimer.timing();
}

//////////////////////////////////////////////////////////////////////
// reset simulation
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::reset()
{
  for (int x = 0; x < _partitions; x++)
    _integrators[x]->reset();
}

//////////////////////////////////////////////////////////////////////
// compute the spring energy of a partition
// add the inter-partition spring forces
//////////////////////////////////////////////////////////////////////
Real SUBSPACE_MULTIBODY_INTEGRATOR::springEnergy(int partition)
{
  UNCONSTRAINED_SUBSPACE_INTEGRATOR* integrator = (UNCONSTRAINED_SUBSPACE_INTEGRATOR*)_integrators[partition];
  UNCONSTRAINED_SUBSPACE_TET_MESH* leftMesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_partitionedMesh->mesh(partition);
  vector<VEC3F>& leftVertices = leftMesh->vertices();
  leftMesh->updateFullMesh();

  // get the updated vars
  VEC3F& translationDotDot = integrator->translationAcceleration();
  VEC3F& angularAcceleration = integrator->angularAcceleration();

  // get the old vars
  VEC3F& translationOld = integrator->translationOld();
  VEC3F& translationDotOld = integrator->translationVelocityOld();
  VEC3F& translationDotDotOld = integrator->translationAccelerationOld();

  QUATERNION& rotationOld = integrator->rotationOld();
  VEC3F& angularVelocityOld = integrator->angularVelocityOld();
  VEC3F& angularAccelerationOld = integrator->angularAccelerationOld();

  VECTOR& qDotOld = integrator->velocityOld();
  VECTOR& qDotDotOld = integrator->accelerationOld();
  VECTOR& qDotDot = integrator->acceleration();

  // compute new angular quantities
  Real* alpha = integrator->accelerationAlpha();
  VEC3F angularDelta = alpha[2] * angularVelocityOld + 
                       alpha[3] * angularAccelerationOld + 
                       alpha[4] * angularAcceleration; 
  VEC3F angularVelocity = angularVelocityOld + 
                          alpha[0] * angularAccelerationOld + 
                          alpha[1] * angularAcceleration;
  QUATERNION update = QUATERNION::fromAxisAngle(angularDelta);
  QUATERNION rotationQuaternion = update * rotationOld;

  // compute new translation quantities
  VEC3F translation = translationOld + alpha[2] * translationDotOld + alpha[3] * translationDotDotOld + alpha[4] * translationDotDot;
  VEC3F translationVelocity = translationDotOld + alpha[0] * translationDotDotOld + alpha[1] * translationDotDot;

  // compute new deformation quantities
  VECTOR qDot = qDotOld + alpha[0] * qDotDotOld + alpha[1] * qDotDot;

  int totalClones = 0;
  VECTOR fullSprings(leftMesh->U().rows());
  VECTOR fullDampers(leftMesh->U().rows());
  VECTOR fullVelocity = leftMesh->U() * qDot;

  Real finalEnergy = 0.0;
  for (int x = 0; x < _partitions; x++)
  {
    vector<pair<int,int> > clones = _partitionedMesh->clonedVertices(partition, x);
    if (clones.size() == 0) continue;

    SUBSPACE_TET_MESH* rightMesh = (SUBSPACE_TET_MESH*)_partitionedMesh->mesh(x);
    vector<VEC3F>& rightVertices = rightMesh->vertices();
    rightMesh->updateFullMesh();

    MATRIX3 rightRotation(MATRIX3::I());
    VEC3F rightTranslation;
    if (_partitionedMesh->unconstrained(x))
    {
      rightTranslation = _partitionedMesh->rigidTranslation(x);
    }
    for (unsigned int y = 0; y < clones.size(); y++)
    {
      int leftIndex = clones[y].first;
      int rightIndex = clones[y].second;

      VEC3F leftVertex = leftVertices[leftIndex];
      VEC3F rightVertex = rightVertices[rightIndex];

      leftVertex = rotationQuaternion.rotates(leftVertex) + translation;
      rightVertex = rightRotation * rightVertex + rightTranslation;

      VEC3F diff = leftVertex - rightVertex;
      pair<int, int> xy(partition, x);
      finalEnergy += 0.5 * _interfaceSpringConsts[xy] * norm2(diff);
    }
  }

  return finalEnergy;
}

//////////////////////////////////////////////////////////////////////
// compute ODE versions of the sandwiches -- the rest poses change
//////////////////////////////////////////////////////////////////////
void SUBSPACE_MULTIBODY_INTEGRATOR::precomputeOdeSandwiches(string dataPath)
{
  // force the rest poses to recompute -- this has to be done
  // irrespective of whether the sandwiches already exist or not
  ((PARTITIONED_SKINNED_SUBSPACE_TET_MESH*)_partitionedMesh)->loadOdeSkeletonFrame(dataPath, 0);

  // check to see if the ODE version has been computed before
  string sandwichFile = _partitionedMesh->partitionPath() + string(".ode.sandwiches");
  cout << " Checking for ODE sandwich cache " << sandwichFile.c_str() << " ... "; flush(cout);
  FILE* file = fopen(sandwichFile.c_str(), "rb");
  if (file != NULL)
  {
    cout << " cache found!" << endl;
    fclose(file);
    return;
  }
  //fclose(file);
  cout << " no cache found, precomputing." << endl;

  // recompute new sandwiches
  precomputeSandwiches();

  // write out the flag file
  sandwichFile = _partitionedMesh->partitionPath() + string(".ode.sandwiches");
  file = fopen(sandwichFile.c_str(), "w");

  fprintf(file, "ODE version computed.\n");

  fclose(file);
}
