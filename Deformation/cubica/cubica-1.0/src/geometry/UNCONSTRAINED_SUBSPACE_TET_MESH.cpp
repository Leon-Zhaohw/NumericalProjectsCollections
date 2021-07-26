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
// UNCONSTRAINED_SUBSPACE_TET_MESH.cpp: implementation of the UNCONSTRAINED_SUBSPACE_TET_MESH class.
//
//////////////////////////////////////////////////////////////////////

#include "UNCONSTRAINED_SUBSPACE_TET_MESH.h"

#if !USING_OSX
#include <GL/glut.h>
#endif

//////////////////////////////////////////////////////////////////////
// Constructor for tetrahedra
//////////////////////////////////////////////////////////////////////
UNCONSTRAINED_SUBSPACE_TET_MESH::UNCONSTRAINED_SUBSPACE_TET_MESH(
                                     const char* filename, 
                                     MATERIAL** materials, 
                                     int totalMaterials,
                                     bool simulateFullspace,
                                     const char* eigenvectors, 
                                     const char* eigenvalues,
                                     bool projectTranslation) :
  SUBSPACE_TET_MESH(filename, materials, totalMaterials, simulateFullspace, eigenvectors, eigenvalues)
{
  // set rigid components to zero
  _translation.clear();
  _translationOld.clear();
  //_rotation = MATRIX3::I();
  //_rotationOld = MATRIX3::I();

  //_rotationDebug = MATRIX3::I();

  _totalRotationsSeen = 0;
  _totalRotationsKept = 0;
  _meanRotationError = 0;

  // Georgii-style rotation vars
  _rotationLambda = 0.0;
  _rotationQuaternion = QUATERNION(MATRIX3::I());
  _rotationQuaternionOld = QUATERNION(MATRIX3::I());

  // print out the center of mass
  computeCenterOfMass();
  computeRestCenterOfMass();

  // set center of mass as a translation
  _translationOld = _centerOfMass;
  _translation = _centerOfMass;

  _originalCenterOfMass = _centerOfMass;

  // make center of mass the origin for all vertices
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _vertices[x] -= _centerOfMass;
  for (unsigned int x = 0; x < _restPose.size(); x++)
    _restPose[x] -= _centerOfMass;

  // if we're not simulating in the subspace, we're training, so don't attach a rigid frame
  // or intiialize any additional subspace quantities
  if (simulateFullspace) return;

  // shouldn't be necessary if MultiZeroEigsLMA was run -- will just introduce
  // some numerical noise
#if 0
  cout << " Projecting translation out of basis ..."; flush(cout);
 
  if (projectTranslation)
  { 
    if (!readTranslationlessBasis())
    {
      cout << " no cache found ... "; flush(cout);

      // subtract translation from the basis entirely
      MATRIX minusTranslations(_UBasis.rows(), _UBasis.cols() + 3);
      
      // copy the old basis into the new basis
      for (int y = 0; y < _UBasis.cols(); y++)
        for (int x = 0; x < _UBasis.rows(); x++)
          minusTranslations(x, y + 3) = _UBasis(x,y);

      // front-load translation into the new basis
      for (int x = 0; x < _UBasis.rows(); x++)
      {
        if (x % 3 == 0) minusTranslations(x, 0) = 1;
        if (x % 3 == 1) minusTranslations(x, 1) = 1;
        if (x % 3 == 2) minusTranslations(x, 2) = 1;
      }
      minusTranslations.orthogonalize();

      // set the basis to the new translationless basis
      MATRIX newBasis(_UBasis.rows(), _UBasis.cols());
      for (int x = 0; x < _UBasis.rows(); x++)
        for (int y = 0; y < _UBasis.cols(); y++)
          newBasis(x,y) = minusTranslations(x, y + 3);
      updateBasis(newBasis);

      // write out the result
      writeTranslationlessBasis();
    }
    else
      cout << " cache found! ... "; flush(cout);

    cout << "done." << endl;
  }
#endif

  /* 
  VECTOR restVector(_restPose.size() * 3);
  for (unsigned int x = 0; x < _restPose.size(); x++)
  {
    restVector[3 * x] = _restPose[x][0];
    restVector[3 * x + 1] = _restPose[x][1];
    restVector[3 * x + 2] = _restPose[x][2];
  }
  */

  // set some defaults to compute rotation stats
  _totalGeorgiiIterations = 0;
  _maxGeorgiiIterations = 0;
  _totalSteps = 0;

  // precompute the center of mass projection
  if (_totalKeyTets != 0)
  {
    cout << " Computing variable mass matrix ..."; flush(cout);
    if (!readVariableMassMatrix())
    {
      cout << " no cache found! ..."; flush(cout);
      cacheSubspaceCenterOfMass();

      // recompute the mass matrix vars since we moved everything
      cacheInertiaVars();
      cacheRotationDefoVars();
      cacheMassMatrixVars();

      writeVariableMassMatrix();
    }
    else
      cout << " cache found! ..."; flush(cout);
    cout << " done." << endl;
  }

  reset();

  updateFullMesh();
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
UNCONSTRAINED_SUBSPACE_TET_MESH::~UNCONSTRAINED_SUBSPACE_TET_MESH()
{
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_TET_MESH::resetRigidTranslation()
{
  /*
  Real xNormalized = _UBasis(0,0);
  Real yNormalized = _UBasis(1,1);
  Real zNormalized = _UBasis(2,2);
 
  // update q
  _q[0] = _centerOfMass[0] / xNormalized;
  _q[1] = _centerOfMass[1] / yNormalized;
  _q[2] = _centerOfMass[2] / zNormalized;
  _qOld[0] = _centerOfMass[0] / xNormalized;
  _qOld[1] = _centerOfMass[1] / yNormalized;
  _qOld[2] = _centerOfMass[2] / zNormalized;

  cout << " translation: " << _q[0] << " " << _q[1] << " " << _q[2] << endl;
  */
}

//////////////////////////////////////////////////////////////////////
// Compute the rotation using shape matching
//////////////////////////////////////////////////////////////////////
bool UNCONSTRAINED_SUBSPACE_TET_MESH::updateRotation(int debug, VECTOR& velocity, VECTOR& acceleration)
{
#if 0
  // backup current rotation
  //_rotationNewtonOld = _rotation;
 
  // backup the current state 
  _qNewtonOld = _q;

  // DEBUG: trying filtering with derivatives as well
  VECTOR velocityNewtonOld = velocity;
  VECTOR accelerationNewtonOld = acceleration;

  /*
  // compute the new rotation
  _rotation = computeShapeMatchingRotation();
  cout << " Shape matching:   " << _rotation << endl;
  cout << "                   " << QUATERNION(_rotation) << endl;

  QUATERNION georgiiRotation = computeGeorgiiRotation();
  cout << " Georgii rotation: " << georgiiRotation.toRotationMatrix() << endl;
  cout << "                   " << georgiiRotation << endl;

  _rotation = georgiiRotation.toRotationMatrix();

  QUATERNION fullGeorgiiRotation = computeFullRankGeorgiiRotation();
  cout << " Full georgii rotation: " << fullGeorgiiRotation.toRotationMatrix() << endl;
  cout << "                        " << fullGeorgiiRotation<< endl;

  static int count = 0;
  if (count == 2) exit(0);
  count++;
  */
  assert(centerOfMassIsZero());

  //QUATERNION georgiiRotation = computeGeorgiiRotation(debug);
  //_rotation = _rotation * georgiiRotation.toRotationMatrix();
  
  //_rotation = computeShapeMatchingRotation();

  // update q
  //MATRIX3 R = _rotation.transpose();
  //R = R * _rotationNewtonOld;

  MATRIX sandwich = _basisSandwich.transform(R);
  VECTOR restSandwich = _restSandwich.vectorTransform(R);
  //MATRIX3 transpose = _rotation.transpose(); 
  MATRIX transSandwich= _translationSandwich.transform(transpose);
  VEC3F diff3 = _translationNewtonOld - _translation;
  VECTOR diff(3); diff[0] = diff3[0]; diff[1] = diff3[1]; diff[2] = diff3[2];

  /*
  // DEBUG
  VEC3F translationBefore;
  translationBefore[0] = _q[0];
  translationBefore[1] = _q[1];
  translationBefore[2] = _q[2];
  VECTOR qBefore = _q;
  if (debug == 4)
  {
    cout << " q before: " << _q << endl;
    cout << " q before norm2: " << _q.norm2() << endl;
    cout << " translation before: " << translationBefore << endl;
  }
  */

  _q =  sandwich * _q;
  _q += restSandwich;
  _q -= _projectedRestPose;

  //cout << " SKIPPING TRANSLATION UPDATE " << endl;
  //_q += transSandwich * diff; 

  /*
  // DEBUG
  VEC3F translationAfter;
  translationAfter[0] = _q[0];
  translationAfter[1] = _q[1];
  translationAfter[2] = _q[2];
  if (debug == 4)
  {
    cout << " q after: " << _q << endl;
    cout << " q after norm2: " << _q.norm2() << endl;
    VECTOR qDiff = qBefore - _q;
    cout << " q diff norm2: " << qDiff.norm2() << endl;
    cout << " translation after: " << translationAfter << endl;
    VEC3F translationDiff = translationBefore - translationAfter;
    cout << " translation diff: " << translationDiff * translationDiff << endl;
  }
  */

  // DEBUG: trying filtering with derivatives as well
  VECTOR velocityNew = sandwich * velocityNewtonOld;
  VECTOR accelerationNew = sandwich * accelerationNewtonOld;

  //Real relativeError              = computeRotationErrorNorm(_q, _rotation, _translation, _qNewtonOld, _rotationNewtonOld, _translationNewtonOld, debug);
  //Real relativeErrorVelocity      = computeRotationErrorNorm(velocityNew, _rotation, _translation, velocityNewtonOld, _rotationNewtonOld, _translationNewtonOld, debug);
  //Real relativeErrorAcceleration  = computeRotationErrorNorm(accelerationNew, _rotation, _translation, accelerationNewtonOld, _rotationNewtonOld, _translationNewtonOld, debug);

  // DEBUG
  if (debug == 4)
  {
    cout << " relative position error: " << relativeError << endl;
    cout << " relative velocity error: " << relativeErrorVelocity << endl;
    cout << " relative acceleration error: " << relativeErrorAcceleration << endl ;

    cout << " Absolute velocity magnitude: " << velocityNew * velocityNew << endl;
    cout << " Absolute acceleration magnitude: " << accelerationNew * accelerationNew << endl;
    cout << endl;
  }

  bool accepted = false;
  //Real relativeThreshold = 1e-3;
  //Real relativeThreshold = 1e-2;
  //Real relativeThreshold = 1e-1;
  Real relativeThreshold = 1;
  //Real relativeThreshold = 0;
  //Real relativeThreshold = 1e-9;

  // if there's too much relative error
  //if (relativeError > 7.5e-3)
  //if (relativeError > 1e-2 || relativeErrorVelocity > 1e-3 || relativeErrorAcceleration > 1e-3)
  if (relativeError > relativeThreshold || 
      relativeErrorVelocity > relativeThreshold || 
      relativeErrorAcceleration > relativeThreshold)
  //if (relativeError > 1e-2)
  //if (relativeError > 1e-4)
  {
    // roll back the rotation -- reject the new one
    _q = _qNewtonOld;
    //_rotation = _rotationNewtonOld;

    //if (debug == 1)
    //  cout << " ROTATION REJECTED " << endl;
  }
  else
  {
    /*
    // the rotation is good --
    // go ahead and commit the qOld now as well
    //
    // SHOULDN'T BE DOING THIS
    _qOld = sandwich * _qOld;
    _qOld += restSandwich;
    _qOld -= _projectedRestPose;
    */

    _totalRotationsKept++;
    accepted = true;
    //if (debug == 1)
    //  cout << " ROTATION ACCEPTED" << endl;
  }

  Real angle;
  VEC3F axis;
  //QUATERNION quaternion(_rotation);
  //quaternion.axisAngle(axis, angle);
  /*
  cout << " translation: " << _translation << endl;
  cout << " center of mass: " << _centerOfMass << endl;
  cout << " rest center of mass: " << _restCenterOfMass << endl;
  cout << " rotation: " << QUATERNION(_rotation) << endl;
  cout << " Axis: " << axis << " Angle: " << angle << endl << endl;
  */

  _totalRotationsSeen++;
  return accepted;

  /*
  if (debug == 1)
  {
    cout << " position error:     " << relativeError << endl;
    cout << " velocity error:     " << relativeErrorVelocity << endl;
    cout << " acceleration error: " << relativeErrorAcceleration << endl;

    //cout << " old rotation: " << _rotationNewtonOld << endl;
    //cout << " new rotation: " << _rotation << endl;
    MATRIX3 diff = _rotationNewtonOld - _rotation;
    //cout << " diff:         " << _rotationNewtonOld - _rotation << endl;
    //cout << " diff norm:  " << diff.absSum() << endl;
    //cout << " error norm: " << relativeError << endl;
  }
  */

  //cout << " Rotations kept: " << ((Real)_totalRotationsKept / _totalRotationsSeen) * 100.0 << "% (" << _totalRotationsKept << " of " << _totalRotationsSeen << ")" << endl;
#endif
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " DEPRECATED. " << endl;
  exit(0);
}

//////////////////////////////////////////////////////////////////////
// compute the 2-norm error of the new rotation
//////////////////////////////////////////////////////////////////////
Real UNCONSTRAINED_SUBSPACE_TET_MESH::computeRotationErrorNorm(VECTOR& q,    MATRIX3& rotation,    VEC3& translation,
                                                               VECTOR& qOld, MATRIX3& rotationOld, VEC3& translationOld, int debug)
{
  /*
  // Brute force computation
  VECTOR xOld = _UBasis * qOld;
  VECTOR xNew = _UBasis * q;
  for (unsigned int x = 0; x < _restPose.size(); x++)
  {
    VEC3F old3;
    VEC3F new3;

    old3[0] = _restPose[x][0] + xOld[3 * x];
    old3[1] = _restPose[x][1] + xOld[3 * x + 1];
    old3[2] = _restPose[x][2] + xOld[3 * x + 2];
    old3 = rotationOld * old3 + translationOld;

    new3[0] = _restPose[x][0] + xNew[3 * x];
    new3[1] = _restPose[x][1] + xNew[3 * x + 1];
    new3[2] = _restPose[x][2] + xNew[3 * x + 2];
    new3 = rotation * new3 + translation;

    xOld[3 * x] = old3[0];
    xOld[3 * x + 1] = old3[1];
    xOld[3 * x + 2] = old3[2];
    
    xNew[3 * x] = new3[0];
    xNew[3 * x + 1] = new3[1];
    xNew[3 * x + 2] = new3[2];
  }

  // term-by-term computation
  VECTOR p0(_restPose.size() * 3);
  VECTOR t0(_restPose.size() * 3);
  VECTOR t1(_restPose.size() * 3);

  VECTOR x0 = _UBasis * qOld;
  VECTOR x1 = _UBasis * q;
  for (unsigned int x = 0; x < _restPose.size(); x++)
  {
    p0[3 * x] = _restPose[x][0];
    p0[3 * x + 1] = _restPose[x][1];
    p0[3 * x + 2] = _restPose[x][2];
    
    t0[3 * x] = translationOld[0];
    t0[3 * x + 1] = translationOld[1];
    t0[3 * x + 2] = translationOld[2];
    
    t1[3 * x] = translation[0];
    t1[3 * x + 1] = translation[1];
    t1[3 * x + 2] = translation[2];
  }
  */

  MATRIX3 R0R1 = rotationOld.transpose() * rotation;
  MATRIX3 R0T  = rotationOld.transpose();
  MATRIX3 R1R0 = rotation.transpose() * rotationOld;
  MATRIX3 R1T  = rotation.transpose();
  MATRIX3& R0  = rotationOld;
  MATRIX3& R1  = rotation;
  Real size = (Real)_restPose.size();

  // vector versions for type reasons;
  VECTOR vTranslation(3);
  VECTOR vTranslationOld(3);
  for (int x = 0; x < 3; x++)
  {
    vTranslation[x] = translation[x];
    vTranslationOld[x] = translationOld[x];
  }

  /*
  // UNREDUCED VERSION
  Real term1 = p0 * p0;
  Real term2 = x0 * x0;
  Real term3 = t0 * t0;
  Real term4 = p0 * p0;
  Real term5 = x1 * x1;
  Real term6 = t1 * t1;
  Real term7_8 = 0;
  for (unsigned int x = 0; x < _restPose.size(); x++)
    term7_8 += -2.0 * _restPose[x] * (R0R1 * _restPose[x]);
  Real term9_10 = 0;
  for (unsigned int x = 0; x < _restPose.size(); x++)
  {
    VEC3F oldX, newX;
    oldX[0] = x0[x * 3]; oldX[1] = x0[x * 3 + 1]; oldX[2] = x0[x * 3 + 2];
    newX[0] = x1[x * 3]; newX[1] = x1[x * 3 + 1]; newX[2] = x1[x * 3 + 2];
    term9_10 += -2.0 * oldX * (R0R1 * newX);
  }
  Real term11_12 = -2.0 * (t0 * t1);
  Real term13_14 = 2.0 * (p0 * x0);
  Real term15_16 = 2.0 * (p0 * t0);
  Real term17_18 = 0.0;
  for (unsigned int x = 0; x < _restPose.size(); x++)
  {
    VEC3F newX;
    newX[0] = x1[x * 3]; newX[1] = x1[x * 3 + 1]; newX[2] = x1[x * 3 + 2];
    term17_18 += -2.0 * _restPose[x] * (R0R1 * newX);
  }
  Real term19_20 = 0.0;
  for (unsigned int x = 0; x < _restPose.size(); x++)
    term19_20 += -2.0 * _restPose[x] * (R0T * translation);
  Real term21_22 = 0.0;
  for (unsigned int x = 0; x < _restPose.size(); x++)
  {
    VEC3F oldX;
    oldX[0] = x0[x * 3]; oldX[1] = x0[x * 3 + 1]; oldX[2] = x0[x * 3 + 2];
    term21_22 += -2.0 * _restPose[x] * (R1R0 * oldX);
  }
  Real term23_24 = 0.0;
  for (unsigned int x = 0; x < _restPose.size(); x++)
    term23_24 += -2.0 * _restPose[x] * (R1T * translationOld);
  Real term25_26 = 2.0 * (p0 * x1);
  Real term27_28 = 0.0;
  for (unsigned int x = 0; x < _restPose.size(); x++)
    term27_28 += 2.0 * _restPose[x] * (R1T * translation);
  Real term29_30 = 0.0;
  for (unsigned int x = 0; x < _restPose.size(); x++)
  {
    VEC3F oldX;
    oldX[0] = x0[x * 3]; oldX[1] = x0[x * 3 + 1]; oldX[2] = x0[x * 3 + 2];
    term29_30 += 2.0 * oldX * (R0T * translationOld);
  }
  Real term31_32 = 0.0;
  for (unsigned int x = 0; x < _restPose.size(); x++)
  {
    VEC3F oldX;
    oldX[0] = x0[x * 3]; oldX[1] = x0[x * 3 + 1]; oldX[2] = x0[x * 3 + 2];
    term31_32 += -2.0 * oldX * (R0T * translation);
  }
  Real term33_34 = 0.0;
  for (unsigned int x = 0; x < _restPose.size(); x++)
  {
    VEC3F newX;
    newX[0] = x1[x * 3]; newX[1] = x1[x * 3 + 1]; newX[2] = x1[x * 3 + 2];
    term33_34 += -2.0 * newX * (R1T * translationOld);
  }
  Real term35_36 = 0.0;
  for (unsigned int x = 0; x < _restPose.size(); x++)
  {
    VEC3F newX;
    newX[0] = x1[x * 3]; newX[1] = x1[x * 3 + 1]; newX[2] = x1[x * 3 + 2];
    term35_36 += 2.0 * newX * (R1T * translation);
  }
  */

  // All reduced version
  MATRIX sandwich9_10  = _UTRU.transform(R0R1);
  VECTOR sandwich17_18 = _UTRp.vectorTransform(R1R0);
  VECTOR sandwich19_20 = _ITRp.vectorTransform(R0);
  VECTOR sandwich21_22 = _UTRp.vectorTransform(R0R1);
  VECTOR sandwich23_24 = _ITRp.vectorTransform(R1);
  VECTOR& sandwich27_28 = sandwich23_24;
  MATRIX sandwich29_30 = _ITRU.transform(R0);
  MATRIX& sandwich31_32 = sandwich29_30;
  MATRIX sandwich33_34 = _ITRU.transform(R1);
  MATRIX& sandwich35_36 = sandwich33_34;

  Real term1 = _pDot;
  Real term2 = qOld * qOld; 
  Real term3 = size * (translationOld * translationOld);
  Real term4 = _pDot;
  Real term5 = q * q;
  Real term6 = size * (translation * translation);
  Real term7_8 = -2.0 * _pTRp.scalarTransform(R0R1);
  Real term9_10 = -2.0 * (qOld * (sandwich9_10 * q));
  Real term11_12 = -2.0 * size * translation * translationOld;
  Real term13_14 = 2.0 * (qOld * _UTp); 
  Real term15_16 = 2.0 * (_ITp * translationOld);
  Real term17_18 = -2.0 * (q * sandwich17_18);
  Real term19_20 = -2.0 * (sandwich19_20 * vTranslation);
  Real term21_22 = -2.0 * (sandwich21_22 * qOld);
  Real term23_24 = -2.0 * (sandwich23_24 * vTranslationOld);
  Real term25_26 = 2.0 * (q * _UTp);
  Real term27_28 = 2.0 * (sandwich27_28 * vTranslation);
  Real term29_30 = 2.0 * (vTranslationOld * (sandwich29_30 * qOld));
  Real term31_32 = -2.0 * (vTranslation * (sandwich31_32 * qOld));
  Real term33_34 = -2.0 * (vTranslationOld * (sandwich33_34 * q));
  Real term35_36 = 2.0 * (vTranslation * (sandwich35_36 * q));

  Real finalSum = term1 + term2 + term3 + term4 + term5 + term6 +
                  term7_8 + term9_10 + term11_12 + term13_14 + term15_16 +
                  term17_18 + term19_20 + term21_22 + term23_24 + term25_26 +
                  term27_28 + term29_30 + term31_32 + term33_34 + term35_36;
 
  Real xOldSum = term1 + term2 + term3 + term13_14 + term29_30 +
                 2.0 * (sandwich19_20 * vTranslationOld);
                 //2.0 * (_UTp * qOld) + 
                 //2.0 * (sandwich19_20 * vTranslationOld) +
                 //2.0 * (vTranslationOld * (sandwich29_30 * qOld)); 

  /*
  VECTOR diffX = xOld - xNew;
  if (debug == 1)
  {
    cout << "***********************************************" << endl;
    cout << " 2 norm of rotation error: " << sqrt(diffX * diffX) << endl;
    cout << " Explicit 2 norm error:    " << sqrt(finalSum) << endl;

    cout << " 2 norm of xOld:          " << sqrt(xOld * xOld) << endl;
    cout << " Explicit 2 norm of xOld: " << sqrt(xOldSum) << endl;

    cout << " Inside sqrt: " << sqrt(finalSum / xOldSum) << endl;
    cout << " Outside sqrt: " << sqrt(finalSum) / sqrt(xOldSum) << endl;
    cout << "***********************************************" << endl;
  }

  Real finalNorm = diffX.norm2() / xOld.norm2();

  return finalNorm;
  */

  /*
  if (debug == 1)
  {
    cout << " finalSum: " << finalSum << endl;
    cout << " xOldSum:  " << xOldSum << endl;
  }
  */

  Real div = finalSum / xOldSum;

  if (div < 1e-8) return div;

  return sqrt(div);
}

//////////////////////////////////////////////////////////////////////
// Compute the rotation using a monolithic cross product
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_TET_MESH::computeCrossProductRotation()
{
  // just do a brute force here to see if it's worth trying in the
  // subspace
  VEC3F translation;
  Real scale = _UBasis(0,0);
  translation[0] = _q[0] * scale;
  translation[1] = _q[1] * scale;
  translation[2] = _q[2] * scale;
  
  int size = _restPose.size();
 
  // compute q cross f
  VEC3 sum;
  for (int x = 0; x < size; x++)
    //sum += _restPose[x] ^ (_vertices[x] - translation);
    sum += (_vertices[x] - translation) ^ _restPose[x];
  sum *= 1.0 / size;

  Real magnitude = norm(sum);
  Real theta = magnitude;
  if (magnitude < 1e-8) magnitude = 1.0;
  VEC3F a = sum / magnitude;
  //_rotation = MATRIX3::rotation(a, theta);
}

//////////////////////////////////////////////////////////////////////
// Compute the full rank rotation using shape matching
//////////////////////////////////////////////////////////////////////
MATRIX3 UNCONSTRAINED_SUBSPACE_TET_MESH::computeFullRankShapeMatchingRotation(bool centerCheck)
{
  int size = _restPose.size();
  TET_MESH::updateFullMesh();

  if (centerCheck)
  {
    // is the center of mass for _vertices being used?
    VEC3F vertexSum;
    for (int x = 0; x < size; x++)
      vertexSum += _vertices[x];
    vertexSum *= 1.0 / size;

    // if this is tripped, then the center of mass aren't centered at zero
    assert(vertexSum * vertexSum < 1e-4);
    //cout << " mean center of vertex positions: " << vertexSum << endl;

    VEC3F restSum;
    for (int x = 0; x < size; x++)
      restSum += _restPose[x];
    restSum *= 1.0 / size;
    
    // if this is tripped, then the center of mass aren't centered at zero
    assert(restSum * restSum < 1e-4);
    //cout << " sum of rest positions: " << restSum << endl;
  }

  // give a known rotation
  MATRIX3 known = MATRIX3::rotation(VEC3F(1,0,0), 90);
    
  // compute q cross f
  MATRIX3 sum;
  for (int x = 0; x < size; x++)
    sum += MATRIX3::outer_product(_restPose[x], _vertices[x]);
  sum *= 1.0 / size;
  sum = sum.transpose();

  // compute the final rotation
  MATRIX Apq(sum);
  MATRIX U(3,3);
  VECTOR S(3);
  MATRIX VT(3,3);
  Apq.SVD(U, S, VT);
  MATRIX R = U * VT;

  MATRIX3 rotation;
  for (int x = 0; x < 3; x++)
    for (int y = 0; y < 3; y++)
      rotation(x,y) = R(x,y);

  return rotation;
}

//////////////////////////////////////////////////////////////////////
// Compute the rotation using shape matching
//////////////////////////////////////////////////////////////////////
MATRIX3 UNCONSTRAINED_SUBSPACE_TET_MESH::computeShapeMatchingRotation()
{
  // unreduced version
#if 0
  int size = _restPose.size();
  /*
  // just do a brute force here to see if it's worth trying in the
  // subspace
  VEC3F translation;
  Real scale = _UBasis(0,0);
  translation[0] = _q[0] * scale;
  translation[1] = _q[1] * scale;
  translation[2] = _q[2] * scale;

  VEC3F centerOfMass;
  for (int x = 0; x < size; x++)
    centerOfMass += _vertices[x];
  centerOfMass *= 1.0 / size;

  //VEC3F centerOfMass = computeSubspaceCenterOfMass();

  cout << " translation:    " << translation << endl;
  cout << " center of mass: " << centerOfMass << endl;
  cout << " rigid translation: " << _translation << endl;
  */

  updateFullMesh();
 
  // is the center of mass for _vertices being used?
  VEC3F vertexSum;
  for (int x = 0; x < size; x++)
    vertexSum += _vertices[x];
  vertexSum *= 1.0 / size;

  // if this is tripped, then the center of mass aren't centered at zero
  assert(vertexSum * vertexSum < 1e-4);
  //cout << " mean center of vertex positions: " << vertexSum << endl;

  VEC3F restSum;
  for (int x = 0; x < size; x++)
    restSum += _restPose[x];
  restSum *= 1.0 / size;
  
  // if this is tripped, then the center of mass aren't centered at zero
  assert(restSum * restSum < 1e-4);
  //cout << " sum of rest positions: " << restSum << endl;

  // give a known rotation
  MATRIX3 known = MATRIX3::rotation(VEC3F(1,0,0), 90);
    
  // compute q cross f
  MATRIX3 sum;
  for (int x = 0; x < size; x++)
  {
    //sum += MATRIX3::outer_product(_restPose[x], _vertices[x]);
    //sum += MATRIX3::outer_product(_restPose[x], _rotation * (_vertices[x] - vertexSum));
    //sum += MATRIX3::outer_product(_restPose[x], _rotationDebug * _vertices[x]);
    sum += MATRIX3::outer_product(_restPose[x], _rotation * _vertices[x]);
    //sum += MATRIX3::outer_product(_restPose[x], _vertices[x]);
    //sum += MATRIX3::outer_product(_restPose[x], known * _restPose[x]);
    //sum += MATRIX3::outer_product(_rotation * _restPose[x], _rotation * (_vertices[x] - vertexSum));
    //sum += MATRIX3::outer_product(_rotationDebug * _restPose[x], _rotationDebug * _vertices[x]);
  }
  sum *= 1.0 / size;
  sum = sum.transpose();
#else

  /* 
  // initialize if this is the first time
  if (_AU.rows() == 0)
    initializeShapeMatchingRotation();

  MATRIX3 sum = _Ap;
  VECTOR displacement = _AU * _q;

  // the MATRIX3 constructor is unfortunately column-first
  sum += MATRIX3(displacement.data()).transpose();
  sum *= 1.0 / _restPose.size();
  */
  /*
  VECTOR pTripleTRp = _pTripleTRp.vectorTransform(_rotation);
  MATRIX pTripleTUR = _pTripleTRU.transform(_rotation);

  VECTOR displacement = pTripleTUR * _q;

  displacement += pTripleTRp;
  MATRIX3 sum(displacement.data());
  sum = sum.transpose();
  sum *= 1.0 / _restPose.size();
  */
#endif

  /*
  MATRIX Apq(sum);
  MATRIX U(3,3);
  VECTOR S(3);
  MATRIX VT(3,3);
  Apq.SVD(U, S, VT);
  MATRIX R = U * VT;

  MATRIX3 rotation;
  for (int x = 0; x < 3; x++)
    for (int y = 0; y < 3; y++)
      rotation(x,y) = R(x,y);
    
  //_rotation = rotationDelta * _rotation;
  //cout << " rotation old: " << rotationOld << endl;
  //cout << " diff: " << _rotationDebug - rotationOld << endl;
  //cout << " Known rotation: " << known << endl;
  //cout << " Diff: " << _rotationDebug - known << endl;
  
  return rotation;
  */
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " DEPRECATED. " << endl;
  exit(0);
}

//////////////////////////////////////////////////////////////////////
// compute shape matching matrices
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_TET_MESH::initializeShapeMatchingRotation()
{
  // compute Ap
  int size = _restPose.size();
  VECTOR Ap(9);
  for (int x = 0; x < size; x++)
  {
    Ap[0] += _restPose[x][0] * _restPose[x][0];
    Ap[1] += _restPose[x][0] * _restPose[x][1];
    Ap[2] += _restPose[x][0] * _restPose[x][2];

    Ap[3] += _restPose[x][1] * _restPose[x][0];
    Ap[4] += _restPose[x][1] * _restPose[x][1];
    Ap[5] += _restPose[x][1] * _restPose[x][2];
    
    Ap[6] += _restPose[x][2] * _restPose[x][0];
    Ap[7] += _restPose[x][2] * _restPose[x][1];
    Ap[8] += _restPose[x][2] * _restPose[x][2];
  }
  _Ap = MATRIX3(Ap.data());

  // compute AU
  _AU.resizeAndWipe(9, _q.size());

  for (int y = 0; y < _q.size(); y++)
    for (int x = 0; x < size; x++)
    {
      _AU(0, y) += _restPose[x][0] * _UBasis(3 * x, y);
      _AU(1, y) += _restPose[x][1] * _UBasis(3 * x, y);
      _AU(2, y) += _restPose[x][2] * _UBasis(3 * x, y);
      
      _AU(3, y) += _restPose[x][0] * _UBasis(3 * x + 1, y);
      _AU(4, y) += _restPose[x][1] * _UBasis(3 * x + 1, y);
      _AU(5, y) += _restPose[x][2] * _UBasis(3 * x + 1, y);
      
      _AU(6, y) += _restPose[x][0] * _UBasis(3 * x + 2, y);
      _AU(7, y) += _restPose[x][1] * _UBasis(3 * x + 2, y);
      _AU(8, y) += _restPose[x][2] * _UBasis(3 * x + 2, y);
    }
}

//////////////////////////////////////////////////////////////////////
// compute _crossMatrix
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_TET_MESH::initializeCrossProductRotation()
{
  /*
  cout << " Initializing cross product matrix ... ";
  flush(cout);
  int fullRank = _x.size();
  int totalVertices = fullRank / 3;
  
  // build big I
  MATRIX bigI(3, fullRank);
  for (int x = 0; x < fullRank; x++)
  {
    if (x % 3 == 0) bigI(0,x) = 1.0 / totalVertices;
    if (x % 3 == 1) bigI(1,x) = 1.0 / totalVertices;
    if (x % 3 == 2) bigI(2,x) = 1.0 / totalVertices;
  }

  // build big C
  MATRIX bigC(fullRank, fullRank);
  for (int x = 0; x < fullRank / 3; x++)
  {
    // build the cross product
    MATRIX3 cross = MATRIX3::cross(_restPose[x]);

    // patch it into bigC
    int x3 = 3 * x;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        bigC(x3 + i, x3 + j) = cross(i,j);
  }

  // build the big cross product matrix
  _crossMatrix = bigI * bigC;
  _crossMatrix = _crossMatrix * _UBasis;

  cout << "done." << endl;
  */
}

//////////////////////////////////////////////////////////////////////
// cache the subspace center of mass matrix
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_TET_MESH::cacheSubspaceCenterOfMass()
{
  SPARSE_MATRIX massI(3, _UBasis.rows());

  for (int x = 0; x < _UBasis.rows() / 3; x++)
  {
    massI(0, 3 * x) = mass(x);
    massI(1, 3 * x + 1) = mass(x);
    massI(2, 3 * x + 2) = mass(x);
  }

  VECTOR restVector(_UBasis.rows());
  //for (unsigned int x = 0; x < _restPose.size(); x++)
  for (int x = 0; x < _unconstrainedSize; x++)
  {
    restVector[3 * x] = _restPose[x][0];
    restVector[3 * x + 1] = _restPose[x][1];
    restVector[3 * x + 2] = _restPose[x][2];
  }

  _subspaceCenterOfMass = massI * _UBasis;
  _centerOfMassRest = massI * restVector;

  /*
  _subspaceCenterOfMass.resizeAndWipe(3, _q.size());

  for (int z = 0; z < _q.size(); z++)
    for (int y = 0; y < 3; y++)
      for (int x = 0; x < _UBasis.rows(); x++)
      {
        if ((x % 3) == y)
          _subspaceCenterOfMass(y,z) += _UBasis(x,z);
      }

  // get the total mass
  Real totalMass = 0;
  for (int x = 0; x < _UBasis.rows(); x++)
    totalMass += mass(x);

  _subspaceCenterOfMass *= 1.0 / totalMass;
  */
}

//////////////////////////////////////////////////////////////////////
// compute the subspace center of mass
//////////////////////////////////////////////////////////////////////
VEC3F UNCONSTRAINED_SUBSPACE_TET_MESH::computeSubspaceCenterOfMass()
{
  /*
  if (_subspaceCenterOfMass.rows() == 0)
  {
    _subspaceCenterOfMass.resizeAndWipe(3, _q.size());

    for (int z = 0; z < _q.size(); z++)
      for (int y = 0; y < 3; y++)
        for (int x = 0; x < _UBasis.rows(); x++)
        {
          if ((x % 3) == y)
            _subspaceCenterOfMass(y,z) += _UBasis(x,z);
        }

    int size = _restPose.size();
    _subspaceCenterOfMass *= 1.0 / size;
  }
  */

  VECTOR center(3); 
  //center = _subspaceCenterOfMass * _q;
  center = _subspaceCenterOfMass * _q + _centerOfMassRest;
  VEC3F toReturn(center[0], center[1], center[2]);
  return toReturn;
}

//////////////////////////////////////////////////////////////////////
// Draw only the surface faces, taking into account rigid components,
// but with no deformation
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_TET_MESH::drawRigidDebug()
{
#if USING_GLVU
  // This version doesn't trust OGL with anything --
  // Screw you and your coordinate system OGL!!!
  glPushMatrix();
  glBegin(GL_TRIANGLES);

  // verify the hard way
  for (unsigned int x = 0; x < _surfaceFaces.size(); x++)
  {
    TRIANGLE triangle = _tets[_surfaceFaces[x].first].face(_surfaceFaces[x].second);
   
    VEC3F v0 = *(restVersion(triangle.vertex(0)));
    VEC3F v1 = *(restVersion(triangle.vertex(1)));
    VEC3F v2 = *(restVersion(triangle.vertex(2)));

    v0 = _rotationDebug * v0;
    v1 = _rotationDebug * v1;
    v2 = _rotationDebug * v2;

    v0 += _translation;
    v1 += _translation;
    v2 += _translation;
    
    VEC3F normal = cross(v1 - v0, v2 - v0);
    normal.normalize();

    // PetSc doesn't like glNormal3fv
    glNormal3d(normal[0], normal[1], normal[2]);
    glVertex3dv(v0);
    glVertex3dv(v1);
    glVertex3dv(v2);
  }
  glEnd();
  glPopMatrix();
#endif
}

//////////////////////////////////////////////////////////////////////
// Draw a sphere at the center of mass
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_TET_MESH::drawCenterOfMass()
{
  // do a brute force recompute of the center of mass to be sure
  updateFullMesh();
  VEC3F vertexSum;
  int size = _vertices.size();
  for (int x = 0; x < size; x++)
    vertexSum += _vertices[x];
  vertexSum *= 1.0 / size;

#if USING_GLVU
  Real bounds[6];
  boundingBox(bounds);
  Real scale = bounds[1] - bounds[0];
  scale = (bounds[3] - bounds[2] > scale) ? bounds[3] - bounds[2] : scale;
  scale = (bounds[5] - bounds[4] > scale) ? bounds[5] - bounds[4] : scale;

  glPushMatrix();
    glTranslatef(_translation[0] + vertexSum[0], 
                 _translation[1] + vertexSum[1], 
                 _translation[2] + vertexSum[2]); 
    scale *= 0.05;
    glScalef(scale, scale, scale);

    glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
    glutSolidSphere(1.0, 10,10);
  glPopMatrix();
#endif
}


//////////////////////////////////////////////////////////////////////
// Draw axes at the center of mass
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_TET_MESH::drawRigidFrame()
{
#if USING_GLVU
  Real bounds[6];
  boundingBox(bounds);
  Real scale = bounds[1] - bounds[0];
  scale = (bounds[3] - bounds[2] > scale) ? bounds[3] - bounds[2] : scale;
  scale = (bounds[5] - bounds[4] > scale) ? bounds[5] - bounds[4] : scale;

  /*
  QUATERNION rotation(_rotation);
  Real angle;
  VEC3F axis;
  rotation.axisAngle(axis, angle);
  */
  Real angle;
  VEC3F axis;
  QUATERNION copy = _rotationQuaternion;
  copy.normalize();
  copy.axisAngle(axis, angle);

  glPushMatrix();
//#if EXCLUDE_TRANSLATIONS
  glTranslatef(_translation[0], _translation[1], _translation[2]);
  /*
#else
  VEC3F translation(_translation);
  Real U00 = _UBasis(0,0);
  Real U11 = _UBasis(1,1);
  Real U22 = _UBasis(2,2);
  translation[0] += _q[0] * U00;
  translation[1] += _q[1] * U11;
  translation[2] += _q[2] * U22;
  glTranslatef(translation[0], translation[1], translation[2]);
#endif
*/
  glRotatef(angle, axis[0], axis[1], axis[2]);
  scale *= 0.1;
  glScalef(scale, scale, scale);

  /*
  // draw a black outline
  glLineWidth(10.0f);
  glBegin(GL_LINES);
    // x axis is red
    glColor4f(0.0f, 0.0f, 0.0f, 0.1f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glColor4f(0.0f, 0.0f, 0.0f, 0.0f);
    glVertex3f(1.0f, 0.0f, 0.0f);

    // y axis is green 
    glColor4f(0.0f, 0.0f, 0.0f, 0.1f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glColor4f(0.0f, 0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 1.0f, 0.0f);
    
    // z axis is blue
    glColor4f(0.0f, 0.0f, 0.0f, 0.1f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glColor4f(0.0f, 0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 1.0f);
  glEnd();
  */

  glLineWidth(4.0f);
  glBegin(GL_LINES);
  /*
    // x axis is red
    glColor4f(10.0f, 0.0f, 0.0f, 10.1f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glColor4f(10.0f, 0.0f, 0.0f, 0.0f);
    glVertex3f(1.0f, 0.0f, 0.0f);

    // y axis is green 
    glColor4f(0.0f, 10.0f, 0.0f, 10.1f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glColor4f(0.0f, 10.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, -1.0f, 0.0f);
    
    // z axis is blue
    glColor4f(0.0f, 0.0f, 10.0f, 10.1f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glColor4f(0.0f, 0.0f, 10.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 1.0f);
    */
  glEnd();

  glPushMatrix();
    glRotatef(30, 0,1,0);
    {
    glPushMatrix();
    glColor4f(1,0,0,1);
    glRotatef(90, 1,0,0);
    glutSolidCone(0.2,2, 10, 10);
    glPopMatrix();

    glPushMatrix();
    glColor4f(0,0,1,1);
    glRotatef(90, 0,1,0);
    glutSolidCone(0.2,2, 10, 10);
    glPopMatrix();

    glPushMatrix();
    glColor4f(0,1,0,1);
    glRotatef(180, 1,0,0);
    glutSolidCone(0.2,2, 10, 10);
    glPopMatrix();
    }
  glPopMatrix();

  glPopMatrix();

  /*
  // draw translation
  glPointSize(100);
  glPushMatrix();
    glTranslatef(_translation[0], _translation[1], _translation[2]);
    cout << " translation: " << _translation << endl;
    glBegin(GL_POINTS);
      glColor4f(1,0,0,0);
      glVertex3f(0,0,0);
    glEnd();
  glPopMatrix();

  // draw original center of mass
  glPushMatrix();
    glTranslatef(_originalCenterOfMass[0], _originalCenterOfMass[1], _originalCenterOfMass[2]);
    cout << " original center of mass: " << _originalCenterOfMass << endl;
    glBegin(GL_POINTS);
      glColor4f(0,1,0,0);
      glVertex3f(0,0,0);
    glEnd();
  glPopMatrix();

  // draw rest center of mass
  glPushMatrix();
    glTranslatef(_restCenterOfMass[0], _restCenterOfMass[1], _restCenterOfMass[2]);

    cout << " rest center of mass: " << _restCenterOfMass << endl;
    glBegin(GL_POINTS);
      glColor4f(0,0,1,0);
      glVertex3f(0,0,0);
    glEnd();
  glPopMatrix();
  */

#endif
}

//////////////////////////////////////////////////////////////////////
// Draw only the surface faces, taking into account rigid components,
// but with no deformation
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_TET_MESH::drawRigidOnly()
{
  /*
  // just go ahead and estimate the most current rotation
  cout << " Estimating rotation and translation per draw" << endl;

  _x = _UBasis * _q;
  cout << " deformation norm: " << _x.norm2() << endl;

  //QUATERNION guess(_rotation); 
  QUATERNION guess(MATRIX3::I()); 
  QUATERNION reducedQuaternion = computeGeorgiiRotation();
  QUATERNION fullQuaternion = computeFullRankGeorgiiRotation(6);
  MATRIX3 shapeMatchingRotation = computeShapeMatchingRotation();
  QUATERNION matchingQuaternion(shapeMatchingRotation);

  cout << " full quaternion:    " << fullQuaternion << endl;
  cout << " reduced quaternion: " << reducedQuaternion << endl;
  cout << " shape matching:     " << matchingQuaternion << endl;

  MATRIX3 rotation = reducedQuaternion.toRotationMatrix();
  //MATRIX3 rotation = fullQuaternion.toRotationMatrix();
  
  //VEC3F translation = computeRigidTranslation();
  //translation += _originalCenterOfMass;
  VEC3F translation = _translation;
  */

#if USING_GLVU
  // This version doesn't trust OGL with anything --
  // Screw you and your coordinate system OGL!!!
  glPushMatrix();
  glBegin(GL_TRIANGLES);

  MATRIX3 rotation = _rotationQuaternion.toExplicitMatrix3x3();

  // verify the hard way
  for (unsigned int x = 0; x < _surfaceFaces.size(); x++)
  {
    TRIANGLE triangle = _tets[_surfaceFaces[x].first].face(_surfaceFaces[x].second);
   
    VEC3F v0 = *(restVersion(triangle.vertex(0)));
    VEC3F v1 = *(restVersion(triangle.vertex(1)));
    VEC3F v2 = *(restVersion(triangle.vertex(2)));

    v0 = rotation * v0;
    v1 = rotation * v1;
    v2 = rotation * v2;
    v0 += _translation;
    v1 += _translation;
    v2 += _translation;

    /*
    v0 = rotation * v0;
    v1 = rotation * v1;
    v2 = rotation * v2;
    v0 += translation;
    v1 += translation;
    v2 += translation;
    */
    
    VEC3F normal = cross(v1 - v0, v2 - v0);
    normal.normalize();

    // PetSc doesn't like glNormal3fv
    glNormal3d(normal[0], normal[1], normal[2]);
    glVertex3dv(v0);
    glVertex3dv(v1);
    glVertex3dv(v2);
  }
  glEnd();
  glPopMatrix();
  
  /*
  // now draw without the rotation
  glColor4f(1,1,1,1);

  // This version doesn't trust OGL with anything --
  // Screw you and your coordinate system OGL!!!
  glPushMatrix();
  glBegin(GL_TRIANGLES);

  // verify the hard way
  for (int x = 0; x < _surfaceFaces.size(); x++)
  {
    TRIANGLE triangle = _tets[_surfaceFaces[x].first].face(_surfaceFaces[x].second);
   
    VEC3F v0 = *(restVersion(triangle.vertex(0)));
    VEC3F v1 = *(restVersion(triangle.vertex(1)));
    VEC3F v2 = *(restVersion(triangle.vertex(2)));

    v0 += translation;
    v1 += translation;
    v2 += translation;
    
    VEC3F normal = cross(v1 - v0, v2 - v0);
    normal.normalize();

    // PetSc doesn't like glNormal3fv
    glNormal3d(normal[0], normal[1], normal[2]);
    glVertex3dv(v0);
    glVertex3dv(v1);
    glVertex3dv(v2);
  }
  glEnd();
  glPopMatrix();
  */

#endif
}

//////////////////////////////////////////////////////////////////////
// Draw all the tets
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_TET_MESH::drawAllTets()
{
  // colors for the different partitions
  float colors[8][4];

  // cache the colors
  colors[0][0] = 1.0f; colors[0][1] = 0.0f; colors[0][2] = 0.0f; colors[0][3] = 1.0f;
  colors[1][0] = 0.0f; colors[1][1] = 1.0f; colors[1][2] = 0.0f; colors[1][3] = 1.0f;
  colors[2][0] = 0.0f; colors[2][1] = 0.0f; colors[2][2] = 1.0f; colors[2][3] = 1.0f;
  colors[3][0] = 1.0f; colors[3][1] = 0.0f; colors[3][2] = 1.0f; colors[3][3] = 1.0f;
  colors[4][0] = 0.0f; colors[4][1] = 1.0f; colors[4][2] = 1.0f; colors[4][3] = 1.0f;
  colors[5][0] = 1.0f; colors[5][1] = 1.0f; colors[5][2] = 0.0f; colors[5][3] = 1.0f;
  colors[6][0] = 1.0f; colors[6][1] = 1.0f; colors[6][2] = 1.0f; colors[6][3] = 1.0f;
  colors[7][0] = 0.5f; colors[7][1] = 0.5f; colors[7][2] = 0.5f; colors[7][3] = 1.0f;
  srand(123456);

#if USING_GLVU
  // This version doesn't trust OGL with anything --
  // Screw you and your coordinate system OGL!!!
  glPushMatrix();
  glBegin(GL_TRIANGLES);

  // verify the hard way
  MATRIX3 rotation = _rotationQuaternion.toExplicitMatrix3x3();
  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    int randColor = 8 * ((Real)rand() / RAND_MAX);
    int i = randColor % 8;
    glColor4f(colors[i][0], colors[i][1], colors[i][2], 1.0);

    //if ((*(_tets[x].vertices[0]))[2] <= 0.0) continue;

    for (int y = 0; y < 4; y++)
    {
      TRIANGLE triangle = _tets[x].face(y);
     
      VEC3F v0 = *(triangle.vertex(0));
      VEC3F v1 = *(triangle.vertex(1));
      VEC3F v2 = *(triangle.vertex(2));
     
      //v0 = _rotation * v0;
      //v1 = _rotation * v1;
      //v2 = _rotation * v2;

      v0 = rotation * v0;
      v1 = rotation * v1;
      v2 = rotation * v2;
    
      v0 += _translation;
      v1 += _translation;
      v2 += _translation;

      VEC3F normal = cross(v1 - v0, v2 - v0);
      normal.normalize();

      // PetSc doesn't like glNormal3fv
      glNormal3d(normal[0], normal[1], normal[2]);
      glVertex3dv(v0);
      glVertex3dv(v1);
      glVertex3dv(v2);
    }
  }
  glEnd();
  glPopMatrix();
#endif
}

//////////////////////////////////////////////////////////////////////
// Draw only the surface faces
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_TET_MESH::drawSurfaceFaces()
{
#if USING_GLVU
  // This version doesn't trust OGL with anything --
  // Screw you and your coordinate system OGL!!!
  glPushMatrix();
  glBegin(GL_TRIANGLES);

  // verify the hard way
  for (unsigned int x = 0; x < _surfaceFaces.size(); x++)
  {
    TRIANGLE triangle = _tets[_surfaceFaces[x].first].face(_surfaceFaces[x].second);
   
    VEC3F v0 = *(triangle.vertex(0));
    VEC3F v1 = *(triangle.vertex(1));
    VEC3F v2 = *(triangle.vertex(2));
   
    //v0 = _rotation * v0;
    //v1 = _rotation * v1;
    //v2 = _rotation * v2;

    MATRIX3 rotation = _rotationQuaternion.toExplicitMatrix3x3();
    v0 = rotation * v0;
    v1 = rotation * v1;
    v2 = rotation * v2;
  
    v0 += _translation;
    v1 += _translation;
    v2 += _translation;
    
    VEC3F normal = cross(v1 - v0, v2 - v0);
    normal.normalize();

    // PetSc doesn't like glNormal3fv
    glNormal3d(normal[0], normal[1], normal[2]);
    glVertex3dv(v0);
    glVertex3dv(v1);
    glVertex3dv(v2);
  }
  glEnd();
  glPopMatrix();
#endif
}

//////////////////////////////////////////////////////////////////////
// Given a position, return the closest surface node
//////////////////////////////////////////////////////////////////////
VEC3F* UNCONSTRAINED_SUBSPACE_TET_MESH::closestSurfaceNode(VEC3F point)
{
  // make sure a surface list was built 
  if (_surfaceVertices.size() == 0) return NULL;

  MATRIX3 rotation = _rotationQuaternion.toExplicitMatrix3x3();

  // tentatively set it to the first in the list
  VEC3F* closest;
  closest = _surfaceVertices[0];
  VEC3F transformed = rotation * (*_surfaceVertices[0]) + _translation;
  Real minDistance = norm2(point - transformed);

  // loop through the rest of the vertices
  for (unsigned int x = 1; x < _surfaceVertices.size(); x++)
  {
    // check if this one is closer
    VEC3F transformed = rotation * (*_surfaceVertices[x]) + _translation;

    Real distance = norm2(point - transformed);
    if (distance < minDistance)
    {
      minDistance = distance;
      closest = _surfaceVertices[x];
    }
  }
  return closest;
}

//////////////////////////////////////////////////////////////////////////////
// Eqn. 4 from the Georgii and Westermann paper
//////////////////////////////////////////////////////////////////////////////
QUATERNION UNCONSTRAINED_SUBSPACE_TET_MESH::fullRotationDistance(QUATERNION& q, VEC3F& rest, VEC3F& displacement)
{
  //QUATERNION rhs(rest - _restCenterOfMass);
  QUATERNION rhs(rest);
  //QUATERNION lhs(rest + displacement - _centerOfMass);

  MATRIX3 rotation = _rotationQuaternion.toExplicitMatrix3x3();
  QUATERNION lhs(rotation * (rest + displacement));
  lhs = q * lhs * q.conjugate();

  return lhs - rhs;
}

//////////////////////////////////////////////////////////////////////////////
// Second derivatives from Eqn. 7 from the Georgii and Westermann paper
//////////////////////////////////////////////////////////////////////////////
QUATERNION UNCONSTRAINED_SUBSPACE_TET_MESH::fullDistanceHessian(QUATERNION& q, QUATERNION& p, int i, int j)
{
  //dxdx =
  if (i == 0 && j == 0)
    return 2 * QUATERNION( p.x(), -p.y(), -p.z(), p.w());

  //dydy = 
  if (i == 1 && j == 1)
    return 2 * QUATERNION( -p.x() , p.y(), -p.z(), p.w());

  //dzdz = 
  if (i == 2 && j == 2)
    return 2 * QUATERNION( -p.x() , -p.y(), p.z(), p.w());

  //dwdw = 
  if (i == 3 && j == 3)
    return 2 * QUATERNION( p.x() , p.y(), p.z(), p.w());

  //dxdy = 
  if ((i == 0 && j == 1) || (i == 1 && j == 0))
    return 2 * QUATERNION( p.y() , p.x(), 0, 0);

  //dxdz = 
  if ((i == 0 && j == 2) || (i == 2 && j == 0))
    return 2 * QUATERNION( p.z(), 0, p.x(), 0);

  //dxdw= 
  if ((i == 0 && j == 3) || (i == 3 && j == 0))
    return 2 * QUATERNION( 0 , -p.z(), p.y(), 0);

  //dydz= 
  if ((i == 1 && j == 2) || (i == 2 && j == 1))
    return 2 * QUATERNION( 0 , p.z(), p.y(), 0);

  //dydw= 
  if ((i == 1 && j == 3) || (i == 3 && j == 1))
    return 2 * QUATERNION( p.z() , 0, -p.x(), 0);

  //dzdw= 
  if ((i == 2 && j == 3) || (i == 3 && j == 2))
    return 2 * QUATERNION( -p.y(), p.x(), 0, 0);

  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " Missed a Hessian case! i = " << i << " j = " << j << endl;
  exit(1);
}

//////////////////////////////////////////////////////////////////////////////
// First derivatives from Eqns. 5 and 7 from the Georgii and Westermann paper
//////////////////////////////////////////////////////////////////////////////
QUATERNION UNCONSTRAINED_SUBSPACE_TET_MESH::fullDistanceGradient(QUATERNION& q, QUATERNION& p, int which)
{
  switch (which)
  {
    case 0:
       return 2 * QUATERNION( q.z()*p.z()+q.y()*p.y()+q.x()*p.x(), 
                             -q.x()*p.y()+q.y()*p.x()-q.w()*p.z(),
                              q.z()*p.x()-q.x()*p.z()+q.w()*p.y(),
                              p.w()*q.x());
        break;

    case 1:
        return 2 * QUATERNION( q.x()*p.y()-q.y()*p.x()+q.w()*p.z(),
                               q.z()*p.z()+q.y()*p.y()+q.x()*p.x(),
                              -q.y()*p.z()+q.z()*p.y()-q.w()*p.x(),
                               p.w()*q.y());
        break;

    case 2:
        return 2 * QUATERNION( -q.z()*p.x()+q.x()*p.z()-q.w()*p.y(),
                                q.y()*p.z()-q.z()*p.y()+q.w()*p.x(),
                                q.z()*p.z()+q.y()*p.y()+q.x()*p.x(),
                                p.w()*q.z());
        break;
    case 3:
        return 2 * QUATERNION( q.y()*p.z()-q.z()*p.y()+q.w()*p.x(),
                               q.z()*p.x()-q.x()*p.z()+q.w()*p.y(),
                               q.x()*p.y()-q.y()*p.x()+q.w()*p.z(),
                               q.w()*p.w());
        break;
  }

  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " distanceGradient has been called incorrectly!" << endl;
  exit(1);
  return QUATERNION();
}

//////////////////////////////////////////////////////////////////////////////
// Eqns. 7 and 8  from the Georgii and Westermann paper
//////////////////////////////////////////////////////////////////////////////
MATRIX UNCONSTRAINED_SUBSPACE_TET_MESH::fullEnergyHessian(QUATERNION& q, Real lambda)
#if 0
{
  MATRIX finalHessian(5,5);

  // build a big rest vector
  VECTOR restVector(3 * _vertices.size());
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    restVector[3 * x]     = _restPose[x][0];
    restVector[3 * x + 1] = _restPose[x][1];
    restVector[3 * x + 2] = _restPose[x][2];
  }
  // build a big center of mass vector
  VECTOR restCenter(3 * _vertices.size());
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    restCenter[3 * x]     = _restCenterOfMass[0];
    restCenter[3 * x + 1] = _restCenterOfMass[1];
    restCenter[3 * x + 2] = _restCenterOfMass[2];
  }
  // build a big center of mass vector
  VECTOR deformedCenter(3 * _vertices.size());
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    deformedCenter[3 * x]     = _centerOfMass[0];
    deformedCenter[3 * x + 1] = _centerOfMass[1];
    deformedCenter[3 * x + 2] = _centerOfMass[2];
  }

  // compute d - this can be computed outside the loop
  //_x = _UBasis * _q;
  VECTOR deformed = (restVector + _x);
  deformed -= deformedCenter;
  VECTOR rest = restVector - restCenter;

  // for each entry in the Hessian
  for (int x = 0; x < 4; x++)
    for (int y = 0; y < 4; y++)
    {
      MATRIX rotation = q.toExplicitMatrix();
      MATRIX paddedRotation(4, 3);
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          paddedRotation(i,j) = rotation(i,j);

      MATRIX paddedI(4, 3);
      paddedI(0,0) = 1;
      paddedI(1,1) = 1;
      paddedI(2,2) = 1;

      BLOCK_MATRIX blockRotation(paddedRotation, _vertices.size());
      BLOCK_MATRIX blockI(paddedI, _vertices.size());

      VECTOR d = (blockRotation * deformed);
      d -= (blockI * rest);
      
      // compute p - this can be computed outside the loop
      VECTOR p = blockI * deformed;

      // computer partialQi and partialQj
      MATRIX qi = gradientMatrix(q, x);
      MATRIX qj = gradientMatrix(q, y);
      
      BLOCK_MATRIX blockQi(qi, _vertices.size());
      BLOCK_MATRIX blockQj(qj, _vertices.size());

      VECTOR partialDpartialQi = blockQi * p;
      VECTOR partialDpartialQj = blockQj * p;

      // compute partialQiQj
      MATRIX qiqj = hessianMatrix(q, x,y);
      BLOCK_MATRIX blockQiQj(qiqj, _vertices.size());
      VECTOR partialDpartialQiQj = blockQiQj * p;

      // combine into the final
      Real final = 0.0; 
      final += partialDpartialQi ^ partialDpartialQj;
      final += partialDpartialQiQj ^ d;

      if (x == y)
        final += 2 * lambda * _vertices.size();

      finalHessian(x,y) = 2 * final;
    }

  // fill out the Lagrange multiplier entries
  finalHessian(0,4) = finalHessian(4,0) = 2 * q[0];
  finalHessian(1,4) = finalHessian(4,1) = 2 * q[1];
  finalHessian(2,4) = finalHessian(4,2) = 2 * q[2];
  finalHessian(3,4) = finalHessian(4,3) = 2 * q[3];
  finalHessian(4,4) = 0.0;

  return finalHessian;
}
#else
{
  MATRIX finalHessian(5,5);

  MATRIX3 rotation = _rotationQuaternion.toExplicitMatrix3x3();

  // for each entry in the Hessian
  for (int x = 0; x < 4; x++)
    for (int y = 0; y < 4; y++)
    {
      // for each vertex
      for (unsigned int z = 0; z < _vertices.size(); z++)
      {
        VEC3F diff = _vertices[z] - _restPose[z];

        QUATERNION d = fullRotationDistance(q, _restPose[z], diff);
        //QUATERNION p(_vertices[z] - _centerOfMass);
        QUATERNION p(rotation * _vertices[z]);

        QUATERNION partialDpartialQi = fullDistanceGradient(q, p, x);
        QUATERNION partialDpartialQj = fullDistanceGradient(q, p, y);
        QUATERNION partialDpartialQiQj = fullDistanceHessian(q, p, x, y);

        Real final = 0.0;
        final += partialDpartialQi ^ partialDpartialQj;
        final += partialDpartialQiQj ^ d;

        if (x == y)
          final += 2 * lambda;

        finalHessian(x,y) += final;
      }

      // take into account the 2 in front of the integral
      finalHessian(x,y) *= 2.0;
    }

  // fill out the Lagrange multiplier entries
  finalHessian(0,4) = finalHessian(4,0) = 2 * q[0];
  finalHessian(1,4) = finalHessian(4,1) = 2 * q[1];
  finalHessian(2,4) = finalHessian(4,2) = 2 * q[2];
  finalHessian(3,4) = finalHessian(4,3) = 2 * q[3];
  finalHessian(4,4) = 0.0;

  return finalHessian;
}
#endif

//////////////////////////////////////////////////////////////////////////////
// Eqn. 5 and 6 from the Georgii and Westermann paper
//////////////////////////////////////////////////////////////////////////////
VECTOR UNCONSTRAINED_SUBSPACE_TET_MESH::projectedEnergyGradient(QUATERNION& q, Real lambda)
{
  QUATERNION finalGradient;

  VECTOR restCenter3(3);
  restCenter3[0] = _restCenterOfMass[0];
  restCenter3[1] = _restCenterOfMass[1];
  restCenter3[2] = _restCenterOfMass[2];
  VECTOR center3(3);
  center3[0] = _centerOfMass[0];
  center3[1] = _centerOfMass[1];
  center3[2] = _centerOfMass[2];
  VECTOR centerDiff = center3 - restCenter3;

  // build quaternion-friendly U
  BLOCK_MATRIX qU;
  buildQU(qU);
  BLOCK_MATRIX blockU;
  buildBlockU(blockU);
  BLOCK_MATRIX qUtranspose = qU.transpose();

  VECTOR xBar = _projectedRestPose;
  VECTOR uBar = _q;
  VECTOR xPlusU = xBar + uBar;

  // for each component of the quaternion
  for (int x = 0; x < 4; x++)
  {
    MATRIX rotation = q.toExplicitMatrix();
    MATRIX paddedRotation(4, 3);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        paddedRotation(i,j) = rotation(i,j);

    MATRIX paddedI(4, 3);
    paddedI(0,0) = 1;
    paddedI(1,1) = 1;
    paddedI(2,2) = 1;

    BLOCK_MATRIX Mq(paddedRotation, _vertices.size());
    BLOCK_MATRIX IB(paddedI, _vertices.size());

    BLOCK_MATRIX IV(_vertices.size(), 1);
    for (unsigned int y = 0; y < _vertices.size(); y++)
      IV.add(paddedI, y, 0);

    MATRIX dq = gradientMatrix(q, x);
    BLOCK_MATRIX Qi(dq, _vertices.size());

    // compute the lhs term
    BLOCK_MATRIX IBU = IB * blockU;
    VECTOR IBUxu = IBU * xPlusU;
    VECTOR IVcd = IV * center3;

    VECTOR term1 = Qi * IBUxu;
    VECTOR term2 = Qi * IVcd;

    // compute the rhs term
    BLOCK_MATRIX MqU = Mq * blockU;
    VECTOR term3 = MqU * xPlusU;
    VECTOR term4 = IV * centerDiff;
    VECTOR term5 = IBU * xBar;

    Real termA = term1 ^ term3;
    Real termB = term1 ^ term4;
    Real termC = term1 ^ term5;

    Real termD = term2 ^ term3;
    Real termE = term2 ^ term4;
    Real termF = term2 ^ term5;

    finalGradient[x] = termA - termB - termC - termD + termE + termF; 
    finalGradient[x] += 2 * lambda * q[x] * _vertices.size();
  }

  // take into account the 2 in front of the integral
  finalGradient *= 2.0;

  VECTOR finalVector(5);
  finalVector[0] = finalGradient[0];
  finalVector[1] = finalGradient[1];
  finalVector[2] = finalGradient[2];
  finalVector[3] = finalGradient[3];
  finalVector[4] = (q ^ q) - 1;
  return finalVector;
}

//////////////////////////////////////////////////////////////////////////////
// Eqn. 5 and 6 from the Georgii and Westermann paper
//////////////////////////////////////////////////////////////////////////////
VECTOR UNCONSTRAINED_SUBSPACE_TET_MESH::fullEnergyGradient(QUATERNION& q, Real lambda)
#if 0
{
  QUATERNION finalGradient;

  // build a big rest vector
  VECTOR restVector(3 * _vertices.size());
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    restVector[3 * x]     = _restPose[x][0];
    restVector[3 * x + 1] = _restPose[x][1];
    restVector[3 * x + 2] = _restPose[x][2];
  }
  // build a big center of mass vector
  VECTOR restCenter(3 * _vertices.size());
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    restCenter[3 * x]     = _restCenterOfMass[0];
    restCenter[3 * x + 1] = _restCenterOfMass[1];
    restCenter[3 * x + 2] = _restCenterOfMass[2];
  }
  // build a big center of mass vector
  VECTOR deformedCenter(3 * _vertices.size());
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    deformedCenter[3 * x]     = _centerOfMass[0];
    deformedCenter[3 * x + 1] = _centerOfMass[1];
    deformedCenter[3 * x + 2] = _centerOfMass[2];
  }
  
  // precompute recentered versions
  //_x = _UBasis * _q;
  VECTOR deformed = (restVector + _x);
  deformed -= deformedCenter;
  VECTOR rest = restVector - restCenter;

  // for each component of the quaternion
  for (int x = 0; x < 4; x++)
  {
    // compute d
    MATRIX rotation = q.toExplicitMatrix();
    MATRIX paddedRotation(4, 3);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        paddedRotation(i,j) = rotation(i,j);

    MATRIX paddedI(4, 3);
    paddedI(0,0) = 1;
    paddedI(1,1) = 1;
    paddedI(2,2) = 1;

    BLOCK_MATRIX blockRotation(paddedRotation, _vertices.size());
    BLOCK_MATRIX blockI(paddedI, _vertices.size());

    VECTOR d = (blockRotation * deformed);
    d -= (blockI * rest);

    // compute p
    VECTOR p = blockI * deformed;

    // compute the partial
    MATRIX dq = gradientMatrix(q, x);
    BLOCK_MATRIX blockDq(dq, _vertices.size());
    VECTOR partialDpartialQ = blockDq * p;

    finalGradient[x] = partialDpartialQ ^ d;
    
    finalGradient[x] += 2 * lambda * q[x] * _vertices.size();
  }

  // take into account the 2 in front of the integral
  finalGradient *= 2.0;

  VECTOR finalVector(5);
  finalVector[0] = finalGradient[0];
  finalVector[1] = finalGradient[1];
  finalVector[2] = finalGradient[2];
  finalVector[3] = finalGradient[3];
  finalVector[4] = (q ^ q) - 1;
  return finalVector;
}
#else
{
  QUATERNION finalGradient;
  //VEC3F centerOfMass = _centerOfMass;
  //VEC3F centerOfMass = _translation;
  VEC3F centerOfMass; // should actually be zero
  vector<VEC3F>& vertices = _vertices;
  vector<VEC3F>& restVertices = _restPose;
  MATRIX3 rotation = _rotationQuaternion.toExplicitMatrix3x3();

  // for each component of the quaternion
  //Real dx = 1.0 / vertices.size();
  for (int x = 0; x < 4; x++)
  {
    // for each vertex
    for (unsigned int y = 0; y < vertices.size(); y++)
    {
      VEC3F diff   = vertices[y] - restVertices[y];
      QUATERNION d = fullRotationDistance(q, restVertices[y], diff);
      //QUATERNION p(vertices[y] - centerOfMass);
      QUATERNION p(rotation * vertices[y]);

      QUATERNION partialDpartialQ = fullDistanceGradient(q, p, x);

      //Real finalComponent = (partialDpartialQ ^ d) * dx + 2 * lambda * q[x];
      Real finalComponent = (partialDpartialQ ^ d);
      //finalComponent *= dx;
      finalComponent += 2 * lambda * q[x];
      finalGradient[x] += finalComponent;
    }
  }

  // take into account the 2 in front of the integral
  finalGradient *= 2.0;

  VECTOR finalVector(5);
  finalVector[0] = finalGradient[0];
  finalVector[1] = finalGradient[1];
  finalVector[2] = finalGradient[2];
  finalVector[3] = finalGradient[3];
  finalVector[4] = (q ^ q) - 1;
  return finalVector;
}
#endif

//////////////////////////////////////////////////////////////////////////////
// Eqn. 4 from the Georgii and Westermann paper
//////////////////////////////////////////////////////////////////////////////
Real UNCONSTRAINED_SUBSPACE_TET_MESH::fullRotationEnergy(QUATERNION& q, Real lambda)
{
  Real final = 0;
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    VEC3F diff = _vertices[x] - _restPose[x];
    QUATERNION distance = fullRotationDistance(q, _restPose[x], diff);
    final += (distance ^ distance) + lambda * ((q ^ q) - 1);
  }
  final *= 1.0 / _vertices.size();
  return final;
}

//////////////////////////////////////////////////////////////////////////////
// compute the rigid rotation according to Georgii and Westermann's
// "Corotated Finite Elements Made Fast and Stable"
//////////////////////////////////////////////////////////////////////////////
QUATERNION UNCONSTRAINED_SUBSPACE_TET_MESH::computeFullRankGeorgiiRotation(unsigned int digits, QUATERNION guess, bool useQ)
{
  /*
  static int counter = 0;
  counter++;
  cout << " =========================== " << endl;
  cout << " Full Georgii " << counter <<endl;
  cout << " =========================== " << endl;
  */
  
  // might as well warm start it with the shape matching solution
  //MATRIX3 guess = computeShapeMatchingRotation();
  //QUATERNION q(guess);

  //MATRIX3 I = MATRIX3::I();
  //QUATERNION q(0,0,0,1);
  Real lambda = 0;
  QUATERNION q(guess);

  // make sure the entire mesh knows about _q
  if (useQ)
  {
    _x = _UBasis * _q;
    updateFullMesh();
  }
  // or, if the deformation is in x, make sure that the center of mass
  // is clamped to the origin for the solve
  else
  {
    TET_MESH::updateFullMesh();
    TET_MESH::centerOfMassToOrigin();
  }

  VECTOR gradient = fullEnergyGradient(q, lambda);
  //VECTOR gradient = projectedEnergyGradient(q, lambda);
  //Real magnitude = (gradient ^ gradient) - gradient[4] * gradient[4];

  int count = 0;
  Real convergence = 1;
  Real threshold = pow(10.0, -(Real)digits);
  while (convergence > threshold && count < 1000)
  {
    MATRIX hessian = fullEnergyHessian(q, lambda);
    VECTOR update = gradient;
    hessian.factorLU();
    hessian.solveLU(update);

    q[0]   -= update[0];
    q[1]   -= update[1];
    q[2]   -= update[2];
    q[3]   -= update[3];
    
    // this seems to fix everything -- otherwise lambda either
    // explodes or the number of iterations necessary becomes large because
    // the update step size needs to be tuned down
    //
    // if quaternion magnitude is in fact 0, then the lambda solution
    // becomes non-unique, so it goes nuts.
    //
    // Removing the lambda from the derivatives entirely ruins the convergence
    // as well though, so the hack stays
    //lambda -= update[4];
    
    gradient = fullEnergyGradient(q, lambda);
    //gradient = projectedEnergyGradient(q, lambda);
    //cout << " gradient: " << gradient << endl;
    //cout << " update: " << update << endl;

    // see what the convergence looks like without lambda
    update[4] = 0.0;
    convergence = update ^ update;
    count++;
  }
  _georgiiIterations = count;
  _totalGeorgiiIterations += count;
  q.normalize();
  q.negateIm();
  //cout << " Full Georgii iterations: " << count << endl;

  return q;
}

//////////////////////////////////////////////////////////////////////////////
// compute the rigid rotation according to Georgii and Westermann's
// "Corotated Finite Elements Made Fast and Stable"
//////////////////////////////////////////////////////////////////////////////
QUATERNION UNCONSTRAINED_SUBSPACE_TET_MESH::computeGeorgiiRotation(int debug)
{
  // might as well warm start it with the shape matching solution
  //MATRIX3 guess = computeShapeMatchingRotation();
  // or, warm start it with the previous rotation (probably better)
  //QUATERNION q(_rotation);
  //QUATERNION q(computeShapeMatchingRotation());
 
  /* 
  MATRIX3 shapeMatching = computeShapeMatchingRotation();
  MATRIX3 guess = _rotation.transpose() * shapeMatching;
  QUATERNION q(guess);
  */

  /*
  static int counter = 0;
  counter++;
  cout << " =========================== " << endl;
  cout << " Reduced Georgii " << counter << endl;
  cout << " =========================== " << endl;
  */

  MATRIX3 I = MATRIX3::I();
  QUATERNION q(I);
  //QUATERNION q(0,0,0,1);
  Real lambda = 0;

  VECTOR gradient = energyGradient(q, lambda);

  int count = 0;
  Real convergence = 1;
  //while (convergence > 1e-8 && count < 100)
  while (convergence > 1e-6 && count < 100)
  //while (convergence > 1e-4 && count < 100)
  //while (convergence > 1e-2 && count < 100)
  {
    MATRIX hessian = energyHessian(q, lambda);
    //cout << " hessian: " << hessian << endl;
    VECTOR update = gradient;
    hessian.factorLU();
    hessian.solveLU(update);

    q[0]   -= update[0];
    q[1]   -= update[1];
    q[2]   -= update[2];
    q[3]   -= update[3];
    
    // this seems to fix everything -- otherwise lambda either
    // explodes or the number of iterations necessary becomes large because
    // the update step size needs to be tuned down
    //
    // if quaternion magnitude is in fact 0, then the lambda solution
    // becomes non-unique, so it goes nuts.
    //
    // Removing the lambda from the derivatives entirely ruins the convergence
    // as well though, so the hack stays
    //lambda -= update[4];
    
    gradient = energyGradient(q, lambda);
    //cout << "gradient: " << gradient << endl;

    // see what the convergence looks like without lambda
    update[4] = 0.0;
    convergence = update ^ update;
    count++;
  }

  _totalSteps++;

  _georgiiIterations = count;
  _totalGeorgiiIterations += _georgiiIterations;
  if (_georgiiIterations > _maxGeorgiiIterations)
    _maxGeorgiiIterations = _georgiiIterations;
  q.normalize();
  //if (debug == 4)
  {
    cout << " reduced Georgii iterations: " << count << " convergence: " << convergence << endl;
    //cout << " Shape matching guess: " << guess << endl;
    //cout << " Reduced Georgii result: " << q.toRotationMatrix() << endl;
    //cout << " Total shape matching solution: " << shapeMatching << endl;
  }
  //cout << " Mean Georgii: " << _totalGeorgiiIterations / (float)_totalSteps ;
  //cout << " \tMax Georgii: " << _maxGeorgiiIterations << endl;

  // undo the transpose
  q.negateIm();

  return q;
}

//////////////////////////////////////////////////////////////////////////////
// Eqns. 7 and 8  from the Georgii and Westermann paper
//////////////////////////////////////////////////////////////////////////////
MATRIX UNCONSTRAINED_SUBSPACE_TET_MESH::energyHessian(QUATERNION& q, Real lambda)
{
  MATRIX finalHessian(5,5);

  VECTOR centerOfMass(3);
  VECTOR restCenterOfMass(3);
  for (int x = 0; x < 3; x++)
  {
    //centerOfMass[x] = _centerOfMass[x];
    //centerOfMass[x] = _translation[x];
    //restCenterOfMass[x] = _restCenterOfMass[x];
  }

  VECTOR centerDiff = centerOfMass - restCenterOfMass;

  // generate Mq, the q full quaternion multiply
  MATRIX tmp = q.toExplicitMatrix();
  MATRIX Mq(4,3);
  for (int x = 0; x < 4; x++)
    for (int y = 0; y < 3; y++)
      Mq(x,y) = tmp(x,y);

  VECTOR& xBar = _projectedRestPose;
  VECTOR& uBar = _q;
  VECTOR xPlusU = xBar + uBar;

  // create a padded rotation
  MATRIX3 rotation = _rotationQuaternion.toExplicitMatrix3x3();
  MATRIX R4x4(4,4);
  for (int x = 0; x < 3; x++)
    for (int y = 0; y < 3; y++)
      R4x4(x,y) = rotation(x,y);

  // precache all the Q gradients
  MATRIX Qi[4];
  for (int x = 0; x < 4; x++)
  {
    Qi[x] = gradientMatrix(q, x);

    // add the rotation this way
    Qi[x] = Qi[x] * R4x4;
  }

  // for each entry in the Hessian
  for (int x = 0; x < 4; x++)
    for (int y = 0; y < 4; y++)
    {
      Real final = 0.0;

      // compute the product of partials
      //MATRIX QiQj = Qi[x] ^ Qi[y];
      MATRIX QiQj = Qi[x] ^ Qi[y];
      MATRIX IBU_IBU = _IBU_IBU.transform(QiQj);
      //MATRIX IBU_IV = _IBU_IV.transform(QiQj);
      //MATRIX IV_IBU = _IV_IBU.transform(QiQj);
      //MATRIX IV_IV = _IV_IV.transform(QiQj);

      final += xPlusU ^ (IBU_IBU * xPlusU);
      //final -= xPlusU ^ (IBU_IV * centerOfMass);
      //final -= centerOfMass ^ (IV_IBU * xPlusU);
      //final += centerOfMass ^ (IV_IV * centerOfMass);

      // compute the second derivative
      MATRIX Qij = hessianMatrix(q, x, y);
      MATRIX QijTranspose = Qij.transpose();
      //MATRIX QijMq = Qij ^ Mq;
      MATRIX QijMq = Qij ^ Mq;

      MATRIX IBU_U = _IBU_U.transform(QijMq);
      //IBU_IV = _IBU_IV.transform(QijTranspose);
      IBU_IBU = _IBU_IBU.transform(QijTranspose);
      //MATRIX IV_U = _IV_U.transform(QijMq);
      //IV_IV = _IV_IV.transform(QijTranspose);
      //IV_IBU = _IV_IBU.transform(QijTranspose);

      final += xPlusU ^ (IBU_U * xPlusU);
      //final -= xPlusU ^ (IBU_IV * centerDiff);
      final -= xPlusU ^ (IBU_IBU * xBar);
      //final -= centerOfMass ^ (IV_U * xPlusU);
      //final += centerOfMass ^ (IV_IV * centerDiff);
      //final += centerOfMass ^ (IV_IBU * xBar);

      if (x == y)
        final += 2 * lambda * _vertices.size();

      finalHessian(x,y) = 2 * final;
    }

  // fill out the Lagrange multiplier entries
  finalHessian(0,4) = finalHessian(4,0) = 2 * q[0];
  finalHessian(1,4) = finalHessian(4,1) = 2 * q[1];
  finalHessian(2,4) = finalHessian(4,2) = 2 * q[2];
  finalHessian(3,4) = finalHessian(4,3) = 2 * q[3];
  finalHessian(4,4) = 0.0;

  return finalHessian;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
VECTOR UNCONSTRAINED_SUBSPACE_TET_MESH::energyGradient(QUATERNION& q, Real lambda)
{
  QUATERNION finalGradient;
  VECTOR centerOfMass(3);
  VECTOR restCenterOfMass(3);
  for (int x = 0; x < 3; x++)
  {
    //centerOfMass[x] = _centerOfMass[x];
    //centerOfMass[x] = _translation[x];
    //restCenterOfMass[x] = _restCenterOfMass[x];
  }
  //VECTOR centerDiff = centerOfMass + restCenterOfMass;
  VECTOR centerDiff = centerOfMass - restCenterOfMass;

  // generate Mq, the q full quaternion multiply
  MATRIX tmp = q.toExplicitMatrix();
  MATRIX Mq(4,3);
  for (int x = 0; x < 4; x++)
    for (int y = 0; y < 3; y++)
      Mq(x,y) = tmp(x,y);

  // create a padded rotation
  MATRIX R4x4(4,4);
  MATRIX3 rotation = _rotationQuaternion.toExplicitMatrix3x3();
  for (int x = 0; x < 3; x++)
    for (int y = 0; y < 3; y++)
      R4x4(x,y) = rotation(x,y);

  VECTOR& xBar = _projectedRestPose;
  VECTOR& uBar = _q;
  VECTOR xPlusU = xBar + uBar;

  // for each component of the quaternion
  for (int x = 0; x < 4; x++)
  {
    MATRIX Qi = gradientMatrix(q, x);
    Qi = Qi * R4x4;
    MATRIX QiTranspose = Qi.transpose();
    MATRIX QiMq = QiTranspose * Mq;

    MATRIX IBU_U = _IBU_U.transform(QiMq);
    //MATRIX IBU_IV = _IBU_IV.transform(QiTranspose);
    MATRIX IBU_IBU = _IBU_IBU.transform(QiTranspose);
    //MATRIX IV_U = _IV_U.transform(QiMq);
    //MATRIX IV_IV = _IV_IV.transform(QiTranspose);
    //MATRIX IV_IBU = _IV_IBU.transform(QiTranspose);

    // if the centers of mass are zero, can skip a few of these
    Real final = 0.0;
    final += xPlusU ^ (IBU_U * xPlusU);
    //final -= xPlusU ^ (IBU_IV * centerDiff);
    final -= xPlusU ^ (IBU_IBU * xBar);
    //final -= centerOfMass ^ (IV_U * xPlusU);
    //final += centerOfMass ^ (IV_IV * centerDiff);
    //final += centerOfMass ^ (IV_IBU * xBar);

    finalGradient[x] = final;
    
    finalGradient[x] += _vertices.size() * 2 * lambda * q[x];
  }

  // take into account the 2 in front of the integral
  finalGradient *= 2.0;

  VECTOR finalVector(5);
  finalVector[0] = finalGradient[0];
  finalVector[1] = finalGradient[1];
  finalVector[2] = finalGradient[2];
  finalVector[3] = finalGradient[3];
  finalVector[4] = (q ^ q) - 1;
  return finalVector;
}

//////////////////////////////////////////////////////////////////////////////
// First derivatives from Eqns. 5 and 7 from the Georgii and Westermann paper,
// but phrased as a matrix multiply
//////////////////////////////////////////////////////////////////////////////
MATRIX UNCONSTRAINED_SUBSPACE_TET_MESH::gradientMatrix(QUATERNION& q, int which)
{
  MATRIX dq(4,4);
  switch (which)
  {
    case 0:
      dq(0,0) =  q[0]; dq(0,1) =  q[1]; dq(0,2) =  q[2];
      dq(1,0) =  q[1]; dq(1,1) = -q[0]; dq(1,2) = -q[3];
      dq(2,0) =  q[2]; dq(2,1) =  q[3]; dq(2,2) = -q[0];
      dq(3,3) =  q[0];
      break;
    case 1:
      dq(0,0) = -q[1]; dq(0,1) =  q[0]; dq(0,2) =  q[3];
      dq(1,0) =  q[0]; dq(1,1) =  q[1]; dq(1,2) =  q[2];
      dq(2,0) = -q[3]; dq(2,1) =  q[2]; dq(2,2) = -q[1];
      dq(3,3) =  q[1];
      break;
    case 2:
      dq(0,0) = -q[2]; dq(0,1) = -q[3]; dq(0,2) =  q[0];
      dq(1,0) =  q[3]; dq(1,1) = -q[2]; dq(1,2) =  q[1];
      dq(2,0) =  q[0]; dq(2,1) =  q[1]; dq(2,2) =  q[2];
      dq(3,3) =  q[2];
      break;
    case 3:
      dq(0,0) =  q[3]; dq(0,1) = -q[2]; dq(0,2) =  q[1];
      dq(1,0) =  q[2]; dq(1,1) =  q[3]; dq(1,2) = -q[0];
      dq(2,0) = -q[1]; dq(2,1) =  q[0]; dq(2,2) =  q[3];
      dq(3,3) =  q[3];
      break;
    default:
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " gradientMatrix incorrectly called! " << endl;
      break;
  }
  dq *= 2;
  return dq;
}

//////////////////////////////////////////////////////////////////////////////
// Second derivatives from Eqn. 7 from the Georgii and Westermann paper,
// but phrased as a matrix multiply
//////////////////////////////////////////////////////////////////////////////
MATRIX UNCONSTRAINED_SUBSPACE_TET_MESH::hessianMatrix(QUATERNION& q, int i, int j)
{
  MATRIX final(4,4);

  //dxdx =
  if (i == 0 && j == 0)
  {
    final(0,0) = 1;
    final(1,1) = -1;
    final(2,2) = -1;
    final(3,3) = 1;
  }

  //dydy = 
  if (i == 1 && j == 1)
  {
    final(0,0) = -1;
    final(1,1) = 1;
    final(2,2) = -1;
    final(3,3) = 1;
  }

  //dzdz = 
  if (i == 2 && j == 2)
  {
    final(0,0) = -1;
    final(1,1) = -1;
    final(2,2) = 1;
    final(3,3) = 1;
  }

  //dwdw = 
  if (i == 3 && j == 3)
  {
    final(0,0) = 1;
    final(1,1) = 1;
    final(2,2) = 1;
    final(3,3) = 1;
  }

  //dxdy = 
  if ((i == 0 && j == 1) || (i == 1 && j == 0))
  {
    final(0,1) = 1;
    final(1,0) = 1;
  }

  //dxdz = 
  if ((i == 0 && j == 2) || (i == 2 && j == 0))
  {
    final(0,2) = 1;
    final(2,0) = 1;
  }

  //dxdw= 
  if ((i == 0 && j == 3) || (i == 3 && j == 0))
  {
    final(1,2) = -1;
    final(2,1) = 1;
  }

  //dydz= 
  if ((i == 1 && j == 2) || (i == 2 && j == 1))
  {
    final(1,2) = 1;
    final(2,1) = 1;
  }

  //dydw= 
  if ((i == 1 && j == 3) || (i == 3 && j == 1))
  {
    final(0,2) = 1;
    final(2,0) = -1;
  }

  //dzdw= 
  if ((i == 2 && j == 3) || (i == 3 && j == 2))
  {
    final(0,1) = -1;
    final(1,0) = 1;
  }

  final *= 2;
  return final;
}

//////////////////////////////////////////////////////////////////////////////
// for debugging only -- unpack, project, expand, and repack a vector
//////////////////////////////////////////////////////////////////////////////
VECTOR UNCONSTRAINED_SUBSPACE_TET_MESH::projectQuaternionVector(VECTOR& d)
{
  VECTOR repack(d.size() / 4 * 3);

  for (int x = 0; x < d.size() / 4; x++)
  {
    repack[3 * x] = d[4 * x];
    repack[3 * x + 1] = d[4 * x + 1];
    repack[3 * x + 2] = d[4 * x + 2];
  }

  VECTOR temp = _UBasis ^ repack;
  repack = _UBasis * temp;

  for (int x = 0; x < d.size() / 4; x++)
  {
    d[4 * x] = repack[3 * x];
    d[4 * x + 1] = repack[3 * x + 1];
    d[4 * x + 2] = repack[3 * x + 2];
  }

  return repack;
}

//////////////////////////////////////////////////////////////////////////////
// for debugging only -- make a padded U for quaternions
//////////////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_TET_MESH::buildQU(BLOCK_MATRIX& qU)
{
  int cols = _UBasis.cols();
  int N = _UBasis.rows() /3;
  qU.resizeAndWipe(N, 1);

  for (int x = 0; x < N; x++)
  {
    MATRIX blockU(4, cols);
    for (int y = 0; y < cols; y++)
    {
      blockU(0,y) = _UBasis(3 * x, y);
      blockU(1,y) = _UBasis(3 * x + 1, y);
      blockU(2,y) = _UBasis(3 * x + 2, y);
    }

    qU.add(blockU, x, 0);
  }
}

//////////////////////////////////////////////////////////////////////////////
// for debugging only -- make a blockes version of U
//////////////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_TET_MESH::buildBlockU(BLOCK_MATRIX& blockU)
{
  int cols = _UBasis.cols();
  int N = _UBasis.rows() /3;
  blockU.resizeAndWipe(N, 1);

  for (int x = 0; x < N; x++)
  {
    MATRIX block(3, cols);
    for (int y = 0; y < cols; y++)
    {
      block(0,y) = _UBasis(3 * x, y);
      block(1,y) = _UBasis(3 * x + 1, y);
      block(2,y) = _UBasis(3 * x + 2, y);
    }

    blockU.add(block, x, 0);
  }
}

/*
//////////////////////////////////////////////////////////////////////////////
// compute the rigid translation without actually altering the state
//////////////////////////////////////////////////////////////////////////////
VEC3F UNCONSTRAINED_SUBSPACE_TET_MESH::computeRigidTranslation()
{
  // extract the translation
  VECTOR& q = _q;

  MATRIX& U = _UBasis;
  VEC3F normalization(U(0,0), U(1,1), U(2,2));
  VEC3F qTranslation;
  qTranslation[0] = q[0] * normalization[0];
  qTranslation[1] = q[1] * normalization[1];
  qTranslation[2] = q[2] * normalization[2];

  // undo the rigid rotation
  MATRIX3 rotation = _rotation;
  qTranslation = rotation * qTranslation;
 
  // zero out q translation
  //qOldCopy[0] -= q[0]; qOldCopy[1] -= q[1]; qOldCopy[2] -= q[2];
  VEC3F qBackup;
  qBackup[0] = q[0];
  qBackup[1] = q[1];
  qBackup[2] = q[2];
  q[0] = q[1] = q[2] = 0.0;

  // get the center of mass
  VEC3F centerOfMass = computeSubspaceCenterOfMass();

  // undo the rigid rotation
  centerOfMass = rotation * centerOfMass;

  // add center of mass to total translation
  qTranslation += centerOfMass;

  // restore q to its previous state
  q[0] = qBackup[0];
  q[1] = qBackup[1];
  q[2] = qBackup[2];
}
*/

//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_TET_MESH::updateRigidTranslation()
{
  /*
  // extract the translation
  VECTOR& q = _q;

  MATRIX& U = _UBasis;
  VEC3F normalization(U(0,0), U(1,1), U(2,2));

  VEC3F centerOfMass = computeSubspaceCenterOfMass();

  q[0] -= centerOfMass[0] / normalization[0];
  q[1] -= centerOfMass[1] / normalization[1];
  q[2] -= centerOfMass[2] / normalization[2];
  //cout << " translation qs after: " << q[0] << " " << q[1] << " " << q[2] << endl;
 
  //cout << " rigid translation before: " << _translation << endl; 
  _translation += centerOfMass;
  //cout << " rigid translation after: " << _translation << endl; 
  */
  
  // deactivate everything for now
  /*
  // extract the translation
  VECTOR& q = _q;

  MATRIX& U = _UBasis;
  VEC3F normalization(U(0,0), U(1,1), U(2,2));
  VEC3F qTranslation;
  qTranslation[0] = q[0] * normalization[0];
  qTranslation[1] = q[1] * normalization[1];
  qTranslation[2] = q[2] * normalization[2];

  // undo the rigid rotation
  MATRIX3 rotation = _rotation;
 
  // add to the total translation
  _translationNewtonOld = _translation;
  _translation += qTranslation;

  // zero out q translation
  q[0] = q[1] = q[2] = 0.0;

  //exit(0);

  // get the center of mass
  VEC3F centerOfMass = computeSubspaceCenterOfMass();

  // subtract center of mass from q translation
  q[0] -= centerOfMass[0] / normalization[0];
  q[1] -= centerOfMass[1] / normalization[1];
  q[2] -= centerOfMass[2] / normalization[2];

  // undo the rigid rotation
  centerOfMass = rotation * centerOfMass;

  // add center of mass to total translation
  _translation += centerOfMass;
  */
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
bool UNCONSTRAINED_SUBSPACE_TET_MESH::centerOfMassIsZero()
{
  updateFullMesh();
 
  // is the center of mass for _vertices being used?
  VEC3F vertexSum;
  Real totalMass = 0;
  int size = _vertices.size();
  for (int x = 0; x < size; x++)
  {
    vertexSum += mass(x) * _vertices[x];
    totalMass += mass(x);
  }
  vertexSum *= 1.0 / totalMass;

  // if this is tripped, then the center of mass aren't centered at zero
  if (vertexSum * vertexSum > 1e-4)
    return false;
  return true;
}

//////////////////////////////////////////////////////////////////////
// reset state to zero
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_TET_MESH::reset()
{
  _rotationLambda = 0.0;
  _rotationQuaternion = QUATERNION(MATRIX3::I());

  /*
  // zero translation version
  _translation *= 0; 
  _translationOld *= 0;
  _centerOfMass *= 0;

  Real U00 = _UBasis(0,0);
  Real U11 = _UBasis(1,1);
  Real U22 = _UBasis(2,2);

  _q *= 0.0;
  _qOld *= 0.0;
  _q[0] = _qOld[0] = _originalCenterOfMass[0] / U00;
  _q[1] = _qOld[1] = _originalCenterOfMass[1] / U11;
  _q[2] = _qOld[2] = _originalCenterOfMass[2] / U22;

  cout << " Reset q: " << _q << endl;
  */

  // with translations version
  _translation = _originalCenterOfMass;
  _translationOld = _originalCenterOfMass;
  _centerOfMass = _originalCenterOfMass;

  _q *= 0.0;
  _qOld *= 0.0;

  updateFullMesh();
}

//////////////////////////////////////////////////////////////////////
// read basis with translation projected out
//////////////////////////////////////////////////////////////////////
bool UNCONSTRAINED_SUBSPACE_TET_MESH::readTranslationlessBasis()
{
  string basisFile = _filename + string(".unconstrained.basis");
  FILE* file = fopen(basisFile.c_str(), "rb");

  if (file == NULL) return false;

  _UBasis.read(file);
  _q.read(file);
  _qOld.read(file);
  _reducedM.read(file);
  _surfaceU.read(file);
  _Kreduced.read(file);

  fclose(file);
  return true;
}

//////////////////////////////////////////////////////////////////////
// write basis with translation projected out
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_TET_MESH::writeTranslationlessBasis()
{
  string basisFile = _filename + string(".unconstrained.basis");
  FILE* file = fopen(basisFile.c_str(), "wb");

  if (file == NULL)
  {
    cout << " Attempt to write cache, " << basisFile.c_str() << " failed! " << endl;
    return;
  }

  _UBasis.write(file);
  _q.write(file);
  _qOld.write(file);
  _reducedM.write(file);
  _surfaceU.write(file);
  _Kreduced.write(file);

  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// read variable mass matrix vars
//////////////////////////////////////////////////////////////////////
bool UNCONSTRAINED_SUBSPACE_TET_MESH::readVariableMassMatrix(string filename)
{
  string massFile = _filename + string(".unconstrained.inertia");
  FILE* file = fopen(massFile.c_str(), "rb");

  if (file == NULL) return false;

  _subspaceCenterOfMass.read(file);
  _centerOfMassRest.read(file);

  // cacheInertiaVars
  fread((void*)&_restXrestX, sizeof(Real), 1, file);
  fread((void*)&_restYrestY, sizeof(Real), 1, file);
  fread((void*)&_restZrestZ, sizeof(Real), 1, file);

  fread((void*)&_restXrestY, sizeof(Real), 1, file);
  fread((void*)&_restXrestZ, sizeof(Real), 1, file);
  fread((void*)&_restYrestZ, sizeof(Real), 1, file);

  _restXUX.read(file);
  _restYUY.read(file);
  _restZUZ.read(file);
  
  _restXUY.read(file);
  _restXUZ.read(file);
  _restYUZ.read(file);

  _restYUX.read(file);
  _restZUX.read(file);
  _restZUY.read(file);

  _UXTUX.read(file);
  _UYTUY.read(file);
  _UZTUZ.read(file);

  _UXTUY.read(file);
  _UXTUZ.read(file);
  _UYTUZ.read(file);

  _restYUZrestZUY.read(file);
  _restXUZrestZUX.read(file);
  _restXUYrestYUX.read(file);

  _workspaceL.read(file);

  // cacheRotationDefoVars
  _restIthetaf.read(file);
  _UijTildeUij.read(file);

  // cacheMassMatrixVars
  _SiBar.read(file);
 
  VECTOR massSummed(3);
  massSummed.read(file);
  _massSummedRestPoses = VEC3F(massSummed);

  fclose(file);
  return true;
}

//////////////////////////////////////////////////////////////////////
// write variable mass matrix vars
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_TET_MESH::writeVariableMassMatrix(string filename)
{
  string cacheName("");
  cacheName = _filename + string(".unconstrained.inertia");
  FILE* file;
  file = fopen(cacheName.c_str(), "wb");

  _subspaceCenterOfMass.write(file);
  _centerOfMassRest.write(file);

  // cacheInertiaVars
  fwrite((void*)&_restXrestX, sizeof(Real), 1, file);
  fwrite((void*)&_restYrestY, sizeof(Real), 1, file);
  fwrite((void*)&_restZrestZ, sizeof(Real), 1, file);

  fwrite((void*)&_restXrestY, sizeof(Real), 1, file);
  fwrite((void*)&_restXrestZ, sizeof(Real), 1, file);
  fwrite((void*)&_restYrestZ, sizeof(Real), 1, file);

  _restXUX.write(file);
  _restYUY.write(file);
  _restZUZ.write(file);
  
  _restXUY.write(file);
  _restXUZ.write(file);
  _restYUZ.write(file);

  _restYUX.write(file);
  _restZUX.write(file);
  _restZUY.write(file);

  _UXTUX.write(file);
  _UYTUY.write(file);
  _UZTUZ.write(file);

  _UXTUY.write(file);
  _UXTUZ.write(file);
  _UYTUZ.write(file);

  _restYUZrestZUY.write(file);
  _restXUZrestZUX.write(file);
  _restXUYrestYUX.write(file);

  _workspaceL.write(file);

  // cacheRotationDefoVars
  _restIthetaf.write(file);
  _UijTildeUij.write(file);

  // cacheMassMatrixVars
  _SiBar.write(file);
 
  VECTOR massSummed = _massSummedRestPoses.toVector(); 
  massSummed.write(file);

  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// Draw the constrained nodes as GL points
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_TET_MESH::drawConstrainedNodes()
{
  MATRIX3 R = _rotationQuaternion.toExplicitMatrix3x3();
  glBegin(GL_POINTS);
  for (unsigned int x = _unconstrainedSize; x < _vertices.size(); x++)
  {
    VEC3F vertex = _vertices[x];
    vertex = R * vertex + _translation;
#ifdef SINGLE_PRECISION
    glVertex3fv(vertex);
#else
    glVertex3dv(vertex);
#endif
  }
  glEnd();
}

//////////////////////////////////////////////////////////////////////
// recenter the mesh after modifying the mass matrix
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_TET_MESH::recenterMesh()
{
  // add the old center of mass back to everything
  _translationOld *= 0;
  _translation *= 0;
  _originalCenterOfMass *= 0;
  // make center of mass the origin for all vertices
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _vertices[x] += _centerOfMass;
  for (unsigned int x = 0; x < _restPose.size(); x++)
    _restPose[x] += _centerOfMass;

  // compute new center of mass
  computeCenterOfMass();
  computeRestCenterOfMass();

  // set center of mass as a translation
  _translationOld = _centerOfMass;
  _translation = _centerOfMass;
  _originalCenterOfMass = _centerOfMass;

  // make center of mass the origin for all vertices
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _vertices[x] -= _centerOfMass;
  for (unsigned int x = 0; x < _restPose.size(); x++)
    _restPose[x] -= _centerOfMass;

  // precompute the center of mass projection
  if (_totalKeyTets != 0)
  {
    cout << " no cache found! ..."; flush(cout);
    cacheSubspaceCenterOfMass();

    // recompute the mass matrix vars since we moved everything
    cacheInertiaVars();
    cacheRotationDefoVars();
    cacheMassMatrixVars();
    cout << " done." << endl;
  }

  reset();
}
