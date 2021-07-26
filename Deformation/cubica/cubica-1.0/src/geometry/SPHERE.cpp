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
// SPHERE.cpp: implementation of the SPHERE class.
//
//////////////////////////////////////////////////////////////////////

#include "SPHERE.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

SPHERE::SPHERE(float radius, float potential) :
  _radius(radius), _potential(potential)
{
  _center[0] = 0.5;
  _center[1] = 0.5;
  _center[2] = 0.5;
  _type.assign("SPHERE");

  //_collisionStiffness = 100;
  //_collisionDamping = 0.1;
}

SPHERE::SPHERE(float radius, float* center) :
  _radius(radius), _potential(0.0f)
{
  _center[0] = center[0];
  _center[1] = center[1];
  _center[2] = center[2];
  _type.assign("SPHERE");

  //_collisionStiffness = 100;
  //_collisionDamping = 0.1;
}

SPHERE::SPHERE(float radius, VEC3F& center) :
  _radius(radius), _potential(0.0f)
{
  _center[0] = center[0];
  _center[1] = center[1];
  _center[2] = center[2];
  _type.assign("SPHERE");

  //_collisionStiffness = 1000;

  // WAS ON THIS 2/11/11
  //_collisionStiffness = 200;
  // 32^3 ARMADILLO WAS STABLE ON THIS - 2/13/11
  _collisionStiffness = 50;
  //_collisionStiffness = 5;
  //_collisionStiffness = 0;
  _collisionDamping = 0.1;
  //_collisionDamping = 0.5;
  //_collisionDamping = 0.0;

  _stickiness = 0.0;
}

SPHERE::SPHERE(const SPHERE& sphere) :
  _radius(sphere._radius), _potential(sphere._potential)
{
  _center[0] = sphere._center[0];
  _center[1] = sphere._center[1];
  _center[2] = sphere._center[2];
  _type.assign("SPHERE");

  //_collisionStiffness = 100;
  //_collisionDamping = 0.1;
}

SPHERE::~SPHERE()
{

}

//////////////////////////////////////////////////////////////////////
// check if inside the sphere
//////////////////////////////////////////////////////////////////////
bool SPHERE::inside(float* point)
{
  float translated[3];
  translated[0] = point[0] - _center[0]; 
  translated[1] = point[1] - _center[1];
  translated[2] = point[2] - _center[2];
  
  float magnitude = translated[0] * translated[0] + 
                    translated[1] * translated[1] + 
                    translated[2] * translated[2];

  //return magnitude <= _radius * _radius;
  return magnitude <= (_radius + _stickiness) * (_radius + _stickiness);
}

//////////////////////////////////////////////////////////////////////
// get distance value in relation to sphere
//////////////////////////////////////////////////////////////////////
float SPHERE::distance(float* point)
{
  float translated[3];
  translated[0] = point[0] - _center[0];
  translated[1] = point[1] - _center[1];
  translated[2] = point[2] - _center[2];
  
  float magnitude = translated[0] * translated[0] + 
                    translated[1] * translated[1] + 
                    translated[2] * translated[2];

  return fabs(sqrt(magnitude) - _radius);
}

//////////////////////////////////////////////////////////////////////
// see if two spheres overlap
//////////////////////////////////////////////////////////////////////
bool SPHERE::intersect(SPHERE& rightSphere)
{
  VEC3F diff = _center - rightSphere.center();
  Real distanceSq = norm2(diff);

  Real combinedRadii = _radius + rightSphere.radius();
  return distanceSq <= (combinedRadii * combinedRadii);
  /*
  Real distance = sqrt(diff * diff);
  Real combinedRadii = _radius + rightSphere.radius();
  return distance <= combinedRadii;
  */
}

//////////////////////////////////////////////////////////////////////
// compute the force at the current point and velocity
//////////////////////////////////////////////////////////////////////
VEC3F SPHERE::force(const VEC3F& collisionPoint, const VEC3F& collisionVelocity) {
  VEC3F direction = collisionPoint - _center;
  Real dot = direction * direction;

  VEC3F normal = direction;
  normal = normal / sqrt(dot);
  Real velocityDot = normal * collisionVelocity;

  VEC3F springForce = _collisionStiffness * ((normal * _radius) - direction);
  VEC3F dampingForce = -_collisionDamping * normal * velocityDot;

  return springForce + dampingForce;
}

//////////////////////////////////////////////////////////////////////
// compute the force jacobian at the current point
//////////////////////////////////////////////////////////////////////
MATRIX SPHERE::springJacobian(const VEC3F& collisionPoint) {
  VEC3F direction = collisionPoint - _center;
  Real dot = direction * direction;
  Real sqrtDot = sqrt(dot);
  //Real invSqrtDot = 1.0 / sqrtDot;
  MATRIX final(3,3);

  Real same = _radius / sqrtDot - 1.0;
  final(0,0) = final(1,1) = final(2,2) = same;

  VEC3F scaledDirection = -pow(dot, -1.5) * _radius * direction;

  //VEC3F column0 = -pow(dot, -1.5) * direction[0] * _radius * direction;
  VEC3F column0 = direction[0] * scaledDirection;
  final(0,0) += column0[0];
  final(1,0) += column0[1];
  final(2,0) += column0[2];

  //VEC3F column1 = -pow(dot, -1.5) * direction[1] * _radius * direction;
  VEC3F column1 = direction[1] * scaledDirection;
  final(0,1) += column1[0];
  final(1,1) += column1[1];
  final(2,1) += column1[2];
  
  //VEC3F column2 = -pow(dot, -1.5) * direction[2] * _radius * direction;
  VEC3F column2 = direction[2] * scaledDirection;
  final(0,2) += column2[0];
  final(1,2) += column2[1];
  final(2,2) += column2[2];

  return -_collisionStiffness * final;
}

//////////////////////////////////////////////////////////////////////
// compute the damping jacobian at the current point
//////////////////////////////////////////////////////////////////////
MATRIX SPHERE::dampingJacobian(const VEC3F& collisionPoint, const VEC3F& collisionVelocity) { 
  MATRIX final(3,3);
  VEC3F direction = collisionPoint - _center;
  Real dot = direction * direction;
  Real invSqrtDot = 1.0 / sqrt(dot);
  VEC3F normal = direction * invSqrtDot;

  VEC3F column0 = direction[0] * invSqrtDot * normal;
  final(0,0) += column0[0];
  final(1,0) += column0[1];
  final(2,0) += column0[2];
  VEC3F column1 = direction[1] * invSqrtDot * normal;
  final(0,1) += column1[0];
  final(1,1) += column1[1];
  final(2,1) += column1[2];
  VEC3F column2 = direction[2] * invSqrtDot * normal;
  final(0,2) += column2[0];
  final(1,2) += column2[1];
  final(2,2) += column2[2];
  /*
  final(0,0) = direction[0] * direction[0];
  final(1,1) = direction[1] * direction[1];
  final(2,2) = direction[2] * direction[2];

  final(0,1) = final(1,0) = direction[1] * direction[0];
  final(0,2) = final(2,0) = direction[2] * direction[0];
  final(1,2) = final(2,1) = direction[1] * direction[2];
  */

  //Real velocityDot = normal * collisionVelocity;

  return -_collisionDamping * final;
}

//////////////////////////////////////////////////////////////////////
// draw a sphere
//////////////////////////////////////////////////////////////////////
void SPHERE::draw()
{
  glPushMatrix();
    glTranslatef(_center[0], _center[1], _center[2]);
    glutSolidSphere(_radius, 100, 100);
  glPopMatrix();
}

//////////////////////////////////////////////////////////////////////
// multibody version -- adds the necessary quantities to 3x3 block 
// matrix 'blockJacobian'
//////////////////////////////////////////////////////////////////////
void SPHERE::springJacobian(const VEC3F& collisionPoint, 
                            const VEC3F& collisionVelocity,
                            const VEC3F& collisionForce,
                            const VEC3F& localPoint,
                            const MATRIX& U, 
                            BLOCK_MATRIX& systemMatrix)
{
  int rank = U.cols();
  MATRIX translationPartialTranslation(3,3);
  MATRIX translationPartialAngular(3,3);
  MATRIX angularPartialTranslation(3,3);
  MATRIX translationPartialDefo(3, rank);
  MATRIX defoPartialAngular(rank, 3);
  MATRIX angularPartialDefo(3, rank);
  MATRIX angularPartialAngular(3,3);
  MATRIX defoPartialDefo(rank,rank);

  Real* accelerationAlpha = _collisionAlpha;
  VEC3F& translation = _collisionTranslation;
  MATRIX3& rotation = _collisionRotation;
  TENSOR3& rotationPartial = _collisionRotationPartial;

  // trans-trans
  translationPartialTranslation += accelerationAlpha[4] * springJacobian(collisionPoint);
  // no damping term -- doesn't appear in the velocity equation
  //translationPartialTranslation += accelerationAlpha[1] * dampingJacobian(collisionPoint, collisionVelocity);

  // trans-defo
  VECTOR direction = (collisionPoint - _center).toVector();
  Real dd = direction * direction;
  Real ddHalf = pow(dd, -0.5);
  Real ddThreeHalves = pow(dd, -1.5);
  MATRIX RU = rotation * U;
  MATRIX RUT = RU.transpose();

  translationPartialDefo = (RU * ddHalf) * _radius - RU;
  VECTOR ddPartialDefo = 2.0 * (RUT * direction);
  ddPartialDefo *= -0.5 * ddThreeHalves;
  MATRIX outerDefo = MATRIX::outerProduct(direction, ddPartialDefo);
  translationPartialDefo += _radius * outerDefo;
  translationPartialDefo *= -_collisionStiffness;
  translationPartialDefo *= accelerationAlpha[4]; 

  // defo-defo spring
  defoPartialDefo = RUT * translationPartialDefo;

  // angular-defo
  MATRIX RT = rotation.transpose();
  MATRIX uBarTildeT = MATRIX::cross(localPoint).transpose();
  MATRIX AuBarT = uBarTildeT * RT;
  angularPartialDefo = -1.0 * AuBarT * translationPartialDefo;

  VECTOR springForce = this->force(collisionPoint, collisionVelocity).toVector();
  TENSOR3 tildeU = TENSOR3::cross(U);
  TENSOR3 tildeUT = tildeU.transpose();

  VECTOR RTf = RT * springForce;
  MATRIX crossTerm = tildeUT.modeOneProduct(RTf);
  angularPartialDefo += accelerationAlpha[4] * crossTerm;

  // trans-angular
  MATRIX partialAngular = rotationPartial.modeOneProduct(localPoint.toVector());
  VECTOR ddPartialAngular = 2.0 * (partialAngular ^ direction);
  VECTOR ddPartialHalfAngular = -0.5 * ddThreeHalves * ddPartialAngular;

  MATRIX outerAngular = MATRIX::outerProduct(direction, ddPartialHalfAngular);

  translationPartialAngular = _radius * (ddHalf * partialAngular + outerAngular) - partialAngular;
  translationPartialAngular *= -_collisionStiffness;

  // angular-angular
  angularPartialAngular = -1.0 * AuBarT * translationPartialAngular;
  TENSOR3 rotationPartialT = rotationPartial.transpose();
  MATRIX rotationPartialTForce = rotationPartialT.modeOneProduct(springForce);
  angularPartialAngular += uBarTildeT * rotationPartialTForce;

  // angular-defo
  defoPartialAngular -= U ^ rotationPartialTForce;
  defoPartialAngular += RUT * translationPartialAngular;

  // angular-trans
  angularPartialTranslation = -1.0 * AuBarT * translationPartialTranslation;

  // trans-defo damping
  VECTOR velocity = collisionVelocity.toVector();
  VECTOR normal = ddHalf * direction;
  MATRIX normalPartialQ = RU * ddHalf + outerDefo;
  Real nDotV = normal * velocity;
 
  MATRIX translationPartialDefoDamping = accelerationAlpha[4] * (normalPartialQ * nDotV + 
                                                                 MATRIX::outerProduct(normal, normalPartialQ ^ velocity));
  translationPartialDefoDamping += accelerationAlpha[1] * MATRIX::outerProduct(normal, RU ^ normal);
  translationPartialDefoDamping *= _collisionDamping;

  // angular-defo damping
  MATRIX angularPartialDefoDamping = -1.0 * AuBarT * translationPartialDefoDamping;

  // defo-defo damping
  MATRIX defoPartialDefoDamping = RUT * translationPartialDefoDamping;

  // trans-angular damping
  VECTOR localVelocity = RT * velocity;
  MATRIX normalPartialAngular = partialAngular * ddHalf + outerAngular;
  MATRIX velocityPartialAngular = rotationPartial.modeOneProduct(localVelocity);

  MATRIX translationPartialAngularDamping = normalPartialAngular * nDotV + 
                                            MATRIX::outerProduct(normal, normalPartialAngular ^ velocity);
  translationPartialAngularDamping += MATRIX::outerProduct(normal, velocityPartialAngular ^ normal);
  translationPartialAngularDamping *= _collisionDamping;

  // angular-angular damping
  MATRIX angularPartialAngularDamping = -1.0 * AuBarT * translationPartialAngularDamping;

  // defo-angular damping
  MATRIX defoPartialAngularDamping = RUT * translationPartialAngularDamping;

  // translation-translation damping
  MATRIX normalPartialTrans(3,3);
  normalPartialTrans.setToIdentity();
  normalPartialTrans *= ddHalf;
  normalPartialTrans -= MATRIX::outerProduct(direction, ddThreeHalves * direction);

  MATRIX translationPartialTranslationDamping = normalPartialTrans * nDotV +
                                                MATRIX::outerProduct(normal, normalPartialTrans ^ velocity);
  translationPartialTranslationDamping *= accelerationAlpha[4] * _collisionDamping;

  // angular-translation damping
  MATRIX angularPartialTranslationDamping = -1.0 * AuBarT * translationPartialTranslationDamping;

  // defo-translation damping
  MATRIX defoPartialTranslationDamping = RUT * translationPartialTranslationDamping;

  MATRIX defoPartialTranslation = translationPartialDefo.transpose();

  systemMatrix.add(translationPartialTranslation, 0,0);
  systemMatrix.add(translationPartialTranslationDamping, 0,0);

  systemMatrix.add(angularPartialAngular, 1,1);
  systemMatrix.add(angularPartialAngularDamping, 1,1);

  systemMatrix.add(defoPartialDefo, 2,2);
  systemMatrix.add(defoPartialDefoDamping, 2,2);

  systemMatrix.add(defoPartialAngular, 2,1);
  systemMatrix.add(defoPartialAngularDamping, 2,1);

  systemMatrix.add(angularPartialDefo, 1,2);
  systemMatrix.add(angularPartialDefoDamping, 1,2);

  systemMatrix.add(translationPartialDefo, 0,2);
  systemMatrix.add(translationPartialDefoDamping, 0,2);

  systemMatrix.add(defoPartialTranslation, 2,0);
  systemMatrix.add(defoPartialTranslationDamping, 2,0);

  systemMatrix.add(translationPartialAngular, 0,1);
  systemMatrix.add(translationPartialAngularDamping, 0,1);

  systemMatrix.add(angularPartialTranslation, 1,0);
  systemMatrix.add(angularPartialTranslationDamping, 1,0);
}

//////////////////////////////////////////////////////////////////////
// multibody version -- adds the necessary quantities to 3x3 block 
// matrix 'blockJacobian'
//////////////////////////////////////////////////////////////////////
void SPHERE::springJacobianDebug(const VEC3F& collisionPoint, 
                           const VEC3F& collisionVelocity,
                           const VEC3F& collisionForce,
                           const VEC3F& localPoint,
                           const MATRIX& U, 
                           BLOCK_MATRIX& systemMatrix,
                           map<string, double>& timingBreakdown)
{
  TIMER transTrans;
  int rank = U.cols();
  MATRIX translationPartialTranslation(3,3);
  MATRIX translationPartialAngular(3,3);
  MATRIX angularPartialTranslation(3,3);
  MATRIX translationPartialDefo(3, rank);
  MATRIX defoPartialAngular(rank, 3);
  MATRIX angularPartialDefo(3, rank);
  MATRIX angularPartialAngular(3,3);
  MATRIX defoPartialDefo(rank,rank);

  Real* accelerationAlpha = _collisionAlpha;
  VEC3F& translation = _collisionTranslation;
  MATRIX3& rotation = _collisionRotation;
  TENSOR3& rotationPartial = _collisionRotationPartial;

  // trans-trans
  translationPartialTranslation += accelerationAlpha[4] * springJacobian(collisionPoint);
  // no damping term -- doesn't appear in the velocity equation
  //translationPartialTranslation += accelerationAlpha[1] * dampingJacobian(collisionPoint, collisionVelocity);
  timingBreakdown["Collision Trans-trans"] += transTrans.timing();

  // trans-defo
  TIMER transDefo;
  VECTOR direction = (collisionPoint - _center).toVector();
  Real dd = direction * direction;
  Real ddHalf = pow(dd, -0.5);
  Real ddThreeHalves = pow(dd, -1.5);
  MATRIX RU = rotation * U;
  MATRIX RUT = RU.transpose();

  translationPartialDefo = (RU * ddHalf) * _radius - RU;
  VECTOR ddPartialDefo = 2.0 * (RUT * direction);
  ddPartialDefo *= -0.5 * ddThreeHalves;
  MATRIX outerDefo = MATRIX::outerProduct(direction, ddPartialDefo);
  translationPartialDefo += _radius * outerDefo;
  translationPartialDefo *= -_collisionStiffness;
  translationPartialDefo *= accelerationAlpha[4]; 
  timingBreakdown["Collision Trans-defo"] += transDefo.timing();

  // defo-defo spring
  TIMER defoDefo;
  defoPartialDefo = RUT * translationPartialDefo;
  timingBreakdown["Sphere collision Defo-defo"] += defoDefo.timing();

  // angular-defo
  TIMER angularDefo;
  MATRIX RT = rotation.transpose();
  MATRIX uBarTildeT = MATRIX::cross(localPoint).transpose();
  MATRIX AuBarT = uBarTildeT * RT;
  angularPartialDefo = -1.0 * AuBarT * translationPartialDefo;

  //VECTOR springForce = this->force(collisionPoint, collisionVelocity).toVector();
  VECTOR springForce = collisionForce.toVector();
  TENSOR3 tildeU = TENSOR3::cross(U);
  TENSOR3 tildeUT = tildeU.transpose();

  VECTOR RTf = RT * springForce;
  //MATRIX crossTerm = tildeUT.modeOneProduct(RTf);
  //angularPartialDefo += accelerationAlpha[4] * crossTerm;
  tildeUT.modeOneAxpy(RTf, accelerationAlpha[4], angularPartialDefo);

  timingBreakdown["Sphere collision Angular-defo"] += angularDefo.timing();

  // trans-angular
  TIMER transAngular;
  MATRIX partialAngular = rotationPartial.modeOneProduct(localPoint.toVector());
  VECTOR ddPartialAngular = 2.0 * (partialAngular ^ direction);
  VECTOR ddPartialHalfAngular = -0.5 * ddThreeHalves * ddPartialAngular;

  MATRIX outerAngular = MATRIX::outerProduct(direction, ddPartialHalfAngular);

  translationPartialAngular = _radius * (ddHalf * partialAngular + outerAngular) - partialAngular;
  translationPartialAngular *= -_collisionStiffness;
  timingBreakdown["Collision Trans-angular"] += transAngular.timing();

  // angular-angular
  TIMER angularAngular;
  angularPartialAngular = -1.0 * AuBarT * translationPartialAngular;
  TENSOR3 rotationPartialT = rotationPartial.transpose();
  MATRIX rotationPartialTForce = rotationPartialT.modeOneProduct(springForce);
  angularPartialAngular += uBarTildeT * rotationPartialTForce;
  timingBreakdown["Sphere collision Angular-angular"] += angularAngular.timing();

  // defo-angualr
  TIMER defoAngular;
  defoPartialAngular -= U ^ rotationPartialTForce;
  defoPartialAngular += RUT * translationPartialAngular;
  timingBreakdown["Collision Defo-angular"] += defoAngular.timing();

  // angular-trans
  TIMER angularTrans;
  angularPartialTranslation = -1.0 * AuBarT * translationPartialTranslation;
  timingBreakdown["Collision angular-trans"] += angularTrans.timing();

  // trans-defo damping
  TIMER dampingTimer;
  VECTOR velocity = collisionVelocity.toVector();
  VECTOR normal = ddHalf * direction;
  MATRIX normalPartialQ = RU * ddHalf + outerDefo;
  Real nDotV = normal * velocity;
 
  MATRIX translationPartialDefoDamping = accelerationAlpha[4] * (normalPartialQ * nDotV + 
                                                                 MATRIX::outerProduct(normal, normalPartialQ ^ velocity));
  translationPartialDefoDamping += accelerationAlpha[1] * MATRIX::outerProduct(normal, RU ^ normal);
  translationPartialDefoDamping *= _collisionDamping;

  // angular-defo damping
  MATRIX angularPartialDefoDamping = -1.0 * AuBarT * translationPartialDefoDamping;

  // defo-defo damping
  MATRIX defoPartialDefoDamping = RUT * translationPartialDefoDamping;

  // trans-angular damping
  VECTOR localVelocity = RT * velocity;
  MATRIX normalPartialAngular = partialAngular * ddHalf + outerAngular;
  MATRIX velocityPartialAngular = rotationPartial.modeOneProduct(localVelocity);

  MATRIX translationPartialAngularDamping = normalPartialAngular * nDotV + 
                                            MATRIX::outerProduct(normal, normalPartialAngular ^ velocity);
  translationPartialAngularDamping += MATRIX::outerProduct(normal, velocityPartialAngular ^ normal);
  translationPartialAngularDamping *= _collisionDamping;

  // angular-angular damping
  MATRIX angularPartialAngularDamping = -1.0 * AuBarT * translationPartialAngularDamping;

  // defo-angular damping
  MATRIX defoPartialAngularDamping = RUT * translationPartialAngularDamping;

  // translation-translation damping
  MATRIX normalPartialTrans(3,3);
  normalPartialTrans.setToIdentity();
  normalPartialTrans *= ddHalf;
  normalPartialTrans -= MATRIX::outerProduct(direction, ddThreeHalves * direction);

  MATRIX translationPartialTranslationDamping = normalPartialTrans * nDotV +
                                                MATRIX::outerProduct(normal, normalPartialTrans ^ velocity);
  translationPartialTranslationDamping *= accelerationAlpha[4] * _collisionDamping;

  // angular-translation damping
  MATRIX angularPartialTranslationDamping = -1.0 * AuBarT * translationPartialTranslationDamping;

  // defo-translation damping
  MATRIX defoPartialTranslationDamping = RUT * translationPartialTranslationDamping;

  MATRIX defoPartialTranslation = translationPartialDefo.transpose();
  timingBreakdown["Collision damping"] += dampingTimer.timing();

  TIMER finalCommit;
  systemMatrix.add(translationPartialTranslation, 0,0);
  systemMatrix.add(translationPartialTranslationDamping, 0,0);

  systemMatrix.add(angularPartialAngular, 1,1);
  systemMatrix.add(angularPartialAngularDamping, 1,1);

  systemMatrix.add(defoPartialDefo, 2,2);
  systemMatrix.add(defoPartialDefoDamping, 2,2);

  systemMatrix.add(defoPartialAngular, 2,1);
  systemMatrix.add(defoPartialAngularDamping, 2,1);

  systemMatrix.add(angularPartialDefo, 1,2);
  systemMatrix.add(angularPartialDefoDamping, 1,2);

  systemMatrix.add(translationPartialDefo, 0,2);
  systemMatrix.add(translationPartialDefoDamping, 0,2);

  systemMatrix.add(defoPartialTranslation, 2,0);
  systemMatrix.add(defoPartialTranslationDamping, 2,0);

  systemMatrix.add(translationPartialAngular, 0,1);
  systemMatrix.add(translationPartialAngularDamping, 0,1);

  systemMatrix.add(angularPartialTranslation, 1,0);
  systemMatrix.add(angularPartialTranslationDamping, 1,0);
  timingBreakdown["Collision commit"] += finalCommit.timing();
}

//////////////////////////////////////////////////////////////////////
// return the bounding box dimensions
//////////////////////////////////////////////////////////////////////
void SPHERE::boundingBox(VEC3F& mins, VEC3F& maxs)
{
  for (int x = 0; x < 3; x++)
  {
    mins[x] = _center[x] - _radius;
    maxs[x] = _center[x] + _radius;
  }
}
