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
// PLANE.cpp: implementation of the PLANE class.
//
//////////////////////////////////////////////////////////////////////

#include "PLANE.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

PLANE::PLANE(VEC3F& point, VEC3F& normal) :
  _point(point), _normal(normal), _springJacobian(3,3), _dampingJacobian(3,3)
{
  _normal.normalize();

  VEC3F xAxis(1,0,0);
  VEC3F zNew = cross(xAxis, _normal);
  zNew.normalize();

  // check for the case where the normal is the x axis
  if (fabs(xAxis * _normal) - 1.0 <=  1e-4)
    zNew = VEC3F(0,0,1);

  VEC3F xNew = cross(_normal, zNew);

  _rotation(0,0) = xNew[0];
  _rotation(1,0) = xNew[1];
  _rotation(2,0) = xNew[2];
  _rotation(0,1) = _normal[0];
  _rotation(1,1) = _normal[1];
  _rotation(2,1) = _normal[2];
  _rotation(0,2) = zNew[0];
  _rotation(1,2) = zNew[1];
  _rotation(2,2) = zNew[2];

  QUATERNION quaternion(_rotation);
  quaternion.axisAngle(_axis, _angle);

  // with stiffness = 100 and damping = 0, the bunny appears stable
  //
  // try pushing up the spring const and seeing when it goes unstable, then push up
  // the damping and see if it stabilizes again
  //
  // seems to get unhappy around 400

  //_collisionStiffness = 10000;
  //_collisionStiffness = 4000;
  //_collisionStiffness = 1000;
  //_collisionStiffness = 800;
  //_collisionStiffness = 400;

  // WAS ON THIS - 2/1/11
  //_collisionStiffness = 200;

  //_collisionStiffness = 100;

  // 32^3 ARMADILLO WAS STABLE ON THIS - 2/13/11
  _collisionStiffness = 50;
  
  //_collisionStiffness = 5;
  //_collisionStiffness = 10;
  //_collisionStiffness = 2;
  //_collisionStiffness = 1;
  //_collisionDamping = 1.0;
  //_collisionDamping = 0.5;
  _collisionDamping = 0.1;
  //_collisionDamping = 0.01;
  //_collisionDamping = 0.0;
  
  //_stickiness = 0.05;
  //_stickiness = 0.03;
  //_stickiness = 0.01;
  //_stickiness = 0.001;
  _stickiness = 0.0;

  /*
  _springJacobian(0,0) = _normal[0] * _normal[0];
  _springJacobian(0,1) = _normal[0] * _normal[1];
  _springJacobian(0,2) = _normal[0] * _normal[2];
  _springJacobian(1,0) = _normal[1] * _normal[0];
  _springJacobian(1,1) = _normal[1] * _normal[1];
  _springJacobian(1,2) = _normal[1] * _normal[2];
  _springJacobian(2,0) = _normal[2] * _normal[0];
  _springJacobian(2,1) = _normal[2] * _normal[1];
  _springJacobian(2,2) = _normal[2] * _normal[2];
  */
  _springJacobian = MATRIX(MATRIX3::outer_product(_normal, _normal));
  _dampingJacobian = _springJacobian; 
  
  _springJacobian *= _collisionStiffness;
  _dampingJacobian *= _collisionDamping;
  _type.assign("PLANE");
}

PLANE::~PLANE()
{
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
bool PLANE::inside(float* point) {
  VEC3F diff;
  diff[0] = point[0] - _point[0];
  diff[1] = point[1] - _point[1];
  diff[2] = point[2] - _point[2];

  Real dot = _normal * diff;
  if (dot < 0.0) return true;

  // for resting contact, make points stick to the plane
  if (dot < _stickiness) return true;

  return false;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
bool PLANE::inside(VEC3F& point) 
{
  VEC3F diff = point - _point;
  
  Real dot = _normal * diff;
  if (dot < 0.0) return true;

  // for resting contact, make points stick to the plane
  if (dot < _stickiness) return true;

  return false;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
float PLANE::potential() { 
  return 0.0f; 
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
float PLANE::distance(float* point) {
  VEC3F diff;
  diff[0] = point[0] - _point[0];
  diff[1] = point[1] - _point[1];
  diff[2] = point[2] - _point[2];

  return (diff * _normal);
}

void PLANE::draw()
{
  VEC3F v0(-1, 0, -1);
  VEC3F v1(-1, 0, 1);
  VEC3F v2(1, 0, 1);
  VEC3F v3(1, 0, -1);

  glPushMatrix();
#ifdef DOUBLE_PRECISION
    glTranslated(_point[0], _point[1], _point[2]);
#else
    glTranslatef(_point[0], _point[1], _point[2]);
#endif
    glRotatef(_angle, _axis[0], _axis[1], _axis[2]);
    glBegin(GL_TRIANGLES);
      glNormal3f(0,1,0);
#ifdef DOUBLE_PRECISION
      glVertex3dv(v0);
      glVertex3dv(v1);
      glVertex3dv(v2);

      glVertex3dv(v0);
      glVertex3dv(v2);
      glVertex3dv(v3);
#else
      glVertex3fv(v0);
      glVertex3fv(v1);
      glVertex3fv(v2);

      glVertex3fv(v0);
      glVertex3fv(v2);
      glVertex3fv(v3);
#endif
    glEnd();
    glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
    glLineWidth(3.0f);
    glBegin(GL_LINES);
      glVertex3f(0.0f, 0.0f, 0.0f);
      glVertex3f(0.0f, 0.5f, 0.0f);
    glEnd();
  glPopMatrix();
}

//////////////////////////////////////////////////////////////////////
// compute the force on a point
//////////////////////////////////////////////////////////////////////
VEC3F PLANE::force(const VEC3F& collisionPoint, const VEC3F& collisionVelocity)
{
  // compute distance to plane
  VEC3F diff = collisionPoint - _point;

  Real distance = _normal * diff;
  //if (distance > 0) return VEC3F();

  Real velocityComponent = _normal * collisionVelocity;

  VEC3F springForce = -_collisionStiffness * (distance * _normal);
  VEC3F dampingForce = -_collisionDamping * (velocityComponent * _normal);

  return springForce + dampingForce;
}

//////////////////////////////////////////////////////////////////////
// compute the force on a point
//////////////////////////////////////////////////////////////////////
void PLANE::verifySpringJacobian(VEC3F& collisionPoint)
{
  VEC3F velocity;
  VEC3F originalForce = force(collisionPoint, velocity);
  VEC3F originalPoint = collisionPoint;

  cout << " original force: " << originalForce << endl;

  Real delta = 1e-6;
  vector<VECTOR> columns;
  for (int x = 0; x < 3; x++)
  {
    VEC3F perturb = originalPoint;
    perturb[x] += delta;

    VEC3F newForce = force(perturb, velocity);
    VEC3F diff = newForce - originalForce;
    diff *= 1.0 / delta;
    columns.push_back(diff.toVector());
  }

  MATRIX finiteDiff(columns);
  cout << " finite diff: " << finiteDiff << endl;
  cout << " computed: " << _springJacobian << endl;
}

//////////////////////////////////////////////////////////////////////
// spring jacobian -- only return if the point is actually inside
// the plane
//////////////////////////////////////////////////////////////////////
MATRIX PLANE::springJacobian(const VEC3F& collisionPoint)
{
  VEC3F diff = collisionPoint - _point;
  Real distance = _normal * diff;
  //if (distance > 0) return MATRIX(3,3);

  return _springJacobian;
}

//////////////////////////////////////////////////////////////////////
// damping jacobian -- only return if the point is actually inside
// the plane
//////////////////////////////////////////////////////////////////////
MATRIX PLANE::dampingJacobian(const VEC3F& collisionPoint, const VEC3F& collisionVelocity)
{
  VEC3F diff = collisionPoint - _point;
  Real distance = _normal * diff;
  //if (distance > 0) return MATRIX(3,3);

  return _dampingJacobian;
}

//////////////////////////////////////////////////////////////////////
// multibody defo-defo jacobian for a plane
//////////////////////////////////////////////////////////////////////
MATRIX PLANE::defoSpringJacobian(const VEC3F& collisionPoint, const MATRIX& U, const MATRIX3& rotation, const VEC3F& translation)
{
  MATRIX rotatedBasis = rotation * U;
  VECTOR normal = _normal.toVector();
  VECTOR dotNormals = rotatedBasis ^ normal;

  MATRIX outer = MATRIX::outerProduct(normal, dotNormals);
  MATRIX rotatedOuter = rotation.transpose() * outer;

  return _collisionStiffness * (U ^ rotatedOuter);
}

//////////////////////////////////////////////////////////////////////
// multibody defo-defo jacobian for a plane
//////////////////////////////////////////////////////////////////////
MATRIX PLANE::defoDampingJacobian(const VEC3F& collisionPoint, const MATRIX& U, const MATRIX3& rotation)
{
  MATRIX rotatedBasis = rotation * U;
  VECTOR normal = _normal.toVector();
  VECTOR dotNormals = rotatedBasis ^ normal;

  MATRIX outer = MATRIX::outerProduct(normal, dotNormals);
  MATRIX rotatedOuter = rotation.transpose() * outer;

  return _collisionDamping * (U ^ rotatedOuter);
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
MATRIX PLANE::angularSpringJacobian(const VEC3F& collisionPoint, const VEC3F& collisionVelocity, const VEC3F& localPoint, const MATRIX3& rotation, const TENSOR3& rotationPartial)
{
  MATRIX uBarTildeT = MATRIX::cross(localPoint).transpose();

  MATRIX rotationT = rotation.transpose();
  TENSOR3 rotationPartialT = rotationPartial.transpose();

  // partial with respect to the outer rotation
  VEC3F springForce = force(collisionPoint, collisionVelocity);
  MATRIX final = rotationPartialT.modeOneProduct(springForce.toVector());
  final = uBarTildeT * final;

  // partial with respect to the inner rotation, spring
  MATRIX partialLocal = rotationPartial.modeOneProduct(localPoint.toVector());
  VECTOR partialLocalNormal = partialLocal.transpose() * _normal.toVector();
  
  MATRIX outer = MATRIX::outerProduct(_normal.toVector(), partialLocalNormal);
  final -= _collisionStiffness * (uBarTildeT * rotationT * outer);

  // partial with respect to the inner rotation, damping
  VECTOR localVelocity = rotationT * collisionVelocity.toVector();
  MATRIX partialVelocity = rotationPartial.modeOneProduct(localVelocity);
  VECTOR partialVelocityNormal = partialVelocity.transpose() * _normal.toVector();
  
  MATRIX outerVelocity = MATRIX::outerProduct(_normal.toVector(), partialVelocityNormal);
  final -= _collisionDamping * (uBarTildeT * rotationT * outerVelocity);

  return final;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
MATRIX PLANE::translationDefoJacobian(const MATRIX& U, const MATRIX3& rotation)
{
  MATRIX RUT = rotation * U;
  RUT = RUT.transpose();

  VECTOR normal = _normal.toVector();
  VECTOR RUTn = RUT * normal;

  MATRIX final = MATRIX::outerProduct(normal, RUTn);
  final *= _collisionStiffness;

  return final;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
MATRIX PLANE::translationDefoDampingJacobian(const MATRIX& U, const MATRIX3& rotation)
{
  MATRIX RUT = rotation * U;
  RUT = RUT.transpose();

  VECTOR normal = _normal.toVector();
  VECTOR RUTn = RUT * normal;

  MATRIX final = MATRIX::outerProduct(normal, RUTn);
  final *= _collisionDamping;

  return final;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
MATRIX PLANE::translationAngularJacobian(const VEC3F& localPoint, const TENSOR3& partialRotation)
{
  MATRIX partialLocal = partialRotation.modeOneProduct(localPoint.toVector());
  VECTOR normal = _normal.toVector();

  VECTOR partialLocalNormal = partialLocal.transpose() * normal;

  MATRIX final = MATRIX::outerProduct(normal, partialLocalNormal);
  final *= _collisionStiffness;

  return final;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
MATRIX PLANE::translationAngularDampingJacobian(const VEC3F& localVelocity, const TENSOR3& partialRotation)
{
  MATRIX partialLocal = partialRotation.modeOneProduct(localVelocity.toVector());
  VECTOR normal = _normal.toVector();

  VECTOR partialLocalNormal = partialLocal.transpose() * normal;

  MATRIX final = MATRIX::outerProduct(normal, partialLocalNormal);
  final *= _collisionDamping;

  return final;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
MATRIX PLANE::angularTranslationJacobian(const VEC3F& localPoint, const MATRIX3& rotation)
{
  MATRIX uBarTildeT = MATRIX::cross(localPoint).transpose();
  MATRIX RT = rotation.transpose();

  return -1.0 * uBarTildeT * RT * _springJacobian;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
MATRIX PLANE::defoAngularJacobian(const VEC3F& collisionPoint, const VEC3F& collisionVelocity, const VEC3F& localPoint, const MATRIX3& rotation, const TENSOR3& rotationPartial, const MATRIX& U)
{
  MATRIX rotationT = rotation.transpose();
  TENSOR3 rotationPartialT = rotationPartial.transpose();

  VEC3F springForce = force(collisionPoint, collisionVelocity);
  MATRIX modeOne = rotationPartialT.modeOneProduct(springForce.toVector());
  MATRIX final = U ^ modeOne;

  VECTOR normal = _normal.toVector();
  MATRIX partialLocal = rotationPartial.modeOneProduct(localPoint.toVector());
  VECTOR partialLocalNormal = partialLocal.transpose() * normal;
  MATRIX outer = MATRIX::outerProduct(normal, partialLocalNormal);
  final -= _collisionStiffness * (U ^ rotationT * outer);

  VECTOR localVelocity = rotationT * collisionVelocity.toVector();
  MATRIX partialLocalVelocity = rotationPartial.modeOneProduct(localVelocity);
  VECTOR partialLocalVelocityNormal = partialLocalVelocity.transpose() * normal;
  MATRIX outerVelocity = MATRIX::outerProduct(normal, partialLocalVelocityNormal);
  final -= _collisionDamping * (U ^ rotationT * outerVelocity);

  return -1.0 * final;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
MATRIX PLANE::angularDefoJacobian(const VEC3F& collisionPoint, const VEC3F& collisionVelocity, const VEC3F& localPoint, const MATRIX& U, const MATRIX3& rotation)
{
  MATRIX defoJacobian = translationDefoJacobian(U, rotation);

  MATRIX uBarTilde = MATRIX::cross(localPoint);
  MATRIX RT = rotation.transpose();
  MATRIX AuBar = rotation * uBarTilde;
  MATRIX AuBarT = AuBar.transpose();

  MATRIX final = -1.0 * AuBarT * defoJacobian;

  VEC3F f = force(collisionPoint, collisionVelocity);
  TENSOR3 tildeU = TENSOR3::cross(U);
  TENSOR3 tildeUT = tildeU.transpose();

  VECTOR RTf = RT * f.toVector();
  MATRIX crossTerm = tildeUT.modeOneProduct(RTf);
  final += crossTerm;

  return final;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
MATRIX PLANE::angularDefoJacobianDamping(const VEC3F& collisionPoint, const VEC3F& collisionVelocity, const VEC3F& localPoint, const MATRIX& U, const MATRIX3& rotation)
{
  MATRIX defoJacobian = translationDefoJacobian(U, rotation);
  // give it the correct coefficient
  defoJacobian *= 1.0 / _collisionStiffness;
  defoJacobian *= _collisionDamping;

  MATRIX uBarTilde = MATRIX::cross(localPoint);
  MATRIX RT = rotation.transpose();
  MATRIX AuBar = rotation * uBarTilde;
  MATRIX AuBarT = AuBar.transpose();

  MATRIX final = -1.0 * AuBarT * defoJacobian;

  return final;
}

//////////////////////////////////////////////////////////////////////
// multibody version -- adds the necessary quantities to 3x3 block 
// matrix 'blockJacobian'
//////////////////////////////////////////////////////////////////////
void PLANE::springJacobian(const VEC3F& collisionPoint, 
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
  MATRIX translationPartialDefoDamping(3, rank);
  MATRIX defoPartialAngular(rank, 3);
  MATRIX angularPartialDefo(3, rank);
  MATRIX angularPartialAngular(3,3);
  MATRIX defoPartialDefo(rank,rank);

  Real* accelerationAlpha = _collisionAlpha;
  VEC3F& translation = _collisionTranslation;
  MATRIX3& rotation = _collisionRotation;
  TENSOR3& rotationPartial = _collisionRotationPartial;
  VECTOR localPointVector = localPoint.toVector();

  // refresh the spring jacobian
  _springJacobian = _collisionStiffness * MATRIX3::outer_product(_normal, _normal);

  // add to easiest term, translation wrt itself
  translationPartialTranslation += accelerationAlpha[4] * _springJacobian;

  //////////////////////////////////////////////////////////////////////////////
  // defo-defo
  //////////////////////////////////////////////////////////////////////////////
  MATRIX rotatedBasis = rotation * U;
  VECTOR normal = _normal.toVector();
  VECTOR dotNormals = rotatedBasis ^ normal;
  MATRIX outer = MATRIX::outerProduct(normal, dotNormals);
  MATRIX UrotatedOuter = rotatedBasis ^ outer;
  defoPartialDefo.axpy((accelerationAlpha[4] * _collisionStiffness +
                       accelerationAlpha[1] * _collisionDamping), UrotatedOuter);

  //////////////////////////////////////////////////////////////////////////////
  // angular-angular
  //////////////////////////////////////////////////////////////////////////////
  MATRIX uBarTilde = MATRIX::cross(localPoint);
  MATRIX rotationT = rotation.transpose();
  MATRIX uBarTildeTrotationT = uBarTilde ^ rotationT;
  TENSOR3 rotationPartialT = rotationPartial.transpose();

  // partial with respect to the outer rotation
  VECTOR springForce = collisionForce.toVector();
  MATRIX final = rotationPartialT.modeOneProduct(springForce);
  final = uBarTilde ^ final;

  // partial with respect to the inner rotation, spring
  MATRIX partialLocal = rotationPartial.modeOneProduct(localPointVector);
  VECTOR partialLocalNormal = partialLocal ^ normal;
  
  MATRIX angularAngularOuter = MATRIX::outerProduct(_collisionStiffness * normal, partialLocalNormal);

  // partial with respect to the inner rotation, damping
  VECTOR localVelocity = rotationT * collisionVelocity.toVector();
  MATRIX partialVelocity = rotationPartial.modeOneProduct(localVelocity);
  VECTOR partialVelocityNormal = partialVelocity ^ normal;
  
  MATRIX outerVelocity = MATRIX::outerProduct(_collisionDamping * normal, partialVelocityNormal);
  final -= (uBarTildeTrotationT * (outerVelocity + angularAngularOuter));
  angularPartialAngular += final;

  //////////////////////////////////////////////////////////////////////////////
  // trans-defo and transpose
  //////////////////////////////////////////////////////////////////////////////
  MATRIX RUT = rotation * U;
  RUT = RUT.transpose();
  VECTOR RUTn = RUT * normal;
  outer = MATRIX::outerProduct(normal, RUTn);

  translationPartialDefo += accelerationAlpha[4] * _collisionStiffness * outer;
  translationPartialDefoDamping += accelerationAlpha[1] * _collisionDamping * outer;

  //////////////////////////////////////////////////////////////////////////////
  // trans-angular and transpose
  //////////////////////////////////////////////////////////////////////////////
  translationPartialAngular += (1.0 + _collisionDamping / _collisionStiffness) * angularAngularOuter;
  angularPartialTranslation -= accelerationAlpha[4] * uBarTildeTrotationT * _springJacobian;

  //////////////////////////////////////////////////////////////////////////////
  // defo-angular
  //////////////////////////////////////////////////////////////////////////////
  MATRIX modeOne = rotationPartialT.modeOneProduct(springForce);
  MATRIX defoAngularFinal = (U ^ (rotationT * (angularAngularOuter + outerVelocity) - modeOne));
  defoAngularFinal *= -1.0;
  defoPartialAngular -= defoAngularFinal;

  //////////////////////////////////////////////////////////////////////////////
  // angular-defo
  //////////////////////////////////////////////////////////////////////////////
  // assumes outer = MATRIX::outerProduct(normal, RUTn);
  MATRIX angularDefoFinal = -_collisionStiffness * uBarTildeTrotationT * outer;

  VECTOR RTf = rotationT * springForce;
  TENSOR3 tildeUT = TENSOR3::crossTranspose(U);
  tildeUT.modeOneAxpy(RTf, 1.0, angularDefoFinal);

  angularPartialDefo.axpy(accelerationAlpha[4], angularDefoFinal);
  
  //////////////////////////////////////////////////////////////////////////////
  // store final quantities
  //////////////////////////////////////////////////////////////////////////////
  //translationPartialDefo += translationPartialDefoDamping;
  MATRIX defoPartialTranslation = translationPartialDefo.transpose();
  
  translationPartialDefo += translationPartialDefoDamping;

  systemMatrix.add(translationPartialTranslation, 0,0);
  systemMatrix.add(angularPartialAngular, 1,1);
  systemMatrix.add(defoPartialDefo, 2,2);
  systemMatrix.add(defoPartialAngular, 2,1);
  systemMatrix.add(angularPartialDefo, 1,2);
  systemMatrix.add(translationPartialDefo, 0,2);
  systemMatrix.add(defoPartialTranslation, 2,0);
  systemMatrix.add(translationPartialAngular, 0,1);
  systemMatrix.add(angularPartialTranslation, 1,0);
}
/*
{
  int rank = U.cols();
  MATRIX translationPartialTranslation(3,3);
  MATRIX translationPartialAngular(3,3);
  MATRIX angularPartialTranslation(3,3);
  MATRIX translationPartialDefo(3, rank);
  MATRIX translationPartialDefoDamping(3, rank);
  MATRIX defoPartialAngular(rank, 3);
  MATRIX angularPartialDefo(3, rank);
  MATRIX angularPartialAngular(3,3);
  MATRIX defoPartialDefo(rank,rank);

  Real* accelerationAlpha = _collisionAlpha;
  VEC3F& translation = _collisionTranslation;
  MATRIX3& rotation = _collisionRotation;
  TENSOR3& rotationPartial = _collisionRotationPartial;

  // add to easiest term, translation wrt itself
  translationPartialTranslation += accelerationAlpha[4] * _springJacobian;

  //////////////////////////////////////////////////////////////////////////////
  // defo-defo
  //////////////////////////////////////////////////////////////////////////////
  MATRIX rotatedBasis = rotation * U;
  VECTOR normal = _normal.toVector();
  VECTOR dotNormals = rotatedBasis ^ normal;
  MATRIX outer = MATRIX::outerProduct(normal, dotNormals);
  MATRIX rotatedOuter = rotation.transpose() * outer;
  MATRIX UrotatedOuter = U ^ rotatedOuter;
  defoPartialDefo += (accelerationAlpha[4] * _collisionStiffness +
                      accelerationAlpha[1] * _collisionDamping)  * UrotatedOuter;

  //////////////////////////////////////////////////////////////////////////////
  // angular-angular
  //////////////////////////////////////////////////////////////////////////////
  MATRIX uBarTildeT = MATRIX::cross(localPoint).transpose();
  MATRIX rotationT = rotation.transpose();
  TENSOR3 rotationPartialT = rotationPartial.transpose();

  // partial with respect to the outer rotation
  VEC3F springForce = force(collisionPoint, collisionVelocity);
  MATRIX final = rotationPartialT.modeOneProduct(springForce.toVector());
  final = uBarTildeT * final;

  // partial with respect to the inner rotation, spring
  MATRIX partialLocal = rotationPartial.modeOneProduct(localPoint.toVector());
  VECTOR partialLocalNormal = partialLocal.transpose() * _normal.toVector();
  
  MATRIX angularAngularOuter = MATRIX::outerProduct(normal, partialLocalNormal);
  final -= _collisionStiffness * (uBarTildeT * rotationT * angularAngularOuter);

  // partial with respect to the inner rotation, damping
  VECTOR localVelocity = rotationT * collisionVelocity.toVector();
  MATRIX partialVelocity = rotationPartial.modeOneProduct(localVelocity);
  VECTOR partialVelocityNormal = partialVelocity.transpose() * normal;
  
  MATRIX outerVelocity = MATRIX::outerProduct(normal, partialVelocityNormal);
  final -= _collisionDamping * (uBarTildeT * rotationT * outerVelocity);
  angularPartialAngular += final;

  //////////////////////////////////////////////////////////////////////////////
  // trans-defo and transpose
  //////////////////////////////////////////////////////////////////////////////
  MATRIX RUT = rotation * U;
  RUT = RUT.transpose();
  VECTOR RUTn = RUT * normal;
  outer = MATRIX::outerProduct(normal, RUTn);

  translationPartialDefo += accelerationAlpha[4] * _collisionStiffness * outer;
  translationPartialDefoDamping += accelerationAlpha[1] * _collisionDamping * outer;

  //////////////////////////////////////////////////////////////////////////////
  // trans-angular and transpose
  //////////////////////////////////////////////////////////////////////////////
  translationPartialAngular += (_collisionStiffness + _collisionDamping) * angularAngularOuter;
  angularPartialTranslation -= accelerationAlpha[4] * uBarTildeT * rotationT * _springJacobian;

  //////////////////////////////////////////////////////////////////////////////
  // defo-angular
  //////////////////////////////////////////////////////////////////////////////
  MATRIX modeOne = rotationPartialT.modeOneProduct(springForce.toVector());
  MATRIX defoAngularFinal = U ^ modeOne;
  defoAngularFinal -= _collisionStiffness * (U ^ rotationT * angularAngularOuter);
  defoAngularFinal -= _collisionDamping * (U ^ rotationT * outerVelocity);
  defoPartialAngular -= defoAngularFinal;

  //////////////////////////////////////////////////////////////////////////////
  // angular-defo
  //////////////////////////////////////////////////////////////////////////////
  // assumes outer = MATRIX::outerProduct(normal, RUTn);
  MATRIX defoJacobian = _collisionStiffness * outer;
  MATRIX AuBarT = uBarTildeT * rotationT;
  MATRIX angularDefoFinal = -1.0 * AuBarT * defoJacobian;

  TENSOR3 tildeU = TENSOR3::cross(U);
  TENSOR3 tildeUT = tildeU.transpose();

  VECTOR RTf = rotationT * springForce.toVector();
  MATRIX crossTerm = tildeUT.modeOneProduct(RTf);
  angularDefoFinal += crossTerm;

  // should this be a -=?
  angularPartialDefo += accelerationAlpha[4] * angularDefoFinal;
  
  //////////////////////////////////////////////////////////////////////////////
  // store final quantities
  //////////////////////////////////////////////////////////////////////////////
  MATRIX defoPartialTranslation = translationPartialDefo.transpose();
  translationPartialDefo += translationPartialDefoDamping;

  systemMatrix.add(translationPartialTranslation, 0,0);
  systemMatrix.add(angularPartialAngular, 1,1);
  systemMatrix.add(defoPartialDefo, 2,2);
  systemMatrix.add(defoPartialAngular, 2,1);
  systemMatrix.add(angularPartialDefo, 1,2);
  systemMatrix.add(translationPartialDefo, 0,2);
  systemMatrix.add(defoPartialTranslation, 2,0);
  systemMatrix.add(translationPartialAngular, 0,1);
  systemMatrix.add(angularPartialTranslation, 1,0);
}
*/

//////////////////////////////////////////////////////////////////////
// multibody version -- adds the necessary quantities to 3x3 block 
// matrix 'blockJacobian'
//////////////////////////////////////////////////////////////////////
void PLANE::springJacobianDebug(const VEC3F& collisionPoint, 
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
  MATRIX translationPartialDefoDamping(3, rank);
  MATRIX defoPartialAngular(rank, 3);
  MATRIX angularPartialDefo(3, rank);
  MATRIX angularPartialAngular(3,3);
  MATRIX defoPartialDefo(rank,rank);

  Real* accelerationAlpha = _collisionAlpha;
  VEC3F& translation = _collisionTranslation;
  MATRIX3& rotation = _collisionRotation;
  TENSOR3& rotationPartial = _collisionRotationPartial;
  VECTOR localPointVector = localPoint.toVector();

  // add to easiest term, translation wrt itself
  translationPartialTranslation += accelerationAlpha[4] * _springJacobian;
  timingBreakdown["Collision Trans-trans"] += transTrans.timing();

  //////////////////////////////////////////////////////////////////////////////
  // defo-defo
  //////////////////////////////////////////////////////////////////////////////
  TIMER defoDefo;
  MATRIX rotatedBasis = rotation * U;
  VECTOR normal = _normal.toVector();
  VECTOR dotNormals = rotatedBasis ^ normal;
  MATRIX outer = MATRIX::outerProduct(normal, dotNormals);
  MATRIX UrotatedOuter = rotatedBasis ^ outer;
  defoPartialDefo.axpy((accelerationAlpha[4] * _collisionStiffness +
                       accelerationAlpha[1] * _collisionDamping), UrotatedOuter);
  timingBreakdown["Plane collision Defo-defo"] += defoDefo.timing();

  //////////////////////////////////////////////////////////////////////////////
  // angular-angular
  //////////////////////////////////////////////////////////////////////////////
  TIMER angularAngular;
  //TIMER angularAngular1;
  MATRIX uBarTilde = MATRIX::cross(localPoint);
  MATRIX rotationT = rotation.transpose();
  MATRIX uBarTildeTrotationT = uBarTilde ^ rotationT;
  TENSOR3 rotationPartialT = rotationPartial.transpose();
  //timingBreakdown["Plane collision Angular-angular 1"] += angularAngular1.timing();

  // partial with respect to the outer rotation
  //TIMER angularAngular2;
  VECTOR springForce = collisionForce.toVector();
  MATRIX final = rotationPartialT.modeOneProduct(springForce);
  final = uBarTilde ^ final;
  //timingBreakdown["Plane collision Angular-angular 2"] += angularAngular2.timing();

  // partial with respect to the inner rotation, spring
  //TIMER angularAngular3;
  MATRIX partialLocal = rotationPartial.modeOneProduct(localPointVector);
  VECTOR partialLocalNormal = partialLocal ^ normal;
  //timingBreakdown["Plane collision Angular-angular 3"] += angularAngular3.timing();
  
  //TIMER angularAngular4;
  MATRIX angularAngularOuter = MATRIX::outerProduct(_collisionStiffness * normal, partialLocalNormal);
  //timingBreakdown["Plane collision Angular-angular 4"] += angularAngular4.timing();

  // partial with respect to the inner rotation, damping
  //TIMER angularAngular5;
  VECTOR localVelocity = rotationT * collisionVelocity.toVector();
  MATRIX partialVelocity = rotationPartial.modeOneProduct(localVelocity);
  VECTOR partialVelocityNormal = partialVelocity ^ normal;
  //timingBreakdown["Plane collision Angular-angular 5"] += angularAngular5.timing();
  
  //TIMER angularAngular6;
  MATRIX outerVelocity = MATRIX::outerProduct(_collisionDamping * normal, partialVelocityNormal);
  final -= (uBarTildeTrotationT * (outerVelocity + angularAngularOuter));
  angularPartialAngular += final;
  //timingBreakdown["Plane collision Angular-angular 6"] += angularAngular6.timing();
  timingBreakdown["Plane collision Angular-angular"] += angularAngular.timing();

  //////////////////////////////////////////////////////////////////////////////
  // trans-defo and transpose
  //////////////////////////////////////////////////////////////////////////////
  TIMER transDefo;
  MATRIX RUT = rotation * U;
  RUT = RUT.transpose();
  VECTOR RUTn = RUT * normal;
  outer = MATRIX::outerProduct(normal, RUTn);

  translationPartialDefo += accelerationAlpha[4] * _collisionStiffness * outer;
  translationPartialDefoDamping += accelerationAlpha[1] * _collisionDamping * outer;
  timingBreakdown["Collision Trans-defo"] += transDefo.timing();

  //////////////////////////////////////////////////////////////////////////////
  // trans-angular and transpose
  //////////////////////////////////////////////////////////////////////////////
  TIMER transAngular;
  translationPartialAngular += (1.0 + _collisionDamping / _collisionStiffness) * angularAngularOuter;
  angularPartialTranslation -= accelerationAlpha[4] * uBarTildeTrotationT * _springJacobian;
  timingBreakdown["Collision Trans-angular"] += transAngular.timing();

  //////////////////////////////////////////////////////////////////////////////
  // defo-angular
  //////////////////////////////////////////////////////////////////////////////
  TIMER defoAngular;
  MATRIX modeOne = rotationPartialT.modeOneProduct(springForce);
  MATRIX defoAngularFinal = (U ^ (rotationT * (angularAngularOuter + outerVelocity) - modeOne));
  defoAngularFinal *= -1.0;
  defoPartialAngular -= defoAngularFinal;
  timingBreakdown["Collision Defo-angular"] += defoAngular.timing();

  //////////////////////////////////////////////////////////////////////////////
  // angular-defo
  //////////////////////////////////////////////////////////////////////////////
  TIMER angularDefo;
  // assumes outer = MATRIX::outerProduct(normal, RUTn);
  MATRIX angularDefoFinal = -_collisionStiffness * uBarTildeTrotationT * outer;
  timingBreakdown["Plane collision Angular-defo"] += angularDefo.timing();

  TIMER angularDefoTensorTimer;
  VECTOR RTf = rotationT * springForce;
  TENSOR3 tildeUT = TENSOR3::crossTranspose(U);
  tildeUT.modeOneAxpy(RTf, 1.0, angularDefoFinal);

  angularPartialDefo.axpy(accelerationAlpha[4], angularDefoFinal);
  timingBreakdown["Plane collision tensor product"] += angularDefoTensorTimer.timing();
  
  //////////////////////////////////////////////////////////////////////////////
  // store final quantities
  //////////////////////////////////////////////////////////////////////////////
  TIMER finalCommit;
  //translationPartialDefo += translationPartialDefoDamping;
  MATRIX defoPartialTranslation = translationPartialDefo.transpose();
  
  translationPartialDefo += translationPartialDefoDamping;

  systemMatrix.add(translationPartialTranslation, 0,0);
  systemMatrix.add(angularPartialAngular, 1,1);
  systemMatrix.add(defoPartialDefo, 2,2);
  systemMatrix.add(defoPartialAngular, 2,1);
  systemMatrix.add(angularPartialDefo, 1,2);
  systemMatrix.add(translationPartialDefo, 0,2);
  systemMatrix.add(defoPartialTranslation, 2,0);
  systemMatrix.add(translationPartialAngular, 0,1);
  systemMatrix.add(angularPartialTranslation, 1,0);
  timingBreakdown["Collision commit"] += finalCommit.timing();
}
