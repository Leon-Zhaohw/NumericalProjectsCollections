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
// UNCONSTRAINED_TET_MESH.h: interface for the UNCONSTRAINED_TET_MESH class.
//
//////////////////////////////////////////////////////////////////////

#include "UNCONSTRAINED_TET_MESH.h"

//////////////////////////////////////////////////////////////////////
// Constructor for tetrahedra
//////////////////////////////////////////////////////////////////////
UNCONSTRAINED_TET_MESH::UNCONSTRAINED_TET_MESH(const char* filename, MATERIAL** materials, int totalMaterials, bool simulate) :
  TET_MESH(filename, materials, totalMaterials, simulate),
  _georgiiIterations(-1)
{
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " Shouldn't be calling this !!! " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  exit(0);

  // set rigid components to zero
  _translation.clear();
  _translationOld.clear();

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

  // make center of mass the origin for all vertices
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _vertices[x] -= _centerOfMass;
  for (unsigned int x = 0; x < _restPose.size(); x++)
    _restPose[x] -= _centerOfMass;
}

//////////////////////////////////////////////////////////////////////
// Draw all the tets, including internal faces, taking into account
// rigid components
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_TET_MESH::drawAllTets()
{
#if USING_GLVU
  glPushMatrix();

  // get the modelview matrix
  glMatrixMode(GL_MODELVIEW);
  float modelMatrix[16];
  glGetFloatv(GL_MODELVIEW_MATRIX, modelMatrix);

  // force the rotation
  MATRIX3 transpose = _rotationQuaternion.toExplicitMatrix3x3().transpose();
  for (int y = 0; y < 3; y++)
    for (int x = 0; x < 3; x++)
      modelMatrix[x + y * 4] = transpose(x,y);

  // add the translation
  modelMatrix[0 + 12] += _translation[0];
  modelMatrix[1 + 12] += _translation[1];
  modelMatrix[2 + 12] += _translation[2];

  // write the matrix back
  glLoadMatrixf(modelMatrix);

  // draw the polys  
  for (unsigned int x = 0; x < _tets.size(); x++)
    _tets[x].drawTriangles();
  glPopMatrix();
#endif
}

//////////////////////////////////////////////////////////////////////
// Draw only the surface faces, taking into account rigid components
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_TET_MESH::drawSurfaceFaces()
{
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
  glEnd();
  glPopMatrix();
#endif
}

//////////////////////////////////////////////////////////////////////
// Draw only the surface faces, not taking into account rigid components
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_TET_MESH::drawFacesWithoutRigidTransform()
{
#if USING_GLVU
  glPushMatrix();
    TET_MESH::drawSurfaceFaces();
  glPopMatrix();
#endif
}

//////////////////////////////////////////////////////////////////////
// Given a position, return the closest surface node
//////////////////////////////////////////////////////////////////////
VEC3F* UNCONSTRAINED_TET_MESH::closestSurfaceNode(VEC3F point)
{
  // make sure a surface list was built 
  if (_surfaceVertices.size() == 0) return NULL;

  // tentatively set it to the first in the list
  VEC3F* closest;
  closest = _surfaceVertices[0];
  VEC3F transformed = _rotationQuaternion.toExplicitMatrix3x3() * (*_surfaceVertices[0]) + _translation;
  Real minDistance = norm2(point - transformed);

  // loop through the rest of the vertices
  for (unsigned int x = 1; x < _surfaceVertices.size(); x++)
  {
    // check if this one is closer
    VEC3F transformed = _rotationQuaternion.toExplicitMatrix3x3() * (*_surfaceVertices[x]) + _translation;

    Real distance = norm2(point - transformed);
    if (distance < minDistance)
    {
      minDistance = distance;
      closest = _surfaceVertices[x];
    }
  }
  return closest;
}

//////////////////////////////////////////////////////////////////////
// Draw only the surface faces of the rest pose
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_TET_MESH::drawRestPose()
{
  VECTOR backup = _x;
  _x *= 0;
  updateFullMesh();

#if USING_GLVU
  glPushMatrix();
  glTranslatef(_translation[0], _translation[1], _translation[2]);
  for (unsigned int x = 0; x < _surfaceFaces.size(); x++)
  {
    TRIANGLE triangle = _tets[_surfaceFaces[x].first].face(_surfaceFaces[x].second);
    triangle.draw();
  }
  glPopMatrix();
#endif
  _x = backup;
  updateFullMesh();
}

//////////////////////////////////////////////////////////////////////
// Compute the full rank rotation using shape matching
//////////////////////////////////////////////////////////////////////
MATRIX3 UNCONSTRAINED_TET_MESH::computeShapeMatchingRotation()
{
  int size = _restPose.size();
  TET_MESH::updateFullMesh();
 
  // is the center of mass for _vertices being used?
  VEC3F vertexSum;
  for (int x = 0; x < size; x++)
    vertexSum += _vertices[x];
  vertexSum *= 1.0 / size;

  // if this is tripped, then the center of mass aren't centered at zero
  assert(vertexSum * vertexSum < 1e-4);

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

//////////////////////////////////////////////////////////////////////////////
// Support functions for computeGeorgiiRotation
// Eqn. 4 from the Georgii and Westermann paper
//////////////////////////////////////////////////////////////////////////////
QUATERNION UNCONSTRAINED_TET_MESH::rotationDistance(QUATERNION& q, VEC3F& rest, VEC3F& displacement)
{
  QUATERNION rhs(rest - _restCenterOfMass);
  QUATERNION lhs(rest + displacement - _centerOfMass);
  lhs = q * lhs * q.conjugate();

  return lhs - rhs;
}

//////////////////////////////////////////////////////////////////////////////
// Support functions for computeGeorgiiRotation
// Second derivatives from Eqn. 7 from the Georgii and Westermann paper,
// but phrased as a matrix multiply
//////////////////////////////////////////////////////////////////////////////
MATRIX UNCONSTRAINED_TET_MESH::hessianMatrix(QUATERNION& q, int i, int j)
{
  MATRIX final(4,4);

  //dxdx =
  if (i == 0 && j == 0)
  {
    final(0,0) = 1;
    final(1,1) = -1;
    final(2,2) = -1;
    final(3,3) = 1;
    //return QUATERNION( p.x(), -p.y(), -p.z(), p.w());
  }

  //dydy = 
  if (i == 1 && j == 1)
  {
    final(0,0) = -1;
    final(1,1) = 1;
    final(2,2) = -1;
    final(3,3) = 1;
    //return QUATERNION( -p.x() , p.y(), -p.z(), p.w());
  }

  //dzdz = 
  if (i == 2 && j == 2)
  {
    final(0,0) = -1;
    final(1,1) = -1;
    final(2,2) = 1;
    final(3,3) = 1;
    //return QUATERNION( -p.x() , -p.y(), p.z(), p.w());
  }

  //dwdw = 
  if (i == 3 && j == 3)
  {
    final(0,0) = 1;
    final(1,1) = 1;
    final(2,2) = 1;
    final(3,3) = 1;
    //return QUATERNION( p.x() , p.y(), p.z(), p.w());
  }

  //dxdy = 
  if ((i == 0 && j == 1) || (i == 1 && j == 0))
  {
    final(0,1) = 1;
    final(1,0) = 1;
    //return QUATERNION( p.y() , p.x(), 0, 0);
  }

  //dxdz = 
  if ((i == 0 && j == 2) || (i == 2 && j == 0))
  {
    final(0,2) = 1;
    final(2,0) = 1;
    //return QUATERNION( p.z(), 0, p.x(), 0);
  }

  //dxdw= 
  if ((i == 0 && j == 3) || (i == 3 && j == 0))
  {
    final(1,2) = -1;
    final(2,1) = 1;
    //return QUATERNION( 0 , -p.z(), p.y(), 0);
  }

  //dydz= 
  if ((i == 1 && j == 2) || (i == 2 && j == 1))
  {
    final(1,2) = 1;
    final(2,1) = 1;
    //return QUATERNION( 0 , p.z(), p.y(), 0);
  }

  //dydw= 
  if ((i == 1 && j == 3) || (i == 3 && j == 1))
  {
    final(0,2) = 1;
    final(2,0) = -1;
    //return QUATERNION( p.z() , 0, -p.x(), 0);
  }

  //dzdw= 
  if ((i == 2 && j == 3) || (i == 3 && j == 2))
  {
    final(0,1) = -1;
    final(1,0) = 1;
    //return QUATERNION( -p.y(), p.x(), 0, 0);
  }

  final *= 2;
  return final;
}

//////////////////////////////////////////////////////////////////////////////
// Support functions for computeGeorgiiRotation
// Second derivatives from Eqn. 7 from the Georgii and Westermann paper
//////////////////////////////////////////////////////////////////////////////
QUATERNION UNCONSTRAINED_TET_MESH::distanceHessian(QUATERNION& q, QUATERNION& p, int i, int j)
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
// Support functions for computeGeorgiiRotation
// First derivatives from Eqns. 5 and 7 from the Georgii and Westermann paper,
// but phrased as a matrix multiply
//////////////////////////////////////////////////////////////////////////////
MATRIX UNCONSTRAINED_TET_MESH::gradientMatrix(QUATERNION& q, int which)
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
// Support functions for computeGeorgiiRotation
// First derivatives from Eqns. 5 and 7 from the Georgii and Westermann paper
//////////////////////////////////////////////////////////////////////////////
QUATERNION UNCONSTRAINED_TET_MESH::distanceGradient(QUATERNION& q, QUATERNION& p, int which)
{
  switch (which)
  {
    case 0:
       return 2 * QUATERNION( q.x() * p.x() + q.y() * p.y() + q.z() * p.z(),
                              q.y() * p.x() - q.x() * p.y() - q.w() * p.z(),
                              q.z() * p.x() + q.w() * p.y() - q.x() * p.z(),
                              q.x() * p.w());
        break;

    case 1:
        return 2 * QUATERNION( -q.y() * p.x() + q.x() * p.y() + q.w() * p.z(),
                                q.x() * p.x() + q.y() * p.y() + q.z() * p.z(), 
                               -q.w() * p.x() + q.z() * p.y() - q.y() * p.z(), 
                                q.y() * p.w());
        break;

    case 2:
        return 2 * QUATERNION( -q.z() * p.x() - q.w() * p.y() + q.x() * p.z(), 
                                q.w() * p.x() - q.z() * p.y() + q.y() * p.z(), 
                                q.x() * p.x() + q.y() * p.y() + q.z() * p.z(),
                                q.z() * p.w());
        break;
    case 3:
        return 2 * QUATERNION( q.w() * p.x() - q.z() * p.y() + q.y() * p.z(), 
                               q.z() * p.x() + q.w() * p.y() - q.x() * p.z(),
                              -q.y() * p.x() + q.x() * p.y() + q.w() * p.z(),
                               q.w() * p.w());
        break;
  }

  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " distanceGradient has been called incorrectly!" << endl;
  exit(1);
  return QUATERNION();
}

//////////////////////////////////////////////////////////////////////////////
// Support functions for computeGeorgiiRotation
// Eqns. 7 and 8  from the Georgii and Westermann paper
//////////////////////////////////////////////////////////////////////////////
MATRIX UNCONSTRAINED_TET_MESH::energyHessian(QUATERNION& q, Real lambda)
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

  // for each entry in the Hessian
  for (int x = 0; x < 4; x++)
    for (int y = 0; y < 4; y++)
    {
      // compute d - this can be computed outside the loop
      VECTOR deformed = (restVector + _x);
      deformed -= deformedCenter;
      VECTOR rest = restVector - restCenter;

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
      Real final = partialDpartialQi ^ partialDpartialQj;
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

//////////////////////////////////////////////////////////////////////////////
// Support functions for computeGeorgiiRotation
// Eqn. 5 and 6 from the Georgii and Westermann paper
//////////////////////////////////////////////////////////////////////////////
VECTOR UNCONSTRAINED_TET_MESH::energyGradient(QUATERNION& q, Real lambda)
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

  // for each component of the quaternion
  for (int x = 0; x < 4; x++)
  {
    // compute d
    VECTOR deformed = (restVector + _x);
    deformed -= deformedCenter;
    VECTOR rest = restVector - restCenter;

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

//////////////////////////////////////////////////////////////////////////////
// Support functions for computeGeorgiiRotation
// Eqn. 4 from the Georgii and Westermann paper
//////////////////////////////////////////////////////////////////////////////
Real UNCONSTRAINED_TET_MESH::rotationEnergy(QUATERNION& q, Real lambda)
{
  Real final = 0;
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    VEC3F diff = _vertices[x] - _restPose[x];
    QUATERNION distance = rotationDistance(q, _restPose[x], diff);
    final += (distance ^ distance) + lambda * ((q ^ q) - 1);
  }
  final *= 1.0 / _vertices.size();
  return final;
}

//////////////////////////////////////////////////////////////////////////////
// compute the rigid rotation according to Georgii and Westermann's
// "Corotated Finite Elements Made Fast and Stable"
//////////////////////////////////////////////////////////////////////////////
QUATERNION UNCONSTRAINED_TET_MESH::computeGeorgiiRotation()
{
  // might as well warm start it with the shape matching solution
  MATRIX3 guess = computeShapeMatchingRotation();
  QUATERNION q(guess);

  //MATRIX3 I = MATRIX3::I();
  //QUATERNION q(0,0,0,1);
  Real lambda = 0;

  VECTOR gradient = energyGradient(q, lambda);
  //Real magnitude = (gradient ^ gradient) - gradient[4] * gradient[4];

  int count = 0;
  Real convergence = 1;
  while (convergence > 1e-6 && count < 100)
  {
    MATRIX hessian = energyHessian(q, lambda);
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

    // see what the convergence looks like without lambda
    update[4] = 0.0;
    convergence = update ^ update;
    count++;
    cout << " count: " << count << endl;
    cout << " energy: " << rotationEnergy(q, 0);
  }
  _georgiiIterations = count;
  q.normalize();

  return q;
}

//////////////////////////////////////////////////////////////////////
// Draw axes at the center of mass
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_TET_MESH::drawRigidFrame()
{
#if USING_GLVU
  Real bounds[6];
  boundingBox(bounds);
  Real scale = bounds[1] - bounds[0];
  scale = (bounds[3] - bounds[2] > scale) ? bounds[3] - bounds[2] : scale;
  scale = (bounds[5] - bounds[4] > scale) ? bounds[5] - bounds[4] : scale;

  Real angle;
  VEC3F axis;
  QUATERNION copy = _rotationQuaternion;
  copy.normalize();
  copy.axisAngle(axis, angle);

  glPushMatrix();
  glTranslatef(_translation[0], _translation[1], _translation[2]);
  glRotatef(angle, axis[0], axis[1], axis[2]);
  scale *= 0.1;
  glScalef(scale, scale, scale);

  glLineWidth(4.0f);
  glBegin(GL_LINES);
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
#endif
}

//////////////////////////////////////////////////////////////////////
// Draw the constrained nodes as GL points
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_TET_MESH::drawConstrainedNodes()
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
// compute coupling term between rotation and deformation
//////////////////////////////////////////////////////////////////////
MATRIX UNCONSTRAINED_TET_MESH::rotationDefoTensor()
{
  MATRIX final(3, _unconstrainedSize * 3);
  updateFullMesh();

  for (int x = 0; x < _unconstrainedSize; x++)
  {
    MATRIX vertexTilde = MATRIX::cross(_vertices[x]);

    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        final(i,j + 3 * x) = vertexTilde(i,j);
  }

  return final;
}
