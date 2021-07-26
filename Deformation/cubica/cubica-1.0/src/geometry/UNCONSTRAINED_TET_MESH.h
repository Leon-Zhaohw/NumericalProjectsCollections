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
// This child class of TET_MESH handles unconstrained meshes by
// tracking the rigid and deformable components separately. The
// original TET_MESH class handles unconstrained meshes just fine,
// but this class is a lead-up to handling reduced unconstrained meshes.
//
// The methodology is from:
//
// "Physically Based Models with Rigid and Deformable Components"
// Terzopoulos and Witkin, IEEE CG&A, 1988
//
// Their notation is followed wherever possible. Unless otherwise
// noted, the equation numbers refer to the numbers from that paper
//
//////////////////////////////////////////////////////////////////////

#ifndef UNCONSTRAINED_TET_MESH_H
#define UNCONSTRAINED_TET_MESH_H

#include "TET_MESH.h"
#include <QUATERNION.h>
#include <BLOCK_MATRIX.h>

//////////////////////////////////////////////////////////////////////
// Unconstrained tetrahedron mesh class that tracks the rigid
// component explicitly
//////////////////////////////////////////////////////////////////////
class UNCONSTRAINED_TET_MESH : public TET_MESH {

public:
  // "filename" is the tet mesh file
  // "materials" is the array of constitutive models
  // "totalMaterials" is the total number of materials in the model
  // "simulate" decides if we want to allocate the gigantic matrices for
  //            direct simulation
  UNCONSTRAINED_TET_MESH(const char* filename, MATERIAL** materials, int totalMaterials, bool simulate = false);
  virtual ~UNCONSTRAINED_TET_MESH() {};

  // accessors
  VEC3F& rigidTranslation() { return _translation; };
  VEC3F& rigidTranslationOld() { return _translationOld; };
  //MATRIX3& rigidRotation() { return _rotation; };
  //MATRIX3& rigidRotationOld() { return _rotationOld; };
  QUATERNION& rotationQuaternion() { return _rotationQuaternion; };
  QUATERNION& rotationQuaternionOld() { return _rotationQuaternionOld; };
  Real& rotationLambda() { return _rotationLambda; };
  
  // drawing routines, taking into account the rigid components
  void drawAllTets();
  void drawSurfaceFaces();
  void drawFacesWithoutRigidTransform();
  void drawRestPose();
  void drawRigidFrame();
  virtual void drawConstrainedNodes();

  // find closest node to an arbitrary point, taking into account
  // the rigid components
  virtual VEC3F* closestSurfaceNode(VEC3F point);

  // use shape matching to extract a rotation from the current deformation
  MATRIX3 computeShapeMatchingRotation();

  // compute the rigid rotation according to Georgii and Westermann's
  // "Corotated Finite Elements Made Fast and Stable"
  QUATERNION computeGeorgiiRotation();

  // Support functions for computeGeorgiiRotation, but it's also
  // useful to know the rotation energy when debugging
  Real rotationEnergy(QUATERNION& q, Real lambda = 0);

  // query how many Newton iterations the last georgii rotation solve took
  const int georgiiIterations() { return _georgiiIterations; };

  // compute coupling term between rotation and deformation
  MATRIX rotationDefoTensor();

protected:
  ///////////////////////////////////////////////////////////////////
  // Hybrid Rigid/Deformable variables and functions
  ///////////////////////////////////////////////////////////////////
  VEC3F _translationOld;    // rigid translation
  VEC3F _translation;       // rigid translation
  //MATRIX3 _rotation;        // rigid rotation;
  //MATRIX3 _rotationOld;
  QUATERNION _rotationQuaternion; // unnormalized rigid rotation quaternion
  QUATERNION _rotationQuaternionOld; // unnormalized rigid rotation quaternion
  Real _rotationLambda;

  ///////////////////////////////////////////////////////////////////
  // Support functions for computeGeorgiiRotation
  ///////////////////////////////////////////////////////////////////
  QUATERNION rotationDistance(QUATERNION& q, VEC3F& rest, VEC3F& displacement);
  QUATERNION distanceHessian(QUATERNION& q, QUATERNION& p, int i, int j);
  QUATERNION distanceGradient(QUATERNION& q, QUATERNION& p, int which);
  MATRIX energyHessian(QUATERNION& q, Real lambda);
  VECTOR energyGradient(QUATERNION& q, Real lambda);
  MATRIX gradientMatrix(QUATERNION& q, int which);
  MATRIX hessianMatrix(QUATERNION& q, int i, int j);

  int _georgiiIterations; // track how many Newton iterations the rotation took
};

#endif
