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
// INVERTIBLE.h: interface for the INVERTIBLE class.
//
//////////////////////////////////////////////////////////////////////

#ifndef INVERTIBLE_H
#define INVERTIBLE_H

#include <MATERIAL.h>
#include <TIMER.h>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Implements invertible finite elements from [Irving et al. 2004]
//
// The approach here is to filter the deformation gradient as it
// comes and pass it with clamped singular values to a concrete 
// MATERIAL subclass for evaluation.
//////////////////////////////////////////////////////////////////////
class INVERTIBLE : public MATERIAL {

public:
  static int WHICHCASE;
  static int NUMINVERTED;

  //INVERTIBLE(MATERIAL* material, Real epsilon = 0.1);
  INVERTIBLE(MATERIAL* material, bool doFind = true);
  ~INVERTIBLE() { delete _material; };

  // return a copy of this material object
  MATERIAL* copy();
  Real& epsilon() { return _epsilon; };

  // return the material being clamped
  MATERIAL* material() { return _material; };

  MATRIX stiffness(TET& tet, bool diagonal = false);
  MATRIX3 firstPiolaKirchhoff(MATRIX3& F, bool diagonal = false);
  MATRIX3 firstPiolaKirchhoff(MATRIX3& U, MATRIX3& Fhat, MATRIX3& V,
                              bool diagonal = false);
  MATRIX3 secondPiolaKirchhoff(MATRIX3& F, bool diagonal = false);
  //void forceDensity(TET& tet, VEC3F* forces);
  void forceDensity(const Real* F, Real* forces, bool diagonal = false);
  void stiffnessDensity(const Real* F, Real* stiffness, bool diagonal = false);
  virtual TENSOR3 hessian(TET& tet) { return _material->hessian(tet); };

  // diagonalized version
  void stiffnessDensity(const MATRIX3& U, 
                        const MATRIX3& Fhat,
                        const MATRIX3& V,
                        Real* stiffness);
  void stiffnessDensityDebug(const MATRIX3& U, 
                        const MATRIX3& Fhat,
                        const MATRIX3& V,
                        Real* stiffness);

  // diagonalize the deformation gradient
  void diagonalizeF(MATRIX3& F, MATRIX3& U, MATRIX3& Fhat, MATRIX3& V);
  void diagonalizeFDebug(MATRIX3& F, MATRIX3& U, MATRIX3& Fhat, MATRIX3& V);

  void debugF(MATRIX3& F);

  // access counter of how many tets were inverted
  static int& inversions() { return _inversions; };

  Real &stiffnessDensityTime() { return _stiffnessDensityTime; }
  Real &materialStiffnessDensityTime() { return _materialStiffnessDensityTime; }

private:
  MATERIAL* _material;
  Real _epsilon;  // the threshold at which to start clamping F
  Real _implicitEpsilon;  // the threshold at which to start clamping the jacobian
 
  // functions to take the SVD of deformation gradient F
  void removeUReflection(MATRIX3& U, MATRIX3& Fhat);
  void orthogonalizeU(MATRIX3& U, MATRIX3& Fhat);
  void svd(MATRIX3& F, MATRIX3& U, MATRIX3& Fhat, MATRIX3& V);

  Real findRoot(Real find);

  void findClamping(bool scanRightLeft = false);

  Real _slope;
  Real _intercept;
  static int _inversions;

  Real _stiffnessDensityTime;
  Real _materialStiffnessDensityTime;

  virtual void hessian0(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian)
  { _material->hessian0(DmInv, b, vertices, hessian); };
  virtual void hessian1(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian)
  { _material->hessian1(DmInv, b, vertices, hessian); };
  virtual void hessian2(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian)
  { _material->hessian2(DmInv, b, vertices, hessian); };
  virtual void hessian3(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian)
  { _material->hessian3(DmInv, b, vertices, hessian); };
  virtual void hessian4(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian)
  { _material->hessian4(DmInv, b, vertices, hessian); };
  virtual void hessian5(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian)
  { _material->hessian5(DmInv, b, vertices, hessian); };
  virtual void hessian6(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian)
  { _material->hessian6(DmInv, b, vertices, hessian); };
  virtual void hessian7(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian)
  { _material->hessian7(DmInv, b, vertices, hessian); };
  virtual void hessian8(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian)
  { _material->hessian8(DmInv, b, vertices, hessian); };
  virtual void hessian9(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian)
  { _material->hessian9(DmInv, b, vertices, hessian); };
  virtual void hessian10(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian)
  { _material->hessian10(DmInv, b, vertices, hessian); };
  virtual void hessian11(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian)
  { _material->hessian11(DmInv, b, vertices, hessian); };

  virtual void hessian0(const Real* F, const VEC3F* b, MATRIX& hessian)
  { _material->hessian0(F, b, hessian); };
  virtual void hessian1(const Real* F, const VEC3F* b, MATRIX& hessian)
  { _material->hessian1(F, b, hessian); };
  virtual void hessian2(const Real* F, const VEC3F* b, MATRIX& hessian)
  { _material->hessian2(F, b, hessian); };
  virtual void hessian3(const Real* F, const VEC3F* b, MATRIX& hessian)
  { _material->hessian3(F, b, hessian); };
  virtual void hessian4(const Real* F, const VEC3F* b, MATRIX& hessian)
  { _material->hessian4(F, b, hessian); };
  virtual void hessian5(const Real* F, const VEC3F* b, MATRIX& hessian)
  { _material->hessian5(F, b, hessian); };
  virtual void hessian6(const Real* F, const VEC3F* b, MATRIX& hessian)
  { _material->hessian6(F, b, hessian); };
  virtual void hessian7(const Real* F, const VEC3F* b, MATRIX& hessian)
  { _material->hessian7(F, b, hessian); };
  virtual void hessian8(const Real* F, const VEC3F* b, MATRIX& hessian)
  { _material->hessian8(F, b, hessian); };

#if 0
  // Workspaces for stiffness density calculation.  Variable names
  // correspond to variable names from the original stiffness
  // density calculation.  When we append _*, it means that we
  // need multiple copies of the matrix with this name, due to
  // matrix mutliplications, assignments, etc. in the calculation
  MATRIX     _workspaceA_1, _workspaceA_2, _workspaceA_3;
  MATRIX     _workspaceB12_1, workspaceB12_2, workspaceB12_3;
  MATRIX     _workspaceB13;
  MATRIX     _workspaceB23;

  MATRIX     _workspaceEigenvectors3x3_1, _workspaceEigenvectors3x3_2;
  VECTOR     _workspaceEigenvalues3x3;
  MATRIX     _workspaceEigenvectors2x2_1, _workspaceEigenvectors2x2_2;
  VECTOR     _workspaceEigenvalues2x2;

  MATRIX     _workspaceValueMatrix;
#endif

};

#endif


