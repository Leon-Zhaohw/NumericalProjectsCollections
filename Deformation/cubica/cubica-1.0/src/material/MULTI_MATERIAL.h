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
// MULTI_MATERIAL.h: interface for the MULTI_MATERIAL class.
//
//////////////////////////////////////////////////////////////////////

#ifndef MULTI_MATERIAL_H
#define MULTI_MATERIAL_H

#include <MATERIAL.h>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Abstract material class
//////////////////////////////////////////////////////////////////////
class MULTI_MATERIAL : public MATERIAL {

public:
  MULTI_MATERIAL();
  virtual ~MULTI_MATERIAL() {};

  virtual MATRIX stiffness(TET& tet, bool diagonal = false);
  virtual MATRIX3 firstPiolaKirchhoff(MATRIX3& F, bool diagonal = false);
  virtual MATRIX3 secondPiolaKirchhoff(MATRIX3& F, bool diagonal = false);
  virtual MATRIX3 cauchyStress(MATRIX3& F, bool diagonal = false);
  virtual void stiffnessDensity(const Real* F, Real* stiffness,
                                bool diagonal = false);

  // this is just a flattened version of first Piola-Kirchoff, and should
  // be renamed at some point
  virtual void forceDensity(const Real* F, Real* forces,
                            bool diagonal = false);

  // return a copy of this material object
  virtual MATERIAL* copy();

  // get the hessian tensor
  virtual TENSOR3 hessian(TET& tet);

  void addMaterial(MATERIAL* material) { _materials.push_back(material); };

public:  
  virtual void hessian0(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian) {};
  virtual void hessian1(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian) {};
  virtual void hessian2(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian) {};
  virtual void hessian3(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian) {};
  virtual void hessian4(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian) {};
  virtual void hessian5(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian) {};
  virtual void hessian6(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian) {};
  virtual void hessian7(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian) {};
  virtual void hessian8(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian) {};
  virtual void hessian9(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian) {};
  virtual void hessian10(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian) {};
  virtual void hessian11(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian) {};

  virtual void hessian0(const Real* F, const VEC3F* b, MATRIX& hessian) {};
  virtual void hessian1(const Real* F, const VEC3F* b, MATRIX& hessian) {};
  virtual void hessian2(const Real* F, const VEC3F* b, MATRIX& hessian) {};
  virtual void hessian3(const Real* F, const VEC3F* b, MATRIX& hessian) {};
  virtual void hessian4(const Real* F, const VEC3F* b, MATRIX& hessian) {};
  virtual void hessian5(const Real* F, const VEC3F* b, MATRIX& hessian) {};
  virtual void hessian6(const Real* F, const VEC3F* b, MATRIX& hessian) {};
  virtual void hessian7(const Real* F, const VEC3F* b, MATRIX& hessian) {};
  virtual void hessian8(const Real* F, const VEC3F* b, MATRIX& hessian) {};

private:
  vector<MATERIAL*> _materials;
};

#endif
