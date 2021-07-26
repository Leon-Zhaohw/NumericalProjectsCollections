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
// STVK.h: interface for the STVK class.
//
//////////////////////////////////////////////////////////////////////

#ifndef STVK_H
#define STVK_H

#include <TET.h>
#include <MATERIAL.h>

using namespace std;

//////////////////////////////////////////////////////////////////////
// St-VK material class
//////////////////////////////////////////////////////////////////////
class STVK : public MATERIAL {

public:
  STVK(Real lambda, Real mu);
  ~STVK() {};

  MATRIX stiffness(TET& tet, bool diagonal = false);
  MATRIX3 secondPiolaKirchhoff(MATRIX3& F, bool diagonal = false);
  void forceDensity(TET& tet, VEC3F* forces, bool diagonal = false);
  void forceDensity(const Real* F, Real* forces, bool diagonal = false);
  void stiffnessDensity(const Real* F, Real* stiffness, bool diagonal = false);
  
  MATERIAL* copy();

private:
  Real _lambda;
  Real _mu;

  // work arrays
  MATRIX _pf_pF;
  MATRIX _pF_pu;

  void computeStresses(TET& tet);
  virtual void hessian0(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian);
  virtual void hessian1(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian);
  virtual void hessian2(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian);
  virtual void hessian3(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian);
  virtual void hessian4(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian);
  virtual void hessian5(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian);
  virtual void hessian6(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian);
  virtual void hessian7(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian);
  virtual void hessian8(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian);
  virtual void hessian9(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian);
  virtual void hessian10(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian);
  virtual void hessian11(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian);

  virtual void hessian0(const Real* F, const VEC3F* b, MATRIX& hessian);
  virtual void hessian1(const Real* F, const VEC3F* b, MATRIX& hessian);
  virtual void hessian2(const Real* F, const VEC3F* b, MATRIX& hessian);
  virtual void hessian3(const Real* F, const VEC3F* b, MATRIX& hessian);
  virtual void hessian4(const Real* F, const VEC3F* b, MATRIX& hessian);
  virtual void hessian5(const Real* F, const VEC3F* b, MATRIX& hessian);
  virtual void hessian6(const Real* F, const VEC3F* b, MATRIX& hessian);
  virtual void hessian7(const Real* F, const VEC3F* b, MATRIX& hessian);
  virtual void hessian8(const Real* F, const VEC3F* b, MATRIX& hessian);
};

#endif
