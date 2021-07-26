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
// MULTI_MATERIAL.cpp: implementation of the MULTI_MATERIAL class.
//
//////////////////////////////////////////////////////////////////////

#include "MULTI_MATERIAL.h"

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
MULTI_MATERIAL::MULTI_MATERIAL() :
  MATERIAL()
{
  _materialName.assign("MULTI_MATERIAL");
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
MATRIX MULTI_MATERIAL::stiffness(TET& tet, bool diagonal)
{
  MATRIX final(12,12);
  for (int x = 0; x < _materials.size(); x++)
    final += _materials[x]->stiffness(tet, diagonal);

  return final;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
MATRIX3 MULTI_MATERIAL::firstPiolaKirchhoff(MATRIX3& F, bool diagonal)
{
  MATRIX3 final;
  for (int x = 0; x < _materials.size(); x++)
    final += _materials[x]->firstPiolaKirchhoff(F, diagonal);

  return final;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
MATRIX3 MULTI_MATERIAL::secondPiolaKirchhoff(MATRIX3& F, bool diagonal)
{ 
  MATRIX3 final;
  for (int x = 0; x < _materials.size(); x++)
    final += _materials[x]->secondPiolaKirchhoff(F, diagonal);

  return final;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
MATRIX3 MULTI_MATERIAL::cauchyStress(MATRIX3& F, bool diagonal)
{
  MATRIX3 final;
  for (int x = 0; x < _materials.size(); x++)
    final += _materials[x]->cauchyStress(F, diagonal);

  return final;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void MULTI_MATERIAL::stiffnessDensity(const Real* F, Real* stiffness,
                                      bool diagonal)
{
  // clear the initial stiffness
  for (int x = 0; x < 81; x++)
    stiffness[x] = 0.0;

  Real data[81];
  for (int x = 0; x < _materials.size(); x++)
  {
    _materials[x]->stiffnessDensity(F, data, diagonal);
    for (int y = 0; y < 81; y++)
      stiffness[y] += data[y];
  }
}

//////////////////////////////////////////////////////////////////////
// this is just a flattened version of first Piola-Kirchoff, and should
// be renamed at some point
//////////////////////////////////////////////////////////////////////
void MULTI_MATERIAL::forceDensity(const Real* F, Real* forces, bool diagonal)
{
  // clear the initial stiffness
  for (int x = 0; x < 81; x++)
    forces[x] = 0.0;

  Real data[81];
  for (int x = 0; x < _materials.size(); x++)
  {
    _materials[x]->forceDensity(F, data, diagonal);
    for (int y = 0; y < 81; y++)
      forces[y] += data[y];
  }
}

//////////////////////////////////////////////////////////////////////
// return a copy of this material object
//////////////////////////////////////////////////////////////////////
MATERIAL* MULTI_MATERIAL::copy()
{
  MULTI_MATERIAL* newMaterial = new MULTI_MATERIAL();
  for (int x = 0; x < _materials.size(); x++)
    newMaterial->addMaterial(_materials[x]->copy());

  return newMaterial;
}

//////////////////////////////////////////////////////////////////////
// get the hessian tensor
//////////////////////////////////////////////////////////////////////
TENSOR3 MULTI_MATERIAL::hessian(TET& tet)
{
  TENSOR3 final;
  for (int x = 0; x < _materials.size(); x++)
  {
    TENSOR3 temp = _materials[x]->hessian(tet);
    final = final + temp;
  }

  return final;
}
