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
// INVERTIBLE.cpp: implementation of the INVERTIBLE class.
//
//////////////////////////////////////////////////////////////////////

#include "INVERTIBLE.h"
#include <float.h>

int INVERTIBLE::_inversions = 0;

//////////////////////////////////////////////////////////////////////
// Constructor for INVERTIBLE
//////////////////////////////////////////////////////////////////////
INVERTIBLE::INVERTIBLE(MATERIAL* material, bool doFind) :
  _material(material),
  //_implicitEpsilon(0.01),
  //_implicitEpsilon(1e-8)
  _implicitEpsilon(0.0)
{
  _materialName = string("INVERTIBLE");
  //_epsilon = 0.1;

  if (doFind)
  {
    cout << " Material name: " << _material->materialName().c_str() << endl;
    /*
    if (_material->materialName().compare("OGDEN") == 0)
    {
      OGDEN* ogden = (OGDEN*)_material;
      ogden->singularity();
      findClamping(true);
      //findClamping();
    }
    else
    */
      findClamping();

    /*
    cout << __FILE__ << " " << __LINE__ << " : " << endl;
    cout << " SETTING INVERSION EPSILON BY HAND" << endl;
    cout << " found epsilon: " << _epsilon << endl;
    char a;
    cin >> a;
    */
  }
  //_epsilon = 0.7;

  /*
  MATRIX3 F = MATRIX3::I();
  F(0,0) = 0.1;
  MATRIX3 PK = _material->firstPiolaKirchhoff(F);
  cout << " PK at 0.1: " << PK(0,0) << endl;
  */

  //cout << " Found root for material: " << _epsilon << endl;
}

//////////////////////////////////////////////////////////////////////
// make a copy
//////////////////////////////////////////////////////////////////////
MATERIAL* INVERTIBLE::copy()
{
  MATERIAL* baseMaterial = _material->copy();
  INVERTIBLE* toReturn = new INVERTIBLE(baseMaterial, false);

  toReturn->_slope = _slope;
  toReturn->_intercept = _intercept;
  toReturn->_epsilon = _epsilon;

  return toReturn;
}

//////////////////////////////////////////////////////////////////////
// stiffness matrix implementation
//
// Does this need to be filtered? The tet being near inversion in the
// rest pose implies it is ill-conditioned, which shouldn't be happening
// with isosurface stuffing.
//////////////////////////////////////////////////////////////////////
MATRIX INVERTIBLE::stiffness(TET& tet, bool diagonal)
{
  MATRIX3 F = tet.F();
  MATRIX3 U;
  MATRIX3 Fhat;
  MATRIX3 V;
  Real stiffnessData[81];
  diagonalizeF(F, U, Fhat, V);
  stiffnessDensity(U, Fhat, V, stiffnessData);
  MATRIX diagonalStiffness(9,9,stiffnessData);

  // compute Bm matrix
  const VEC3F* b = tet.b();
  
  // get PFPu from the tet - only depends on rest state,
  // so don't need to update tet state at all
  MATRIX pFpu(9,12);
  computePFPu(tet, pFpu);
  MATRIX stiffnessFinal(12,12);

  // compute each column of the matrix
  for (int y = 0; y < 12; y++)
  {
    // extract a column from pFpu (ie pretend 'x' is just a 
    // Kronecker delta)
    VECTOR deltaF(9);
    for (int z = 0; z < 9; z++)
      deltaF(z) = pFpu(z, y);

    // rotate deltaF
    MATRIX3 rotated = U.transpose() * TET::repackF(deltaF) * V;
    deltaF = TET::flattenF(rotated);

    VECTOR contraction = diagonalStiffness  * deltaF;
    MATRIX3 deltaP = TET::repackF(contraction);

    // rotate deltaP back
    deltaP = U * deltaP * V.transpose();
    VEC3F forceVecs[4];
    forceVecs[0] = deltaP * b[0];
    forceVecs[1] = deltaP * b[1];
    forceVecs[2] = deltaP * b[2];
    forceVecs[3] = deltaP * b[3];
   
    // copy result into stiffness column
    for (int z = 0; z < 4; z++)
      for (int a = 0; a < 3; a++)
        stiffnessFinal(z * 3 + a, y) = forceVecs[z][a];
  }
  stiffnessFinal *= -1.0;

  return stiffnessFinal;
}

//////////////////////////////////////////////////////////////////////
// implementation of first PK stress tensor
//////////////////////////////////////////////////////////////////////
MATRIX3 INVERTIBLE::firstPiolaKirchhoff(MATRIX3& F, bool diagonal)
{
  // diagonalize F
  MATRIX3 U;
  MATRIX3 Fhat;
  MATRIX3 V;
  svd(F, U, Fhat, V);

  // create a filtered version, just for Phat
  MATRIX3 FhatFiltered = Fhat;
  bool filtered = false;
  for (int x = 0; x < 3; x++)
    if (FhatFiltered(x,x) < _epsilon)
    {
      FhatFiltered(x,x) = _epsilon;
      x = 4;
      filtered = true;
    }

  // this will create a clamping at the compression limit
  MATRIX3 Phat = _material->firstPiolaKirchhoff(FhatFiltered);

  // linearize by explicitly using the stiffness density matrix
  if (filtered)
  {
    // get the stiffness density at the linearization point
    Real stiffnessData[81];
    stiffnessDensity(U, FhatFiltered, V, stiffnessData);
    MATRIX stiffness(9, 9, stiffnessData);

    // recenter the origin
    MATRIX3 recenteredF = Fhat - FhatFiltered;
    
    // multiply by the original unclamped F to linearize
    VECTOR flatF = TET::flattenF(recenteredF);
    VECTOR linearization = stiffness * flatF;
    linearization *= -1.0;
    MATRIX3 repackedP = TET::repackF(linearization);

    Phat += repackedP;
  }

  /*
  // linearize past the compression limit
  bool epsFound = false;
  for (int x = 0; x < 3; x++)
  {
    if (Fhat(x,x) < _epsilon && !epsFound)
    {
      Phat(x,x) = _intercept + _slope * Fhat(x,x);
      epsFound = true;
    }
  }
  */
  
  return U * Phat * V.transpose();
}

//////////////////////////////////////////////////////////////////////
// Print out various stats about F
//////////////////////////////////////////////////////////////////////
void INVERTIBLE::debugF(MATRIX3& F)
{
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " P without filtering: " << _material->firstPiolaKirchhoff(F) << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;

  // diagonalize F
  MATRIX3 U;
  MATRIX3 Fhat;
  MATRIX3 V;
  svd(F, U, Fhat, V);

  cout << " Decomposition: " << endl;
  cout << " U: " << U << endl;
  cout << " Fhat: " << Fhat << endl;
  cout << " V: " << V << endl;
 
  // filter Fhat
  for (int x = 0; x < 3; x++)
  {
    if (Fhat(x,x) < _epsilon)
    {
      Fhat(x,x) = _epsilon;
      x = 4;
    }
  }
  cout << " Fhat Filtered: " << Fhat << endl;
  
  MATRIX3 Phat = _material->firstPiolaKirchhoff(Fhat);

  // linearize past compression limit
  for (int x = 0; x < 3; x++)
  {
    if (Fhat(x,x) < 0.0)
      Phat(x,x) = _intercept + _slope * Fhat(x,x);
  }
  
  cout << " Phat Filtered: " << Phat << endl;
  cout << " Rotated P: " << U * Phat * V.transpose() << endl;
}

//////////////////////////////////////////////////////////////////////
// implementation of first PK stress tensor
//////////////////////////////////////////////////////////////////////
MATRIX3 INVERTIBLE::firstPiolaKirchhoff(MATRIX3& U, MATRIX3& Fhat, MATRIX3& V, bool diagonal)
{
  // create a filtered version, just for Phat
  MATRIX3 FhatFiltered = Fhat;
  bool filtered = false;
  for (int x = 0; x < 3; x++)
    if (FhatFiltered(x,x) < _epsilon)
    {
      FhatFiltered(x,x) = _epsilon;
      x = 4;
      filtered = true;
    }

  // this will create a clamping at the compression limit
  MATRIX3 Phat = _material->firstPiolaKirchhoff(FhatFiltered);

  // linearize by explicitly using the stiffness density matrix
  if (filtered)
  {
    // get the stiffness density at the linearization point
    Real stiffnessData[81];
    stiffnessDensity(U, FhatFiltered, V, stiffnessData);
    MATRIX stiffness(9, 9, stiffnessData);

    // recenter the origin
    MATRIX3 recenteredF = Fhat - FhatFiltered;
    
    // multiply by the original unclamped F to linearize
    VECTOR flatF = TET::flattenF(recenteredF);
    VECTOR linearization = stiffness * flatF;
    linearization *= -1.0;
    MATRIX3 repackedP = TET::repackF(linearization);

    Phat += repackedP;
  } 

  /*
  // linearize past the compression limit
  bool epsFound = false;
  for (int x = 0; x < 3; x++)
  {
    if (Fhat(x,x) < _epsilon && !epsFound)
    {
      Phat(x,x) = _intercept + _slope * Fhat(x,x);
      epsFound = true;
    }
  }
  */

  return U * Phat * V.transpose();
}

//////////////////////////////////////////////////////////////////////
// diagonalize the deformation gradient
//////////////////////////////////////////////////////////////////////
void INVERTIBLE::diagonalizeF(MATRIX3& F, MATRIX3& U, MATRIX3& Fhat, MATRIX3& V)
{
  // diagonalize F
  svd(F, U, Fhat, V);

  // filter Fhat
  //
  // Still filter F here in lieu of the first PK because it is later used
  // for the stiffnessDensity call in the implicit solve, which does not
  // call the firstPK directly. Some force limiting needs to happen somewhere.
  for (int x = 0; x < 3; x++)
  {
    if (Fhat(x,x) < _epsilon)
    {
      Fhat(x,x) = _epsilon;
      x = 4;
      _inversions++;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// diagonalize the deformation gradient
//////////////////////////////////////////////////////////////////////
void INVERTIBLE::diagonalizeFDebug(MATRIX3& F, MATRIX3& U, MATRIX3& Fhat, MATRIX3& V)
{
  // diagonalize F
  svd(F, U, Fhat, V);
}

//////////////////////////////////////////////////////////////////////
// implementation of second PK stress tensor
//////////////////////////////////////////////////////////////////////
MATRIX3 INVERTIBLE::secondPiolaKirchhoff(MATRIX3& F, bool diagonal)
{
  return F.inverse() * firstPiolaKirchhoff(F);
}

//////////////////////////////////////////////////////////////////////
// Do the computationally expensive thing and compute the diagonalization
// from scratch. If this becomes the bottleneck, the calling object
// should pass the diagonalization explicitly to the other version of
// "stiffnessDensity".
//////////////////////////////////////////////////////////////////////
void INVERTIBLE::stiffnessDensity(const Real* F, Real* stiffness, bool diaongal)
{
  MATRIX3 F3 = TET::repackF(F);
  MATRIX3 U;
  MATRIX3 Fhat;
  MATRIX3 V;
  svd(F3, U, Fhat, V);
  
  stiffnessDensity(U, Fhat, V, stiffness);
  //_material->stiffnessDensity(F, stiffness);
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void INVERTIBLE::stiffnessDensity(const MATRIX3& U, 
                                  const MATRIX3& Fhat,
                                  const MATRIX3& V,
                                  Real* stiffness)
{
  MATRIX3 FhatFiltered = Fhat;
  for (int x = 0; x < 3; x++)
    if (FhatFiltered(x,x) < _epsilon)
    {
      FhatFiltered(x,x) = _epsilon;
      x = 4;
    }

  Real FFiltered[9];
  for (int x = 0; x < 9; x++)
    FFiltered[x] = 0.0;
  FFiltered[0] = FhatFiltered(0,0);
  FFiltered[4] = FhatFiltered(1,1);
  FFiltered[8] = FhatFiltered(2,2);
  _material->stiffnessDensity(FFiltered, stiffness);

  // gather 
  MATRIX A(3,3);
  A(0,0) = -stiffness[0];
  A(0,1) = -stiffness[4];
  A(0,2) = -stiffness[8];
  A(1,0) = -stiffness[0 + 4 * 9];
  A(1,1) = -stiffness[4 + 4 * 9];
  A(1,2) = -stiffness[8 + 4 * 9];
  A(2,0) = -stiffness[0 + 8 * 9];
  A(2,1) = -stiffness[4 + 8 * 9];
  A(2,2) = -stiffness[8 + 8 * 9];

  MATRIX B12(2,2);
  B12(0,0) = -stiffness[1 + 1 * 9];
  B12(0,1) = -stiffness[3 + 1 * 9];
  B12(1,0) = -stiffness[1 + 3 * 9];
  B12(1,1) = -stiffness[3 + 3 * 9];

  MATRIX B13(2,2);
  B13(0,0) = -stiffness[2 + 2 * 9];
  B13(0,1) = -stiffness[6 + 2 * 9];
  B13(1,0) = -stiffness[2 + 6 * 9];
  B13(1,1) = -stiffness[6 + 6 * 9];
  
  MATRIX B23(2,2);
  B23(0,0) = -stiffness[5 + 5 * 9];
  B23(0,1) = -stiffness[7 + 5 * 9];
  B23(1,0) = -stiffness[5 + 7 * 9];
  B23(1,1) = -stiffness[7 + 7 * 9];

  MATRIX eigenvectors3x3(3,3);
  VECTOR eigenvalues3x3(3);
  MATRIX eigenvectors2x2(2,2);
  VECTOR eigenvalues2x2(2);
 
  // clamp and rebuild A
  int clamped = 0;
#if USING_MKL
  A.eigensystem(eigenvalues3x3, eigenvectors3x3);
#elif USING_OSX
  A.eigensystem(eigenvalues3x3, eigenvectors3x3);
#else
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " WARNING: USING A LOW PRECISION EIGENSOLVER" << endl;
  A.eigensystem3x3(eigenvalues3x3, eigenvectors3x3);
#endif
  for (int x = 0; x < 3; x++)
    if (eigenvalues3x3(x) < _implicitEpsilon) 
    {
      //cout << " Bad eigenvalue: " << eigenvalues3x3(x) << endl;
      eigenvalues3x3(x) = _implicitEpsilon;
      clamped++;
    }
  // GCC doesn't like this form
  //A = eigenvectors3x3 * valueMatrix * eigenvectors3x3.transpose();
  MATRIX valueMatrix(eigenvalues3x3);
  A = eigenvectors3x3 * valueMatrix;
  eigenvectors3x3 = eigenvectors3x3.transpose();
  A = A * eigenvectors3x3;
 
  // clamp and rebuild B12
  B12.eigensystem2x2(eigenvalues2x2, eigenvectors2x2);
  if (eigenvalues2x2(0) < _implicitEpsilon) 
  {
    //cout << " Bad eigenvalue: " << eigenvalues2x2(0) << endl;
    eigenvalues2x2(0) = _implicitEpsilon;
    clamped++;
  }
  if (eigenvalues2x2(1) < _implicitEpsilon) 
  {
    //cout << " Bad eigenvalue: " << eigenvalues2x2(1) << endl;
    eigenvalues2x2(1) = _implicitEpsilon;
    clamped++;
  }
  // GCC doesn't like this form
  //B12 = eigenvectors2x2 * MATRIX(eigenvalues2x2) * eigenvectors2x2.transpose();
  valueMatrix = MATRIX(eigenvalues2x2);
  B12 = eigenvectors2x2 * valueMatrix;
  eigenvectors2x2 = eigenvectors2x2.transpose();
  B12 = B12 * eigenvectors2x2;

  // clamp and rebuild B13
  B13.eigensystem2x2(eigenvalues2x2, eigenvectors2x2);
  if (eigenvalues2x2(0) < _implicitEpsilon) 
  {
    //cout << " Bad eigenvalue: " << eigenvalues2x2(0) << endl;
    eigenvalues2x2(0) = _implicitEpsilon;
    clamped++;
  }
  if (eigenvalues2x2(1) < _implicitEpsilon) 
  {
    //cout << " Bad eigenvalue: " << eigenvalues2x2(1) << endl;
    eigenvalues2x2(1) = _implicitEpsilon;
    clamped++;
  }
  // GCC doesn't like this form
  //B13 = eigenvectors2x2 * MATRIX(eigenvalues2x2) * eigenvectors2x2.transpose();
  valueMatrix = MATRIX(eigenvalues2x2);
  B13 = eigenvectors2x2 * valueMatrix;
  eigenvectors2x2 = eigenvectors2x2.transpose();
  B13 = B13 * eigenvectors2x2;

  // clamp and rebuild B23
  B23.eigensystem2x2(eigenvalues2x2, eigenvectors2x2);
  if (eigenvalues2x2(0) < _implicitEpsilon) 
  {
    //cout << " Bad eigenvalue: " << eigenvalues2x2(0) << endl;
    eigenvalues2x2(0) = _implicitEpsilon;
    clamped++;
  }
  if (eigenvalues2x2(1) < _implicitEpsilon) 
  {
    //cout << " Bad eigenvalue: " << eigenvalues2x2(1) << endl;
    eigenvalues2x2(1) = _implicitEpsilon;
    clamped++;
  }
  // GCC doesn't like this form
  //B23 = eigenvectors2x2 * MATRIX(eigenvalues2x2) * eigenvectors2x2.transpose();
  valueMatrix = MATRIX(eigenvalues2x2);
  B23 = eigenvectors2x2 * valueMatrix;
  eigenvectors2x2 = eigenvectors2x2.transpose();
  B23 = B23 * eigenvectors2x2;
  
  // scatter
  stiffness[0] = -A(0,0);
  stiffness[4] = -A(0,1);
  stiffness[8] = -A(0,2);
  stiffness[0 + 4 * 9] = -A(1,0);
  stiffness[4 + 4 * 9] = -A(1,1);
  stiffness[8 + 4 * 9] = -A(1,2);
  stiffness[0 + 8 * 9] = -A(2,0);
  stiffness[4 + 8 * 9] = -A(2,1);
  stiffness[8 + 8 * 9] = -A(2,2);
  stiffness[1 + 1 * 9] = -B12(0,0);
  stiffness[3 + 1 * 9] = -B12(0,1);
  stiffness[1 + 3 * 9] = -B12(1,0);
  stiffness[3 + 3 * 9] = -B12(1,1);
  stiffness[2 + 2 * 9] = -B13(0,0);
  stiffness[6 + 2 * 9] = -B13(0,1);
  stiffness[2 + 6 * 9] = -B13(1,0);
  stiffness[6 + 6 * 9] = -B13(1,1);
  stiffness[5 + 5 * 9] = -B23(0,0);
  stiffness[7 + 5 * 9] = -B23(0,1);
  stiffness[5 + 7 * 9] = -B23(1,0);
  stiffness[7 + 7 * 9] = -B23(1,1);

  if (clamped > 0)
    _inversions += clamped;

  /*
  // DEBUG -- sanity check
  MATRIX stiffnessSanity(9,9, stiffness);
  MATRIX eigenvectors9x9(9,9);
  VECTOR eigenvalues9x9(9);
  stiffnessSanity.eigensystem(eigenvalues9x9, eigenvectors9x9);

  Real smallest = 1000000.0;
  for (int x = 0; x < 9; x++)
    if (fabs(eigenvalues9x9(x)) < smallest)
      smallest = fabs(eigenvalues9x9(x));
 
  if (smallest < 0.01)
  {
    // clamp and rebuild A
    cout << " stiffness eigenvalues: " << eigenvalues9x9 << endl;
  }
  */
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void INVERTIBLE::stiffnessDensityDebug(const MATRIX3& U, 
                                  const MATRIX3& Fhat,
                                  const MATRIX3& V,
                                  Real* stiffness)
{
  Real FFiltered[9];
  for (int x = 0; x < 9; x++)
    FFiltered[x] = 0.0;
  FFiltered[0] = Fhat(0,0);
  FFiltered[4] = Fhat(1,1);
  FFiltered[8] = Fhat(2,2);
  _material->stiffnessDensity(FFiltered, stiffness);

  /*
  // gather 
  MATRIX A(3,3);
  A(0,0) = -stiffness[0];
  A(0,1) = -stiffness[4];
  A(0,2) = -stiffness[8];
  A(1,0) = -stiffness[0 + 4 * 9];
  A(1,1) = -stiffness[4 + 4 * 9];
  A(1,2) = -stiffness[8 + 4 * 9];
  A(2,0) = -stiffness[0 + 8 * 9];
  A(2,1) = -stiffness[4 + 8 * 9];
  A(2,2) = -stiffness[8 + 8 * 9];

  MATRIX B12(2,2);
  B12(0,0) = -stiffness[1 + 1 * 9];
  B12(0,1) = -stiffness[3 + 1 * 9];
  B12(1,0) = -stiffness[1 + 3 * 9];
  B12(1,1) = -stiffness[3 + 3 * 9];

  MATRIX B13(2,2);
  B13(0,0) = -stiffness[2 + 2 * 9];
  B13(0,1) = -stiffness[6 + 2 * 9];
  B13(1,0) = -stiffness[2 + 6 * 9];
  B13(1,1) = -stiffness[6 + 6 * 9];
  
  MATRIX B23(2,2);
  B23(0,0) = -stiffness[5 + 5 * 9];
  B23(0,1) = -stiffness[7 + 5 * 9];
  B23(1,0) = -stiffness[5 + 7 * 9];
  B23(1,1) = -stiffness[7 + 7 * 9];

  MATRIX eigenvectors3x3(3,3);
  VECTOR eigenvalues3x3(3);
  MATRIX eigenvectors2x2(2,2);
  VECTOR eigenvalues2x2(2);
 
  // clamp and rebuild A
  A.eigensystem(eigenvalues3x3, eigenvectors3x3);
  for (int x = 0; x < 3; x++)
    if (eigenvalues3x3(x) < 0.0f) eigenvalues3x3(x) = 0.0f;
  A = eigenvectors3x3 * MATRIX(eigenvalues3x3) * eigenvectors3x3.transpose();
 
  // clamp and rebuild B12
  B12.eigensystem2x2(eigenvalues2x2, eigenvectors2x2);
  if (eigenvalues2x2(0) < 0.0) eigenvalues2x2(0) = 0.0;
  if (eigenvalues2x2(1) < 0.0) eigenvalues2x2(1) = 0.0;
  B12 = eigenvectors2x2 * MATRIX(eigenvalues2x2) * eigenvectors2x2.transpose();

  // clamp and rebuild B13
  B13.eigensystem2x2(eigenvalues2x2, eigenvectors2x2);
  if (eigenvalues2x2(0) < 0.0) eigenvalues2x2(0) = 0.0;
  if (eigenvalues2x2(1) < 0.0) eigenvalues2x2(1) = 0.0;
  B13 = eigenvectors2x2 * MATRIX(eigenvalues2x2) * eigenvectors2x2.transpose();

  // clamp and rebuild B23
  B23.eigensystem2x2(eigenvalues2x2, eigenvectors2x2);
  if (eigenvalues2x2(0) < 0.0) eigenvalues2x2(0) = 0.0;
  if (eigenvalues2x2(1) < 0.0) eigenvalues2x2(1) = 0.0;
  B23 = eigenvectors2x2 * MATRIX(eigenvalues2x2) * eigenvectors2x2.transpose();
  
  // scatter
  stiffness[0] = -A(0,0);
  stiffness[4] = -A(0,1);
  stiffness[8] = -A(0,2);
  stiffness[0 + 4 * 9] = -A(1,0);
  stiffness[4 + 4 * 9] = -A(1,1);
  stiffness[8 + 4 * 9] = -A(1,2);
  stiffness[0 + 8 * 9] = -A(2,0);
  stiffness[4 + 8 * 9] = -A(2,1);
  stiffness[8 + 8 * 9] = -A(2,2);
  stiffness[1 + 1 * 9] = -B12(0,0);
  stiffness[3 + 1 * 9] = -B12(0,1);
  stiffness[1 + 3 * 9] = -B12(1,0);
  stiffness[3 + 3 * 9] = -B12(1,1);
  stiffness[2 + 2 * 9] = -B13(0,0);
  stiffness[6 + 2 * 9] = -B13(0,1);
  stiffness[2 + 6 * 9] = -B13(1,0);
  stiffness[6 + 6 * 9] = -B13(1,1);
  stiffness[5 + 5 * 9] = -B23(0,0);
  stiffness[7 + 5 * 9] = -B23(0,1);
  stiffness[5 + 7 * 9] = -B23(1,0);
  stiffness[7 + 7 * 9] = -B23(1,1);
  */
}

//////////////////////////////////////////////////////////////////////
// force density with filtered F
//
// The copies are nasty, and are due to the fact that MATRIX3 does not
// store data in a flat array. Unless this becomes a bottleneck though,
// optimization is premature.
//////////////////////////////////////////////////////////////////////
void INVERTIBLE::forceDensity(const Real* F, Real* forces, bool diagonal)
{
#if 0
  // gold standard computation
  // do the pass-through
  _material->forceDensity(F, forces);

  // copy forces into a 3x3 matrix
  MATRIX3 forces3;
  forces3(0,0) = forces[0]; forces3(0,1) = forces[1]; forces3(0,2) = forces[2];
  forces3(1,0) = forces[3]; forces3(1,1) = forces[4]; forces3(1,2) = forces[5];
  forces3(2,0) = forces[6]; forces3(2,1) = forces[7]; forces3(2,2) = forces[8];
  cout << " Correct forceDensity: " << forces3 << endl;

  // try computing it through the stiffness density
  Real stiffnessData[81];
  for (int x = 0; x < 81; x++)
    stiffnessData[x] = 0.0;
  _material->stiffnessDensity(F, stiffnessData);

  VECTOR flatF(9);
  /*
  for (int x = 0; x < 9; x++)
    flatF[x] = F[x];
    */
  // multiply by a "circulated" F
  flatF[0] = F[1]; flatF[1] = F[2]; flatF[2] = F[0];
  flatF[3] = F[4]; flatF[4] = F[5]; flatF[5] = F[3];
  flatF[6] = F[7]; flatF[7] = F[8]; flatF[8] = F[6];

  flatF[0] = F[3]; flatF[1] = F[4]; flatF[2] = F[5];
  flatF[3] = F[6]; flatF[4] = F[7]; flatF[5] = F[8];
  flatF[6] = F[0]; flatF[7] = F[1]; flatF[8] = F[2];

  MATRIX stiffnessMatrix(9, 9, stiffnessData);
  VECTOR forceDensity = stiffnessMatrix * flatF;
  /*
  forces3(0,0) = forceDensity[0]; forces3(0,1) = forceDensity[1]; forces3(0,2) = forceDensity[2];
  forces3(1,0) = forceDensity[3]; forces3(1,1) = forceDensity[4]; forces3(1,2) = forceDensity[5];
  forces3(2,0) = forceDensity[6]; forces3(2,1) = forceDensity[7]; forces3(2,2) = forceDensity[8];
  */
  forces3(0,0) = forceDensity[2]; forces3(0,1) = forceDensity[0]; forces3(0,2) = forceDensity[1];
  forces3(1,0) = forceDensity[5]; forces3(1,1) = forceDensity[3]; forces3(1,2) = forceDensity[4];
  forces3(2,0) = forceDensity[8]; forces3(2,1) = forceDensity[6]; forces3(2,2) = forceDensity[7];

  forces3(0,0) = forceDensity[6]; forces3(0,1) = forceDensity[7]; forces3(0,2) = forceDensity[8];
  forces3(1,0) = forceDensity[0]; forces3(1,1) = forceDensity[1]; forces3(1,2) = forceDensity[2];
  forces3(2,0) = forceDensity[3]; forces3(2,1) = forceDensity[4]; forces3(2,2) = forceDensity[5];
  cout << " flatF: " << flatF << endl;
  cout << " stiffnessDensity: " << stiffnessMatrix << endl;
  cout << " forceDensity via stiffnessDensity: " << forces3 << endl;
#endif

#if 0
  MATRIX3 F3;
  F3(0,0) = F[0]; F3(0,1) = F[1]; F3(0,2) = F[2];
  F3(1,0) = F[3]; F3(1,1) = F[4]; F3(1,2) = F[5];
  F3(2,0) = F[6]; F3(2,1) = F[7]; F3(2,2) = F[8];

  // diagonalize F
  MATRIX3 U;
  MATRIX3 Fhat3;
  MATRIX3 V;
  svd(F3, U, Fhat3, V);

  // filter Fhat here because this is called by cubature training,
  // which does not rely on the firstPK, and force needs to be limited
  // somewhere.
  for (int x = 0; x < 3; x++)
  {
    if (Fhat3(x,x) < _epsilon)
    {
      Fhat3(x,x) = _epsilon;
      x = 4;
    }
  }

  // copy Fhat to a flat array
  Real Fhat[9];
  Fhat[0] = Fhat3(0,0); Fhat[1] = Fhat3(0,1); Fhat[2] = Fhat3(0,2);
  Fhat[3] = Fhat3(1,0); Fhat[4] = Fhat3(1,1); Fhat[5] = Fhat3(1,2);
  Fhat[6] = Fhat3(2,0); Fhat[7] = Fhat3(2,1); Fhat[8] = Fhat3(2,2);
 
  // do the pass-through
  _material->forceDensity(Fhat, forces);

  // copy forces into a 3x3 matrix
  MATRIX3 forces3;
  forces3(0,0) = forces[0]; forces3(0,1) = forces[1]; forces3(0,2) = forces[2];
  forces3(1,0) = forces[3]; forces3(1,1) = forces[4]; forces3(1,2) = forces[5];
  forces3(2,0) = forces[6]; forces3(2,1) = forces[7]; forces3(2,2) = forces[8];

  // apply the SVD rotations
  MATRIX3 final = U * forces3 * V.transpose();

  // copy the final result back to the flat array
  forces[0] = final(0,0); forces[1] = final(0,1); forces[2] = final(0,2);
  forces[3] = final(1,0); forces[4] = final(1,1); forces[5] = final(1,2);
  forces[6] = final(2,0); forces[7] = final(2,1); forces[8] = final(2,2);
#endif

  MATRIX3 packedF = TET::repackF(F);
  MATRIX3 P = firstPiolaKirchhoff(packedF);
  P *= -1.0;

  forces[0] = P(0,0);
  forces[1] = P(1,0);
  forces[2] = P(2,0);
  forces[3] = P(0,1);
  forces[4] = P(1,1);
  forces[5] = P(2,1);
  forces[6] = P(0,2);
  forces[7] = P(1,2);
  forces[8] = P(2,2);
  //_material->forceDensity(F, forces);
}

//////////////////////////////////////////////////////////////////////////////
// Correct both Fhat and U if U contains a reflection
//////////////////////////////////////////////////////////////////////////////
void INVERTIBLE::removeUReflection(MATRIX3& U, MATRIX3& Fhat)
{
  // find the smallest singular value
  int smallest = (Fhat(0,0) < Fhat(1,1)) ? 0 : 1;
  smallest = (Fhat(smallest, smallest) < Fhat(2,2)) ? smallest : 2;
 
  // negate it, and the corresponding column in U
  Fhat(smallest, smallest) *= -1.0;
  U(0, smallest) *= -1.0;
  U(1, smallest) *= -1.0;
  U(2, smallest) *= -1.0;
}

//////////////////////////////////////////////////////////////////////////////
// Orthogonalize U if one of the singular values is near zero
//////////////////////////////////////////////////////////////////////////////
void INVERTIBLE::orthogonalizeU(MATRIX3& U, MATRIX3& Fhat)
{
  // record the sign of U (ie whether it is a reflection) so that
  // the sign can be preserved.
  //
  // If it contains a reflection, it will be corrected in
  // removeUReflection, so make sure we don't double-handle it here.
  Real signU = (det(U) >= 0.0) ? 1.0 : -1.0;
  
  Real smallEigenvalue = 1e-4;

  // see if an eigenvalue is really small
  bool invalid[] = {false, false, false};
  int totalInvalid = 0;
  for (int x = 0; x < 3; x++)
    if (Fhat(x,x) < smallEigenvalue)
    {
      invalid[x] = true;
      totalInvalid++;
    }

  // if none are bad, do nothing
  if (totalInvalid == 0) return;

  // correct U according to the number of small eigenvalues
  if (totalInvalid == 1)
  {
    // figure out which column is bad
    int whichInvalid = 0;
    if (invalid[1]) whichInvalid = 1;
    if (invalid[2]) whichInvalid = 2;
    
    // get the other two columns
    VEC3F col0 = U.col((whichInvalid + 1) % 3);
    VEC3F col1 = U.col((whichInvalid + 2) % 3);

    // compute an orthogonal vector
    VEC3F col2 = col0 ^ col1;
    col2.normalize();

    // store it in U
    U(0, whichInvalid) = col2[0];
    U(1, whichInvalid) = col2[1];
    U(2, whichInvalid) = col2[2];

    // see if we changed any reflections
    if (signU * det(U) < 0.0)
    {
      U(0, whichInvalid) = -col2[0];
      U(1, whichInvalid) = -col2[1];
      U(2, whichInvalid) = -col2[2];
    }
    
    return;
  }
  else if (totalInvalid == 2)
  {
    // find the good vector
    int whichValid = 0;
    if (!invalid[1]) whichValid = 1;
    if (!invalid[2]) whichValid = 2;
    VEC3F validVec = U.col(whichValid);
   
    // find the smallest entry in the good vector
    int smallest = (fabs(validVec[0]) < fabs(validVec[1])) ? 0 : 1;
    smallest = (fabs(validVec[2]) < fabs(validVec[smallest])) ? 2 : smallest;
    
    // compute something arbitrarily orthogonal to
    // the good vector
    VEC3F secondValid;
    int next = (smallest + 1) % 3;
    int nextnext = (smallest + 2) % 3;
    secondValid[smallest] = 0.0;
    secondValid[next] = -validVec[nextnext];
    secondValid[nextnext] = validVec[next];

    // compute a third good vector
    VEC3F thirdValid = validVec ^ secondValid;
    thirdValid.normalize();

    // copy the good vectors into U
    next = (whichValid + 1) % 3;
    nextnext = (whichValid + 2) % 3;
    U(0, next) = secondValid[0];
    U(1, next) = secondValid[1];
    U(2, next) = secondValid[2];
    U(0, nextnext) = thirdValid[0];
    U(1, nextnext) = thirdValid[1];
    U(2, nextnext) = thirdValid[2];

    // see if we changed any reflections
    if (signU * det(U) < 0.0)
    {
      // negate the one with the smallest F, which is the one
      // removeUReflection will be negating again
      if (Fhat(next,next) < Fhat(nextnext, nextnext))
      {
        U(0, next) = -secondValid[0];
        U(1, next) = -secondValid[1];
        U(2, next) = -secondValid[2];
      }
      else
      {
        U(0, nextnext) = -thirdValid[0];
        U(1, nextnext) = -thirdValid[1];
        U(2, nextnext) = -thirdValid[2];
      }
    }
    return;
  }

  // all three are bad, just use identity
  U = MATRIX3::I();
}

//////////////////////////////////////////////////////////////////////////////
// Get the SVD of F
//////////////////////////////////////////////////////////////////////////////
void INVERTIBLE::svd(MATRIX3& F, MATRIX3& U, MATRIX3& Fhat, MATRIX3& V)
{
  // compute the SVD
  MATRIX3 Fnormal3 = F.transpose() * F;
  MATRIX Fnormal(Fnormal3);

  // There is a bug in MKL that if you have an identity matrix with
  // near-zero off-diagonal entries, the eigensolver returns total garbage.
  //
  // Very frustrating!!!!
  Real* data = Fnormal.data();
  for (int x = 0; x < Fnormal.cols() * Fnormal.rows(); x++)
    if (fabs(data[x]) < 1e-12)
      data[x] = 0.0;
 
  VECTOR eigenvalues(3);
  MATRIX eigenvectors(3,3);
#if USING_MKL
  Fnormal.eigensystem(eigenvalues, eigenvectors);
#elif USING_OSX
  Fnormal.eigensystem(eigenvalues, eigenvectors);
#else
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " WARNING: USING A LOW PRECISION EIGENSOLVER" << endl;
  Fnormal.eigensystem3x3(eigenvalues, eigenvectors);
#endif

  // populate the diagonal
  Fhat.clear();
  Fhat(0,0) = sqrt(eigenvalues(0));
  Fhat(1,1) = sqrt(eigenvalues(1));
  Fhat(2,2) = sqrt(eigenvalues(2));

#ifdef _WIN32  
  if (_isnan(Fhat(0,0))) Fhat(0,0) = 0.0;
  if (_isnan(Fhat(1,1))) Fhat(1,1) = 0.0;
  if (_isnan(Fhat(2,2))) Fhat(2,2) = 0.0;
#else
  if (isnan(Fhat(0,0))) Fhat(0,0) = 0.0;
  if (isnan(Fhat(1,1))) Fhat(1,1) = 0.0;
  if (isnan(Fhat(2,2))) Fhat(2,2) = 0.0;
#endif

  // populate V
  eigenvectors.copiesInto(V);
  
  // if V is a reflection, multiply a column by -1
  // to ensure it is a rotation
  Real detV = det(V);
  if (detV < 0.0)
  {
    V(0,0) *= -1.0;
    V(1,0) *= -1.0;
    V(2,0) *= -1.0;
  }

  // compute U
  U.clear();
  U(0,0) = (Fhat(0,0) > 0.0f) ? 1.0f / Fhat(0,0) : 0.0f;
  U(1,1) = (Fhat(1,1) > 0.0f) ? 1.0f / Fhat(1,1) : 0.0f;
  U(2,2) = (Fhat(2,2) > 0.0f) ? 1.0f / Fhat(2,2) : 0.0f;
  U = F * V * U;
  orthogonalizeU(U, Fhat);
 
  // correct any reflections in U to ensure it is a rotation
  if (det(F) < 0.0)
    removeUReflection(U, Fhat);

  /*
  ///////////////////////////////////////////////////////////////
  // All error checking starting here. Fast version would macro
  // this out
  ///////////////////////////////////////////////////////////////
  
  // basic error checking
  if (fabs(det(U) - 1.0) > 1e-4 || 
      fabs(det(V) - 1.0) > 1e-4)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " U and V of SVD have drifted from pure rotation!" << endl;
    cout << " det(U): " << det(U) << endl;
    cout << " det(V): " << det(V) << endl;
  }
  
  // sanity check that this decomposition does hand back F when combined
  MATRIX3 sanity = U * Fhat * V.transpose();
  for (int x = 0; x < 3; x++)
    for (int y = 0; y < 3; y++)
    {
      Real diff = sanity(x,y) - F(x,y);
      if (F(x,y) > 0.0)
        diff = diff / F(x,y);
      if (diff > 1e-4)
      {
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
        cout << "SVD and original matrix differ!" << endl;
        cout << " relative diff: " << diff << " sanity: " << sanity(x,y) << " original: " << F(x,y) << endl;
        cout << " SVD product: " << sanity << endl;
        cout << " original: " << F << endl;
        cout << " normal F: " << Fnormal3 << endl;
        cout << " eigenvalues: " << eigenvalues << endl;
        cout << " eigenvectors: " << eigenvectors << endl;
        cout << " eigenproduct: " << eigenvectors * MATRIX(eigenvalues) * eigenvectors.transpose();
        x = 3;
        y = 3;
      }
    }
    */
}

//////////////////////////////////////////////////////////////////////////////
// Perform midpoint method to find the root of the material
//////////////////////////////////////////////////////////////////////////////
Real INVERTIBLE::findRoot(Real find)
{
  MATRIX3 F = MATRIX3::I();
  MATRIX3 PK;
  Real leftEval, rightEval, middleEval;
  Real tolerance = 1e-6;
  Real left = tolerance;
  Real right = 1.0 + tolerance;
  Real middle = (left + right) * 0.5;

  // Reasonable for StVK
  //Real find = -0.1;

  F(0,0) = left;
  PK = _material->firstPiolaKirchhoff(F);
  leftEval = PK(0,0) - find;

  F(0,0) = right;
  PK = _material->firstPiolaKirchhoff(F);
  rightEval = PK(0,0) - find;

  F(0,0) = middle;
  PK = _material->firstPiolaKirchhoff(F);
  middleEval = PK(0,0) - find;

  int iterations = 0;
  while (fabs(middleEval) > 1e-6 && iterations < 1000)
  {
    if (middleEval * leftEval < 0.0)
    {
      right = middle;
      rightEval = middleEval;
    }
    else
    {
      left = middle;
      leftEval = middleEval;
    }
    middle = (right + left) * 0.5;
    F(0,0) = middle;
    PK = _material->firstPiolaKirchhoff(F);
    middleEval = PK(0,0) - find;

    iterations++;
  }
  if (iterations == 1000)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " No root found! You need to set the epsilon by hand. " << endl;
  }

  _slope = (rightEval - middleEval) / (right - middle);
  _intercept = middle - _slope * middle;

  return middle;
}

//////////////////////////////////////////////////////////////////////////////
// Use stress estimation to find a good clamping point
//////////////////////////////////////////////////////////////////////////////
void INVERTIBLE::findClamping(bool scanRightLeft)
{
  int samplePoints = 100000;
  //int samplePoints = 1000;
  VECTOR xSamples(samplePoints);
  VECTOR ySamples(samplePoints);
  VECTOR stiffnessSamples(samplePoints);
  MATRIX3 F = MATRIX3::I();
  
  // probe along [1,10]
  //Real rightLimit = (scanRightLeft) ? 2.0 : 10.0;
  Real rightLimit = 2.0;
  for (int x = 0; x < xSamples.size(); x++)
    xSamples[x] = 1.0 + (rightLimit - 1.0) * x / (xSamples.size() - 1);

  Real maxFound = 0.0;
  Real xFound = 0.0;
  for (int x = 0; x < xSamples.size(); x++)
  {
    F(0,0) = xSamples(x);
    VECTOR flatF = TET::flattenF(F);
    Real stiffnessData[81];
    _material->stiffnessDensity(flatF.data(), stiffnessData);
    stiffnessSamples[x] = -stiffnessData[0];

    if (x == 0 || stiffnessSamples[x] > maxFound)
    {
      maxFound = stiffnessSamples[x];
      xFound = xSamples(x);
    }
  }

  // find the max slope to search for
  //Real probeFor = 1000.0 * maxFound;
  //Real probeFor = 100.0 * maxFound;
  // OLD? -- works much better
  Real probeFor = 10.0 * maxFound;
  //Real probeFor = 5.0 * maxFound;
  // NEW? -- makes the Newton problem much harder
  //Real probeFor = 2.0 * maxFound;
  cout << " Max found: " << maxFound << endl;
  cout << " Probing for: " << probeFor << endl;
  cout << " Found probing value at: " << xFound << endl;

  // probe along [0,1]
  for (int x = 0; x < xSamples.size(); x++)
    if (scanRightLeft)
      xSamples[x] = 1.0 - (Real)x / xSamples.size();
    else
      xSamples[x] = (Real)x / (xSamples.size() - 1);

  Real foundX = -1;
  Real foundY = 0;
  Real foundStiffness = 0;
  for (int x = 1; x < xSamples.size(); x++)
  {
    F(0,0) = xSamples(x);
    VECTOR flatF = TET::flattenF(F);
    Real stiffnessData[81];
    _material->stiffnessDensity(flatF.data(), stiffnessData);
    Real sample = -stiffnessData[0];

    // if the stiffness is positive and it is the first one encountered
    // that is below the needed threshold
    bool probing = (scanRightLeft) ? sample > probeFor : sample < probeFor;
    if (sample > 0 && probing && foundX < 0)
    {
      foundX = xSamples[x];
      MATRIX3 leftSample = _material->firstPiolaKirchhoff(F);
      foundY = leftSample(0,0);
      foundStiffness = sample;
    }
  }

  /*
  // DEBUG
  cout << " xSamples = [";
  for (int x = 1; x < xSamples.size(); x++)
    cout << xSamples[x] << " ";
  cout << "];" << endl;
  cout << " stiffness = [";
  for (int x = 1; x < xSamples.size(); x++)
  {
    F(0,0) = xSamples(x);
    VECTOR flatF = TET::flattenF(F);
    Real stiffnessData[81];
    _material->stiffnessDensity(flatF.data(), stiffnessData);
    Real sample = -stiffnessData[0];
    cout << sample << " ";
  }
  cout << "];" << endl;
  */

  /*
  // DEBUG force inversion to some value
  //foundX = 0.8;  // blowup
  //foundX = 0.7268;  // blowup
  //foundX = 0.71; // blowup
  //foundX = 0.696; // stable
  foundX = 0.694; // stable
  //foundX = 0.68; // stable
  //foundX = 0.6; // blowup
  F(0,0) = foundX;
  VECTOR flatF = TET::flattenF(F);
  Real stiffnessData[81];
  _material->stiffnessDensity(flatF.data(), stiffnessData);
  Real sample = -stiffnessData[0];
  MATRIX3 leftSample = _material->firstPiolaKirchhoff(F);
  foundY = leftSample(0,0);
  foundStiffness = sample;
  */

  cout << " found x: " << foundX << endl;
  cout << " found y: " << foundY << endl;
  cout << " found stiffness: " << foundStiffness << endl;

  // sanity check
  if (foundX < 0.0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " INVERSION CLAMPING FAILED!!!! " << endl;
  }

  // set the filtering consts
  _epsilon = foundX;
  _slope = foundStiffness;
  _intercept = foundY - foundStiffness * foundX;
}
