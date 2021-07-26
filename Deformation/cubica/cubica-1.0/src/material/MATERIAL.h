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
// MATERIAL.h: interface for the MATERIAL class.
//
//////////////////////////////////////////////////////////////////////

#ifndef MATERIAL_H
#define MATERIAL_H

#include <TET.h>
#include <SPARSE_MATRIX.h>
#include <MATRIX.h>
#include <TENSOR3.h>

using namespace std;

// FIXME: I still need to add diagonal versions of the
// incompressibility stuff.

//////////////////////////////////////////////////////////////////////
// Abstract material class
//////////////////////////////////////////////////////////////////////
class MATERIAL {

public:
  MATERIAL()
    : _materialName(""), stiffnessDensityTime( 0.0 )
  {}
  virtual ~MATERIAL() {};

  // Set the parameter diagonal = true in cases when the deformation
  // gradient is explicitly known to be diagonal (cases in which
  // we diagonalize it ourselves).
  virtual MATRIX stiffness(TET& tet, bool diagonal = false) = 0;
  virtual MATRIX3 firstPiolaKirchhoff(MATRIX3& F, bool diagonal = false)
  {
    return F * secondPiolaKirchhoff(F, diagonal);
  }
  virtual MATRIX3 secondPiolaKirchhoff(MATRIX3& F, bool diagonal = false)
  {
    return F.inverse() * firstPiolaKirchhoff(F, diagonal);
  }
  virtual MATRIX3 cauchyStress(MATRIX3& F, bool diagonal = false);
  virtual void stiffnessDensity(const Real* F, Real* stiffness,
                                bool diagonal = false) = 0;

  // this is just a flattened version of first Piola-Kirchoff, and should
  // be renamed at some point
  //virtual void forceDensity(TET& tet, VEC3F* forces) = 0;
  virtual void forceDensity(const Real* F, Real* forces, bool diagonal = false) = 0;

  // compute the invariant derivatives, 
  // \frac{\partial \Psi}{\partial I}
  // \frac{\partial \Psi}{\partial II}
  // \frac{\partial \Psi}{\partial III}
  // These are stacked in that order in the returned VEC3F
  VEC3F invariantPartials(const Real* F);

  void computePFPu(TET& tet, MATRIX& pFpu);

  // return a copy of this material object
  virtual MATERIAL* copy() = 0;
  string& materialName() { return _materialName; };

  // get the hessian tensor
  virtual TENSOR3 hessian(TET& tet);
 
  // get the firstPK corresponding to the penalty
  MATRIX3 penaltyFirstPK(Real kBulkMod, MATRIX3& F) {
    MATRIX3 C = F.transpose() * F;
    return F * penaltySecondPK(kBulkMod, C);
  };

  enum PenaltyType {
    LOG_QUADRATIC = 0,    // W = K log(J)^2
    QUADRATIC,            // W = (J - 1)^2
    POLY_TWO_TERM,
    POLY_THREE_TERM,      // Of the form sum_{i=1}^{N}( (1/Di)(J - 1)^2i )
    POLY_FOUR_TERM,
    LOG                   // W = 2 K log(J)
  };

  MATRIX3 penaltyFirstPK(Real kBulkMod, MATRIX3& F, PenaltyType type) {
    MATRIX3 C = F.transpose() * F;
    return F * penaltySecondPK(kBulkMod, C, type);
  };

  // Copies upper triangular part of a symmetric stiffness
  // density matrix in to its lower triangular part
  void copySymmetricStiffness( Real *stiffness )
  {
    for ( int i = 0; i < 9; i++ )
    for ( int j = i + 1; j < 9; j++ )
    {
      stiffness[ j * 9 + i ] = stiffness[ i * 9 + j ];
    }
  }

protected:
  string _materialName;

  void addPenalty(Real kBulkmod, TET &tet, MATRIX &accum,
                  PenaltyType type)
  {
    switch( type )
    {
      case LOG_QUADRATIC:
      default:
      {
        addPenalty(kBulkmod, tet, accum);
        break;
      }
      case QUADRATIC:
      {
        addPenaltyQuadratic(kBulkmod, tet, accum);
        break;
      }
      case LOG:
      {
        cout << "Error: not implemented" << endl;
        break;
      }
    }
  }

  MATRIX3 penaltySecondPK(Real kBulkmod, MATRIX3 &C, PenaltyType type)
  {
    switch( type )
    {
      case LOG_QUADRATIC:
      default:
      {
        return penaltySecondPK(kBulkmod, C);
      }
      case QUADRATIC:
      {
        return penaltySecondPKQuadratic(kBulkmod, C);
      }
      case LOG:
      {
        cout << "Error: not implemented" << endl;
        break;
      }
    }

    return MATRIX3();
  }

  void penaltyForceDensity(Real kBulkmod, const Real *F, Real *forces,
                           PenaltyType type)
  {
    switch( type )
    {
      case LOG_QUADRATIC:
      default:
      {
        penaltyForceDensity(kBulkmod, F, forces);
        break;
      }
      case QUADRATIC:
      {
        penaltyForceDensityQuadratic(kBulkmod, F, forces);
        break;
      }
      case LOG:
      {
        cout << "Error: not implemented" << endl;
        break;
      }
    }
  }

  void penaltyStiffnessDensity(Real kBulkmod, const Real *F, Real *stiffness,
                               PenaltyType type)
  {
    switch( type )
    {
      case LOG_QUADRATIC:
      default:
      {
        penaltyStiffnessDensity(kBulkmod, F, stiffness);
        break;
      }
      case QUADRATIC:
      {
        penaltyStiffnessDensityQuadratic(kBulkmod, F, stiffness);

        copySymmetricStiffness( stiffness );

        break;
      }
      case LOG:
      {
        cout << "Error: not implemented" << endl;
        break;
      }
    }
  }

  // incompressibility functions for Mooney-Rivlin and Arruda-Boyce.
  // StVK doesn't use it, but since they don't add any member variables,
  // creating a separate INCOMPRESSIBLE_MATERIAL class is overkill
  void addPenalty(Real kBulkmod, TET& tet, MATRIX& accum);
  MATRIX3 penaltySecondPK(Real kBulkmod, MATRIX3& C);
  void penaltyForceDensity(Real kBulkmod, TET& tet, VEC3F* forces);
  void penaltyForceDensity(Real kBulkmod, const Real* F, Real* forces);
  void penaltyStiffnessDensity(Real kBulkmod, const Real* F, Real* stiffness);

  // Quadratic (nonsingular) incompressibility terms
  void addPenaltyQuadratic(Real kBulkMod, TET &tet, MATRIX &accum);
  void addPenaltyQuadratic(Real kBulkMod, const Real *F, TET &tet,
                           MATRIX &accum);
  MATRIX3 penaltySecondPKQuadratic(Real kBulkMod, MATRIX3 &C);
  void penaltyForceDensityQuadratic(Real kBulkMod, const Real *F, Real *forces);
  void penaltyStiffnessDensityQuadratic(Real kBulkMod, const Real *F,
                                        Real *stiffness);

  // Nonsingular polynomial incompressibility terms
  void addPenaltyPolynomial2term(Real D1, Real D2,
                                 TET &tet, MATRIX &accum);
  void addPenaltyPolynomial2term(Real D1, Real D2,
                                 const Real *F, TET &tet, MATRIX &accum);
  MATRIX3 penaltySecondPKPolynomial2term(Real D1, Real D2,
                                         MATRIX3 &C);
  void penaltyForceDensityPolynomial2term(Real D1, Real D2,
                                          const Real *F, Real *forces);
  void penaltyStiffnessDensityPolynomial2term(Real D1, Real D2,
                                              const Real *F, Real *stiffness);

  void addPenaltyPolynomial3term(Real D1, Real D2, Real D3,
                                 TET &tet, MATRIX &accum);
  void addPenaltyPolynomial3term(Real D1, Real D2, Real D3,
                                 const Real *F, TET &tet, MATRIX &accum);
  MATRIX3 penaltySecondPKPolynomial3term(Real D1, Real D2, Real D3,
                                         MATRIX3 &C);
  void penaltyForceDensityPolynomial3term(Real D1, Real D2, Real D3,
                                          const Real *F, Real *forces);
  void penaltyStiffnessDensityPolynomial3term(Real D1, Real D2, Real D3,
                                              const Real *F, Real *stiffness);

  void addPenaltyPolynomial4term(Real D1, Real D2, Real D3, Real D4,
                                 TET &tet, MATRIX &accum);
  void addPenaltyPolynomial4term(Real D1, Real D2, Real D3, Real D4,
                                 const Real *F, TET &tet, MATRIX &accum);
  MATRIX3 penaltySecondPKPolynomial4term(Real D1, Real D2, Real D3, Real D4,
                                         MATRIX3 &C);
  void penaltyForceDensityPolynomial4term(Real D1, Real D2, Real D3, Real D4,
                                          const Real *F, Real *forces);
  void penaltyStiffnessDensityPolynomial4term(Real D1, Real D2, Real D3, Real D4,
                                              const Real *F, Real *stiffness);

  Real stiffnessDensityTime;

public:  
  virtual void hessian0(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian) = 0;
  virtual void hessian1(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian) = 0;
  virtual void hessian2(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian) = 0;
  virtual void hessian3(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian) = 0;
  virtual void hessian4(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian) = 0;
  virtual void hessian5(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian) = 0;
  virtual void hessian6(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian) = 0;
  virtual void hessian7(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian) = 0;
  virtual void hessian8(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian) = 0;
  virtual void hessian9(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian) = 0;
  virtual void hessian10(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian) = 0;
  virtual void hessian11(MATRIX3& DmInv, const VEC3F* b, const VEC3F* vertices, MATRIX& hessian) = 0;

  virtual void hessian0(const Real* F, const VEC3F* b, MATRIX& hessian) = 0;
  virtual void hessian1(const Real* F, const VEC3F* b, MATRIX& hessian) = 0;
  virtual void hessian2(const Real* F, const VEC3F* b, MATRIX& hessian) = 0;
  virtual void hessian3(const Real* F, const VEC3F* b, MATRIX& hessian) = 0;
  virtual void hessian4(const Real* F, const VEC3F* b, MATRIX& hessian) = 0;
  virtual void hessian5(const Real* F, const VEC3F* b, MATRIX& hessian) = 0;
  virtual void hessian6(const Real* F, const VEC3F* b, MATRIX& hessian) = 0;
  virtual void hessian7(const Real* F, const VEC3F* b, MATRIX& hessian) = 0;
  virtual void hessian8(const Real* F, const VEC3F* b, MATRIX& hessian) = 0;

  // incompressibility hessian
  void hessian0(const Real k, const Real* F, const VEC3F* b, MATRIX& hessian);
  void hessian1(const Real k, const Real* F, const VEC3F* b, MATRIX& hessian);
  void hessian2(const Real k, const Real* F, const VEC3F* b, MATRIX& hessian);
  void hessian3(const Real k, const Real* F, const VEC3F* b, MATRIX& hessian);
  void hessian4(const Real k, const Real* F, const VEC3F* b, MATRIX& hessian);
  void hessian5(const Real k, const Real* F, const VEC3F* b, MATRIX& hessian);
  void hessian6(const Real k, const Real* F, const VEC3F* b, MATRIX& hessian);
  void hessian7(const Real k, const Real* F, const VEC3F* b, MATRIX& hessian);
  void hessian8(const Real k, const Real* F, const VEC3F* b, MATRIX& hessian);
};

/*
Much of the code in this class and its subclasses was generated using Matlab and Maple.
For completeness, the generating code will be listed when possible. Some search-and-replace
is necessary to get the final code used in these classes, but the renamings necessary
should be self-evident.

The following are Matlab scripts used by all the material subclasses:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup_symbols_wrt_u.m
% setup symbols for when you want derivatives wrt u.
% This is more expensive than doing wrt just F
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = sym('a'); b = sym('b'); c = sym('c'); d = sym('d');
e = sym('e'); f = sym('f'); g = sym('g'); h = sym('h');
i = sym('i'); j = sym('j'); k = sym('k'); l = sym('l');    
m = sym('m'); n = sym('n'); o = sym('o'); p = sym('p');
q = sym('q'); r = sym('r'); s = sym('s'); t = sym('t');
v = sym('u');    
n0 = [a b c];
n1 = [d e f];
n2 = [g h i];
n3 = [j k l];
u = mytrans([a b c d e f g h i j k l]);
x = mytrans([
    n1-n0
    n2-n0
    n3-n0
    ]);
Xinv = [
    m n o
    p q r
    s t v
    ];
F = x * Xinv;
C = mytrans(F)*F;
I_C = trace(C);
II_C = dblcon(C,C);
E = 0.5 * (C - eye(3));
J = det(F);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup_symbols_wrt_F.m
% setup symbols for when you only wants things wrt just F, not u
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms F0 F1 F2 F3 F4 F5 F6 F7 F8;
F = [F0 F3 F6 
     F1 F4 F7 
     F2 F5 F8];
C = mytrans(F)*F;
I_C = trace(C);
II_C = dblcon(C,C);
E = 0.5 * (C - eye(3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dblcon.m 
% symbolic double contraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = dblcon(A, B)
    [m,n] = size(A);
    s = 0;
    for i = 1:m
        for j = 1:n
            s = s + A(i,j)*B(i,j);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mytrans.m
% symbolic transpose
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N = mytrans(M)
    [m,n] = size(M);
    N = sym('a')*eye(n,m);
    for i = 1:m
        for j = 1:n
            N(j,i) = M(i,j);
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% codegen_density.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [forceDensity, stiffnessDensity] = codegen_density(strain_energy, F)
    Psi = strain_energy;
    F = reshape(F, 9, 1);
    display '-- computing stress force density --';
    Rd_sym = sym('DUMMY_TEMP')*ones(1,9);
    for i = 1:9
        Rd_sym(1,i) = - diff( Psi, F(i));
    end
    display '-- computing stiffness density --';
    Kd_sym = sym('DUMMY_TEMP')*ones(9*9,1);
    for i = 1:9
        for j = 1:9
            Kd_sym((i-1)*9+j) = diff( Rd_sym(1,i), F(j) );
        end
    end
    display '-- codegenning --';
    forceDensity = maple('codegen[C]', Rd_sym, 'optimized');
    stiffnessDensity = maple('codegen[C]', Kd_sym, 'optimized');
end
*/

#endif
