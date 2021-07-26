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
// MOONEY_RIVLIN.cpp: implementation of the MOONEY_RIVLIN class.
//
//////////////////////////////////////////////////////////////////////

#include "MOONEY_RIVLIN.h"

//////////////////////////////////////////////////////////////////////
// Constructor for MOONEY_RIVLIN
//////////////////////////////////////////////////////////////////////
MOONEY_RIVLIN::MOONEY_RIVLIN(Real mu01, Real mu10, Real k) :
  _mu01(mu01), _mu10(mu10), _k(k),
  // work arrays
  _pf_pF(12,9),
  _pF_pu(9,12)
{
  _materialName.assign("MOONEY_RIVLIN");
}

//////////////////////////////////////////////////////////////////////
// make a copy
//////////////////////////////////////////////////////////////////////
MATERIAL* MOONEY_RIVLIN::copy()
{
  MATERIAL* material = new MOONEY_RIVLIN(_mu01, _mu10, _k);
  return material;
}

//////////////////////////////////////////////////////////////////////
// stiffness matrix implementation
//////////////////////////////////////////////////////////////////////
MATRIX MOONEY_RIVLIN::stiffness(TET& tet, bool diagonal)
{
  computeStresses(tet);
  addPenalty(_k, tet, _pf_pF);
  computePFPu(tet, _pF_pu);
  MATRIX product = _pf_pF * _pF_pu;
  return product;
}

//////////////////////////////////////////////////////////////////////
// implementation of second PK stress tensor
//////////////////////////////////////////////////////////////////////
MATRIX3 MOONEY_RIVLIN::secondPiolaKirchhoff(MATRIX3& F, bool diagonal)
{
  MATRIX3 C = F.transpose() * F;
  MATRIX3 I = MATRIX3::I();
  MATRIX3 Sp = penaltySecondPK(_k, C);
  MATRIX3 S = 2*( _mu10*I + _mu01*trace(C)*I - _mu01*C) + Sp;
  return S;
}

//////////////////////////////////////////////////////////////////////
// autogenerated as "forceDensity" from:
//
/*
 setup_symbols_wrt_u;
 syms mu01 mu10;
 strain_energy = mu10*( I_C - 3 ) + mu01/2*( I_C^2 - II_C - 6 );
 Rd_sym = sym('DUMMY_TEMP')*ones(1,12);
 for i = 1:12
    Rd_sym(1,i) = -diff( strain_energy, u(i));
 end
 forceDensity = maple('codegen[C]', Rd_sym, 'optimized')
*/
//
// see "MATERIAL.h" for more details
//////////////////////////////////////////////////////////////////////
/*
void MOONEY_RIVLIN::forceDensity(TET& tet, VEC3F* forces)
{
  for (int x = 0; x < 4; x++)
    forces[x].clear();

  penaltyForceDensity(_k, tet, forces);
  VEC3F** nodes = tet.vertices;
  const MATRIX3 matInv = tet.DmInv();
  const double mu10 = _mu10;
  const double mu01 = _mu01;
  
  const double a = (*nodes[0])[0];
  const double b = (*nodes[0])[1];
  const double c = (*nodes[0])[2];
  const double d = (*nodes[1])[0];
  const double e = (*nodes[1])[1];
  const double f = (*nodes[1])[2];
  const double g = (*nodes[2])[0];
  const double h = (*nodes[2])[1];
  const double i = (*nodes[2])[2];
  const double j = (*nodes[3])[0];
  const double k = (*nodes[3])[1];
  const double l = (*nodes[3])[2];

  const double m = matInv(0,0);
  const double n = matInv(0,1);
  const double o = matInv(0,2);
  const double p = matInv(1,0);
  const double q = matInv(1,1);
  const double r = matInv(1,2);
  const double s = matInv(2,0);
  const double t = matInv(2,1);
  const double u = matInv(2,2);

  const double t1 = d-a;
  const double t3 = g-a;
  const double t5 = j-a;
  const double t7 = t1*m+t3*p+t5*s;
  const double t8 = -m-p-s;
  const double t13 = t1*n+t3*q+t5*t;
  const double t14 = -n-q-t;
  const double t19 = t1*o+t3*r+t5*u;
  const double t20 = -o-r-u;
  const double t22 = t7*t8+t13*t14+t19*t20;
  const double t24 = t7*t7;
  const double t25 = e-b;
  const double t27 = h-b;
  const double t29 = k-b;
  const double t31 = t25*m+t27*p+t29*s;
  const double t32 = t31*t31;
  const double t33 = f-c;
  const double t35 = i-c;
  const double t37 = l-c;
  const double t39 = t33*m+t35*p+t37*s;
  const double t40 = t39*t39;
  const double t41 = t13*t13;
  const double t45 = t25*n+t27*q+t29*t;
  const double t46 = t45*t45;
  const double t50 = t33*n+t35*q+t37*t;
  const double t51 = t50*t50;
  const double t52 = t19*t19;
  const double t56 = t25*o+t27*r+t29*u;
  const double t57 = t56*t56;
  const double t61 = t33*o+t35*r+t37*u;
  const double t62 = t61*t61;
  const double t63 = t24+t32+t40+t41+t46+t51+t52+t57+t62;
  const double t66 = t24+t32+t40;
  const double t67 = t66*t7;
  const double t73 = t7*t13+t31*t45+t39*t50;
  const double t82 = t7*t19+t31*t56+t39*t61;
  const double t88 = t41+t46+t51;
  const double t89 = t88*t13;
  const double t95 = t13*t19+t45*t56+t50*t61;
  const double t101 = t52+t57+t62;
  const double t102 = t101*t19;
  const double t112 = t31*t8+t45*t14+t56*t20;
  const double t116 = t66*t31;
  const double t129 = t88*t45;
  const double t137 = t101*t56;
  const double t147 = t39*t8+t50*t14+t61*t20;
  const double t151 = t66*t39;
  const double t164 = t88*t50;
  const double t172 = t101*t61;
  const double t182 = t7*m+t13*n+t19*o;
  const double t214 = t31*m+t45*n+t56*o;
  const double t246 = t39*m+t50*n+t61*o;
  const double t278 = t7*p+t13*q+t19*r;
  const double t310 = t31*p+t45*q+t56*r;
  const double t342 = t39*p+t50*q+t61*r;
  const double t374 = t7*s+t13*t+t19*u;
  const double t406 = t31*s+t45*t+t56*u;
  const double t438 = t39*s+t50*t+t61*u;
  forces[0][0] += -2.0*mu10*t22-mu01*(4.0*t63*t22-4.0*t67*t8-4.0*t73*(t8*t13+t7*t14)-4.0*t82*(t8*t19+t7*t20)-4.0*t89*t14-4.0*t95*(t14*t19+t13*t20)-4.0*t102*t20)/2.0;
  forces[0][1] += -2.0*mu10*t112-mu01*(4.0*t63*t112-4.0*t116*t8-4.0*t73*(t8*t45+t31*t14)-4.0*t82*(t8*t56+t31*t20)-4.0*t129*t14-4.0*t95*(t14*t56+t45*t20)-4.0*t137*t20)/2.0;
  forces[0][2] += -2.0*mu10*t147-mu01*(4.0*t63*t147-4.0*t151*t8-4.0*t73*(t8*t50+t39*t14)-4.0*t82*(t8*t61+t39*t20)-4.0*t164*t14-4.0*t95*(t14*t61+t50*t20)-4.0*t172*t20)/2.0;
  forces[1][0] += -2.0*mu10*t182-mu01*(4.0*t63*t182-4.0*t67*m-4.0*t73*(m*t13+t7*n)-4.0*t82*(m*t19+t7*o)-4.0*t89*n-4.0*t95*(n*t19+t13*o)-4.0*t102*o)/2.0;
  forces[1][1] += -2.0*mu10*t214-mu01*(4.0*t63*t214-4.0*t116*m-4.0*t73*(m*t45+t31*n)-4.0*t82*(m*t56+t31*o)-4.0*t129*n-4.0*t95*(n*t56+t45*o)-4.0*t137*o)/2.0;
  forces[1][2] += -2.0*mu10*t246-mu01*(4.0*t63*t246-4.0*t151*m-4.0*t73*(m*t50+t39*n)-4.0*t82*(m*t61+t39*o)-4.0*t164*n-4.0*t95*(n*t61+t50*o)-4.0*t172*o)/2.0;
  forces[2][0] += -2.0*mu10*t278-mu01*(4.0*t63*t278-4.0*t67*p-4.0*t73*(p*t13+t7*q)-4.0*t82*(p*t19+t7*r)-4.0*t89*q-4.0*t95*(q*t19+t13*r)-4.0*t102*r)/2.0;
  forces[2][1] += -2.0*mu10*t310-mu01*(4.0*t63*t310-4.0*t116*p-4.0*t73*(p*t45+t31*q)-4.0*t82*(p*t56+t31*r)-4.0*t129*q-4.0*t95*(q*t56+t45*r)-4.0*t137*r)/2.0;
  forces[2][2] += -2.0*mu10*t342-mu01*(4.0*t63*t342-4.0*t151*p-4.0*t73*(p*t50+t39*q)-4.0*t82*(p*t61+t39*r)-4.0*t164*q-4.0*t95*(q*t61+t50*r)-4.0*t172*r)/2.0;
  forces[3][0] += -2.0*mu10*t374-mu01*(4.0*t63*t374-4.0*t67*s-4.0*t73*(s*t13+t7*t)-4.0*t82*(s*t19+t7*u)-4.0*t89*t-4.0*t95*(t*t19+t13*u)-4.0*t102*u)/2.0;
  forces[3][1] += -2.0*mu10*t406-mu01*(4.0*t63*t406-4.0*t116*s-4.0*t73*(s*t45+t31*t)-4.0*t82*(s*t56+t31*u)-4.0*t129*t-4.0*t95*(t*t56+t45*u)-4.0*t137*u)/2.0;
  forces[3][2] += -2.0*mu10*t438-mu01*(4.0*t63*t438-4.0*t151*s-4.0*t73*(s*t50+t39*t)-4.0*t82*(s*t61+t39*u)-4.0*t164*t-4.0*t95*(t*t61+t50*u)-4.0*t172*u)/2.0;
}
*/

//////////////////////////////////////////////////////////////////////
// autogenerated as "stiffnessDensity" from:
//
/*
 setup_symbols_wrt_F;
 syms mu01 mu10;
 strain_energy = mu10*( I_C - 3 ) + mu01/2*( I_C^2 - II_C - 6 );
 [forceDensity stiffnessDensity] = codegen_density(strain_energy,F)
*/
//
// see "MATERIAL.h" for more details
//////////////////////////////////////////////////////////////////////
void MOONEY_RIVLIN::stiffnessDensity(const Real* F, Real* stiffness,
                                     bool diagonal)
{
  const double mu10 = _mu10;
  const double mu01 = _mu01;

  const double t1 = 2.0*mu10;
  const double t2 = F[4]*F[4];
  const double t3 = F[5]*F[5];
  const double t4 = F[7]*F[7];
  const double t5 = F[8]*F[8];
  const double t10 = F[4]*F[3];
  const double t11 = F[7]*F[6];
  const double t14 = 2.0*mu01*(-t10-t11);
  const double t15 = F[5]*F[3];
  const double t16 = F[8]*F[6];
  const double t19 = 2.0*mu01*(-t15-t16);
  const double t20 = F[1]*F[4];
  const double t21 = F[2]*F[5];
  const double t24 = 2.0*mu01*(-t20-t21);
  const double t25 = F[0]*F[4];
  const double t27 = F[1]*F[3];
  const double t31 = mu01*(8.0*t25-4.0*t27)/2.0;
  const double t32 = F[0]*F[5];
  const double t34 = F[2]*F[3];
  const double t38 = mu01*(8.0*t32-4.0*t34)/2.0;
  const double t39 = F[1]*F[7];
  const double t40 = F[2]*F[8];
  const double t43 = 2.0*mu01*(-t39-t40);
  const double t44 = F[0]*F[7];
  const double t46 = F[1]*F[6];
  const double t50 = mu01*(8.0*t44-4.0*t46)/2.0;
  const double t51 = F[0]*F[8];
  const double t53 = F[2]*F[6];
  const double t57 = mu01*(8.0*t51-4.0*t53)/2.0;
  const double t58 = F[3]*F[3];
  const double t59 = F[6]*F[6];
  const double t64 = F[5]*F[4];
  const double t65 = F[8]*F[7];
  const double t68 = 2.0*mu01*(-t64-t65);
  const double t73 = mu01*(8.0*t27-4.0*t25)/2.0;
  const double t74 = F[0]*F[3];
  const double t77 = 2.0*mu01*(-t74-t21);
  const double t78 = F[1]*F[5];
  const double t80 = F[2]*F[4];
  const double t84 = mu01*(8.0*t78-4.0*t80)/2.0;
  const double t89 = mu01*(8.0*t46-4.0*t44)/2.0;
  const double t90 = F[0]*F[6];
  const double t93 = 2.0*mu01*(-t90-t40);
  const double t94 = F[1]*F[8];
  const double t96 = F[2]*F[7];
  const double t100 = mu01*(8.0*t94-4.0*t96)/2.0;
  const double t109 = mu01*(8.0*t34-4.0*t32)/2.0;
  const double t114 = mu01*(8.0*t80-4.0*t78)/2.0;
  const double t117 = 2.0*mu01*(-t74-t20);
  const double t122 = mu01*(8.0*t53-4.0*t51)/2.0;
  const double t127 = mu01*(8.0*t96-4.0*t94)/2.0;
  const double t130 = 2.0*mu01*(-t90-t39);
  const double t131 = F[1]*F[1];
  const double t132 = F[2]*F[2];
  const double t137 = F[1]*F[0];
  const double t140 = 2.0*mu01*(-t137-t11);
  const double t141 = F[2]*F[0];
  const double t144 = 2.0*mu01*(-t141-t16);
  const double t145 = F[4]*F[7];
  const double t146 = F[5]*F[8];
  const double t149 = 2.0*mu01*(-t145-t146);
  const double t150 = F[3]*F[7];
  const double t152 = F[6]*F[4];
  const double t156 = mu01*(8.0*t150-4.0*t152)/2.0;
  const double t157 = F[3]*F[8];
  const double t159 = F[6]*F[5];
  const double t163 = mu01*(8.0*t157-4.0*t159)/2.0;
  const double t164 = F[0]*F[0];
  const double t169 = F[2]*F[1];
  const double t172 = 2.0*mu01*(-t169-t65);
  const double t177 = mu01*(8.0*t152-4.0*t150)/2.0;
  const double t178 = F[3]*F[6];
  const double t181 = 2.0*mu01*(-t178-t146);
  const double t182 = F[4]*F[8];
  const double t184 = F[7]*F[5];
  const double t188 = mu01*(8.0*t182-4.0*t184)/2.0;
  const double t197 = mu01*(8.0*t159-4.0*t157)/2.0;
  const double t202 = mu01*(8.0*t184-4.0*t182)/2.0;
  const double t205 = 2.0*mu01*(-t178-t145);
  const double t212 = 2.0*mu01*(-t137-t10);
  const double t215 = 2.0*mu01*(-t141-t15);
  const double t222 = 2.0*mu01*(-t169-t64);
  stiffness[0] = -t1-2.0*mu01*(t2+t3+t4+t5);
  stiffness[1] = -t14;
  stiffness[2] = -t19;
  stiffness[3] = -t24;
  stiffness[4] = -t31;
  stiffness[5] = -t38;
  stiffness[6] = -t43;
  stiffness[7] = -t50;
  stiffness[8] = -t57;
  stiffness[9] = -t14;
  stiffness[10] = -t1-2.0*mu01*(t58+t3+t59+t5);
  stiffness[11] = -t68;
  stiffness[12] = -t73;
  stiffness[13] = -t77;
  stiffness[14] = -t84;
  stiffness[15] = -t89;
  stiffness[16] = -t93;
  stiffness[17] = -t100;
  stiffness[18] = -t19;
  stiffness[19] = -t68;
  stiffness[20] = -t1-2.0*mu01*(t58+t2+t59+t4);
  stiffness[21] = -t109;
  stiffness[22] = -t114;
  stiffness[23] = -t117;
  stiffness[24] = -t122;
  stiffness[25] = -t127;
  stiffness[26] = -t130;
  stiffness[27] = -t24;
  stiffness[28] = -t73;
  stiffness[29] = -t109;
  stiffness[30] = -t1-2.0*mu01*(t131+t132+t4+t5);
  stiffness[31] = -t140;
  stiffness[32] = -t144;
  stiffness[33] = -t149;
  stiffness[34] = -t156;
  stiffness[35] = -t163;
  stiffness[36] = -t31;
  stiffness[37] = -t77;
  stiffness[38] = -t114;
  stiffness[39] = -t140;
  stiffness[40] = -t1-2.0*mu01*(t164+t132+t59+t5);
  stiffness[41] = -t172;
  stiffness[42] = -t177;
  stiffness[43] = -t181;
  stiffness[44] = -t188;
  stiffness[45] = -t38;
  stiffness[46] = -t84;
  stiffness[47] = -t117;
  stiffness[48] = -t144;
  stiffness[49] = -t172;
  stiffness[50] = -t1-2.0*mu01*(t164+t131+t59+t4);
  stiffness[51] = -t197;
  stiffness[52] = -t202;
  stiffness[53] = -t205;
  stiffness[54] = -t43;
  stiffness[55] = -t89;
  stiffness[56] = -t122;
  stiffness[57] = -t149;
  stiffness[58] = -t177;
  stiffness[59] = -t197;
  stiffness[60] = -t1-2.0*mu01*(t131+t132+t2+t3);
  stiffness[61] = -t212;
  stiffness[62] = -t215;
  stiffness[63] = -t50;
  stiffness[64] = -t93;
  stiffness[65] = -t127;
  stiffness[66] = -t156;
  stiffness[67] = -t181;
  stiffness[68] = -t202;
  stiffness[69] = -t212;
  stiffness[70] = -t1-2.0*mu01*(t164+t132+t58+t3);
  stiffness[71] = -t222;
  stiffness[72] = -t57;
  stiffness[73] = -t100;
  stiffness[74] = -t130;
  stiffness[75] = -t163;
  stiffness[76] = -t188;
  stiffness[77] = -t205;
  stiffness[78] = -t215;
  stiffness[79] = -t222;
  stiffness[80] = -t1-2.0*mu01*(t164+t131+t58+t2);

  penaltyStiffnessDensity(_k, F, stiffness);
}

//////////////////////////////////////////////////////////////////////
// autogenerated as "forceDensity" from:
//
/*
 setup_symbols_wrt_F;
 syms mu01 mu10;
 strain_energy = mu10*( I_C - 3 ) + mu01/2*( I_C^2 - II_C - 6 );
 [forceDensity stiffnessDensity] = codegen_density(strain_energy,F)
*/
//
// see "MATERIAL.h" for more details
//////////////////////////////////////////////////////////////////////
void MOONEY_RIVLIN::forceDensity(const Real* F, Real* forces, bool diagonal)
{
  const double mu10 = _mu10;
  const double mu01 = _mu01;

  const double t3 = F[0]*F[0];
  const double t4 = F[1]*F[1];
  const double t5 = F[2]*F[2];
  const double t6 = F[3]*F[3];
  const double t7 = F[4]*F[4];
  const double t8 = F[5]*F[5];
  const double t9 = F[6]*F[6];
  const double t10 = F[7]*F[7];
  const double t11 = F[8]*F[8];
  const double t12 = t3+t4+t5+t6+t7+t8+t9+t10+t11;
  const double t14 = t3+t4+t5;
  const double t19 = F[0]*F[3]+F[1]*F[4]+F[2]*F[5];
  const double t24 = F[0]*F[6]+F[1]*F[7]+F[2]*F[8];
  const double t54 = t6+t7+t8;
  const double t59 = F[3]*F[6]+F[4]*F[7]+F[5]*F[8];
  const double t90 = t9+t10+t11;

  forces[0] = -2.0*mu10*F[0]-2.0*mu01*(t12*F[0]-t14*F[0]-t19*F[3]-t24*F[6]);
  forces[1] = -2.0*mu10*F[1]-2.0*mu01*(t12*F[1]-t14*F[1]-t19*F[4]-t24*F[7]);
  forces[2] = -2.0*mu10*F[2]-2.0*mu01*(t12*F[2]-t14*F[2]-t19*F[5]-t24*F[8]);
  forces[3] = -2.0*mu10*F[3]-2.0*mu01*(t12*F[3]-t19*F[0]-t54*F[3]-t59*F[6]);
  forces[4] = -2.0*mu10*F[4]-2.0*mu01*(t12*F[4]-t19*F[1]-t54*F[4]-t59*F[7]);
  forces[5] = -2.0*mu10*F[5]-2.0*mu01*(t12*F[5]-t19*F[2]-t54*F[5]-t59*F[8]);
  forces[6] = -2.0*mu10*F[6]-2.0*mu01*(t12*F[6]-t24*F[0]-t59*F[3]-t90*F[6]);
  forces[7] = -2.0*mu10*F[7]-2.0*mu01*(t12*F[7]-t24*F[1]-t59*F[4]-t90*F[7]);
  forces[8] = -2.0*mu10*F[8]-2.0*mu01*(t12*F[8]-t24*F[2]-t59*F[5]-t90*F[8]);

  //----------------------------------------
  //  Teran03 penalty force term
  //----------------------------------------
  penaltyForceDensity(_k, F, forces);
}

//////////////////////////////////////////////////////////////////////
// stiffness matrix implementation
//
// The Matlab code to generate the code is below. This computes
// stiffness wrt F with the expectation that it will be multiplied
// by DF/Du later to get the stiffness wrt u
//
// see "MATERIAL.h" for more details
/*
function code = codegen_mooney_rivlin_stiffness_wrt_F()
    syms f00 f01 f02 f10 f11 f12 f20 f21 f22;
    F = [
        f00 f01 f02
        f10 f11 f12
        f20 f21 f22
        ];    
    
    syms b00 b01 b02 b10 b11 b12 b20 b21 b22 b30 b31 b32;
    b = [
        b00 b10 b20 b30
        b01 b11 b21 b31
        b02 b12 b22 b32
        ];
    
    display '-- calculating 2nd PK --';    
    mu01 = sym('mu01');
    mu10 = sym('mu10');
    C = mytrans(F)*F;
    I = eye(3);
    S = 2*( mu10*I + mu01*trace(C)*I - mu01*C );
    P = F * S;  % 1st PK
    
    forces = P * b;
    forces = reshape(forces, 12, 1);
    F = reshape(F, 9, 1);
    
    display '-- differentiating --';
    df_dF = mu01 * ones(12,9);
    for i = 1:12
        for j = 1:9
            df_dF(i,j) = diff( forces(i), F(j) );
        end
    end
    display '-- code gen''ing --';   
    code = maple('codegen[C]', df_dF, 'optimized');
end
*/
//////////////////////////////////////////////////////////////////////
void MOONEY_RIVLIN::computeStresses(TET& tet)
{
  MATRIX3 F = tet.F();
  const VEC3F* areaVecs = tet.b();
  _pf_pF.clear();

  const double f00 = F(0,0);
  const double f01 = F(0,1);
  const double f02 = F(0,2);
  const double f10 = F(1,0);
  const double f11 = F(1,1);
  const double f12 = F(1,2);
  const double f20 = F(2,0);
  const double f21 = F(2,1);
  const double f22 = F(2,2);

  const double b00 = areaVecs[0][0];
  const double b01 = areaVecs[0][1];
  const double b02 = areaVecs[0][2];
  const double b10 = areaVecs[1][0];
  const double b11 = areaVecs[1][1];
  const double b12 = areaVecs[1][2];
  const double b20 = areaVecs[2][0];
  const double b21 = areaVecs[2][1];
  const double b22 = areaVecs[2][2];
  const double b30 = areaVecs[3][0];
  const double b31 = areaVecs[3][1];
  const double b32 = areaVecs[3][2];

  const double t1 = f00*f00;
  const double t2 = f10*f10;
  const double t3 = f20*f20;
  const double t4 = f01*f01;
  const double t5 = f11*f11;
  const double t6 = f21*f21;
  const double t7 = f02*f02;
  const double t8 = f12*f12;
  const double t9 = f22*f22;
  const double t11 = _mu01*(t1+t2+t3+t4+t5+t6+t7+t8+t9);
  const double t13 = _mu01*(t1+t2+t3);
  const double t14 = t4*_mu01;
  const double t15 = t7*_mu01;
  const double t16 = _mu10+t11-t13-t14-t15;
  const double t22 = _mu01*(f00*f01+f10*f11+f20*f21);
  const double t23 = _mu01*f00;
  const double t25 = -t22+t23*f01;
  const double t31 = _mu01*(f00*f02+f10*f12+f20*f22);
  const double t33 = -t31+t23*f02;
  const double t36 = _mu01*f01;
  const double t37 = t36*f11;
  const double t38 = _mu01*f02;
  const double t39 = t38*f12;
  const double t40 = -t37-t39;
  const double t41 = 2.0*t40*b00;
  const double t42 = t23*f11;
  const double t44 = t36*f10;
  const double t46 = -2.0*t42+4.0*t44;
  const double t48 = t23*f12;
  const double t50 = t38*f10;
  const double t52 = -2.0*t48+4.0*t50;
  const double t55 = t36*f21;
  const double t56 = t38*f22;
  const double t57 = -t55-t56;
  const double t58 = 2.0*t57*b00;
  const double t59 = t23*f21;
  const double t61 = t36*f20;
  const double t63 = -2.0*t59+4.0*t61;
  const double t65 = t23*f22;
  const double t67 = t38*f20;
  const double t69 = -2.0*t65+4.0*t67;
  const double t73 = t1*_mu01;
  const double t75 = _mu01*(t4+t5+t6);
  const double t76 = -t73+_mu10+t11-t75-t15;
  const double t82 = _mu01*(f01*f02+f11*f12+f21*f22);
  const double t84 = -t82+t36*f02;
  const double t89 = 4.0*t42-2.0*t44;
  const double t91 = t23*f10;
  const double t92 = -t91-t39;
  const double t93 = 2.0*t92*b01;
  const double t94 = t36*f12;
  const double t96 = t38*f11;
  const double t98 = -2.0*t94+4.0*t96;
  const double t103 = 4.0*t59-2.0*t61;
  const double t105 = t23*f20;
  const double t106 = -t105-t56;
  const double t107 = 2.0*t106*b01;
  const double t108 = t36*f22;
  const double t110 = t38*f21;
  const double t112 = -2.0*t108+4.0*t110;
  const double t118 = _mu01*(t7+t8+t9);
  const double t119 = -t73-t14+_mu10+t11-t118;
  const double t124 = 4.0*t48-2.0*t50;
  const double t128 = 4.0*t94-2.0*t96;
  const double t130 = -t91-t37;
  const double t131 = 2.0*t130*b02;
  const double t135 = 4.0*t65-2.0*t67;
  const double t139 = 4.0*t108-2.0*t110;
  const double t141 = -t105-t55;
  const double t142 = 2.0*t141*b02;
  const double t147 = t5*_mu01;
  const double t148 = t8*_mu01;
  const double t149 = _mu10+t11-t13-t147-t148;
  const double t151 = _mu01*f10;
  const double t153 = -t22+t151*f11;
  const double t156 = -t31+t151*f12;
  const double t159 = _mu01*f11;
  const double t160 = t159*f21;
  const double t161 = _mu01*f12;
  const double t162 = t161*f22;
  const double t163 = -t160-t162;
  const double t164 = 2.0*t163*b00;
  const double t165 = t151*f21;
  const double t167 = t159*f20;
  const double t169 = -2.0*t165+4.0*t167;
  const double t171 = t151*f22;
  const double t173 = t161*f20;
  const double t175 = -2.0*t171+4.0*t173;
  const double t182 = t2*_mu01;
  const double t183 = -t182+_mu10+t11-t75-t148;
  const double t186 = -t82+t159*f12;
  const double t191 = 4.0*t165-2.0*t167;
  const double t193 = t151*f20;
  const double t194 = -t193-t162;
  const double t195 = 2.0*t194*b01;
  const double t196 = t159*f22;
  const double t198 = t161*f21;
  const double t200 = -2.0*t196+4.0*t198;
  const double t208 = -t182-t147+_mu10+t11-t118;
  const double t213 = 4.0*t171-2.0*t173;
  const double t217 = 4.0*t196-2.0*t198;
  const double t219 = -t193-t160;
  const double t220 = 2.0*t219*b02;
  const double t228 = t6*_mu01;
  const double t229 = t9*_mu01;
  const double t230 = _mu10+t11-t13-t228-t229;
  const double t232 = _mu01*f20;
  const double t234 = -t22+t232*f21;
  const double t237 = -t31+t232*f22;
  const double t247 = t3*_mu01;
  const double t248 = -t247+_mu10+t11-t75-t229;
  const double t252 = -t82+f21*_mu01*f22;
  const double t263 = -t247-t228+_mu10+t11-t118;
  const double t270 = 2.0*t40*b10;
  const double t274 = 2.0*t57*b10;
  const double t283 = 2.0*t92*b11;
  const double t287 = 2.0*t106*b11;
  const double t296 = 2.0*t130*b12;
  const double t300 = 2.0*t141*b12;
  const double t309 = 2.0*t163*b10;
  const double t321 = 2.0*t194*b11;
  const double t333 = 2.0*t219*b12;
  const double t369 = 2.0*t40*b20;
  const double t373 = 2.0*t57*b20;
  const double t382 = 2.0*t92*b21;
  const double t386 = 2.0*t106*b21;
  const double t395 = 2.0*t130*b22;
  const double t399 = 2.0*t141*b22;
  const double t408 = 2.0*t163*b20;
  const double t420 = 2.0*t194*b21;
  const double t432 = 2.0*t219*b22;
  const double t468 = 2.0*t40*b30;
  const double t472 = 2.0*t57*b30;
  const double t481 = 2.0*t92*b31;
  const double t485 = 2.0*t106*b31;
  const double t494 = 2.0*t130*b32;
  const double t498 = 2.0*t141*b32;
  const double t507 = 2.0*t163*b30;
  const double t519 = 2.0*t194*b31;
  const double t531 = 2.0*t219*b32;
  _pf_pF(0,0) += 2.0*t16*b00+2.0*t25*b01+2.0*t33*b02;
  _pf_pF(0,1) += t41+t46*b01+t52*b02;
  _pf_pF(0,2) += t58+t63*b01+t69*b02;
  _pf_pF(0,3) += 2.0*t25*b00+2.0*t76*b01+2.0*t84*b02;
  _pf_pF(0,4) += t89*b00+t93+t98*b02;
  _pf_pF(0,5) += t103*b00+t107+t112*b02;
  _pf_pF(0,6) += 2.0*t33*b00+2.0*t84*b01+2.0*t119*b02;
  _pf_pF(0,7) += t124*b00+t128*b01+t131;
  _pf_pF(0,8) += t135*b00+t139*b01+t142;
  _pf_pF(1,0) += t41+t89*b01+t124*b02;
  _pf_pF(1,1) += 2.0*t149*b00+2.0*t153*b01+2.0*t156*b02;
  _pf_pF(1,2) += t164+t169*b01+t175*b02;
  _pf_pF(1,3) += t46*b00+t93+t128*b02;
  _pf_pF(1,4) += 2.0*t153*b00+2.0*t183*b01+2.0*t186*b02;
  _pf_pF(1,5) += t191*b00+t195+t200*b02;
  _pf_pF(1,6) += t52*b00+t98*b01+t131;
  _pf_pF(1,7) += 2.0*t156*b00+2.0*t186*b01+2.0*t208*b02;
  _pf_pF(1,8) += t213*b00+t217*b01+t220;
  _pf_pF(2,0) += t58+t103*b01+t135*b02;
  _pf_pF(2,1) += t164+t191*b01+t213*b02;
  _pf_pF(2,2) += 2.0*t230*b00+2.0*t234*b01+2.0*t237*b02;
  _pf_pF(2,3) += t63*b00+t107+t139*b02;
  _pf_pF(2,4) += t169*b00+t195+t217*b02;
  _pf_pF(2,5) += 2.0*t234*b00+2.0*t248*b01+2.0*t252*b02;
  _pf_pF(2,6) += t69*b00+t112*b01+t142;
  _pf_pF(2,7) += t175*b00+t200*b01+t220;
  _pf_pF(2,8) += 2.0*t237*b00+2.0*t252*b01+2.0*t263*b02;
  _pf_pF(3,0) += 2.0*t16*b10+2.0*t25*b11+2.0*t33*b12;
  _pf_pF(3,1) += t270+t46*b11+t52*b12;
  _pf_pF(3,2) += t274+t63*b11+t69*b12;
  _pf_pF(3,3) += 2.0*t25*b10+2.0*t76*b11+2.0*t84*b12;
  _pf_pF(3,4) += t89*b10+t283+t98*b12;
  _pf_pF(3,5) += t103*b10+t287+t112*b12;
  _pf_pF(3,6) += 2.0*t33*b10+2.0*t84*b11+2.0*t119*b12;
  _pf_pF(3,7) += t124*b10+t128*b11+t296;
  _pf_pF(3,8) += t135*b10+t139*b11+t300;
  _pf_pF(4,0) += t270+t89*b11+t124*b12;
  _pf_pF(4,1) += 2.0*t149*b10+2.0*t153*b11+2.0*t156*b12;
  _pf_pF(4,2) += t309+t169*b11+t175*b12;
  _pf_pF(4,3) += t46*b10+t283+t128*b12;
  _pf_pF(4,4) += 2.0*t153*b10+2.0*t183*b11+2.0*t186*b12;
  _pf_pF(4,5) += t191*b10+t321+t200*b12;
  _pf_pF(4,6) += t52*b10+t98*b11+t296;
  _pf_pF(4,7) += 2.0*t156*b10+2.0*t186*b11+2.0*t208*b12;
  _pf_pF(4,8) += t213*b10+t217*b11+t333;
  _pf_pF(5,0) += t274+t103*b11+t135*b12;
  _pf_pF(5,1) += t309+t191*b11+t213*b12;
  _pf_pF(5,2) += 2.0*t230*b10+2.0*t234*b11+2.0*t237*b12;
  _pf_pF(5,3) += t63*b10+t287+t139*b12;
  _pf_pF(5,4) += t169*b10+t321+t217*b12;
  _pf_pF(5,5) += 2.0*t234*b10+2.0*t248*b11+2.0*t252*b12;
  _pf_pF(5,6) += t69*b10+t112*b11+t300;
  _pf_pF(5,7) += t175*b10+t200*b11+t333;
  _pf_pF(5,8) += 2.0*t237*b10+2.0*t252*b11+2.0*t263*b12;
  _pf_pF(6,0) += 2.0*t16*b20+2.0*t25*b21+2.0*t33*b22;
  _pf_pF(6,1) += t369+t46*b21+t52*b22;
  _pf_pF(6,2) += t373+t63*b21+t69*b22;
  _pf_pF(6,3) += 2.0*t25*b20+2.0*t76*b21+2.0*t84*b22;
  _pf_pF(6,4) += t89*b20+t382+t98*b22;
  _pf_pF(6,5) += t103*b20+t386+t112*b22;
  _pf_pF(6,6) += 2.0*t33*b20+2.0*t84*b21+2.0*t119*b22;
  _pf_pF(6,7) += t124*b20+t128*b21+t395;
  _pf_pF(6,8) += t135*b20+t139*b21+t399;
  _pf_pF(7,0) += t369+t89*b21+t124*b22;
  _pf_pF(7,1) += 2.0*t149*b20+2.0*t153*b21+2.0*t156*b22;
  _pf_pF(7,2) += t408+t169*b21+t175*b22;
  _pf_pF(7,3) += t46*b20+t382+t128*b22;
  _pf_pF(7,4) += 2.0*t153*b20+2.0*t183*b21+2.0*t186*b22;
  _pf_pF(7,5) += t191*b20+t420+t200*b22;
  _pf_pF(7,6) += t52*b20+t98*b21+t395;
  _pf_pF(7,7) += 2.0*t156*b20+2.0*t186*b21+2.0*t208*b22;
  _pf_pF(7,8) += t213*b20+t217*b21+t432;
  _pf_pF(8,0) += t373+t103*b21+t135*b22;
  _pf_pF(8,1) += t408+t191*b21+t213*b22;
  _pf_pF(8,2) += 2.0*t230*b20+2.0*t234*b21+2.0*t237*b22;
  _pf_pF(8,3) += t63*b20+t386+t139*b22;
  _pf_pF(8,4) += t169*b20+t420+t217*b22;
  _pf_pF(8,5) += 2.0*t234*b20+2.0*t248*b21+2.0*t252*b22;
  _pf_pF(8,6) += t69*b20+t112*b21+t399;
  _pf_pF(8,7) += t175*b20+t200*b21+t432;
  _pf_pF(8,8) += 2.0*t237*b20+2.0*t252*b21+2.0*t263*b22;
  _pf_pF(9,0) += 2.0*t16*b30+2.0*t25*b31+2.0*t33*b32;
  _pf_pF(9,1) += t468+t46*b31+t52*b32;
  _pf_pF(9,2) += t472+t63*b31+t69*b32;
  _pf_pF(9,3) += 2.0*t25*b30+2.0*t76*b31+2.0*t84*b32;
  _pf_pF(9,4) += t89*b30+t481+t98*b32;
  _pf_pF(9,5) += t103*b30+t485+t112*b32;
  _pf_pF(9,6) += 2.0*t33*b30+2.0*t84*b31+2.0*t119*b32;
  _pf_pF(9,7) += t124*b30+t128*b31+t494;
  _pf_pF(9,8) += t135*b30+t139*b31+t498;
  _pf_pF(10,0) += t468+t89*b31+t124*b32;
  _pf_pF(10,1) += 2.0*t149*b30+2.0*t153*b31+2.0*t156*b32;
  _pf_pF(10,2) += t507+t169*b31+t175*b32;
  _pf_pF(10,3) += t46*b30+t481+t128*b32;
  _pf_pF(10,4) += 2.0*t153*b30+2.0*t183*b31+2.0*t186*b32;
  _pf_pF(10,5) += t191*b30+t519+t200*b32;
  _pf_pF(10,6) += t52*b30+t98*b31+t494;
  _pf_pF(10,7) += 2.0*t156*b30+2.0*t186*b31+2.0*t208*b32;
  _pf_pF(10,8) += t213*b30+t217*b31+t531;
  _pf_pF(11,0) += t472+t103*b31+t135*b32;
  _pf_pF(11,1) += t507+t191*b31+t213*b32;
  _pf_pF(11,2) += 2.0*t230*b30+2.0*t234*b31+2.0*t237*b32;
  _pf_pF(11,3) += t63*b30+t485+t139*b32;
  _pf_pF(11,4) += t169*b30+t519+t217*b32;
  _pf_pF(11,5) += 2.0*t234*b30+2.0*t248*b31+2.0*t252*b32;
  _pf_pF(11,6) += t69*b30+t112*b31+t498;
  _pf_pF(11,7) += t175*b30+t200*b31+t531;
  _pf_pF(11,8) += 2.0*t237*b30+2.0*t252*b31+2.0*t263*b32;
}
