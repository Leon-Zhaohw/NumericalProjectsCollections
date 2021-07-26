//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include "svd-test.h"
#include <cstdio>

random_type rnd;

double eps=1e-6;

typedef double T;
typedef vec<double,3> V;
typedef mat<double,3> M;

T compute_diag(const V& sig,V& dE,M& ddE,V& div)
{
    T a=sum(sig);
    T b=sig.dot(sig);
    T E=cos(a)+sin(b);

    V da=ones<T,3>();
    V db=sig*2;
    dE=-sin(a)*da+cos(b)*db;

    M ddb=M::id()*2;
    ddE=-cos(a)*outer(da,da)-sin(b)*outer(db,db)+cos(b)*ddb;

    for(int i=0;i<3;i++)
    {
        int j=(i+1)%3;
        int k=(i+2)%3;
        div[i]=(dE[j]-dE[k])/(sig[j]-sig[k]);
    }

    return E;
}

// ddE = Hessian : dF;
T compute(const V& sig, const M& UU, const M& VV, const M& dF, M& dE, M& ddE)
{
    diff_helper<T,3> h;
    V diag_dE;
    T E = compute_diag(sig,diag_dE,h.H,h.A);
    h.B = sum_over_sum(diag_dE,sig);
    h.flip=false;
    dE = UU * diag(diag_dE) * VV.t();
    ddE = UU * h.apply(UU.t() * dF * VV) * VV.t();
    return E;
}

int main(int argc, char* argv[])
{
    random_type rand;
    M UU=random_rotation(rand,M_PI),VV=random_rotation(rand,M_PI);
    V Sig;
    fill_random(rand,Sig,-1.,1.);

    M dUU=random_rotation(rand,eps),dVV=random_rotation(rand,eps);
    V dSig;
    fill_random(rand,dSig,-eps,eps);

    M UU1=UU*dUU;
    M VV1=VV*dVV;
    V Sig1=Sig+dSig;

    M F=UU*diag(Sig)*VV.t();
    M F1=UU1*diag(Sig1)*VV1.t();
    M dF=F1-F;

    M dE0,dE1;
    M ddE0,ddE1;

    T E0=compute(Sig,UU,VV,dF,dE0,ddE0);
    T E1=compute(Sig1,UU1,VV1,dF,dE1,ddE1);
    T d0=(E1-E0)/eps;
    T d1=(dE0+dE1).double_contract(dF)/(2*eps);
    printf("%g %g %g\n",d0,d1,fabs(d0-d1)/std::max(std::max(fabs(d0),fabs(d1)),1e-30));

    M u=dE1-dE0,v=(ddE0+ddE1)/2;

    T e0=u.frobenius_norm()/eps;
    T e1=v.frobenius_norm()/eps;
    T e2=(v-u).frobenius_norm()/eps;
    printf("%g %g %g\n",e0,e1,e2/std::max(std::max(e0,e1),1e-30));

    return 0;
}
