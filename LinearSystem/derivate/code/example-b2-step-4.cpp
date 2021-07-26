//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include "la.h"
#include <cstdio>

random_type rnd;

double eps=1e-6;

typedef double T;
typedef vec<double,3> V;
typedef mat<double,3> M;

T compute(const V x[3],V c,T a,T b,T InputArea,V dE[3],M ddE[3][3])
{
    V z=x[0]-c;
    V u=x[0]-x[2];
    V v=x[1]-x[2];
    V N=u.cross(v);
    N=x[2];
    T d=z.dot(N),d2=d*d;
    d=x[1][0];d2=d*d;
    T A2=N.mag2();
    T A=sqrt(A2);
    T g=1/d2;
    T f=b/InputArea;
    T C=2*a*A*A2+f*A2*A2;
    T E=g*C;

    T dC_dA2=3*a*A+2*f*A2;
    T dC_dA2_dA2=1.5*a/A+2*f;
    T dE_dg=C;
    T dE_dA2=g*dC_dA2;
    T dE_dA2_dg=dC_dA2;
    T dE_dA2_dA2=g*dC_dA2_dA2;

    T dg_dd=-2*g*g*d;
    T dg_dd_dd=6*g*g;
    T dE_dd=dE_dg*dg_dd;
    T dE_dd_dd=dE_dg*dg_dd_dd;
    T dE_dA2_dd=dE_dA2_dg*dg_dd;

    V dA2_dN=N*2;
    M dA2_dN_dN=M::id()*2;
    V dE_dN=dE_dA2*dA2_dN;
    M dE_dN_dN=dE_dA2_dA2*outer(dA2_dN,dA2_dN)+dE_dA2*dA2_dN_dN;
    V dE_dN_dd=dE_dA2_dd*dA2_dN;

    dE[1][0]=dE_dd;
    dE[2]=dE_dN;
    ddE[1][1](0,0)=dE_dd_dd;
    ddE[2][1].set_column(dE_dN_dd,0);
    ddE[1][2].set_row(dE_dN_dd,0);
    ddE[2][2]=dE_dN_dN;

    return E;
}

int main(int argc, char* argv[])
{
    random_type rand;
    V c;
    double a,b,InputArea;
    fill_random(rand,c,-1.0,1.0);
    fill_random(rand,a,.1,1.0);
    fill_random(rand,b,.1,1.0);
    fill_random(rand,InputArea,.1,1.0);

    V x[3],dx[3],x1[3];
    for(int i=0;i<3;i++)
    {
        fill_random(rand,x[i],-1.0,1.0);
        fill_random(rand,dx[i],-eps,eps);
        x1[i]=x[i]+dx[i];
    }

    V dE0[3],dE1[3];
    M ddE0[3][3],ddE1[3][3];
    for(int i=0;i<3;i++)
    {
        dE0[i].make_zero();
        dE1[i].make_zero();
        for(int j=0;j<3;j++)
        {
            ddE0[i][j].make_zero();
            ddE1[i][j].make_zero();
        }
    }

    T E0=compute(x,c,a,b,InputArea,dE0,ddE0);
    T E1=compute(x1,c,a,b,InputArea,dE1,ddE1);
    T d0=(E1-E0)/eps;
    T d1=0;
    for(int i=0;i<3;i++) d1+=(dE0[i]+dE1[i]).dot(dx[i]);
    d1/=2*eps;
    printf("%g %g %g\n",d0,d1,fabs(d0-d1)/std::max(std::max(fabs(d0),fabs(d1)),1e-30));

    T e0=0,e1=0,e2=0;
    for(int i=0;i<3;i++)
    {
        V u=dE1[i]-dE0[i],v;
        v.make_zero();
        for(int j=0;j<3;j++)
            v+=(ddE0[i][j]+ddE1[i][j])*dx[j]/2;
        e0+=u.mag2();
        e1+=v.mag2();
        e2+=(v-u).mag2();
    }
    e0=sqrt(e0)/eps;
    e1=sqrt(e1)/eps;
    e2=sqrt(e2)/eps;
    printf("%g %g %g\n",e0,e1,e2/std::max(std::max(e0,e1),1e-30));

    return 0;
}
