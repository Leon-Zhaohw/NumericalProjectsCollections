//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include "la.h"
#include <cstdio>

random_type rnd;

double eps=1e-6;

typedef double T;
typedef vec<double,3> V3;
typedef vec<double,2> V2;
typedef mat<double,2> M2;
typedef mat<double,3> M3;
typedef mat<double,3,2> M;

T compute(const V3 x[3],const M2& R_X_inv,T A_X,T lambda,T mu,V3 dE[3],M3 ddE[3][3])
{
    M D;
    D.set_column(x[0]-x[2],0);
    D.set_column(x[1]-x[2],1);
    M F_bar=D*R_X_inv;
    M2 E_bar=0.5*(F_bar.t()*F_bar-1);
    T B=E_bar.tr();
    T W=A_X*(lambda/2*B*B + mu*(E_bar.double_contract(E_bar.t())));
    M2 E=R_X_inv*E_bar;
    M F=F_bar*R_X_inv.t();
    M3 G=F_bar*F_bar.t();
    M dW=A_X*(lambda*B*F + 2*mu*F_bar*E.t());
    M2 H=R_X_inv*R_X_inv.t();
    M2 C=lambda*B*H + 2*mu*R_X_inv*E.t();
    dE[0]=dW.column(0);
    dE[1]=dW.column(1);
    dE[2]=-dE[0]-dE[1];
    for(int s=0;s<2;s++)
        for(int v=0;v<2;v++)
            ddE[s][v]=A_X*(mu*H(s,v)*G+C(s,v)
                +outer(lambda*F.column(s),F.column(v))
                +outer(mu*F.column(v),F.column(s)));
    for(int s=0;s<2;s++)
    {
        ddE[s][2]=-ddE[s][0]-ddE[s][1];
        ddE[2][s]=-ddE[0][s]-ddE[1][s];
    }
    ddE[2][2]=-ddE[0][2]-ddE[1][2];

    return W;
}

int main(int argc, char* argv[])
{
    random_type rand;
    M2 R_X_inv;
    double A_X,lambda,mu;
    fill_random(rand,R_X_inv,-1.0,1.0);
    fill_random(rand,A_X,.1,1.0);
    fill_random(rand,lambda,.1,1.0);
    fill_random(rand,mu,.1,1.0);

    V3 x[3],dx[3],x1[3];
    for(int i=0;i<3;i++)
    {
        fill_random(rand,x[i],-1.0,1.0);
        fill_random(rand,dx[i],-eps,eps);
        x1[i]=x[i]+dx[i];
    }

    V3 dE0[3],dE1[3];
    M3 ddE0[3][3],ddE1[3][3];
    T E0=compute(x,R_X_inv,A_X,lambda,mu,dE0,ddE0);
    T E1=compute(x1,R_X_inv,A_X,lambda,mu,dE1,ddE1);
    T d0=(E1-E0)/eps;
    T d1=0;
    for(int i=0;i<3;i++) d1+=(dE0[i]+dE1[i]).dot(dx[i]);
    d1/=2*eps;
    printf("%g %g %g\n",d0,d1,fabs(d0-d1)/std::max(std::max(fabs(d0),fabs(d1)),1e-30));

    T e0=0,e1=0,e2=0;
    for(int i=0;i<3;i++)
    {
        V3 u=dE1[i]-dE0[i],v;
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
