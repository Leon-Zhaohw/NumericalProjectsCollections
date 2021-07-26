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
    T f=b/InputArea;
    V z=x[0]-c;
    V u=x[0]-x[2];
    V v=x[1]-x[2];
    V N=u.cross(v);
    T A2=N.mag2();
    T A=sqrt(A2);
    T d=z.dot(N),d2=d*d;
    T g=1/d2;
    T C=(2*a*A+f*A2)*A2;
    T E=g*C;
    T h=6*a*A+4*f*A2;
    T s=2*g*g,s1=s*d,s2=s1*C,s3=g*h;
    V dE_dN=-s2*z+s3*N;
    V dE_dz=-s2*N;

    V dN_dx[3]={-v,u,v-u};

    for(int i=0;i<3;i++)
        dE[i]=-dN_dx[i].cross(dE_dN);
    dE[0]+=dE_dz;

    T m=s*(4*g*d2-1),s4=m*C;
    V p=s4*z-s1*h*N;
    M ddE_dN_dN=outer(p,z)+outer(g*((6*a/A+8*f)*N-2*g*d*h*z),N)+s3;
    M ddE_dN_dz=outer(p,N)-s2;
    M ddE_dz_dz=s4*outer(N,N);

    int ddN_dx_dx[3][3]={{0,1,-1},{-1,0,1},{1,-1,0}};
    M cp[3],q[3],r[3];
    for(int i=0;i<3;i++)
    {
        cp[i]=cpm(dN_dx[i]);
        q[i]=-cp[i]*ddE_dN_dz;
        r[i]=cp[i]*ddE_dN_dN;
    }

    for(int i=0;i<3;i++)
        for(int j=i;j<3;j++)
            ddE[i][j]=-r[i]*cp[j]-cpm(dE_dN)*ddN_dx_dx[i][j];
    for(int i=0;i<3;i++)
        ddE[0][i]+=q[i].t();
    ddE[0][0]+=ddE_dz_dz+q[0];
    for(int i=0;i<3;i++)
        for(int j=0;j<i;j++)
            ddE[i][j]=ddE[j][i].t();

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
