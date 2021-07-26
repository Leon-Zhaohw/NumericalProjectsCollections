//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include "la.h"
#include <cstdio>
#include <gsl/gsl_poly.h>

random_type rnd;

double eps=1e-6;

typedef double T;
typedef vec<double,3> V;
typedef mat<double,3> M;

T psi(V x,V v,V& dE,M& ddE,M& v_dddE,M& vv_ddddE)
{
    T a=x.dot(x);
    T b=x.dot(v);
    T d=v.dot(v);
    T c=cos(a);
    T s=sin(a);
    dE=2*c*x;
    ddE=-s*4*outer(x,x)+c*2;
    v_dddE=-8*c*b*outer(x,x)-4*s*outer(v,x)-4*s*outer(x,v)-4*s*b;
    vv_ddddE=8*(2*s*b*b-c*d)*outer(x,x)-16*c*b*outer(v,x)-16*c*b*outer(x,v)-s*8*outer(v,v)-4*(2*c*b*b+s*d);

    return s;
}

T compute(V x,V r,T a,V& dE,M& ddE)
{
    V v=x-r;

    T C0=r[0]*r[1]*r[2]-a;
    T C1=v[0]*r[1]*r[2]+r[0]*v[1]*r[2]+r[0]*r[1]*v[2];
    T C2=v[0]*v[1]*r[2]+r[0]*v[1]*v[2]+v[0]*r[1]*v[2];
    T C3=v[0]*v[1]*v[2];

    T x0=-1,x1=-1,x2=-1;
    int n=gsl_poly_solve_cubic(C2/C3,C1/C3,C0/C3,&x0,&x1,&x2);
    assert(n>0);
    if(x0<0 || x0>1)
    {
        assert(n>1);
        x0=x1;
        if(x0<0 || x0>1)
        {
            assert(n>2);
            x0=x2;
            assert(!(x0<0 || x0>1));
        }
    }
    T s=x0,t=1-s,t2h=t*t/2;
    V q=r+v*s;
    M Mq;
    Mq.x[0]=Mq.x[4]=Mq.x[8]=0;
    Mq.x[1]=Mq.x[3]=q[2];
    Mq.x[2]=Mq.x[6]=q[1];
    Mq.x[5]=Mq.x[7]=q[0];

    V Mqq=Mq*q,Mqv=Mq*v;

    V Tc=Mqq*s;
    T Ts=Mqq.dot(v);
    T Tss=Mqv.dot(v)*2;
    V Tcs=Mqq+Mqv*(2*s);
    M Tcc=2*s*s*Mq;

    V ds=-Tc/Ts;
    M dds=-(Tss*outer(ds,ds)+outer(Tcs,ds)+outer(ds,Tcs)+Tcc)/Ts;
    M dq=outer(v,ds)+s;

    V dZ;
    M ddZ,v_dddZ,vv_ddddZ;
    T Z=psi(q,v,dZ,ddZ,v_dddZ,vv_ddddZ);
    V vddZ=ddZ*v;
    V vvdddZ=v_dddZ*v;
    M A=outer(v,ds)-1.5*t;

    T E=Z+dZ.dot(v)*t+vddZ.dot(v)*t2h;
    dE=dZ+vddZ*t+dq.t()*vvdddZ*t2h;
    ddE=ddZ-t*A.t()*v_dddZ*A+(2.5*t2h+1)*t*v_dddZ+vvdddZ.dot(v)*t2h*dds+t2h*dq.t()*vv_ddddZ*dq;

    return E;
}

int main(int argc, char* argv[])
{
    random_type rand;
    V c;
    double a;
    fill_random(rand,a,0.8,0.9);

    V x,r,dx,x1;
    fill_random(rand,r,1.0,1.1);
    fill_random(rand,x,-1.0,0.9);
    fill_random(rand,dx,-eps,eps);
    x1=x+dx;

    V dE0,dE1;
    M ddE0,ddE1;
    T E0=compute(x,r,a,dE0,ddE0);
    T E1=compute(x1,r,a,dE1,ddE1);
    T d0=(E1-E0)/eps;
    T d1=(dE0+dE1).dot(dx)/(2*eps);
    printf("%g %g %g\n",d0,d1,fabs(d0-d1)/std::max(std::max(fabs(d0),fabs(d1)),1e-30));

    V u=dE1-dE0,v=(ddE0+ddE1)*dx/2;
    T e0=u.mag()/eps;
    T e1=v.mag()/eps;
    T e2=(v-u).mag()/eps;
    printf("%g %g %g\n",e0,e1,e2/std::max(std::max(e0,e1),1e-30));

    return 0;
}
