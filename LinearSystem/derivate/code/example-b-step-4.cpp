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

T compute_orig(const V x[3],V c,T a,T b,T InputArea,V dE[3],M ddE[3][3])
{
    T d_min = (x[0] - c).dot((x[2]-x[0]).cross(x[1]-x[0]))/(((x[2]-x[0]).cross(x[1]-x[0])).mag());
    T cot_alpha = ((x[2]-x[0]).dot(x[1]-x[0]))/(((x[2]-x[0]).cross(x[1]-x[0])).mag());
    T cot_beta = ((x[0]-x[1]).dot(x[2]-x[1]))/(((x[0]-x[1]).cross(x[2]-x[1])).mag());
    T cot_gamma = ((x[1]-x[2]).dot(x[0]-x[2]))/(((x[1]-x[2]).cross(x[0]-x[2])).mag());
    T E_Dirichlet = 1/(d_min*d_min)*(cot_alpha*(x[1]-x[2]).mag2() + cot_beta*(x[0]-x[2]).mag2() + cot_gamma*(x[0]-x[1]).mag2());
    T E_area = 1/(d_min*d_min)*(((x[0]-x[2]).cross(x[1]-x[2])).mag2())/InputArea;
    T E = a*E_Dirichlet + b*E_area;
    return E;
}

T compute(const V x[3],V c,T a,T b,T InputArea,V dE[3],M ddE[3][3])
{
    V z=x[0]-c;
    V u=x[0]-x[2];
    V v=x[1]-x[2];
    V N=u.cross(v);
    T A2=N.mag2();
    T A=sqrt(A2);
    T d_min = z.dot(N);
    T E_Dirichlet = 1/(d_min*d_min)*2*A2*A;
    T E_area = 1/(d_min*d_min)*A2*A2/InputArea;
    T E = a*E_Dirichlet + b*E_area;
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

    V x[3],dE[3];
    M ddE[3][3];
    for(int i=0;i<3;i++)
    {
        fill_random(rand,x[i],-1.0,1.0);
    }

    T E0=compute(x,c,a,b,InputArea,dE,ddE);
    T E1=compute_orig(x,c,a,b,InputArea,dE,ddE);

    printf("%.16g %.16g %.16g\n",E0,E1,E1-E0);

    return 0;
}
