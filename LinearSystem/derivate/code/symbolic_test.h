//#####################################################################
// Copyright 2015, Daniel Ram, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef __symbolic_test__
#define __symbolic_test__
#include "symbolic_tensor.h"

extern double eps;

// Derivative computed explicitly.
void Test(const char* name,const symbolic_tensor& dx,const symbolic_tensor& a,
    const symbolic_tensor& b,const symbolic_tensor& da,const symbolic_tensor& db)
{
    std::string norm_ind;
    std::string diff_ind;
    for(size_t i=0;i<a.size.size();i++)
        norm_ind+='a'+i;
    for(size_t i=0;i<dx.size.size();i++)
        diff_ind+='r'+i;
    std::string full_ind=norm_ind+diff_ind;
    symbolic_tensor E(norm_ind,b(norm_ind)-a(norm_ind));
    symbolic_tensor F(norm_ind,(da(full_ind)+db(full_ind))*dx(diff_ind)/2);
    symbolic_tensor G(norm_ind,E(norm_ind)-F(norm_ind));
    double e=norm(E.x);
    double f=norm(F.x);
    double g=norm(G.x);
    printf("DIFF %s: %.16g %.16g -> %.16g\n",name,e/eps,f/eps,g/max(max(abs(e),abs(f)),1e-30));
}
#define TEST(x) Test(#x,dx,z0.x,z1.x,z0.d##x,z1.d##x);

// Derivative applied to dx.
void Test_Hess(const char* name,const symbolic_tensor& dx,const symbolic_tensor& a,
    const symbolic_tensor& b,const symbolic_tensor& da,const symbolic_tensor& db)
{
    std::string norm_ind;
    for(size_t i=0;i<a.size.size();i++)
        norm_ind+='a'+i;
    symbolic_tensor E(norm_ind,b(norm_ind)-a(norm_ind));
    symbolic_tensor F(norm_ind,(da(norm_ind)+db(norm_ind))/2);
    symbolic_tensor G(norm_ind,E(norm_ind)-F(norm_ind));
    double e=norm(E.x);
    double f=norm(F.x);
    double g=norm(G.x);
    printf("DIFF %s: %.16g %.16g -> %.16g\n",name,e/eps,f/eps,g/max(max(abs(e),abs(f)),1e-30));
}
#define TEST_H(x) Test_Hess(#x,dx,z0.x,z1.x,z0.d##x,z1.d##x);

// Derivative computed explicitly.
void Test(const char* name,const symbolic_tensor dx[],const symbolic_tensor& a,
    const symbolic_tensor& b,const symbolic_tensor da[],const symbolic_tensor db[],int n)
{
    std::string norm_ind;
    std::string diff_ind;
    for(size_t i=0;i<a.size.size();i++)
        norm_ind+='a'+i;
    for(size_t i=0;i<dx[0].size.size();i++)
        diff_ind+='r'+i;
    std::string full_ind=norm_ind+diff_ind;
    symbolic_tensor E(norm_ind,b(norm_ind)-a(norm_ind));
    symbolic_tensor F(norm_ind,(da[0](full_ind)+db[0](full_ind))*dx[0](diff_ind)/2);
    for(int i=1;i<n;i++)
        F.set(norm_ind,F(norm_ind)+(da[i](full_ind)+db[i](full_ind))*dx[i](diff_ind)/2);
    symbolic_tensor G(norm_ind,E(norm_ind)-F(norm_ind));
    double e=norm(E.x);
    double f=norm(F.x);
    double g=norm(G.x);
    printf("DIFF %s: %.16g %.16g -> %.16g\n",name,e/eps,f/eps,g/max(max(abs(e),abs(f)),1e-30));
}
#define TEST_N(x,n) Test(#x,dx,z0.x,z1.x,z0.d##x,z1.d##x,n);

#define PR(z_){printf("%s: ",#z_);for(auto q_:z_.x) printf("%.16g ",q_);printf("\n");}

#endif
