//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include "symbolic_test.h"

random_type rnd;

double eps=1e-6;

struct HELPER
{
    symbolic_tensor x[3],dx[3];
    symbolic_tensor z,u,v;
    symbolic_tensor N,dN[3],ddN[3][3];
    symbolic_tensor A2;
    symbolic_tensor A;
    symbolic_tensor d,dd_dN,dd_dz;
    symbolic_tensor d2;
    symbolic_tensor g;
    symbolic_tensor C;
    symbolic_tensor E,dE_dN,dE_dz,ddE_dN_dN,ddE_dN_dz,ddE_dz_dz;
    symbolic_tensor dE[3],ddE[3][3];

    symbolic_tensor c,InputArea,a,b;

    void init()
    {
        for(int i=0;i<3;i++)
        {
            x[i].set_zero({3});
            x[i].random(rnd,-1,1);
        }

        c.set_zero({3});
        c.random(rnd,0.1,1);

        a.set_zero({});
        a.random(rnd,0.1,1);

        b.set_zero({});
        b.random(rnd,0.1,1);

        InputArea.set_zero({});
        InputArea.random(rnd,0.1,1);
    }

    void fill()
    {
        symbolic_tensor id,e;
        id.set_id(3);
        e.set_perm();
        symbolic_tensor f=b/InputArea;

        z.set("i",x[0]("i")-c("i"));
        u.set("i",x[0]("i")-x[2]("i"));
        v.set("i",x[1]("i")-x[2]("i"));
        N.set("i",e("ijk")*u("j")*v("k"));
        A2.set("",N("i")*N("i"));
        A.set("",sqrt(A2));
        d.set("",z("i")*N("i"));
        d2.set("",d*d);
        g.set("",1/d2);
        C.set("",(a*2*A+f*A2)*A2);
        E.set("",g*C);

        symbolic_tensor h=6*a*A+4*f*A2;
        dE_dN.set("r",-g*g*2*d*z("r")*C+g*h*N("r"));
        dE_dz.set("r",-g*g*2*d*C*N("r"));

        symbolic_tensor dN_dx[3],dz_dx[3];
        dN_dx[0].set("ir",e("irk")*v("k"));
        dN_dx[1].set("ir",e("ijr")*u("j"));
        dN_dx[2].set("ir",-e("irk")*v("k")-e("ijr")*u("j"));
        dz_dx[0].set("ir",id("ir"));
        dz_dx[1].set_zero({3,3});
        dz_dx[2].set_zero({3,3});

        for(int i=0;i<3;i++) dE[i].set("r",dE_dN("i")*dN_dx[i]("ir")+dE_dz("i")*dz_dx[i]("ir"));

        symbolic_tensor m=2*g*g*(4*g*d*d-1);
        ddE_dN_dN.set("rs",m*z("s")*z("r")*C+g*h*id("rs")+g*(6*a/A+8*f)*N("s")*N("r")-g*g*2*d*z("r")*h*N("s")-g*g*2*d*z("s")*h*N("r"));
        ddE_dN_dz.set("rs",m*N("s")*z("r")*C-g*g*2*d*id("rs")*C-g*g*2*d*N("s")*h*N("r"));
        ddE_dz_dz.set("rs",m*N("s")*N("r")*C);

        int ddN_dx_dx[3][3]={{0,1,-1},{-1,0,1},{1,-1,0}};

        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                ddE[i][j].set("rs",
                    ddE_dN_dN("ij")*dN_dx[i]("ir")*dN_dx[j]("js")+
                    ddE_dN_dz("ij")*dN_dx[i]("ir")*dz_dx[j]("js")+
                    ddE_dN_dz("ji")*dz_dx[i]("ir")*dN_dx[j]("js")+
                    ddE_dz_dz("ij")*dz_dx[i]("ir")*dz_dx[j]("js")+
                    ddN_dx_dx[i][j]*dE_dN("i")*e("irs"));
            }
        }

        // E_{,r} = E_{,a} A_{a,r} + E_{,b} B_{b,r}
        // E_{,rs} = E_{,aA} A_{a,r} A_{A,s} + E_{,bA} B_{b,r} A_{A,s} + E_{,aB} A_{a,r} B_{B,s} + E_{,bB} B_{b,r} B_{B,s} + E_{,a} A_{a,rs} + E_{,b} B_{b,rs}
    }
};

#define TEST_dN(x){symbolic_tensor dN("i",z1.N("i")-z0.N("i"));Test(#x,dN,z0.x,z1.x,z0.d##x##_dN,z1.d##x##_dN);}

void Test(const char* name,const symbolic_tensor& dx,const symbolic_tensor& dy,
    const symbolic_tensor& a,const symbolic_tensor& b,
    const symbolic_tensor& da_dx,const symbolic_tensor& da_dy,
    const symbolic_tensor& db_dx,const symbolic_tensor& db_dy)
{
    std::string norm_ind;
    std::string dx_ind;
    std::string dy_ind;
    for(size_t i=0;i<a.size.size();i++)
        norm_ind+='a'+i;
    for(size_t i=0;i<dx.size.size();i++)
        dx_ind+='r'+i;
    for(size_t i=0;i<dy.size.size();i++)
        dy_ind+='r'+i;
    std::string full_dx_ind=norm_ind+dx_ind;
    std::string full_dy_ind=norm_ind+dy_ind;
    symbolic_tensor E(norm_ind,b(norm_ind)-a(norm_ind));
    symbolic_tensor F(norm_ind,(da_dx(full_dx_ind)+db_dx(full_dx_ind))*dx(dx_ind)/2+(da_dy(full_dy_ind)+db_dy(full_dy_ind))*dy(dy_ind)/2);
    symbolic_tensor G(norm_ind,E(norm_ind)-F(norm_ind));
    double e=norm(E.x);
    double f=norm(F.x);
    double g=norm(G.x);
    printf("DIFF %s: %.16g %.16g -> %.16g\n",name,e/eps,f/eps,g/max(max(abs(e),abs(f)),1e-30));
}

#define TEST_d_(x){symbolic_tensor dN("i",z1.N("i")-z0.N("i")),dz("i",z1.z("i")-z0.z("i"));Test(#x,dN,dz,z0.x,z1.x,z0.d##x##_dN,z0.d##x##_dz,z1.d##x##_dN,z1.d##x##_dz);}

int main(int argc, char* argv[])
{
    HELPER z0;
    z0.init();

    symbolic_tensor dx[3];
    for(int i=0;i<3;i++)
    {
        dx[i].set_zero({3});
        dx[i].random(rnd,-eps,eps);
    }

    HELPER z1=z0;
    for(int i=0;i<3;i++) z1.x[i].set("r",z0.x[i]("r")+dx[i]("r"));

    z0.fill();
    z1.fill();

    TEST_d_(E);
    TEST_N(E,3);

    TEST_d_(dE_dN);
    TEST_N(dE[0],3);
    TEST_N(dE[1],3);
    TEST_N(dE[2],3);

    return 0;
}
