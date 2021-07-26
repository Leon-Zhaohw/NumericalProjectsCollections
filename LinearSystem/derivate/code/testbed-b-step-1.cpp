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
    symbolic_tensor A2,dA2_dN,ddA2_dN_dN;
    symbolic_tensor A,dA_dN,ddA_dN_dN;
    symbolic_tensor d,dd_dN,dd_dz,ddd_dN_dN,ddd_dN_dz,ddd_dz_dz,ddd_dz_dN;
    symbolic_tensor d2,dd2_dN,dd2_dz,ddd2_dN_dN,ddd2_dN_dz,ddd2_dz_dN,ddd2_dz_dz;
    symbolic_tensor g,dg_dN,dg_dz,ddg_dN_dN,ddg_dN_dz,ddg_dz_dN,ddg_dz_dz;
    symbolic_tensor C,dC_dN,ddC_dN_dN;
    symbolic_tensor E,dE_dN,dE_dz,ddE_dN_dN,ddE_dN_dz,ddE_dz_dN,ddE_dz_dz;
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
        
        z.set("i",x[0]("i")-c("i"));
        u.set("i",x[0]("i")-x[2]("i"));
        v.set("i",x[1]("i")-x[2]("i"));
        N.set("i",e("ijk")*u("j")*v("k"));
        A2.set("",N("i")*N("i"));
        A.set("",sqrt(A2));
        d.set("",z("i")*N("i"));
        d2.set("",d*d);
        g.set("",1/d2);
        C.set("",a*2*A*A2+b/InputArea*A2*A2);
        E.set("",g*C);

        dA2_dN.set("r",2*N("r"));
        dA_dN.set("r",N("r")/A);
        dd_dN.set("r",z("r"));
        dd_dz.set("r",N("r"));
        dd2_dN.set("r",2*d*dd_dN("r"));
        dd2_dz.set("r",2*d*dd_dz("r"));
        dg_dN.set("r",-g*g*dd2_dN("r"));
        dg_dz.set("r",-g*g*dd2_dz("r"));
        dC_dN.set("r",(3*a*A+2*b/InputArea*A2)*dA2_dN("r"));
        dE_dN.set("r",dg_dN("r")*C+g*dC_dN("r"));
        dE_dz.set("r",dg_dz("r")*C);

        symbolic_tensor dN_dx[3],dz_dx[3];
        dN_dx[0].set("ir",e("irk")*v("k"));
        dN_dx[1].set("ir",e("ijr")*u("j"));
        dN_dx[2].set("ir",-e("irk")*v("k")-e("ijr")*u("j"));
        dz_dx[0].set("ir",id("ir"));
        dz_dx[1].set_zero({3,3});
        dz_dx[2].set_zero({3,3});

        for(int i=0;i<3;i++) dE[i].set("r",dE_dN("i")*dN_dx[i]("ir")+dE_dz("i")*dz_dx[i]("ir"));

        ddA2_dN_dN.set("rs",2*id("rs"));
        ddA_dN_dN.set("rs",id("rs")/A-N("r")*dA_dN("s")/A2);
        ddd_dN_dN.set_zero({3,3});
        ddd_dz_dN.set("rs",id("rs"));
        ddd_dN_dz.set("rs",id("rs"));
        ddd_dz_dz.set_zero({3,3});
        ddd2_dN_dN.set("rs",2*dd_dN("s")*dd_dN("r"));
        ddd2_dz_dN.set("rs",2*d*ddd_dz_dN("rs")+2*dd_dN("s")*dd_dz("r"));
        ddd2_dN_dz.set("rs",2*d*ddd_dN_dz("rs")+2*dd_dz("s")*dd_dN("r"));
        ddd2_dz_dz.set("rs",2*dd_dz("s")*dd_dz("r"));
        ddg_dN_dN.set("rs",-2*g*dg_dN("s")*dd2_dN("r")-g*g*ddd2_dN_dN("rs"));
        ddg_dz_dN.set("rs",-2*g*dg_dN("s")*dd2_dz("r")-g*g*ddd2_dz_dN("rs"));
        ddg_dN_dz.set("rs",-2*g*dg_dz("s")*dd2_dN("r")-g*g*ddd2_dN_dz("rs"));
        ddg_dz_dz.set("rs",-2*g*dg_dz("s")*dd2_dz("r")-g*g*ddd2_dz_dz("rs"));
        ddC_dN_dN.set("rs",(3*a*A+2*b/InputArea*A2)*ddA2_dN_dN("rs")+(3*a*dA_dN("s")+2*b/InputArea*dA2_dN("s"))*dA2_dN("r"));
        ddE_dN_dN.set("rs",ddg_dN_dN("rs")*C+g*ddC_dN_dN("rs")+dg_dN("r")*dC_dN("s")+dg_dN("s")*dC_dN("r"));
        ddE_dz_dN.set("rs",ddg_dz_dN("rs")*C+dg_dz("r")*dC_dN("s"));
        ddE_dN_dz.set("rs",ddg_dN_dz("rs")*C+dg_dz("s")*dC_dN("r"));
        ddE_dz_dz.set("rs",ddg_dz_dz("rs")*C);

        symbolic_tensor ddN_dx_dx[3][3];
        ddN_dx_dx[0][0].set_zero({3,3,3});
        ddN_dx_dx[1][0].set("irs",e("isr"));
        ddN_dx_dx[2][0].set("irs",-e("isr"));

        ddN_dx_dx[0][1].set("irs",e("irs"));
        ddN_dx_dx[1][1].set_zero({3,3,3});
        ddN_dx_dx[2][1].set("irs",-e("irs"));

        ddN_dx_dx[0][2].set("irs",-e("irs"));
        ddN_dx_dx[1][2].set("irs",-e("isr"));
        ddN_dx_dx[2][2].set("irs",e("irs")+e("isr"));

        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                ddE[i][j].set("rs",
                    ddE_dN_dN("ij")*dN_dx[i]("ir")*dN_dx[j]("js")+
                    ddE_dN_dz("ij")*dN_dx[i]("ir")*dz_dx[j]("js")+
                    ddE_dz_dN("ij")*dz_dx[i]("ir")*dN_dx[j]("js")+
                    ddE_dz_dz("ij")*dz_dx[i]("ir")*dz_dx[j]("js")+
                    dE_dN("i")*ddN_dx_dx[i][j]("irs"));
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

    TEST_dN(A2);
    TEST_dN(A);
    TEST_d_(d);
    TEST_d_(d2);
    TEST_d_(g);
    TEST_dN(C);
    TEST_d_(E);
    TEST_N(E,3);

    TEST_dN(dA2_dN);
    TEST_dN(dA_dN);
    TEST_d_(dd_dN);
    TEST_d_(dd_dz);
    TEST_d_(dd2_dN);
    TEST_d_(dd2_dz);
    TEST_d_(dg_dN);
    TEST_d_(dg_dz);
    TEST_dN(dC_dN);
    TEST_d_(dE_dN);
    TEST_d_(dE_dz);
    TEST_N(dE[0],3);
    TEST_N(dE[1],3);
    TEST_N(dE[2],3);

    return 0;
}
