//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include "symbolic_test.h"
#include <gsl/gsl_poly.h>

random_type rnd;

double eps=1e-6;

void Psi(const symbolic_tensor& x,symbolic_tensor& E,symbolic_tensor& dE,
    symbolic_tensor& ddE,symbolic_tensor& dddE,symbolic_tensor& ddddE)
{
    symbolic_tensor id;
    id.set_id(3);

    symbolic_tensor a("",x("i")*x("i")),c,s;
    c.set_zero({});
    s.set_zero({});
    E.set_zero({});
    c.x[0]=cos(a.x[0]);
    s.x[0]=sin(a.x[0]);

    E.set("",s);
    dE.set("r",c*2*x("r"));
    ddE.set("rs",-s*4*x("r")*x("s")+c*2*id("rs"));
    dddE.set("rst",-c*8*x("r")*x("s")*x("t")-s*4*id("rt")*x("s")-s*4*x("r")*id("st")-s*4*id("rs")*x("t"));
    ddddE.set("rstu",
        s*16*x("r")*x("s")*x("t")*x("u")-c*8*id("ru")*x("s")*x("t")-c*8*x("r")*id("su")*x("t")-c*8*x("r")*x("s")*id("tu")
        -c*8*id("rt")*x("s")*x("u")-s*4*id("rt")*id("su")
        -c*8*x("r")*id("st")*x("u")-s*4*id("ru")*id("st")
        -c*8*id("rs")*x("t")*x("u")-s*4*id("rs")*id("tu"));
}


struct HELPER
{
    symbolic_tensor x,dx;
    symbolic_tensor r,a;
    symbolic_tensor v,dv;
    symbolic_tensor nv2,dnv2,ddnv2;
    symbolic_tensor nv,dnv,ddnv;
    symbolic_tensor inv_nv,dinv_nv,ddinv_nv;
    symbolic_tensor u,du,ddu;
    symbolic_tensor C0,dC0,ddC0;
    symbolic_tensor C1,dC1,ddC1;
    symbolic_tensor C2,dC2,ddC2;
    symbolic_tensor C3,dC3,ddC3;
    symbolic_tensor s,ds,dds;
    symbolic_tensor q,dq,ddq;
    symbolic_tensor h,dh,ddh;
    symbolic_tensor phi,g,H,A,T;
    symbolic_tensor psi,dpsi,ddpsi;
    symbolic_tensor Z,dZ,ddZ,dddZ,ddddZ;

    void init()
    {
        x.set_zero({3});
        x.random(rnd,-1,0.9);

        r.set_zero({3});
        r.random(rnd,1.0,1.1);

        a.set_zero({});
        a.random(rnd,0.8,0.9);
    }

    void fill()
    {
        symbolic_tensor id;
        id.set_id(3);

        v.set("i",x("i")-r("i"));
        nv2.set("",v("i")*v("i"));
        nv.set("",sqrt(nv2));
        inv_nv.set("",1/nv);
        u.set("i",inv_nv*v("i"));

        dv.set("ir",id("ir"));
        dnv2.set("r",2*v("i")*dv("ir"));
        dnv.set("r",dnv2("r")*inv_nv/2);
        dinv_nv.set("r",-inv_nv*inv_nv*dnv("r"));
        du.set("ir",inv_nv*dv("ir")+dinv_nv("r")*v("i"));

        ddnv2.set("rs",2*dv("is")*dv("ir"));
        ddnv.set("rs",ddnv2("rs")*inv_nv/2+dnv2("r")*dinv_nv("s")/2);
        ddinv_nv.set("rs",-2*dinv_nv("s")*inv_nv*dnv("r")-inv_nv*inv_nv*ddnv("rs"));
        ddu.set("irs",dinv_nv("s")*dv("ir")+ddinv_nv("rs")*v("i")+dinv_nv("r")*dv("is"));

        C0.set_zero({});
        C1.set_zero({});
        C2.set_zero({});
        C3.set_zero({});
        C0.x[0]=r(0)*r(1)*r(2)-a.x[0];
        C1.x[0]=r(0)*r(1)*v(2)+r(0)*v(1)*r(2)+v(0)*r(1)*r(2);
        C2.x[0]=r(0)*v(1)*v(2)+v(0)*r(1)*v(2)+v(0)*v(1)*r(2);
        C3.x[0]=v(0)*v(1)*v(2);

        dC0.set_zero({3});
        dC1.set_zero({3});
        dC2.set_zero({3});
        dC3.set_zero({3});
        dC1(0)=r(1)*r(2);
        dC1(1)=r(0)*r(2);
        dC1(2)=r(0)*r(1);
        dC2(0)=v(1)*r(2)+r(1)*v(2);
        dC2(1)=v(0)*r(2)+r(0)*v(2);
        dC2(2)=v(0)*r(1)+r(0)*v(1);
        dC3(0)=v(1)*v(2);
        dC3(1)=v(0)*v(2);
        dC3(2)=v(0)*v(1);

        ddC0.set_zero({3,3});
        ddC1.set_zero({3,3});
        ddC2.set_zero({3,3});
        ddC3.set_zero({3,3});
        ddC2(0,1)=r(2);
        ddC2(0,2)=r(1);
        ddC2(1,0)=r(2);
        ddC2(1,2)=r(0);
        ddC2(2,0)=r(1);
        ddC2(2,1)=r(0);
        ddC3(0,1)=v(2);
        ddC3(0,2)=v(1);
        ddC3(1,0)=v(2);
        ddC3(1,2)=v(0);
        ddC3(2,0)=v(1);
        ddC3(2,1)=v(0);

        double x0=-1,x1=-1,x2=-1;
        int n=gsl_poly_solve_cubic(C2.x[0]/C3.x[0],C1.x[0]/C3.x[0],C0.x[0]/C3.x[0],&x0,&x1,&x2);
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
        s.set_zero({});
        s.x[0]=x0;

        symbolic_tensor Tc("r",(((dC3("r")*s+dC2("r"))*s+dC1("r"))*s+dC0("r")));
        symbolic_tensor Ts("",(3*C3*s+2*C2)*s+C1);
        symbolic_tensor Tss("",6*C3*s+2*C2);
        symbolic_tensor Tcs("r",(3*dC3("r")*s+2*dC2("r"))*s+dC1("r"));
        symbolic_tensor Tcc("rs",(((ddC3("rs")*s+ddC2("rs"))*s+ddC1("rs"))*s+ddC0("rs")));

        // T = 0
        // Tc + Ts*ds = 0
        // Tcc + Tcs*ds + Tcs*ds + Tss*ds*ds + Ts*dds = 0
        ds.set("r",-Tc("r")/Ts);
        dds.set("rs",(-Tss*ds("r")*ds("s")-Tcs("r")*ds("s")-Tcs("s")*ds("r")-Tcc("rs"))/Ts);

        q.set("i",r("i")+v("i")*s);
        dq.set("ir",dv("ir")*s+v("i")*ds("r"));
        ddq.set("irs",dv("ir")*ds("s")+dv("is")*ds("r")+v("i")*dds("rs"));

        h.set("",(x("i")-q("i"))*u("i"));
        dh.set("r",(id("ir")-dq("ir"))*u("i")+(x("i")-q("i"))*du("ir"));
        ddh.set("rs",(-ddq("irs"))*u("i")+(id("ir")-dq("ir"))*du("is")+(id("is")-dq("is"))*du("ir")+(x("i")-q("i"))*ddu("irs"));

        Psi(x,Z,dZ,ddZ,dddZ,ddddZ);

// \phi &= \Psi(\qq) 
// g_i &= \Psi_{,i}(\qq) 
// H_{ij} &= \Psi_{,ij}(\qq) 
// \psi &= \phi + h \gg\cdot\uu + \frac{1}{2} h^2 \uu^T\HH\uu
    }
};

int main(int argc, char* argv[])
{
    HELPER z0;
    z0.init();

    symbolic_tensor dx;
    dx.set_zero({3});
    dx.random(rnd,-eps,eps);

    HELPER z1=z0;
    z1.x.set("i",z0.x("i")+dx("i"));

    z0.fill();
    z1.fill();

    TEST(v);
    TEST(nv2);
    TEST(nv);
    TEST(inv_nv);
    TEST(u);
    TEST(C0);
    TEST(C1);
    TEST(C2);
    TEST(C3);
    TEST(s);
    TEST(q);
    TEST(h);

    TEST(dnv2);
    TEST(dnv);
    TEST(dinv_nv);
    TEST(du);
    TEST(dC0);
    TEST(dC1);
    TEST(dC2);
    TEST(dC3);
    TEST(ds);
    TEST(dq);
    TEST(dh);

    TEST(Z);
    TEST(dZ);
    TEST(ddZ);
    TEST(dddZ);

    return 0;
}
