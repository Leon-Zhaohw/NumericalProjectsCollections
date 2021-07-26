//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include "symbolic_test.h"

random_type rnd;

double eps=1e-6;

struct HELPER
{
    symbolic_tensor x,dx;
    symbolic_tensor r,a;
    symbolic_tensor v,dv,ddv;
    symbolic_tensor nv2,dnv2,ddnv2;
    symbolic_tensor nv,dnv,ddnv;
    symbolic_tensor inv_nv,dinv_nv,ddinv_nv;
    symbolic_tensor u,du,ddu;
    symbolic_tensor s,ds,dds;
    symbolic_tensor phi,g,H,A,T;
    symbolic_tensor psi,dpsi,ddpsi;

    void init()
    {
        x.set_zero({3});
        x.random(rnd,-1,0.9);

        r.set_zero({3});
        r.random(rnd,0.1,0.8);

        a.set_zero({});
        a.random(rnd,0.8,0.9);
    }

    void fill()
    {
        symbolic_tensor id,e;
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

// \uu &= \frac{\xx - \rr}{\|\xx - \rr\|} 
// \qq &= \rr + (\xx - \rr) s 
// h &= (\xx - \qq) \cdot \uu 
// a &= \prod_\alpha q_\alpha = \prod_\alpha (r_\alpha + (\sigma_\alpha - r_\alpha) s) 
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

    return 0;
}
