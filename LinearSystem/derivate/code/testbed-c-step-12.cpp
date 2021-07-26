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
    symbolic_tensor C0;
    symbolic_tensor C1,dC1;
    symbolic_tensor C2,dC2,ddC2;
    symbolic_tensor C3,dC3,ddC3;
    symbolic_tensor s,ds,dds;
    symbolic_tensor q,dq,ddq;
    symbolic_tensor h,dh,ddh;
    symbolic_tensor psi,dpsi,ddpsi;
    symbolic_tensor Z,dZ,ddZ,dddZ,ddddZ;
    symbolic_tensor e,de,dde;
    symbolic_tensor f,df,ddf;
    symbolic_tensor g,dg,ddg;
    symbolic_tensor k,dk,ddk;
    symbolic_tensor qZ,dqZ,ddqZ;
    symbolic_tensor qdZ,dqdZ,ddqdZ;
    symbolic_tensor qddZ,dqddZ,ddqddZ;
    symbolic_tensor qdddZ,dqdddZ;

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

        symbolic_tensor M({3,3,3});
        M(0,1,2)=1;
        M(1,2,0)=1;
        M(2,0,1)=1;
        M(0,2,1)=1;
        M(1,0,2)=1;
        M(2,1,0)=1;

        v.set("i",x("i")-r("i"));
        nv2.set("",v("i")*v("i"));
        nv.set("",sqrt(nv2));
        inv_nv.set("",1/nv);

        dnv2.set("r",2*v("i")*id("ir"));
        dnv.set("r",dnv2("r")*inv_nv/2);
        dinv_nv.set("r",-inv_nv*inv_nv*dnv("r"));

        ddnv2.set("rs",2*id("is")*id("ir"));
        ddnv.set("rs",ddnv2("rs")*inv_nv/2+dnv2("r")*dinv_nv("s")/2);
        ddinv_nv.set("rs",-2*dinv_nv("s")*inv_nv*dnv("r")-inv_nv*inv_nv*ddnv("rs"));

        C0.set("",M("ijk")*r("i")*r("j")*r("k")/6-a);
        C1.set("",M("ijk")*r("i")*r("j")*v("k")/2);
        C2.set("",M("ijk")*r("i")*v("j")*v("k")/2);
        C3.set("",M("ijk")*v("i")*v("j")*v("k")/6);
        dC1.set("r",M("ijr")*r("i")*r("j")/2);
        dC2.set("r",M("ijr")*r("i")*v("j"));
        dC3.set("r",M("ijr")*v("i")*v("j")/2);
        ddC2.set("rs",M("irs")*r("i"));
        ddC3.set("rs",M("irs")*v("i"));

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

        q.set("i",r("i")+v("i")*s);

        symbolic_tensor Tc("r",M("ijr")*q("i")*q("j")*s/2);
        symbolic_tensor Ts("",M("ijk")*v("k")*q("i")*q("j")/2);
        symbolic_tensor Tss("",M("ijk")*q("i")*v("j")*v("k"));
        symbolic_tensor Tcs("r",M("ijr")*q("i")*(3*v("j")*s+r("j"))/2);
        symbolic_tensor Tcc("rs",M("irs")*q("i")*s*s);

        // T = 0
        // Tc + Ts*ds = 0
        // Tcc + Tcs*ds + Tcs*ds + Tss*ds*ds + Ts*dds = 0
        ds.set("r",-Tc("r")/Ts);
        dds.set("rs",(-Tss*ds("r")*ds("s")-Tcs("r")*ds("s")-Tcs("s")*ds("r")-Tcc("rs"))/Ts);

        dq.set("ir",id("ir")*s+v("i")*ds("r"));
        ddq.set("irs",id("ir")*ds("s")+id("is")*ds("r")+v("i")*dds("rs"));

        h.set("",nv*(1-s));
        dh.set("r",dnv("r")*(1-s)-nv*ds("r"));
        ddh.set("rs",ddnv("rs")*(1-s)-dnv("r")*ds("s")-dnv("s")*ds("r")-nv*dds("rs"));

        Psi(q,Z,dZ,ddZ,dddZ,ddddZ);

        k.set("i",(1-s)*v("i"));
        dk.set("ir",-ds("r")*v("i")+(1-s)*id("ir"));
        ddk.set("irs",-dds("rs")*v("i")-ds("r")*id("is")-ds("s")*id("ir"));

        symbolic_tensor kdZ("",dZ("i")*k("i"));
        symbolic_tensor kddZ("m",ddZ("mi")*k("i"));
        symbolic_tensor kdddZ("mn",dddZ("mni")*k("i"));
        symbolic_tensor kddddZ("mnp",ddddZ("mnpi")*k("i"));
        symbolic_tensor kkddZ("",ddZ("ij")*k("i")*k("j"));
        symbolic_tensor kkdddZ("m",dddZ("mij")*k("i")*k("j"));
        symbolic_tensor kkddddZ("mn",ddddZ("mnij")*k("i")*k("j"));
        symbolic_tensor kkkdddZ("",dddZ("ijk")*k("i")*k("j")*k("k"));
        symbolic_tensor kkkddddZ("m",ddddZ("mijk")*k("i")*k("j")*k("k"));
        symbolic_tensor dkkkddddZ("",ddddZ("ijkl")*k("i")*k("j")*k("k")*k("l"));

        symbolic_tensor vdZ("",dZ("i")*v("i"));
        symbolic_tensor vddZ("m",ddZ("mi")*v("i"));
        symbolic_tensor vdddZ("mn",dddZ("mni")*v("i"));
        symbolic_tensor vddddZ("mnp",ddddZ("mnpi")*v("i"));
        symbolic_tensor vvddZ("",ddZ("ij")*v("i")*v("j"));
        symbolic_tensor vvdddZ("m",dddZ("mij")*v("i")*v("j"));
        symbolic_tensor vvddddZ("mn",ddddZ("mnij")*v("i")*v("j"));
        symbolic_tensor vvvdddZ("",dddZ("ijk")*v("i")*v("j")*v("k"));
        symbolic_tensor vvvddddZ("m",ddddZ("mijk")*v("i")*v("j")*v("k"));
        symbolic_tensor dvvvddddZ("",ddddZ("ijkl")*v("i")*v("j")*v("k")*v("l"));

        qZ.set("",Z);
        qdZ.set("m",dZ("m"));
        qddZ.set("mn",ddZ("mn"));
        qdddZ.set("mnp",dddZ("mnp"));

        dqZ.set("r",qdZ("i")*dq("ir"));
        dqdZ.set("mr",qddZ("mi")*dq("ir"));
        dqddZ.set("mnr",qdddZ("mni")*dq("ir"));
        dqdddZ.set("mnpr",ddddZ("mnpi")*dq("ir"));

        ddqZ.set("rs",dqdZ("is")*dq("ir")+qdZ("i")*ddq("irs"));
        ddqdZ.set("mrs",dqddZ("mis")*dq("ir")+qddZ("mi")*ddq("irs"));
        ddqddZ.set("mnrs",dqdddZ("mnis")*dq("ir")+qdddZ("mni")*ddq("irs"));

        dqZ.set("r",dZ("r")*s+vdZ*ds("r"));
        ddqZ.set("rs",ddZ("rs")*s*s+(vddZ("r")*s+vvddZ*ds("r")+dZ("r"))*ds("s")+vddZ("s")*ds("r")*s+dZ("s")*ds("r")+vdZ*dds("rs"));


        e.set("",qdZ("i")*k("i"));
        de.set("r",dqdZ("ir")*k("i")+qdZ("i")*dk("ir"));
        dde.set("rs",ddqdZ("irs")*k("i")+dqdZ("ir")*dk("is")+dqdZ("is")*dk("ir")+qdZ("i")*ddk("irs"));

        f.set("",qddZ("ij")*k("i")*k("j")/2);
        df.set("r",dqddZ("ijr")*k("i")*k("j")/2+qddZ("ij")*dk("ir")*k("j")/2+qddZ("ij")*k("i")*dk("jr")/2);
        ddf.set("rs",
            ddqddZ("ijrs")*k("i")*k("j")/2
            +dqddZ("ijr")*dk("is")*k("j")/2
            +dqddZ("ijr")*k("i")*dk("js")/2
            +dqddZ("ijs")*dk("ir")*k("j")/2
            +qddZ("ij")*ddk("irs")*k("j")/2
            +qddZ("ij")*dk("ir")*dk("js")/2
            +dqddZ("ijs")*k("i")*dk("jr")/2
            +qddZ("ij")*dk("is")*dk("jr")/2
            +qddZ("ij")*k("i")*ddk("jrs")/2);

        psi.set("",qZ+e+f);

        dpsi.set("r",dZ("r")+vddZ("r")*(1-s)+vvdddZ("r")*(1-s)*(1-s)*s/2+(1-s)*(1-s)*vvvdddZ*ds("r")/2);

        ddpsi.set("rs",
            ddZ("rs")*(s*s-s+1)
            +vddZ("r")*ds("s")*s
            +vddZ("s")*ds("r")*(s-1)
            +vvddZ*ds("r")*ds("s")
            +vvddddZ("ij")*dq("js")*dq("ir")*(1-s)*(1-s)/2
            +vvdddZ("r")*ds("s")*(1-s)*(1-s)/2
            +vvdddZ("s")*ds("r")*(1-s)*(1-s)/2
            +vvvdddZ*dds("rs")*(1-s)*(1-s)/2
            +vdddZ("ij")*dq("jr")*(-ds("s")*v("i")+(1-s)*id("is"))*(1-s)
            +vdddZ("ir")*dq("is")*(1-s)
            +ddZ("ij")*dq("jr")*(-ds("s")*v("i")+(1-s)*id("is"))


//        dq.set("ir",id("ir")*s+v("i")*ds("r"));
//        dk.set("ir",-ds("r")*v("i")+(1-s)*id("ir"));
//            dq("ir")+dk("ir")=id("ir")


);

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

    TEST(nv2);
    TEST(nv);
    TEST(inv_nv);
    TEST(C1);
    TEST(C2);
    TEST(C3);
    TEST(s);
    TEST(q);
    TEST(h);
    TEST(e);
    TEST(f);

    TEST(qZ);
    TEST(dqZ);
    TEST(qdZ);
    TEST(dqdZ);
    TEST(qddZ);
    TEST(dqddZ);

    TEST(psi);

    TEST(dnv2);
    TEST(dnv);
    TEST(dinv_nv);
    TEST(dC2);
    TEST(dC3);
    TEST(ds);
    TEST(dq);
    TEST(dh);
    TEST(de);
    TEST(df);

    TEST(dpsi);

    return 0;
}
