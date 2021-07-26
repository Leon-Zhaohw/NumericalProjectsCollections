//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include "symbolic_test.h"

random_type rnd;

double eps=1e-6;

struct HELPER
{
    symbolic_tensor D,dD;
    symbolic_tensor F_bar,dF_bar;
    symbolic_tensor E_bar,dE_bar,ddE_bar;
    symbolic_tensor W,dW,ddW;

    symbolic_tensor mu,lambda,A_X,R_X_inv;

    void init()
    {
        D.set_zero({3,2});
        D.random(rnd,-1,1);

        mu.set_zero({});
        mu.random(rnd,0.1,1);

        lambda.set_zero({});
        lambda.random(rnd,0.1,1);

        A_X.set_zero({});
        A_X.random(rnd,0.1,1);

        R_X_inv.set_zero({2,2});
        R_X_inv.random(rnd,-1,1);
    }

    void fill()
    {
        symbolic_tensor id2;
        id2.set_id(2);
        symbolic_tensor id3;
        id3.set_id(3);

        F_bar.set("ij",D("ik")*R_X_inv("kj"));
        E_bar.set("ij",0.5*(F_bar("ki")*F_bar("kj") - id2("ij")));
        W.set("",A_X*(lambda/2*E_bar("ii")*E_bar("jj") + mu*E_bar("ij")*E_bar("ji")));
        W.set("",A_X*(lambda/2*E_bar("ii")*E_bar("jj") + mu*E_bar("ij")*E_bar("ji")));

        dD.set("ijrs",id3("ir")*id2("js"));
        dF_bar.set("ijrs",dD("ikrs")*R_X_inv("kj"));
        dE_bar.set("ijrs",0.5*(dF_bar("kirs")*F_bar("kj")+F_bar("ki")*dF_bar("kjrs")));
        dW.set("rs",A_X*(lambda*dE_bar("iirs")*E_bar("jj") + 2*mu*dE_bar("ijrs")*E_bar("ji")));

        ddE_bar.set("ijrsuv",0.5*(dF_bar("kirs")*dF_bar("kjuv")+dF_bar("kiuv")*dF_bar("kjrs")));
        ddW.set("rsuv",A_X*(lambda*ddE_bar("iirsuv")*E_bar("jj") + lambda*dE_bar("iirs")*dE_bar("jjuv") + 2*mu*ddE_bar("ijrsuv")*E_bar("ji") + 2*mu*dE_bar("ijrs")*dE_bar("jiuv")));
    }
};

int main(int argc, char* argv[])
{
    HELPER z0;
    z0.init();

    symbolic_tensor dx({3,2});
    dx.random(rnd,-eps,eps);

    HELPER z1=z0;
    z1.D.set("rs",z0.D("rs")+dx("rs"));

    z0.fill();
    z1.fill();

    TEST(D);
    TEST(F_bar);
    TEST(E_bar);
    TEST(W);
    TEST(dE_bar);
    TEST(dW);

    return 0;
}
