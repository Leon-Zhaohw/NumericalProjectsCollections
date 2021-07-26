#include <cstdio>
#include <cmath>
#include <algorithm>
using std::max;

struct helper
{
    double a,da,dda;
    double b,g,db,dg,ddb;
    double c,dc,ddc;
    double d,dd,ddd;
    double e,de,dde;
    double f,h,df,dh,ddf;
    double z,dz,ddz;

    void compute_hand(double x)
    {
        a = cos(x);
        da = -sin(x);
        dda = -a;
        b = 1 / a;
        g = b * b;
        db = -g * da;
        dg = 2 * b * db;
        ddb = -dg * da-g * dda;
        c = x * b;
        dc = x * db + b;
        ddc = x * ddb + 2 * db;
        d = 1.5 * c;
        dd = 1.5 * dc;
        ddd = 1.5 * ddc;
        e = d * d - c;
        de = 2 * d * dd - dc;
        dde = 2 * (dd * dd + d * ddd) - ddc;
        f = sqrt(e);
        h = 1 / f;
        df = 0.5 * de * h;
        dh = -h * h * df;
        ddf = 0.5 * (dde * h + de * dh);
        z = d + f;
        dz = dd + df;
        ddz = ddd + ddf;
    }
};


void test(const char* name,double z0,double z1,double dz0,double dz1,double dx)
{
    double a=(z1-z0)/dx;
    double b=dz0+dz1;
    printf("%2s %12g %12g %12g\n",name,a,b,fabs(a-b)/max(max(fabs(a),fabs(b)),1e-30));
}

#define TEST(z) test(#z,A.z,B.z,A.d##z,B.d##z,dx)

int main()
{
    double x=1.341,dx=1e-6;

    helper A,B;
    A.compute_hand(x-dx);
    B.compute_hand(x+dx);

    TEST(a);
    TEST(da);
    TEST(b);
    TEST(g);
    TEST(db);
    TEST(c);
    TEST(dc);
    TEST(d);
    TEST(dd);
    TEST(e);
    TEST(de);
    TEST(f);
    TEST(h);
    TEST(df);
    TEST(z);
    TEST(dz);

    return 0;
}
