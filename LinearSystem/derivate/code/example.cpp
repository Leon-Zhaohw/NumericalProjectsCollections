#include <cstdio>
#include <cmath>
#include <algorithm>
using std::max;

void compute_hand(double x, double& z, double& dz, double& ddz)
{
    double a = cos(x);
    double da = -sin(x);
    double dda = -a;
    double b = 1 / a;
    double g = b * b;
    double db = -g * da;
    double dg = 2 * b * db;
    double ddb = -dg * da-g * dda;
    double c = x * b;
    double dc = x * db + b;
    double ddc = x * ddb + 2 * db;
    double d = 1.5 * c;
    double dd = 1.5 * dc;
    double ddd = 1.5 * ddc;
    double e = d * d - c;
    double de = 2 * d * dd - dc;
    double dde = 2 * (dd * dd + d * ddd) - ddc;
    double f = sqrt(e);
    double h = 1 / f;
    double df = 0.5 * de * h;
    double dh = -h * h * df;
    double ddf = 0.5 * (dde * h + de * dh);
    z = d + f;
    dz = dd + df;
    ddz = ddd + ddf;
}

void compute_maple(double x, double& z, double& dz, double& ddz)
{
    double a = cos(x);
    double b = 4 * a - 9 * x;
    double c = sqrt(-b * x);
    double d = sin(x);
    double e = c * d;
    double f = x * x;
    double g = a * a;
    double h = 1 / c;
    double j = g * g;
    double k = g * a;
    double m = f * x;
    double n = c * f;
    double p = f * f;
    double q = c * m;
    double r = m * a;
    double s = 108 * f * d * g - 8 * x * d * k - 54 * e * f * a
        + 24 * e * x * g - 162 * d * r - 12 * f * g + 4 * f * j
        + 24 * a * n + 81 * p * g + 27 * g * q - 54 * m * k
        - 12 * k * n + 4 * j - 162 * p - 54 * q + 108 * r;
    z=1 / a * (3 * x + c) / 2;
    dz=-1 / g * h * (2 * x * d * a - 9 * f * d - 3 * x * e
        - 3 * a * c - 9 * x * a + 2 * g) / 2;
    ddz=1 / k * h / b / x * s / 2;
}
int main()
{
    double x=1.341,dx=1e-6;

    {
        double z0,dz0,ddz0;
        double z1,dz1,ddz1;
        compute_hand(x-dx,z0,dz0,ddz0);
        compute_hand(x+dx,z1,dz1,ddz1);
        double a=(z1-z0)/dx;
        double b=dz0+dz1;
        double da=(dz1-dz0)/dx;
        double db=ddz0+ddz1;
        printf("%g %g %g\n",a,b,fabs(a-b)/max(max(fabs(a),fabs(b)),1e-30));
        printf("%g %g %g\n",da,db,fabs(da-db)/max(max(fabs(da),fabs(db)),1e-30));
    }
    {
        double z0,dz0,ddz0;
        double z1,dz1,ddz1;
        compute_maple(x-dx,z0,dz0,ddz0);
        compute_maple(x+dx,z1,dz1,ddz1);
        double a=(z1-z0)/dx;
        double b=dz0+dz1;
        double da=(dz1-dz0)/dx;
        double db=ddz0+ddz1;
        printf("%g %g %g\n",a,b,fabs(a-b)/max(max(fabs(a),fabs(b)),1e-30));
        printf("%g %g %g\n",da,db,fabs(da-db)/max(max(fabs(da),fabs(db)),1e-30));
    }

    return 0;
}
