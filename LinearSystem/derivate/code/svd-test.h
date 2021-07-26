#ifndef __SVN_TEST__
#define __SVN_TEST__
#include "la.h"

template<class T,int d>
struct diff_helper
{
    mat<T,d> H;
    vec<T,d> A;
    vec<T,d> B;
    bool flip;
    // T_iikk = H_ik
    // T_ikik = (A_ik+B_ik)/2
    // T_ikki = (A_ik-B_ik)/2
    // If flip, then B_ik=-B_ki; otherwise, B_ik=B_ki

    mat<T,d> apply(const mat<T,d>& M) const
    {
        mat<T,d> R;
        R.make_zero();
        for(int i=0;i<3;i++)
            for(int k=0;k<3;k++)
                R(i,i)+=H(i,k)*M(k,k);

        for(int i=0;i<3;i++)
        {
            int j=(i+1)%3;
            int k=(i+2)%3;
            T a=A[i]/2,b=B[i]/2,c=flip?-b:b;
            R(j,k)+=(a+b)*M(j,k);
            R(j,k)+=(a-b)*M(k,j);
            R(k,j)+=(a+c)*M(k,j);
            R(k,j)+=(a-c)*M(j,k);
        }

        return R;
    }
};

template<class T>
mat<T,3> random_rotation(random_type& rand,T max_angle)
{
    vec<T,3> u;
    fill_random(rand,u,-max_angle,max_angle);
    T t=u.mag();
    mat<T,3> R=cpm(u/t);    
    return ((1-cos(t))*R+sin(t))*R+1;
}

template<class T,int d>
T sum(const vec<T,d>& v)
{
    T s=0;
    for(int i=0;i<d;i++)
        s+=v[i];
    return s;
}

template<class T,int d>
vec<T,d> ones()
{
    vec<T,d> v;
    for(int i=0;i<d;i++)
        v[i]=1;
    return v;
}

template<class T,int d>
vec<T,d> sum_over_sum(const vec<T,d>& A,const vec<T,d>& B)
{
    vec<T,d> R;
    for(int i=0;i<3;i++)
    {
        int j=(i+1)%3;
        int k=(i+2)%3;
        R[i]=(A[j]+A[k])/(B[j]+B[k]);
    }
    return R;
}

template<class T,int d>
vec<T,d> diff_over_sum(const vec<T,d>& A,const vec<T,d>& B)
{
    vec<T,d> R;
    for(int i=0;i<3;i++)
    {
        int j=(i+1)%3;
        int k=(i+2)%3;
        R[i]=(A[j]-A[k])/(B[j]+B[k]);
    }
    return R;
}

template<class T,int d>
vec<T,d> diff_over_diff(const vec<T,d>& A,const vec<T,d>& B)
{
    vec<T,d> R;
    for(int i=0;i<3;i++)
    {
        int j=(i+1)%3;
        int k=(i+2)%3;
        R[i]=(A[j]-A[k])/(B[j]-B[k]);
    }
    return R;
}

template<class T,int d,int n>
mat<T,n,d> diff_over_diff(const mat<T,n,d>& A,const vec<T,d>& B)
{
    mat<T,n,d> R;
    for(int m=0;m<n;m++)
        for(int i=0;i<3;i++)
        {
            int j=(i+1)%3;
            int k=(i+2)%3;
            R(m,i)=(A(m,j)-A(m,k))/(B[j]-B[k]);
        }
    return R;
}

template<class T,int d>
vec<T,d> exp(const vec<T,d>& v)
{
    vec<T,d> u;
    for(int i=0;i<d;i++)
        u[i]=exp(v[i]);
    return u;
}

template<class T,int d>
vec<T,d> ln(const vec<T,d>& v)
{
    vec<T,d> u;
    for(int i=0;i<d;i++)
        u[i]=ln(v[i]);
    return u;
}

template<class T,int d>
vec<T,d> cos(const vec<T,d>& v)
{
    vec<T,d> u;
    for(int i=0;i<d;i++)
        u[i]=cos(v[i]);
    return u;
}

template<class T,int d>
vec<T,d> sin(const vec<T,d>& v)
{
    vec<T,d> u;
    for(int i=0;i<d;i++)
        u[i]=sin(v[i]);
    return u;
}

template<class T,int d>
vec<T,d> divided_diff_exp(const vec<T,d>& v)
{
    vec<T,d> u;
    for(int i=0;i<d;i++)
        u[i]=exp(v[i]);
    return u;
}

template<class T,int d>
vec<T,d> divided_diff_ln(const vec<T,d>& v)
{
    vec<T,d> u;
    for(int i=0;i<d;i++)
        u[i]=ln(v[i]);
    return u;
}

template<class T,int d>
vec<T,d> divided_diff_cos(const vec<T,d>& v)
{
    vec<T,d> u;
    for(int i=0;i<d;i++)
        u[i]=cos(v[i]);
    return u;
}

template<class T,int d>
vec<T,d> divided_diff_sin(const vec<T,d>& v)
{
    vec<T,d> u;
    for(int i=0;i<d;i++)
        u[i]=sin(v[i]);
    return u;
}

#endif
