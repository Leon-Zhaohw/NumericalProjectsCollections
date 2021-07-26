#ifndef __LA__
#define __LA__
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <random>

template<class T, int nn>
struct vec
{
    static const int n=nn;
    T x[n+(n==0)];

    vec() {}
    explicit vec(T a) {assert(n==1);x[0]=a;}
    vec(T a, T b) {assert(n==2);x[0]=a;x[1]=b;}
    vec(T a, T b, T c) {assert(n==3);x[0]=a;x[1]=b;x[2]=c;}
    vec(T a, T b, T c, T d) {assert(n==4);x[0]=a;x[1]=b;x[2]=c;x[3]=d;}

    T& operator[](int i) {return x[i];}
    const T& operator[](int i) const {return x[i];}

    void make_zero() {for(int i=0; i<n; i++) x[i] = 0;}
    static vec zero() {vec r;r.make_zero();return r;}

    vec operator + (const vec& V) const {vec r;for(int i=0; i<n; i++) r.x[i] = x[i] + V.x[i]; return r;}
    vec operator - (const vec& V) const {vec r;for(int i=0; i<n; i++) r.x[i] = x[i] - V.x[i]; return r;}

    vec operator + () const {return *this;}
    vec operator - () const {vec r;for(int i=0; i<n; i++) r.x[i] = -x[i]; return r;}

    vec& operator += (const vec& V) {for(int i=0; i<n; i++) x[i] += V.x[i]; return *this;}
    vec& operator -= (const vec& V) {for(int i=0; i<n; i++) x[i] -= V.x[i]; return *this;}

    vec operator* (T a) const {vec v; for(int i=0; i<n; i++) v.x[i] = x[i] * a; return v;}
    vec operator/ (T a) const {vec v; for(int i=0; i<n; i++) v.x[i] = x[i] / a; return v;}

    vec operator * (const vec& V) const {vec r;for(int i=0; i<n; i++) r.x[i] = x[i] * V.x[i]; return r;}

    T dot(const vec& V) const {T r = 0; for(int i=0; i<n; i++) r += x[i] * V.x[i]; return r;}
    T mag2() const {return dot(*this);}
    T mag() const {return sqrt(mag2());}
    vec unit() const {return *this/mag();}
    vec cross(const vec& V) const {assert(nn==3); return vec(x[1]*V.x[2]-x[2]*V.x[1],x[2]*V.x[0]-x[0]*V.x[2],x[0]*V.x[1]-x[1]*V.x[0]);}

    void print(const char* s="",FILE* F=stdout) const {printf("[ ");for(int i=0; i<n; i++) printf("%.16g ",x[i]);printf("]%s",s);}
};

template<class T, int nn>
vec<T,nn> operator* (T a, const vec<T,nn>& v)
{return v*a;}

template<class T, int mm, int nn=mm>
struct mat
{
    static const int n=nn;
    static const int m=mm;
    T x[m*n+(m*n==0)];

    mat() {}
    explicit mat(T a) {assert(m*n==1);x[0]=a;}
    mat(T a, T b) {assert(m*n==2);x[0]=a;x[1]=b;}
    mat(T a, T b, T c) {assert(m*n==3);x[0]=a;x[1]=b;x[2]=c;}
    mat(T a, T b, T c, T d) {assert(m*n==4);x[0]=a;x[1]=b;x[2]=c;x[3]=d;}

    T& operator()(int i, int j) {return x[i*n+j];}
    const T& operator()(int i, int j) const {return x[i*n+j];}

    void make_zero() {int k=m*n;for(int i=0; i<k; i++) x[i] = 0;}
    static mat zero() {mat r;r.make_zero();return r;}
    static mat id() {assert(m==n);mat r;r.make_zero();int k=m*n;for(int i=0; i<k; i+=m+1) r.x[i]=1; return r;}

    mat operator + (const mat& M) const {mat r;int k=m*n;for(int i=0; i<k; i++) r.x[i] = x[i] + M.x[i]; return r;}
    mat operator - (const mat& M) const {mat r;int k=m*n;for(int i=0; i<k; i++) r.x[i] = x[i] - M.x[i]; return r;}

    mat operator * (T a) const {mat r;int k=m*n;for(int i=0; i<k; i++) r.x[i] = x[i] * a; return r;}
    mat operator / (T a) const {mat r;int k=m*n;for(int i=0; i<k; i++) r.x[i] = x[i] / a; return r;}

    mat operator + (T a) const {assert(m==n);mat r(*this);int k=m*n;for(int i=0; i<k; i+=m+1) r.x[i] += a; return r;}
    mat operator - (T a) const {assert(m==n);mat r(*this);int k=m*n;for(int i=0; i<k; i+=m+1) r.x[i] -= a; return r;}

    mat operator + () const {return *this;}
    mat operator - () const {mat r;int k=m*n;for(int i=0; i<k; i++) r.x[i] = -x[i]; return r;}

    template<int p>
    mat<T,m,p> operator * (const mat<T,n,p>& M) const
    {mat<T,m,p> r;r.make_zero();for(int i=0; i<m; i++) for(int j=0; j<n; j++) for(int k=0; k<p; k++) r(i,k) += (*this)(i,j) * M(j,k); return r;}

    mat& operator += (const mat& M) {int k=m*n;for(int i=0; i<k; i++) x[i] += M.x[i]; return *this;}
    mat& operator -= (const mat& M) {int k=m*n;for(int i=0; i<k; i++) x[i] -= M.x[i]; return *this;}

    vec<T,m> operator * (const vec<T,n>& v) const
    {vec<T,m> r;r.make_zero();for(int i=0; i<m; i++) for(int j=0; j<n; j++) r[i] += (*this)(i,j) * v[j]; return r;}

    void set_column(const vec<T,m>& v, int k) {for(int i=0; i<m; i++) (*this)(i,k) = v[i];}
    void get_column(vec<T,m>& v, int k) const {for(int i=0; i<m; i++) v[i] = (*this)(i,k);}
    vec<T,m> column(int k) const {vec<T,m> v;get_column(v,k);return v;}
    void set_row(const vec<T,n>& v, int k) {for(int i=0; i<n; i++) (*this)(k,i) = v[i];}
    void get_row(vec<T,n>& v, int k) const {for(int i=0; i<n; i++) v[i] = (*this)(k,i);}
    vec<T,n> row(int k) const {vec<T,n> v;get_row(v,k);return v;}

    mat<T,n,m> t() const
    {mat<T,n,m> M;for(int i=0;i<m;i++) for(int j=0;j<n;j++) M(j,i)=(*this)(i,j);return M;}

    T tr() const
    {T t=0;int k=m*n;for(int i=0; i<k; i+=m+1) t+=x[i];return t;}

    T double_contract(const mat& M) const
    {T t=0;int k=m*n;for(int i=0; i<k; i++) t+=x[i]*M.x[i];return t;}

    void print(const char* s="",FILE* F=stdout) const
    {
        printf("[ ");
        for(int i=0; i<m; i++)
        {
            if(i) printf("; ");
            for(int j=0; j<n; j++) printf("%.16g ",(*this)(i,j));
        }
        printf("]%s", s);
    }

    T frobenius_norm() const
    {
        return sqrt(double_contract(*this));
    }
};

template<class T, int m, int n>
mat<T,m,n> outer(const vec<T,m>& u, const vec<T,n>& v)
{
    mat<T,m,n> M;
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
            M(i,j)=u[i]*v[j];
    return M;
}

template<class T, int m, int n>
mat<T,m,n> operator* (T a, const mat<T,m,n>& M)
{return M*a;}

template<class T>
T det(const mat<T,0>& M) {return 0;}

template<class T>
T det(const mat<T,1>& M) {return M(0,0);}

template<class T>
T det(const mat<T,2>& M) {return M(0,0)*M(1,1)-M(1,0)*M(0,1);}

template<class T, int n>
T det(const mat<T,n>& M) {exit(1);return 0;}

template<class T>
void solve(const mat<T,0>& M, const vec<T,0>& v, vec<T,0>& x) {}

template<class T>
void solve(const mat<T,1>& M, const vec<T,1>& v, vec<T,1>& x) {x(0) = v(0) / M(0,0);}

template<class T>
void solve(const mat<T,2>& M, const vec<T,2>& v, vec<T,2>& x) {vec<T,2> t(M(1,1)*v(0) - M(0,1)*v(1), M(0,0) * v(1) - M(1,0) * v(0)); x = t / det(M);}

template<class T, int n>
void solve(const mat<T,n>& M, const vec<T,n>& v, vec<T,n>& x) {exit(1);}


typedef std::mt19937 random_type;

template<class T, int n>
void fill_random(random_type& r,vec<T,n>& v,T lo,T hi)
{
    std::uniform_real_distribution<> dist(lo,hi);
    for(int i=0;i<n;i++) v[i]=dist(r);
}

template<class T, int m, int n>
void fill_random(random_type& r,mat<T,m,n>& M,T lo,T hi)
{
    std::uniform_real_distribution<> dist(lo,hi);
    for(int i=0;i<m;i++) for(int j=0;j<n;j++) M(i,j)=dist(r);
}

template<class T>
void fill_random(random_type& r,T& x,T lo,T hi)
{
    std::uniform_real_distribution<> dist(lo,hi);
    x=dist(r);
}

template<class T>
mat<T,3> cpm(const vec<T,3>& v)
{
    mat<T,3> M;
    M(0,0)=0;
    M(0,1)=-v[2];
    M(0,2)=v[1];
    M(1,0)=v[2];
    M(1,1)=0;
    M(1,2)=-v[0];
    M(2,0)=-v[1];
    M(2,1)=v[0];
    M(2,2)=0;
    return M;
}

template<class T,int d>
mat<T,d> diag(const vec<T,d>& v)
{
    mat<T,d> M;
    M.make_zero();
    for(int i=0;i<d;i++)
        M(i,i)=v[i];
    return M;
}

#define PR(z_){printf("%s: ",#z_);for(auto q_:z_.x) printf("%.16g ",q_);printf("\n");}
#define PS(z_) printf("%s: %.16g\n",#z_,(z_));

#endif
