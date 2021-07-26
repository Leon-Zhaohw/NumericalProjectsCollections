//#####################################################################
// Copyright 2015, Daniel Ram, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class symbolic_tensor
//##################################################################### 
#ifndef __symbolic_tensor__
#define __symbolic_tensor__

#include <cassert>
#include <cmath>
#include <functional>
#include <random>
#include <string>
#include <vector>
using namespace std;

typedef mt19937 random_type;

inline double get_uniform(random_type& rand,double lo,double hi)
{return uniform_real_distribution<>(lo,hi)(rand);}

inline double norm(const vector<double>& v)
{
    double r=0;
    for(auto it:v) r+=it*it;
    return sqrt(r);
}

class symbolic_tensor
{
public:
    vector<double> x;
    vector<size_t> size,strides;
    string indices;

    explicit symbolic_tensor();
    explicit symbolic_tensor(const vector<size_t>& size);
    symbolic_tensor(const string& indices,const symbolic_tensor& a);

    void compute_strides();
    symbolic_tensor operator()(const string& new_indices) const;

    void random(random_type& rand,double lo,double hi);

    template<class... I>
    size_t index_helper(size_t k,size_t i,I&&... j) const
    {return index_helper(k+1,j...)+i*strides[k];}

    template<class... I>
    size_t index_helper(size_t k) const
    {return 0;}

    template<class... I>
    size_t index(I&&... j) const
    {return index_helper(0,j...);}

    template<class... I>
    double& operator()(size_t i,I&&... j)
    {assert(size.size()==sizeof...(I)+1);return x[index(i,j...)];}

    template<class... I>
    double operator()(size_t i,I&&... j) const
    {assert(size.size()==sizeof...(I)+1);return x[index(i,j...)];}

    void set_zero(const vector<size_t>& new_size);
    void set(const string& ind,const symbolic_tensor& a);
    void set_id(size_t new_size);
    void set_perm();
};

inline symbolic_tensor operator+ (const symbolic_tensor& a){return a;}
symbolic_tensor operator+ (const symbolic_tensor& a,double b);
inline symbolic_tensor operator+ (double b,const symbolic_tensor& a){return a+b;}
symbolic_tensor operator+ (const symbolic_tensor& a,const symbolic_tensor& b);
inline symbolic_tensor operator- (const symbolic_tensor& a,double b){return a+-b;}
symbolic_tensor operator- (const symbolic_tensor& a,const symbolic_tensor& b);
symbolic_tensor operator* (const symbolic_tensor& a,double b);
inline symbolic_tensor operator* (double b,const symbolic_tensor& a){return a*b;}
inline symbolic_tensor operator- (const symbolic_tensor& a){return a*-1;}
inline symbolic_tensor operator- (double b,const symbolic_tensor& a){return -a+b;}
symbolic_tensor operator* (const symbolic_tensor& a,const symbolic_tensor& b);
inline symbolic_tensor operator/ (const symbolic_tensor& a,double b){return a*(1./b);}
symbolic_tensor operator/ (double b,const symbolic_tensor& a);
symbolic_tensor operator/ (const symbolic_tensor& a,const symbolic_tensor& b);
symbolic_tensor scalar_op(const symbolic_tensor& a,function<double(double)> f);
inline symbolic_tensor sqrt(const symbolic_tensor& a){return scalar_op(a,(double(*)(double))sqrt);}
#endif
