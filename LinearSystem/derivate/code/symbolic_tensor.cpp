//#####################################################################
// Copyright 2015, Daniel Ram, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <algorithm>
#include <cctype>
#include "symbolic_tensor.h"
#include <map>
bool same_indices(string a,string b)
{
    sort(a.begin(),a.end());
    sort(b.begin(),b.end());
    return a==b;
}

void assert_unindexed(const symbolic_tensor& a)
{
    assert(a.indices.size() == 0);
}

void assert_indexed(const symbolic_tensor& a)
{
    assert(a.indices.size() == a.size.size());
}

symbolic_tensor::symbolic_tensor()
{
    assert_unindexed(*this);
}
symbolic_tensor::symbolic_tensor(const vector<size_t>& size)
{
    assert_unindexed(*this);
    set_zero(size);
}
symbolic_tensor::symbolic_tensor(const std::string& indices,const symbolic_tensor& a)
{
    assert_unindexed(*this);
    assert_indexed(a);
    set(indices,a);
}
void symbolic_tensor::set_zero(const vector<size_t>& new_size)
{
    assert_unindexed(*this);
    indices="";
    size=new_size;
    compute_strides();
}
static void reordered_strides(vector<size_t>& new_strides,const std::string& new_indices,const vector<size_t>& old_strides,const std::string& old_indices)
{
    new_strides.resize(new_indices.size());
    for(size_t i=0;i<new_indices.size();i++){
        size_t f=old_indices.find(new_indices[i]);
        new_strides[i]=f!=string::npos?old_strides[f]:0;}
}
static void recurse_assign(symbolic_tensor& r,const symbolic_tensor& a,const vector<size_t>& as,size_t i,size_t ir,size_t ia)
{
    if(i>=r.size.size()){
        r.x[ir]=a.x[ia];
        return;}
    for(size_t s=0;s<r.size[i];s++)
        recurse_assign(r,a,as,i+1,ir+s*r.strides[i],ia+s*as[i]);
}
void symbolic_tensor::set(const std::string& ind,const symbolic_tensor& a)
{
    assert_unindexed(*this);
    assert_indexed(a);
    indices="";
    assert(same_indices(ind,a.indices));
    size.resize(ind.size());
    for(size_t i=0;i<ind.size();i++)
        size[i]=a.size[a.indices.find(ind[i])];
    compute_strides();
    vector<size_t> as;
    reordered_strides(as,ind,a.strides,a.indices);
    recurse_assign(*this,a,as,0,0,0);
}
void symbolic_tensor::set_id(size_t new_size)
{
    assert_unindexed(*this);
    size={new_size,new_size};
    compute_strides();
    for(size_t i=0;i<new_size;i++)
        (*this)(i,i)=1;
}
void symbolic_tensor::set_perm()
{
    assert_unindexed(*this);
    size={3,3,3};
    compute_strides();
    (*this)(0,1,2)=1;
    (*this)(2,0,1)=1;
    (*this)(1,2,0)=1;
    (*this)(0,2,1)=-1;
    (*this)(2,1,0)=-1;
    (*this)(1,0,2)=-1;
}
void symbolic_tensor::compute_strides()
{
    if(size.size()){
        strides.resize(size.size());
        strides.back()=1;
        for(int i=strides.size()-1;i>0;i--)
            strides[i-1]=strides[i]*size[i];
        x.resize(strides[0]*size[0]);}
    else x.resize(1);
    fill(x.begin(),x.end(),0);
}
template<class f> static void recurse_op(symbolic_tensor& r,const symbolic_tensor& a,const symbolic_tensor& b,
    const vector<size_t>& as,const vector<size_t>& bs,f func,size_t i,size_t ir,size_t ia,size_t ib)
{
    if(i>=r.size.size()){
        r.x[ir]=func(a.x[ia],b.x[ib]);
        return;}
    for(size_t s=0;s<r.size[i];s++)
        recurse_op(r,a,b,as,bs,func,i+1,ir+s*r.strides[i],ia+s*as[i],ib+s*bs[i]);
}
symbolic_tensor operator+ (const symbolic_tensor& a,const symbolic_tensor& b)
{
    assert_indexed(a);
    assert_indexed(b);
    symbolic_tensor r(a);
    vector<size_t> bs;
    reordered_strides(bs,r.indices,b.strides,b.indices);
    recurse_op(r,a,b,a.strides,bs,[](double u,double v){return u+v;},0,0,0,0);
    return r;
}
symbolic_tensor operator- (const symbolic_tensor& a,const symbolic_tensor& b)
{
    assert_indexed(a);
    assert_indexed(b);
    symbolic_tensor r(a);
    vector<size_t> bs;
    reordered_strides(bs,r.indices,b.strides,b.indices);
    recurse_op(r,a,b,a.strides,bs,[](double u,double v){return u-v;},0,0,0,0);
    return r;
}
symbolic_tensor operator* (const symbolic_tensor& a,double b)
{
    assert_indexed(a);
    symbolic_tensor r(a);
    for(auto& it:r.x) it*=b;
    return r;
}
static double recurse_mul(const symbolic_tensor& a,const symbolic_tensor& b,
    const std::string& sum,const vector<size_t>& sum_size,const vector<size_t>& as_sum,
    const vector<size_t>& bs_sum,size_t i,size_t ia,size_t ib)
{
    if(i>=sum.size()) return a.x[ia]*b.x[ib];
    double x=0;
    for(size_t s=0;s<sum_size[i];s++)
        x+=recurse_mul(a,b,sum,sum_size,as_sum,bs_sum,i+1,ia+s*as_sum[i],ib+s*bs_sum[i]);
    return x;
}
static void recurse_mul(symbolic_tensor& r,const symbolic_tensor& a,const symbolic_tensor& b,
    const std::string& sum,const vector<size_t>& sum_size,const vector<size_t>& as,const vector<size_t>& bs,
    const vector<size_t>& as_sum,const vector<size_t>& bs_sum,size_t i,size_t ir,size_t ia,size_t ib)
{
    if(i>=r.size.size()){
        r.x[ir]=recurse_mul(a,b,sum,sum_size,as_sum,bs_sum,0,ia,ib);
        return;}
    for(size_t s=0;s<r.size[i];s++)
        recurse_mul(r,a,b,sum,sum_size,as,bs,as_sum,bs_sum,i+1,ir+s*r.strides[i],ia+s*as[i],ib+s*bs[i]);
}
symbolic_tensor operator* (const symbolic_tensor& a,const symbolic_tensor& b)
{
    assert_indexed(a);
    assert_indexed(b);
    symbolic_tensor r;
    std::map<char,size_t> mp;
    size_t dup=-1;
    std::string sum;
    vector<size_t> sum_size;
    for(size_t i=0;i<a.indices.size();i++)
        mp[a.indices[i]]=i;
    for(size_t i=0;i<b.indices.size();i++){
        auto it=mp.find(b.indices[i]);
        if(it==mp.end()) mp[b.indices[i]]=i;
        else{
            assert(it->second!=dup);
            it->second=dup;
            sum+=b.indices[i];
            sum_size.push_back(b.size[i]);}}
    for(size_t i=0;i<a.indices.size();i++){
        size_t& x=mp[a.indices[i]];
        if(x==dup) continue;
        r.indices+=a.indices[i];
        r.size.push_back(a.size[i]);}
    for(size_t i=0;i<b.indices.size();i++){
        size_t& x=mp[b.indices[i]];
        if(x==dup) continue;
        r.indices+=b.indices[i];
        r.size.push_back(b.size[i]);}
    r.compute_strides();
    vector<size_t> as,bs,as_sum,bs_sum;
    reordered_strides(as,r.indices,a.strides,a.indices);
    reordered_strides(bs,r.indices,b.strides,b.indices);
    reordered_strides(as_sum,sum,a.strides,a.indices);
    reordered_strides(bs_sum,sum,b.strides,b.indices);
    recurse_mul(r,a,b,sum,sum_size,as,bs,as_sum,bs_sum,0,0,0,0);
    return r;
}
symbolic_tensor operator+ (const symbolic_tensor& a,double b)
{
    assert_indexed(a);
    assert(a.size.size()==0 && a.x.size()==1);
    symbolic_tensor r(a);
    r.x[0]+=b;
    return r;
}
symbolic_tensor operator/ (double b,const symbolic_tensor& a)
{
    assert_indexed(a);
    assert(a.size.size()==0 && a.x.size()==1);
    symbolic_tensor r(a);
    r.x[0]=b/r.x[0];
    return r;
}
static double recurse_contract(const symbolic_tensor& a,const std::string& sum,const vector<size_t>& sum_size,
    const vector<size_t>& as_sum,size_t i,size_t ia)
{
    if(i>=sum.size()) return a.x[ia];
    double x=0;
    for(size_t s=0;s<sum_size[i];s++)
        x+=recurse_contract(a,sum,sum_size,as_sum,i+1,ia+s*as_sum[i]);
    return x;
}
static void recurse_contract(symbolic_tensor& r,const symbolic_tensor& a,const std::string& sum,
    const vector<size_t>& sum_size,const vector<size_t>& as,const vector<size_t>& as_sum,size_t i,size_t ir,size_t ia)
{
    if(i>=r.size.size()){
        r.x[ir]=recurse_contract(a,sum,sum_size,as_sum,0,ia);
        return;}
    for(size_t s=0;s<r.size[i];s++)
        recurse_contract(r,a,sum,sum_size,as,as_sum,i+1,ir+s*r.strides[i],ia+s*as[i]);
}
symbolic_tensor symbolic_tensor::operator()(const std::string& new_indices) const
{
    assert_unindexed(*this);
    std::map<char,size_t> mp;
    size_t dup=-1;
    std::string sum,reduced_indices;
    vector<size_t> sum_size,sum_strides,new_size,as,as_sum;
    for(size_t i=0;i<new_indices.size();i++){
        assert(islower(new_indices[i]));
        auto it=mp.find(new_indices[i]);
        if(it==mp.end()) mp[new_indices[i]]=i;
        else{
            assert(it->second+1);
            assert(size[i]==size[it->second]);
            sum+=new_indices[i];
            sum_size.push_back(size[i]);
            sum_strides.push_back(strides[i]+strides[it->second]);
            it->second=dup;}}
    for(size_t i=0;i<new_indices.size();i++){
        size_t& x=mp[new_indices[i]];
        if(x==dup) continue;
        reduced_indices+=new_indices[i];
        new_size.push_back(size[i]);}
    if(reduced_indices==new_indices){
        symbolic_tensor r(*this);
        r.indices=new_indices;
        return r;}
    symbolic_tensor r(new_size);
    r.indices=reduced_indices;
    reordered_strides(as,r.indices,strides,new_indices);
    recurse_contract(r,*this,sum,sum_size,as,sum_strides,0,0,0);
    return r;
}
void symbolic_tensor::random(random_type& rand,double lo,double hi)
{
    assert_unindexed(*this);
    uniform_real_distribution<> dist(lo,hi);
    for(size_t i=0;i<x.size();i++)
        x[i]=dist(rand);
}
symbolic_tensor operator/ (const symbolic_tensor& a,const symbolic_tensor& b)
{
    assert_indexed(a);
    assert(b.size.size()==0);
    return a/b.x[0];
}
symbolic_tensor scalar_op(const symbolic_tensor& a,function<double(double)> f)
{
    assert(a.size.size()==0);
    symbolic_tensor r(a);
    r.x[0]=f(a.x[0]);
    return r;
}
