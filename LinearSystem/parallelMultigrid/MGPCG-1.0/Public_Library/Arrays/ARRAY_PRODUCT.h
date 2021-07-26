//#####################################################################
// Copyright 2008, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_PRODUCT
//#####################################################################
#ifndef __ARRAY_PRODUCT__
#define __ARRAY_PRODUCT__

#include <Arrays/ARRAY_EXPRESSION.h>
#include <Matrices_And_Vectors/ARITHMETIC_POLICY.h>
#include <cassert>
namespace PhysBAM{

template<class T_ARRAY1,class T_ARRAY2> class ARRAY_PRODUCT;
template<class T_ARRAY1,class T_ARRAY2> struct IS_ARRAY<ARRAY_PRODUCT<T_ARRAY1,T_ARRAY2> >:public mpl::true_{};
template<class T_ARRAY1,class T_ARRAY2> struct IS_ARRAY_VIEW<ARRAY_PRODUCT<T_ARRAY1,T_ARRAY2> >:public mpl::true_{};

template<class T_ARRAY1,class T_ARRAY2>
class ARRAY_PRODUCT:public ARRAY_EXPRESSION<typename PRODUCT<typename T_ARRAY1::ELEMENT,typename T_ARRAY2::ELEMENT>::TYPE,ARRAY_PRODUCT<T_ARRAY1,T_ARRAY2> >
{
    typedef typename T_ARRAY1::ELEMENT T1;typedef typename T_ARRAY2::ELEMENT T2;
    typedef typename IF<IS_ARRAY_VIEW<T_ARRAY1>::value,const T_ARRAY1,const T_ARRAY1&>::TYPE T_ARRAY1_VIEW; // if it's an array view we can copy it, otherwise store a reference
    typedef typename IF<IS_ARRAY_VIEW<T_ARRAY2>::value,const T_ARRAY2,const T_ARRAY2&>::TYPE T_ARRAY2_VIEW;
    typedef typename PRODUCT<T1,T2>::TYPE T_PRODUCT;
public:
    typedef T_PRODUCT ELEMENT;

    T_ARRAY1_VIEW array1;
    T_ARRAY2_VIEW array2;

    ARRAY_PRODUCT(const T_ARRAY1& array1,const T_ARRAY2& array2)
        :array1(array1),array2(array2)
    {}

    int Size() const
    {int size=array1.Size();assert(size==array2.Size());return size;}

    const T_PRODUCT operator()(const int i) const
    {return array1(i)*array2(i);}

//#####################################################################
};

template<class T1,class T2,class T_ARRAY1,class T_ARRAY2> ARRAY_PRODUCT<T_ARRAY1,T_ARRAY2>
operator*(const ARRAY_BASE<T1,T_ARRAY1>& array1,const ARRAY_BASE<T2,T_ARRAY2>& array2)
{return ARRAY_PRODUCT<T_ARRAY1,T_ARRAY2>(array1.Derived(),array2.Derived());}

}
#endif
