//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_NEGATION
//#####################################################################
#ifndef __ARRAY_NEGATION__
#define __ARRAY_NEGATION__

#include <Arrays/ARRAY_EXPRESSION.h>
#include <Matrices_And_Vectors/ARITHMETIC_POLICY.h>
namespace PhysBAM{

template<class T_ARRAY> class ARRAY_NEGATION;
template<class T_ARRAY> struct IS_ARRAY<ARRAY_NEGATION<T_ARRAY> >:public mpl::true_{};
template<class T_ARRAY> struct IS_ARRAY_VIEW<ARRAY_NEGATION<T_ARRAY> >:public mpl::true_{};

template<class T_ARRAY>
class ARRAY_NEGATION:public ARRAY_EXPRESSION<typename T_ARRAY::ELEMENT,ARRAY_NEGATION<T_ARRAY>,typename T_ARRAY::INDEX>
{
    typedef typename T_ARRAY::ELEMENT T;
    typedef typename IF<IS_ARRAY_VIEW<T_ARRAY>::value,const T_ARRAY,const T_ARRAY&>::TYPE T_ARRAY_VIEW; // if it's an array view we can copy it, otherwise store a reference
public:
    typedef T ELEMENT;typedef typename T_ARRAY::INDEX INDEX;

    T_ARRAY_VIEW array;

    explicit ARRAY_NEGATION(const T_ARRAY& array)
        :array(array)
    {}

    INDEX Size() const
    {return array.Size();}

    const T operator()(const INDEX i) const
    {return -array(i);}

//#####################################################################
};

template<class T,class T_ARRAY> ARRAY_NEGATION<T_ARRAY>
operator-(const ARRAY_BASE<T,T_ARRAY,typename T_ARRAY::INDEX>& array)
{return ARRAY_NEGATION<T_ARRAY>(array.Derived());}

}
#endif
