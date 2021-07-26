//#####################################################################
// Copyright 2006-2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SCALAR_POLICY
//#####################################################################
#ifndef __SCALAR_POLICY__
#define __SCALAR_POLICY__

#include <Arrays/ARRAYS_FORWARD.h>
#include <Matrices_And_Vectors/MATRIX_FORWARD.h>
#include <Matrices_And_Vectors/VECTOR_FORWARD.h>
#include <Utilities/TYPE_UTILITIES.h>
#include <boost/utility/enable_if.hpp>
namespace PhysBAM{

template<class T> struct IS_SCALAR:public mpl::false_{};
template<> struct IS_SCALAR<int>:public mpl::true_{};
template<> struct IS_SCALAR<unsigned int>:public mpl::true_{};
template<> struct IS_SCALAR<float>:public mpl::true_{};
template<> struct IS_SCALAR<double>:public mpl::true_{};

template<class T> struct IS_SCALAR_BLOCK:public IS_SCALAR<T>{}; // true if memory layout is contiguous array of scalars
template<class T> struct IS_SCALAR_VECTOR_SPACE:public IS_SCALAR<T>{}; // true if we can compute vector space operations on the underlying array of scalars

template<class T,class ENABLER=void> struct SCALAR_POLICY{typedef struct UNUSABLE{} TYPE;};
template<class T> struct SCALAR_POLICY<T,typename boost::enable_if<IS_SCALAR<T> >::type>{typedef T TYPE;};
template<class T> struct SCALAR_POLICY<T,typename boost::enable_if<typename FIRST<mpl::true_,typename T::SCALAR>::TYPE>::type>{typedef typename T::SCALAR TYPE;};

}
#endif
