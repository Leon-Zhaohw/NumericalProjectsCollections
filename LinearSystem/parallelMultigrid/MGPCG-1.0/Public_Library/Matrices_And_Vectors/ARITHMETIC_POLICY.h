//#####################################################################
// Copyright 2007, Geoffrey Irving, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARITHMETIC_POLICY
//#####################################################################
#ifndef __ARITHMETIC_POLICY__
#define __ARITHMETIC_POLICY__

#include <Matrices_And_Vectors/SCALAR_POLICY.h>
#include <Matrices_And_Vectors/VECTOR_FORWARD.h>
#include <Utilities/TYPE_UTILITIES.h>
#include <boost/mpl/logical.hpp>
#include <boost/utility/enable_if.hpp>
namespace PhysBAM{

template<class T1,class T2,class ENABLER=void> struct SUM;
template<class T> struct SUM<T,T>{typedef T TYPE;};
template<class T,class TV> struct SUM<T,TV,typename boost::enable_if<boost::mpl::and_<IS_SCALAR<T>,boost::mpl::not_<IS_SCALAR<TV> > > >::type>{typedef TV TYPE;};
template<class T,class TV> struct SUM<TV,T,typename boost::enable_if<boost::mpl::and_<IS_SCALAR<T>,boost::mpl::not_<IS_SCALAR<TV> > > >::type>{typedef TV TYPE;};
template<> struct SUM<float,double>{typedef double TYPE;};
template<> struct SUM<int,double>{typedef double TYPE;};
template<> struct SUM<double,float>{typedef double TYPE;};
template<> struct SUM<double,int>{typedef double TYPE;};
template<> struct SUM<int,float>{typedef float TYPE;};
template<> struct SUM<float,int>{typedef float TYPE;};

template<class T1,class T2,class ENABLER=void> struct PRODUCT;
template<class T,class TV> struct PRODUCT<T,TV,typename boost::enable_if<IS_SCALAR<T> >::type>{typedef TV TYPE;};

template<class T,int m,int n> struct PRODUCT<MATRIX<T,m,n>,VECTOR<T,n> >{typedef VECTOR<T,m> TYPE;};
template<class T,int d> struct PRODUCT<DIAGONAL_MATRIX<T,d>,VECTOR<T,d> >{typedef VECTOR<T,d> TYPE;};
template<class T,int d> struct PRODUCT<SYMMETRIC_MATRIX<T,d>,VECTOR<T,d> >{typedef VECTOR<T,d> TYPE;};
template<class T,int d> struct PRODUCT<UPPER_TRIANGULAR_MATRIX<T,d>,VECTOR<T,d> >{typedef VECTOR<T,d> TYPE;};
template<class T,int d> struct PRODUCT<MATRIX_MXN<T>,VECTOR<T,d> >{typedef VECTOR_ND<T> TYPE;};

template<class T,int m,int n> struct PRODUCT<MATRIX<T,m,n>,VECTOR_ND<T> >{typedef VECTOR_ND<T> TYPE;};
template<class T,int d> struct PRODUCT<DIAGONAL_MATRIX<T,d>,VECTOR_ND<T> >{typedef VECTOR_ND<T> TYPE;};
template<class T,int d> struct PRODUCT<SYMMETRIC_MATRIX<T,d>,VECTOR_ND<T> >{typedef VECTOR_ND<T> TYPE;};
template<class T,int d> struct PRODUCT<UPPER_TRIANGULAR_MATRIX<T,d>,VECTOR_ND<T> >{typedef VECTOR_ND<T> TYPE;};
template<class T> struct PRODUCT<MATRIX_MXN<T>,VECTOR_ND<T> >{typedef VECTOR_ND<T> TYPE;};

template<class T,int m,int n,int k> struct PRODUCT<MATRIX<T,m,k>,MATRIX<T,k,n> >{typedef MATRIX<T,m,n> TYPE;};
template<class T,int m,int n> struct PRODUCT<MATRIX<T,m,n>,DIAGONAL_MATRIX<T,n> >{typedef MATRIX<T,m,n> TYPE;};
template<class T,int m,int n> struct PRODUCT<MATRIX<T,m,n>,SYMMETRIC_MATRIX<T,n> >{typedef MATRIX<T,m,n> TYPE;};
template<class T,int m,int n> struct PRODUCT<MATRIX<T,m,n>,UPPER_TRIANGULAR_MATRIX<T,n> >{typedef MATRIX<T,m,n> TYPE;};
template<class T,int m,int n> struct PRODUCT<MATRIX<T,m,n>,MATRIX_MXN<T> >{typedef MATRIX_MXN<T> TYPE;};

template<class T,int m,int n> struct PRODUCT<DIAGONAL_MATRIX<T,m>,MATRIX<T,m,n> >{typedef MATRIX<T,m,n> TYPE;};
template<class T,int d> struct PRODUCT<DIAGONAL_MATRIX<T,d>,DIAGONAL_MATRIX<T,d> >{typedef DIAGONAL_MATRIX<T,d> TYPE;};
template<class T,int d> struct PRODUCT<DIAGONAL_MATRIX<T,d>,SYMMETRIC_MATRIX<T,d> >{typedef MATRIX<T,d> TYPE;};
template<class T,int d> struct PRODUCT<DIAGONAL_MATRIX<T,d>,UPPER_TRIANGULAR_MATRIX<T,d> >{typedef UPPER_TRIANGULAR_MATRIX<T,d> TYPE;};
template<class T,int d> struct PRODUCT<DIAGONAL_MATRIX<T,d>,MATRIX_MXN<T> >{typedef MATRIX_MXN<T> TYPE;};

template<class T,int m,int n> struct PRODUCT<SYMMETRIC_MATRIX<T,m>,MATRIX<T,m,n> >{typedef MATRIX<T,m,n> TYPE;};
template<class T,int d> struct PRODUCT<SYMMETRIC_MATRIX<T,d>,DIAGONAL_MATRIX<T,d> >{typedef MATRIX<T,d> TYPE;};
template<class T,int d> struct PRODUCT<SYMMETRIC_MATRIX<T,d>,SYMMETRIC_MATRIX<T,d> >{typedef MATRIX<T,d> TYPE;};
template<class T,int d> struct PRODUCT<SYMMETRIC_MATRIX<T,d>,UPPER_TRIANGULAR_MATRIX<T,d> >{typedef MATRIX<T,d> TYPE;};
template<class T,int d> struct PRODUCT<SYMMETRIC_MATRIX<T,d>,MATRIX_MXN<T> >{typedef MATRIX_MXN<T> TYPE;};

template<class T,int m,int n> struct PRODUCT<UPPER_TRIANGULAR_MATRIX<T,m>,MATRIX<T,m,n> >{typedef MATRIX<T,m,n> TYPE;};
template<class T,int d> struct PRODUCT<UPPER_TRIANGULAR_MATRIX<T,d>,DIAGONAL_MATRIX<T,d> >{typedef UPPER_TRIANGULAR_MATRIX<T,d> TYPE;};
template<class T,int d> struct PRODUCT<UPPER_TRIANGULAR_MATRIX<T,d>,SYMMETRIC_MATRIX<T,d> >{typedef MATRIX<T,d> TYPE;};
template<class T,int d> struct PRODUCT<UPPER_TRIANGULAR_MATRIX<T,d>,UPPER_TRIANGULAR_MATRIX<T,d> >{typedef UPPER_TRIANGULAR_MATRIX<T,d> TYPE;};
template<class T,int d> struct PRODUCT<UPPER_TRIANGULAR_MATRIX<T,d>,MATRIX_MXN<T> >{typedef MATRIX_MXN<T> TYPE;};

template<class T,int m,int n> struct PRODUCT<MATRIX_MXN<T>,MATRIX<T,m,n> >{typedef MATRIX_MXN<T> TYPE;};
template<class T,int d> struct PRODUCT<MATRIX_MXN<T>,DIAGONAL_MATRIX<T,d> >{typedef MATRIX_MXN<T> TYPE;};
template<class T,int d> struct PRODUCT<MATRIX_MXN<T>,SYMMETRIC_MATRIX<T,d> >{typedef MATRIX_MXN<T> TYPE;};
template<class T,int d> struct PRODUCT<MATRIX_MXN<T>,UPPER_TRIANGULAR_MATRIX<T,d> >{typedef MATRIX_MXN<T> TYPE;};
template<class T> struct PRODUCT<MATRIX_MXN<T>,MATRIX_MXN<T> >{typedef MATRIX_MXN<T> TYPE;};

template<class TV> struct PRODUCT<FRAME<TV>,TV>{typedef TV TYPE;};

template<class TM> struct TRANSPOSE;
template<class T,int m,int n> struct TRANSPOSE<MATRIX<T,m,n> >{typedef MATRIX<T,n,m> TYPE;};
template<class T,int d> struct TRANSPOSE<SYMMETRIC_MATRIX<T,d> >{typedef SYMMETRIC_MATRIX<T,d> TYPE;};
template<class T,int d> struct TRANSPOSE<DIAGONAL_MATRIX<T,d> >{typedef DIAGONAL_MATRIX<T,d> TYPE;};
template<class T,int d> struct TRANSPOSE<UPPER_TRIANGULAR_MATRIX<T,d> >{typedef MATRIX<T,d> TYPE;};
template<class T> struct TRANSPOSE<MATRIX_MXN<T> >{typedef MATRIX_MXN<T> TYPE;};

template<class T1,class T2> struct PRODUCT_TRANSPOSE{typedef typename PRODUCT<T1,typename TRANSPOSE<T2>::TYPE>::TYPE TYPE;};
template<class T1,class T2> struct TRANSPOSE_PRODUCT{typedef typename PRODUCT<typename TRANSPOSE<T1>::TYPE,T2>::TYPE TYPE;};

}
#endif
