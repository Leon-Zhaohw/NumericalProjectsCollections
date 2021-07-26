//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TYPE_UTILITIES
//#####################################################################
#ifndef __TYPE_UTILITIES__
#define __TYPE_UTILITIES__

#include <boost/mpl/and.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/less_equal.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <boost/type_traits/is_class.hpp>
#include <boost/type_traits/is_empty.hpp>
#include <boost/type_traits/is_enum.hpp>
#include <boost/type_traits/is_fundamental.hpp>
#include <Arrays/ARRAYS_FORWARD.h>
namespace PhysBAM{

namespace mpl=boost::mpl;

template<bool b,class T1,class T2> struct IF{typedef T1 TYPE;};
template<class T1,class T2> struct IF<false,T1,T2>{typedef T2 TYPE;};

template<bool b> struct NOT:public mpl::false_{};
template<> struct NOT<false>:public mpl::true_{};

template<bool b1,bool b2,bool b3=true,bool b4=true> struct AND:public mpl::false_{};
template<> struct AND<true,true,true,true>:public mpl::true_{};

template<bool b1,bool b2,bool b3=false,bool b4=false> struct OR:public mpl::true_{};
template<> struct OR<false,false,false,false>:public mpl::false_{};

template<class T1,class T2,class T3=mpl::true_,class T4=mpl::true_> struct LAZY_AND:public boost::mpl::and_<T1,T2,T3,T4>{};
template<class T1,class T2,class T3=mpl::false_,class T4=mpl::false_> struct LAZY_OR:public boost::mpl::or_<T1,T2,T3,T4>{};

template<int i,int j> struct INTS_EQUAL:public mpl::false_{};
template<int i> struct INTS_EQUAL<i,i>:public mpl::true_{};

template<class T1,class T2> struct IS_SAME:public mpl::false_{};
template<class T> struct IS_SAME<T,T>:public mpl::true_{};

template<class T1,class T2> struct ASSERT_SAME_HELPER;
template<class T> struct ASSERT_SAME_HELPER<T,T>{};

#define STATIC_ASSERT_SAME(T1,T2) \
   typedef ::boost::static_assert_test< \
      sizeof(::PhysBAM::ASSERT_SAME_HELPER<T1,T2>)> \
         BOOST_JOIN(boost_static_assert_typedef_, __LINE__)

template<class T> struct IS_CONST:public mpl::false_{};
template<class T> struct IS_CONST<const T>:public mpl::true_{};

template<class T> struct REMOVE_CONST{typedef T TYPE;};
template<class T> struct REMOVE_CONST<const T>{typedef T TYPE;};

template<class T> struct IS_POINTER:public mpl::false_{};
template<class T> struct IS_POINTER<T*>:public mpl::true_{};

template<class T> struct REMOVE_POINTER;
template<class T> struct REMOVE_POINTER<T*>{typedef T TYPE;};

template<class T> struct IS_MEMBER_POINTER:public mpl::false_{};
template<class T,class TF> struct IS_MEMBER_POINTER<TF T::*>:public mpl::true_{};

template<class T> struct IS_REFERENCE:public mpl::false_{};
template<class T> struct IS_REFERENCE<T&>:public mpl::true_{};

template<class T> struct ADD_REFERENCE{typedef T& TYPE;};
template<class T> struct ADD_REFERENCE<T&>{typedef T& TYPE;};

template<class T> struct REMOVE_REFERENCE{typedef T TYPE;};
template<class T> struct REMOVE_REFERENCE<T&>{typedef T TYPE;};

template<class T1,class T2=void,class T3=void,class T4=void> struct FIRST{typedef T1 TYPE;};

template<class T> struct IS_FLOAT_OR_DOUBLE:public LAZY_OR<IS_SAME<T,float>,IS_SAME<T,double> >{};

// defer to boost for complicated traits
template<class T> struct IS_CLASS:public boost::is_class<T>{};
template<class T> struct IS_EMPTY:public boost::is_empty<T>{};
template<class T> struct IS_ENUM:public boost::is_enum<T>{};
template<class T> struct IS_INTEGRAL:public boost::is_integral<T>{};

template<class T_ARRAY,class ENABLER=void> struct IS_ARRAY:public mpl::false_{};
template<class T_ARRAY> struct IS_ARRAY<const T_ARRAY>:public IS_ARRAY<T_ARRAY>{};

template<class T_ARRAY,class ENABLER=void> struct IS_ARRAY_VIEW{};
template<class T_ARRAY> struct IS_ARRAY_VIEW<const T_ARRAY>:public IS_ARRAY_VIEW<T_ARRAY>{};

template<class T> struct HAS_CHEAP_COPY:public LAZY_OR<boost::is_fundamental<T>,IS_ENUM<T>,IS_ARRAY_VIEW<T> >{};

// boost::is_base_of<T,T> is false_ in boost 1.33 and true_ in 1.34, so we need this wrapper to ensure the 1.34 behavior
template<class B,class D> struct IS_BASE_OF:public mpl::or_<IS_SAME<B,D>,boost::is_base_of<B,D> >{};

}
#endif
